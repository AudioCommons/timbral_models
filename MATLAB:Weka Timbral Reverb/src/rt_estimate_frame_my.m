function [ RT, par, RT_pre ] = rt_estimate_frame_my( frame, par )
% performs an efficient blind estimation of the reverberation time (RT) for frame-wise
% processing based on Laplacian distribution.
%
% INPUT
% frame : (time-domain) segment with reverberant speech
% par   : struct with all parameters and buffers created by the function
%         init_binaural_speech_enhancement_e.m
%
% OUTPUT
% RT    : estimated RT
% par   : struct with updated buffers to enable a frame-wise processing
% RT_pre: raw RT estimate (for debugging and analysis of the algorithm)
%
% Reference:  Löllmann, H. W., Jeub, M., Yilmaz, E., and Vary, P.:
%            An Improved Algorithm for Blind Reverberation Time Estimation,a
%            International Workshop on Acoustic Echo and Noise Control (IWAENC), Tel Aviv, Israel, Aug. 2010.

%            Tariqullah Jan and Wenwu Wang:
%            Blind reverberation time estimation based on Laplacian distribution
%            European Signal Processing Conference (EUSIPCO), 2012.
%
%
% The codes were adapted based on the original codes by Heinrich Loellmann, IND, RWTH Aachen
%
% Authors: Tariqullah Jan, moderated by Wenwu Wang, University of Surrey (2012)


[ M , N ] = size( frame );

if max( [ M, N ] ) < par.N 
    
    error('Incorrect input frame length for RT estimation!')

elseif min( [ M, N ] ) ~= 1
    
    error('No matrix for RT estimation allowed!')
end

if M > N 
    
    frame = frame.'; % ensures a column vector
end

cnt = 0;     % sub-frame counter for pre-selection of possible sound decay\
RTml = -1;   % default RT estimate (-1 indicates no new RT estimate)


% calculate variance, minimum and maximum of first sub-frame
seg = frame( 1 : par.N_sub );

var_pre = var( seg );
min_pre = min( seg );
max_pre = max( seg );

for k = 2 : par.nos_max,
    
    % calculate variance, minimum and maximum of succeding sub-frame
    seg = frame( 1+(k-1)*par.N_sub : k*par.N_sub );
    var_cur = var( seg );
    max_cur = max( seg );
    min_cur = min( seg );
    
    % -- Pre-Selection of suitable speech decays --------------------
    
    if (var_pre > var_cur) && (max_pre > max_cur) && (min_pre < min_cur)
        % if variance, maximum decraease and minimum increase
        % => possible sound decay detected
        
        cnt = cnt + 1;
        
        % current values becomes previous values
        var_pre = var_cur;
        max_pre = max_cur;
        min_pre = min_cur;
        
    else
        
        if cnt >= par.nos_min % minimum length for assumed sound decay achieved?
            
            % -- Maximum Likelihood (ML) Estimation of the RT
            RTml = max_loglf( frame(1:cnt*par.N_sub), par.a, par.Tquant );
            
        end
        
        break
        
    end
    
    if k == par.nos_max % maximum frame length achieved?
        
        RTml = max_loglf( frame(1:cnt*par.N_sub), par.a, par.Tquant ); 
        
    end
    
end % eof sub-frame loop


if RTml >= 0  % new ML estimate calculated
    
    % apply order statistics to reduce outliers
    par.hist_counter = par.hist_counter+1;
   
    for i = 1: par.no_bins,
        
        % find index corresponding to the ML estimate
        if  ( RTml >= par.hist_limits(i) ) && ( RTml <= par.hist_limits(i+1) )
            
            index = i;
            break
        end        
    end
    
    % update histogram with ML estimates for the RT
    par.hist_rt( index ) = par.hist_rt( index ) + 1;
    
    if par.hist_counter > par.buffer_size +1
        % remove old values from histogram
        par.hist_rt( par.buffer( 1 ) ) = par.hist_rt( par.buffer( 1 ) ) - 1;
    end
    
    par.buffer = [ par.buffer(2:end), index ]; % update buffer with indices
    [ dummy, idx ] = max( par.hist_rt );       % find index for maximum of the histogram
     
    par.RT_raw = par.Tquant( idx );   % map index to RT value
    
end

% final RT estimate obtained by recursive smoothing
RT = par.alpha * par.RT_last + (1-par.alpha) * par.RT_raw;
par.RT_last = RT;
    
RT_pre = RTml;    % intermediate ML estimate for later analysis
     
return


function [ ML, ll ] = max_loglf( h, a, Tquant )

% [ ML, ll ] = max_loglf( h, a, Tquant )
%
% returns the maximum of the log-likelihood (LL) function and the LL
% function itself for a finite set of decay rates
%
% INPUT
% h: input frame
% a: finite set of values for which the max. should be found
% T: corresponding RT values for vector a
%
% OUTPUT
% ML : ML estimate for the RT
% ll : underlying LL-function

N = length(h);
n = (0:N-1); % indices for input vector
ll = zeros(length(a),1);

%h_square = (h.^2).';
h_square = (h).';
% for i=1:length(a),
%     for j=1:N
%         if (a(i)^j*h_square(j))>0
%             temp1(j) = a(i)^j*h_square(j);
%         elseif (a(i)^j*h_square(j))<0
%             temp2(j) = a(i)^j*h_square(j);
%         end
%     end
%     sum1 = sum(temp1);
%     sum2 = sum(temp2);
%     
% sigma( i ) = sqrt(1/N*(sum1-sum2));
% %ll(i) = (2^(-N/2)/sigma(i)^N)*exp(-N*sigma(i));
% ll(i) = (-N/2)*log(2)-N*log(sigma(i))-N*sigma(i);
% end


for i=1:length(a),
%     for j=1:N
%         temp(j) = abs(a(i)^j*h_square(j));
%     end
   % sum1 = sum(abs(temp));
   % sum1 = sum(temp);
    sum1  = ( a(i).^(-n) )*abs(h_square) ;
    sum2 = sum(abs(h_square));
sigma = (1/N)*sum1;
ll(i) = -N*log(2)-N*log(sigma)-sum(log(a(i).^n))-(1/sigma)*sum1;
%ll( i ) = -N/2*( (N-1)*log( a(i) ) + log( 2*pi/N * sigma ) + 1 );
%ll(i) = (2^(-N/2)/sigma(i)^N)*exp(-N*sigma(i));
%ll(i) = (-N/2)*log(2)-N*log(sum1)-N*sigma;
end



% for i=1:length(a),
%     
%     Sum  = ( a(i).^(-2*n) ) * h_square ;
%     
%     if Sum < 1e-12
%         ll( i ) = -inf;
%     else
%         ll( i ) = -N/2*( (N-1)*log( a(i) ) + log( 2*pi/N * Sum ) + 1 );
%     end
%     
% end

[ dummy, idx ] = max( ll ); % maximum of the log-likelihood function
ML = Tquant( idx );        % corresponding ML estimate for the RT

return
