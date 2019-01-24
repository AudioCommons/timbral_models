function [ rt_est, par ]  = RT_estimation_my( y, fs )

% performs blind RT estimation
%
% INPUT
% y  : reverberant speech
% fs : sampling frequency
%
% OUTPUT
% rt_est: estimated RT
% par   : struct with parameters used to execute the function
%         rt_estimate_frame_my.m
%
% Codes were adapted from the original codes by Heinrich Loellmann, IND, RWTH Aachen
%
% Authors: Tariqullah Jan, moderated by Wenwu Wang, University of Surrey (2012)


% Initialization
% ---------------------------------------------

par = init_rt_estimate_e( fs ); % struct with all parameters and buffers for frame-wise processing
BL = par.N*par.down;   % to simplify notation

Ly = length(y);
rt_est = zeros( 1, round( Ly/par.N_shift ) );

% frame-wise processing in the time-domain
% ---------------------------------------------

k = 0;
for n = 1 : par.N_shift: Ly-BL+1,
    
    k = k+1;          % frame counter
    ind = n:n+BL-1;   % indices of current frame
        
    % Actual RT estimation    
    [ RT, par, finalrt ] = rt_estimate_frame_my( y( ind(1:par.down:end) ), par );
    
    rt_est(k) = RT;   % store estimated value
    RT_final(k) = finalrt;
end
RT_final(find(RT_final==-1))=0;
aaa = RT_final(find(RT_final~=0));
NN = size(aaa,2);
 for i=2:NN
%      if i==1
%      RT_temp = (1-0.995) * max(aaa);
%      else
     RT_temp_new(i) = 0.49*aaa(i-1) + (1-0.49)*max(aaa);
     RT_temp=RT_temp_new(i);
    % end
 end
 
% NN = size(aaa,2);
% for i=1:NN
%     temp(i)=(1/NN)*aaa(i);
% end
% RTfinal_value = sum(temp);
RTfinal_value = min(aaa);
%RTfinal_value = mean(aaa);
rt_est = rt_est(1:k); 
