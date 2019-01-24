classdef spectralFeaturesProc < Processor
%SPECTRALFEATURESPROC Spectral features processor.
%  This processor computes the following 14 spectral features that summarize 
%  the spectral content of the ratemap representation across auditory filters for
%  individual time frames.
%             'centroid'     : Spectral centroid [1]
%             'crest'        : Spectral crest measure [1]
%             'spread'       : Spectral spread 
%             'entropy'      : Spectral entropy [2]
%             'brightness'   : Spectral brightness [1]
%             'hfc'          : Spectral high-frequency content [3]
%             'decrease'     : Spectral decrease [1]
%             'flatness'     : Spectral flatness [1]
%             'flux'         : Spectral flux [4]
%             'kurtosis'     : Spectral kurtosis [4]
%             'skewness'     : Spectral skewness [4]
%             'irregularity' : Spectral irregularity [3]
%             'rolloff'      : Spectral rolloff [1]
%             'variation'    : Spectral variation [1]
%
%   SPECTRALFEATURESPROC properties:
%        requests        - Cell array of requested spectral features
%        cfHz            - Row vector of audio center frequencies
%
%   See also: Processor, ratemapProc
%
%   Reference:
%   [1] Peeters, G., Giordano, B. L., Susini, P., Misdariis, N., and 
%       McAdams, S. (2011), "The timbre toolbox: Extracting audio descriptors 
%       from musical signals." Journal of the Acoustical Society of America 
%       130(5), pp. 2902?2916.
%   [2] Misra, H., Ikbal, S., Bourlard, H., and Hermansky, H. (2004), 
%       "Spectral entropy based feature for robust ASR," in Proceedings of 
%       the IEEE International Conference on Acoustics, Speech and Signal 
%       Processing (ICASSP), pp. 193?196.
%   [3] Jensen, K. and Andersen, T. H. (2004), "Real-time beat estimation 
%       using feature extraction," in Computer Music Modeling and Retrieval, 
%       edited by U. K. Wiil, Springer, Berlin?Heidelberg, Lecture Notes in 
%       Computer Science, pp. 13?22.
%   [4] Lerch, A. (2012), An Introduction to Audio Content Analysis: 
%       Applications in Signal Processing and Music Informatics, 
%       John Wiley & Sons, Hoboken, NJ, USA.
    
    properties (Dependent = true)
        requests        % Cell array of requested spectral features
    end
    
    properties (SetAccess = private)
        cfHz            % Row vector of audio center frequencies 
    end
    
    properties (GetAccess = private)
        eps             % Factor used to prevent division by 0 (hard-coded)
        brCF           % Cutoff frequency for brightness feature
        flux_buffer     % Buffered last frame of previous chunk for spectral flux
        var_buffer      % Buffered last frame for spectral variation
        ro_eps          % Epsilon value for spectral rolloff (hard-coded)
        roPerc         % threshold value for spectral rolloff
        bUseInterp      % Flag indicating use of interpolation for spectral rolloff (hard-coded)
    end
    
    methods
        function pObj = spectralFeaturesProc(fs,parObj)
		%spectralFeaturesProc   Construct a spectral features extraction processor
        %
        % USAGE:
        %   pObj = spectralFeaturesProc(fs, parObj)
        %
        % INPUT ARGUMENTS:
        %     fs : Input sampling frequency (Hz)
        % parObj : Parameter object instance
        %
        % OUTPUT ARGUMENTS:
        %   pObj : Processor instance
        %
        % NOTE: Parameter object instance, parObj, can be generated using genParStruct.m
        % User-controllable parameters for this processor and their default values can be
        % found by browsing the script parameterHelper.m
        %
        % See also: genParStruct, parameterHelper, Processor
            
            % Checking input parameter
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end
            
            % Call superconstructor
            pObj = pObj@Processor(fs,fs,'spectralFeaturesProc',parObj);
            
            if nargin>0
                
                pObj.flux_buffer = [];
                pObj.var_buffer = [];
            
            end
        end
            
        function out = processChunk(pObj,in)
            %processChunk       Apply the processor to a new chunk of input
            %                   signal
            %
            %USAGE
            %   out = pObj.processChunk(in)
            %
            %INPUT ARGUMENT
            %    in : New chunk of input data
            %
            %OUTPUT ARGUMENT
            %   out : Corresponding output
            %
            %TODO: Spectral spread is dependent on centroid, don't compute
            %it twice
            
            % Number of spectral features to extract
            n_feat = size(pObj.requests,2);
            
            % Output initialization
            out = zeros(size(in,1),n_feat);

            if isempty( in ) %whyever
                return;
            end
            
            % Size of input chunk
            [nFrames, nFreq] = size(in);
                       
            % Main loop on all the features
            for ii = 1:n_feat
                
                % Switch among the requested features
                switch pObj.requests{ii}
                    
                    case 'centroid'     % Spectral centroid
                        % Spectral center of gravity of the spectrum
                        out(:,ii) = sum(repmat(pObj.cfHz,[nFrames 1]).*in,2)./(sum(in,2)+pObj.eps);
                        
                        % Normalize centroid to "nyquist" frequency channel
                        out(:,ii) = out(:,ii) / pObj.cfHz(end);
                        
                    case 'crest'        % Spectral crest
                        % Ratio of maximum to average in every frame
                        out(:,ii) = max(in,[],2)./(mean(in,2)+pObj.eps);
                        
                    case 'decrease'     % Spectral decrease
                        
                        % Vector of inverse frequency bin index (corrected for 0)
                        kinv = 1./[1 ;(1:nFreq-1)'];
                        
                        % Spectral decrease
                        out(:,ii) = ((in-repmat(in(:,1),[1 nFreq]))*kinv)./(sum(in(:,2:end),2)+pObj.eps);
                        
                    case 'spread'       % Spectral spread (bandwidth)
                        % Average deviation to the centroid weigthed by
                        % amplitude
                        
                        % Dependent on the spectral centroid
                        centroid = sum(repmat(pObj.cfHz,[nFrames 1]).*in,2)./(sum(in,2)+pObj.eps);
                    
                        % Temporary nominator
                        nom = (repmat(pObj.cfHz,[nFrames 1]) - repmat(centroid,[1 nFreq])).^2 .* in;
                        
                        % Spectrum bandwidth
                        out(:,ii) = sqrt(sum(nom,2)./(sum(in,2)+pObj.eps));
                        
                        % Normalize spread to "nyquist" frequency channel
                        out(:,ii) = out(:,ii) / pObj.cfHz(end);
                        
                    case 'brightness'   % Spectral brightness
                        % Ratio of energy above cutoff to total energy in
                        % each frame
                        out(:,ii) = sum(in(:,pObj.cfHz>pObj.brCF),2)./(sum(in,2)+pObj.eps);
                        
                    case 'hfc'          % Spectral high frequency content
                        % Average channel amplitude weighted by squared
                        % channel center frequency across channels
                        out(:,ii) = sum(repmat(pObj.cfHz.^2,[nFrames 1]).*in,2)./(sum(in,2)+pObj.eps);
                        
                    case 'entropy'      % Spectral entropy
                        
                        % Normalized spectrum
                        specN = in./(repmat(sum(in,2),[1 nFreq])+pObj.eps);
                        % Entropy
                        out(:,ii) = -sum(specN .* log(specN+pObj.eps),2)./log(nFreq);
                        
                    case 'flatness'     % Spectral flatness (SFM)
                        % Ratio of geometric mean to arithmetic mean across
                        % frequencies
                        out(:,ii) = exp(mean(log(in),2))./(mean(in,2)+pObj.eps);
                        
                    case 'flux'         % Spectral flux
                        
                        % Compressed power spectrum
                        pSpec = 10*log10(in+pObj.eps);
                        
                        % If the buffer is empty, use the first frame of
                        % the current input
                        if isempty(pObj.flux_buffer)
                            pObj.flux_buffer = pSpec(1,:);
                        end
                        
                        % Compute delta across frames, including buffer
                        deltaSpec = diff([pObj.flux_buffer; pSpec],1,1);
                        
                        % Take the norm across frequency
                        out(:,ii) = sqrt(mean(power(deltaSpec,2),2));
                        
                        % Update the buffer
                        if ~isempty( pSpec )
                            pObj.flux_buffer = pSpec(end,:);
                        end
                        
                    case 'kurtosis'
                        
                        % Mean and standard deviation across frequency
                        mu_x  = mean(abs(in),2);
                        std_x = std(abs(in),0,2);
                        
                        % Remove mean from input
                        X = in - repmat(mu_x,[1 nFreq]);
                        
                        % Excess kurtosis
                        out(:,ii) = mean((X.^4)./(repmat(std_x + pObj.eps, [1 nFreq]).^4),2)-3;
                        
                    case 'skewness'
                        
                        % Mean and standard deviation across frequency
                        mu_x  = mean(abs(in),2);
                        std_x = std(abs(in),0,2);
                        
                        % Remove mean from input
                        X = in - repmat(mu_x,[1 nFreq]);
                        
                        % Kurtosis
                        out(:,ii) = mean((X.^3)./(repmat(std_x + pObj.eps, [1 nFreq]).^3),2);
                        
                    case 'irregularity'
                        
                        % Compressed spectrum
                        pSpec = 10*log10(in+pObj.eps);
                        
                        % 2-Norm of spectrum difference across frequency
                        out(:,ii) = sqrt(mean(power(diff(pSpec,1,2),2),2));
                        
                    case 'rolloff'
                        % Extrapolated frequency threshold at each frame
                        % for which pObj.ro_thresh % of the energy is at 
                        % lower frequencies and the remaining above.
                        
                        % Spectral energy across frequencies multiplied by threshold parameter
                        spec_sum_thres = pObj.roPerc * sum(in,2);
                        % Cumulative sum (+ epsilon ensure that cumsum increases monotonically)
                        spec_cumsum = cumsum(in + pObj.ro_eps,2);
                        
                        % Loop over number of frames
                        for jj = 1 : nFrames

                            % Use interpolation
                            if pObj.bUseInterp
                                if spec_sum_thres(jj) > 0
                                    % Detect spectral roll-off
                                    out(jj,ii) = interp1(spec_cumsum(jj,:),pObj.cfHz,spec_sum_thres(jj),'linear','extrap');
                                end
                            else
                                % The feature range of this code is limited to the vector fHz.

                                % Detect spectral roll-off
                                r = find(spec_cumsum(jj,:) > spec_sum_thres(jj),1);

                                % If roll-off is found ...
                                if ~isempty(r)
                                    % Get frequency bin
                                    out(jj,ii) = pObj.cfHz(r(1));
                                end
                            end
                        end
                        
                        % Normalize rolloff to "nyquist" frequency channel
                        out(:,ii) = out(:,ii) / pObj.cfHz(end);
                        
                    case 'variation'
                        
                        % Initialize the buffer if empty
                        if isempty(pObj.var_buffer)
                            pObj.var_buffer = in(1,:);
                        end
                        
                        % Spectrum "shifted" one frame in the past
                        past_spec = [pObj.var_buffer;in(1:end-1,:)];
                        
                        % Cross-product
                        xprod = sum(in .* past_spec,2);
                        % Auto-product
                        aprod = sqrt(sum(in.^2,2)) .* sqrt(sum(past_spec.^2,2));
                        
                        % Noralized cross-correlation
                        out(:,ii) = 1-(xprod./(aprod+pObj.eps));
                        
                        % Update the buffer
                        pObj.var_buffer = in(end,:);
                        
                    otherwise
                        % This should NEVER be reached in a practical case
                        error('Invalid request name')
                end
                
            end
            
            
        end
            
        function reset(pObj)
            %reset      Reset the internal states of the processor
            %
            %USAGE
            %   pObj.reset
            %
            %INPUT ARGUMENT
            % pObj : Spectral features processor instance
            
            % Reset the two buffers
            pObj.flux_buffer = [];
            pObj.var_buffer = [];
            
        end
        
        function output = instantiateOutput(pObj,dObj)
            %INSTANTIATEOUTPUT  Instantiate the output signal for this processor
            %
            %NB: This method is overloaded here from the master Processor class, as
            %feature signals need additional input argument to construct
            
            featureNames = pObj.requests;
            
            sig = feval(pObj.getProcessorInfo.outputType, ...
                        pObj, ...
                        dObj.bufferSize_s, ...
                        pObj.Channel,...
                        [],...
                        featureNames);
            
            dObj.addSignal(sig);
            
            output = {sig};
            
        end
    end
    
    methods (Access=protected)
       
        function verifyParameters(pObj)
            
            % Check request validity...
            
            available = {'all' 'centroid' 'crest' 'spread' 'entropy' ...
                'brightness' 'hfc' 'decrease' 'flatness' 'flux' ... 
                'kurtosis' 'skewness' 'irregularity' 'rolloff' ...
                'variation'};
            
            requests = pObj.parameters.map('sf_requests');
            
            % Check if a single request was given but not as a cell array
            if ischar(requests)
                requests = {requests};
            end
            
            % Check that requests is a cell array of strings
            if ~iscellstr(requests)
                error('Requests for spectral features processor should be provided as a cell array of strings.')
            end
            
            % Check for typos/incorrect feature name
            if ~isequal(union(available,requests,'stable'),available)
                error(['Incorrect request name. Valid names are as follow: '...
                    strjoin(available,', ')])
            end
            
            % Change the requests to actual names if 'all' was requested
            if ismember('all',requests)
                requests = setdiff(available,{'all'},'stable');
            end
            
            % Check if provided brightness cutoff frequency is in a valid
            % range
            if ~isempty(pObj.cfHz)
                if (pObj.parameters.map('sf_brCF') < pObj.cfHz(1) || ...
                        pObj.parameters.map('sf_brCF') > pObj.cfHz(end)) && ...
                        ismember('brightness',requests)
                    error('Brightness cutoff frequency should be in Nyquist range')
                end
            end
            
            pObj.parameters.map('sf_requests') = requests;
            
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
            
            % Hard-coded properties (for the moment)
            pObj.eps = 1E-15;
            pObj.ro_eps = 1E-10;
            pObj.bUseInterp = true;

            % Access center frequencies for computing statistics
            pObj.cfHz = pObj.getDependentParameter('fb_cfHz');
            
            % Read private properties from parameters
            pObj.brCF = pObj.parameters.map('sf_brCF');
            pObj.roPerc = pObj.parameters.map('sf_roPerc');
            
        end
        
    end
    
    % "Getter" methods
    methods
        function requests = get.requests(pObj)
            requests = pObj.parameters.map('sf_requests');
        end
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'ratemap';
        end
        
        function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ...
            %                           gammatoneProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'sf_requests',...
                    'sf_brCF',...
                    'sf_roPerc'};
            
            descriptions = {['List (cell array) of requested spectral features, type ' ...
                '''help SpectralFeaturesProc'' for a list'],...
                    'Cutoff frequency for brightness computation',...
                    'Threshold (re. 1) for spectral rolloff computation'};
            
            defaultValues = {'all',...
                            1500,...
                            0.8};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Spectral features extractor';
            pInfo.label = 'Spectral features';
            pInfo.requestName = 'spectralFeatures';
            pInfo.requestLabel = 'Spectral features';
            pInfo.outputType = 'FeatureSignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
end