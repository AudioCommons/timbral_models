classdef drnlProc < Processor
%DRNLPROC Dual-Resonance Non-Linear auditory filterbank processor.
%   The DRNL filterbank models the frequency selectivity of the peripheral auditory
%   system incorporating basilar membrane nonlinearity, in attempts to 
%   more closely follow the human physiological findings [1,2]. 
%   It operates on a time-domain signal and returns a
%   time-frequency representation of the signal. For each frequency
%   channel, the signal passes through a linear and nonlinear processing
%   paths which consist of combinations of gain, gammatone filters, nonlinear
%   compression, and/or low-pass filters. The Medial Olivo-Cochlear (MOC)
%   efferent feedback path is also realised as the gain at the nonlinear
%   path [1], which can be adjusted by the user for simulations.
%
%   DRNLPROC properties:
%       cfHz       - Characteristic frequencies (Hz)
%       mocIpsi    - Ipsilateral MOC feedback (as nonlinear gain factor)
%       mocContra  - Contralateral MOC feedback (as nonlinear gain factor)
%       model      - DRNL implementation model (based on CASP [2] or MAP [1])
%
%   There are three different ways of setting up a vector of characteristic frequencies
%   (cfHz) when instantiating this processor:
%       1- By providing the lower and upper characteristic frequencies (lowFreqHz and highFreqHz),
%          and the distance between neighboring filters (nERBs).
%       2- By providing the lower and upper characteristic frequencies (lowFreqHz and highFreqHz),
%          and the number of channels that the representation should have.
%       3- By directly providing a vector of characteristic frequencies (cfHz).
%   In case of conflicting arguments, cfHz is generated from one of the three method above
%   with priority order 3 > 2 > 1.
%
%   See also: Processor, gammatoneProc
%
%   Reference:
%   [1] Clark, N. R., Brown, G. J., J?rgens, T., & Meddis, R. (2012). 
%    A frequency-selective feedback model of auditory efferent suppression 
%    and its implications for the recognition of speech in noise. 
%   The Journal of the Acoustical Society of America, 132(3), 1535?41. 
%   [2] Jepsen, M. L., Ewert, S. D., & Dau, T. (2008). 
%   A computational model of human auditory signal processing and perception.
%   The Journal of the Acoustical Society of America, 124(1), 422?38. 

    properties (Dependent = true)
        cfHz                % Characteristic Frequencies
        % NOTE: parameter cfHz here is DIFFERENT FROM cfHz as used in
        % gammatoneProc! - cfHz in gammatoneProc means CENTER FREQUENCY
        % cfHz is used just to follow the framework convention - e.g., some
        % processors after the gammatone filterbank stage expect 'cfHz' for
        % internal operations, and this should be the 'BM characteristic
        % frequency' when the gammatone filterbank is replaced by DRNL 
        % filterbank. See the constructor for the difference in the naming
        % convention: Characteristic Frequencies are denoted by cf, Center
        % Frequencies are denoted by fc
        
        mocIpsi             % ipsilateral MOC factors (nonlinear path gain)
        mocContra           % contralateral MOC factors (nonlinear path gain)
        model               % DRNL model (to be extended)
        
        
    end

    properties
        %cfHz                % Characteristic Frequencies 
        % NOTE: parameter cfHz here is DIFFERENT FROM cfHz as used in
        % gammatoneProc! - cfHz in gammatoneProc means CENTER FREQUENCY
        % cfHz is used just to follow the framework convention - e.g., some
        % processors after the gammatone filterbank stage expect 'cfHz' for
        % internal operations, and this should be the 'BM characteristic
        % frequency' when the gammatone filterbank is replaced by DRNL 
        % filterbank. See the constructor for the difference in the naming
        % convention: Characteristic Frequencies are denoted by cf, Center
        % Frequencies are denoted by fc
        
%         mocIpsi             % ipsilateral MOC factors (nonlinear path gain)
%         mocContra           % contralateral MOC factors (nonlinear path gain)
%         model               % DRNL model (to be extended)       
        
        % Parameters of DRNL filterbank blocks (- could be made private but
        % put here for testing/checking after calculation)
        
        % linear path
        gainLinearPath  
        fcLinPathGammatoneFilter            % linear path GT filter centre frequency (Hz)
        nCascadeLinPathGammatoneFilter      % linear path GT filter # of cascades
        bwLinPathGammatoneFilter            % linear path GT filter bandwidth                   
        cutoffLinPathLowPassFilter = [];    % linear path LP filter cutoff frequency (Hz)
        nCascadeLinPathLowPassFilter = 0;   % linear path LP filter # of cascades
        
        % nonlinear path: has two [cascaded] GT filter stages before and
        % after nonlinearity       
        fcNonlinPathGammatoneFilter         % nonlinear path GT filter centre frequency
        nCascadeNonlinPathGammatoneFilter   % nonlinear path GT filter # of cascades
        bwNonlinPathGammatoneFilter         % nonlinear path GT filter bandwidth
        % nonlinearity section - the parameters a, b, and c may have 
        % different definitions depending on the implementation model (CASP? MAP?)
        aNonlinPath                         % nonlinear path parameter 'a'
        bNonlinPath                         % parameter 'b'(CASP) or 'ctBMdB'(MAP)
        cNonlinPath                         % parameter 'c'
        % 
        nCascadeNonlinPathGammatoneFilter2  % nonlinear path GT filter AFTER BROKEN STICK STAGE, # of cascades
        cutoffNonlinPathLowPassFilter = []; % nonlinear path LPF cutoff
        nCascadeNonlinPathLowPassFilter = 0;% nonlinear path LPF # of cascades
        
        % MAP1_14h-specific parameters (only used with MAP model)
        % highpass stapes filter (1st order HP filter)
        mapStapesHPcutoff = 1000;
        % set scalar. NB Huber gives 2e-9 m at 80 dB, 1 kHz. (==2e-13 at 0 dB SPL)
        mapStapesScalar = 45e-9;

    end
    
    properties (GetAccess = private)
        
        mapTMLowpassFilter = [];    % Tympanic Membrane(TM) low pass filter    
        mapMEHighpassFilter = [];   % Middle Ear high pass filter to simulate stapes inertia
        
        GTFilters_lin               % GT filters for linear path
        GTFilters_nlin              % GT filters for nonlinear path 
        GTFilters_nlin2             % GT filters for nonlinear path, AFTER BROKEN STICK STAGE
        LPFilters_lin = [];         % Low Pass Filters for linear path
        LPFilters_nlin = [];        % Low Pass Filters for nonlinear path
        
    end
        
    methods

        function pObj = drnlProc(fs,parObj)
        %drnlProc   Construct a dual-resonnance, non-linear filterbank processor
        %
        % USAGE:
        %   pObj = drnlProc(fs, parObj)
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
            
            % TODO: Restore the ability to pick parameters according to a given model.
            % Former code is left here commented out for this purpose
            %
            % drnlProc      Construct a DRNL filterbank inheriting
            %                   the "processor" class
            %
            % USAGE
            % - Minimal usage, either (in order of priority)
            %    pObj = drnlProc(fs,[],[],[],[],cfHz)
            %    pObj = drnlProc(fs,flow,fhigh,[],nChan)
            %    pObj = drnlProc(fs,flow,fhigh,nERBs)
            %    pObj = drnlProc(fs,flow,fhigh)
            %
            % - Additional arguments:
            %    pObj = drnlProc(fs,..., cfHz, mocIpsi, mocContra, model)
            %
            % INPUT ARGUMENTS
            %       fs : Sampling frequency (Hz)
            %     flow : Lowest characteristic frequency for the filterbank (in Hz)
            %    fhigh : Highest characteristic frequency for the filterbank (in Hz)
            %    nERBs : Distance in ERBS between neighboring characteristic
            %               frequencies (default: nERBS = 1)
            %    nChan : Number of characteristic frequency channels
            %     cfHz : Vector of channel characteristic frequencies
            %   mocIpsi: Ipsilateral MOC feedback factor as a nonlinear gain
            %               Can be given as a scalar (across all freq.
            %               channels) or a vector (individual gain per
            %               freq. channel)
            % mocContra: Contralateral MOC feedback factor, format same as mocIpsi    
            %     model: implementation model
            %               'CASP' (default) is based on Jepsen et al. 2008 
            %               JASA paper and 'MAP' is based on
            %               MAP1_14h (Ray Meddis, paper to follow)
            %
            % OUTPUT ARGUMENTS
            %    pObj: Processor object
            
%             if nargin>0
            % Checking input arguments - should be between 3 and 9,
            % otherwise error so "if nargin>0" not necessary
%             narginchk(3, 9)
%             
%             if nargin < 9 || isempty(model); model = 'CASP'; end
%             if nargin < 8 || isempty(mocContra); mocContra = 1; end
%             if nargin < 7 || isempty(mocIpsi); mocIpsi = 1; end
            
            if nargin<2||isempty(parObj); parObj = Parameters; end
            if nargin<1; fs = []; end

            % Call super-constructor
            pObj = pObj@Processor(fs,fs,'drnlProc',parObj);
            
            if nargin>0 && ~isempty(fs)
                
            end
            
            
        end
      
        function out = processChunk(pObj,in)
            %processChunk       Passes an input signal through the
            %                   DRNL filterbank
            %
            %USAGE
            %       out = processChunk(pObj,in)
            %       out = pObj.processChunk(in)
            %
            %INPUT ARGUMENTS
            %      pObj : DRNL filterbank object
            %        in : One-dimensional array containing the input signal
            %
            %OUTPUT ARGUMENTS
            %       out : Multi-dimensional array containing the filterbank
            %             outputs
            %
            %SEE ALSO:
            %       drnlProc.m
            
            % TO DO: Indicate that this function is not buit to deal with
            % multiple channels. Multiple channels should be treated with
            % multiple instances of the filterbank.
            
            % Check inputs
            if min(size(in))>1
                error('The input should be a one-dimensional array')
            end
            
            % Turn input into column vector
            in = in(:);
            
            % Get number of channels (CFs)
            nFilter = length(pObj.cfHz);

            % Pre-allocate memory
            out_lin = zeros(length(in),nFilter);
            out_nlin = zeros(length(in),nFilter);

            switch pObj.model
                case 'CASP'
                    % Assuming level scaling and middle ear filtering have already
                    % been taken care of (Preprocessing)
                    
                    % The current implementation of the DRNL works with
                    % dboffset=100, so we must change to this setting.
                    % The output is always the same, so there is no need for changing back.

%                     % Obtain the dboffset currently used
%                     dboffset=dbspl(1);
% 
%                     % Switch signal to the correct scaling
%                     in=gaindb(in, dboffset-100);

                    % Loop through the CF channels (places on BM)
                    % depending on the number of CF elements, all the parameters
                    % (a, b, g, BW, etc.) can be single values or vectors
                    for ii = 1:nFilter
                        % linear path
                        % apply linear gain
                        out_lin(:, ii) = in.*pObj.gainLinearPath(ii);
                        % linear path GT filtering - cascaded "nCascadeLinPathGammatoneFilter" times
                        % already (when the filter objects were initiated)
                            out_lin(:, ii) = ...
                                pObj.GTFilters_lin(ii).filter(out_lin(:, ii));
                        % linear path LP filtering - cascaded "nCascadeLinPathLowPassFilter" times
                            out_lin(:, ii) = ...
                                pObj.LPFilters_lin(ii).filter(out_lin(:, ii));

                        % nonlinear path
                        % MOC attenuation applied (as gain factor)
                        out_nlin(:, ii) = in.*pObj.mocIpsi(ii).*pObj.mocContra(ii);
                        % nonlinear path GT filtering - cascaded "nCascadeNonlinPathGammatoneFilter"
                        % times
                            out_nlin(:, ii) = ...
                                pObj.GTFilters_nlin(ii).filter(out_nlin(:, ii));
                        % broken stick nonlinearity
                        % refer to (Lopez-Poveda and Meddis, 2001) 
                        % note that out_nlin(:, ii) is a COLUMN vector!
                        y_decide = [pObj.aNonlinPath(ii).*abs(out_nlin(:, ii)) ...
                            pObj.bNonlinPath(ii).*abs(out_nlin(:, ii)).^pObj.cNonlinPath];
                        out_nlin(:, ii) = sign(out_nlin(:, ii)).*min(y_decide, [], 2);
                        % nonlinear path GT filtering again afterwards - cascaded
                        % "nCascadeNonlinPathGammatoneFilter2" times
                            out_nlin(:, ii) = ...
                                pObj.GTFilters_nlin2(ii).filter(out_nlin(:, ii));
                        % nonlinear path LP filtering - cascaded "nCascadeNonlinPathLowPassFilter" times
                            out_nlin(:, ii) = ...
                                pObj.LPFilters_nlin(ii).filter(out_nlin(:, ii));

                    end
                case 'MAP'
                    % In this case the original input must be represented
                    % in pascals and transformed into stapes DISPLACEMENT
                    % through a dedicated middle ear filtering process
                    % Assume level scaling is done but middle ear filtering
                    % (at preprocessing) is not
                    
                    % Convert signal into pascals
                    % Obtain the dboffset currently used
                    dboffset=dbspl(1);    
                    % Represent in pascals
                    in = in * 20e-5 * 10^(dboffset/20);

                    % Middle Ear filtering 1: convert input pressure (velocity) to
                    %  tympanic membrane(TM) displacement using low pass filter
                %     [TMdisplacementSegment  OME_TMdisplacementBndry] = ...
                %         filter(TMdisp_b,TMdisp_a,in, ...
                %         OME_TMdisplacementBndry);
                    TMdisplacement = pObj.mapTMLowpassFilter.filter(in);

                    % ME filtering 2: middle ear high pass filter simulate stapes inertia
                    stapesDisplacement = pObj.mapMEHighpassFilter.filter(TMdisplacement);
                %     [stapesDisplacement  OMEhighPassBndry] = ...
                %         filter(stapesDisp_b,stapesDisp_a,TMdisplacementSegment, ...
                %         OMEhighPassBndry);

                    % ME stage 3:  apply stapes scala
                    % in now becomes diplacement (m)
                    in = stapesDisplacement*pObj.mapStapesScalar;
                    
                    % set nonlinear compression-related constants
                    referenceDisplacement = 1e-9;
                    ctBM = referenceDisplacement*10^(pObj.bNonlinPath/20);
                    % CtBM is the displacement knee point (m)
                    % CtS is computed here to avoid repeated division by a later
                    % a==0 means no nonlinear path active
                    if pObj.aNonlinPath>0
                        CtS = ctBM./pObj.aNonlinPath; 
                    else CtS = inf(length(nFilter),1); 
                    end
                    
                    % Loop through the CF channels (places on BM)
                    % depending on the number of CF elements, all the parameters
                    % (a, b, g, BW, etc.) can be single values or vectors
                    for ii = 1:nFilter
                        % linear path
                        % apply linear gain
                        out_lin(:, ii) = in.*pObj.gainLinearPath;
                        % linear path GT filtering - cascaded "nCascadeLinPathGammatoneFilter" times
                        % already (when the filter objects were initiated)
                            out_lin(:, ii) = ...
                                pObj.GTFilters_lin(ii).filter(out_lin(:, ii));

                        % nonlinear path
                        % MOC attenuation applied (as gain factor)
                        out_nlin(:, ii) = in.*pObj.mocIpsi(ii).*pObj.mocContra(ii);
                        % nonlinear path GT filtering - cascaded "nCascadeNonlinPathGammatoneFilter"
                        % times
                            out_nlin(:, ii) = ...
                                pObj.GTFilters_nlin(ii).filter(out_nlin(:, ii));
                        % Nick Clark's compression algorithm
                        % note that out_nlin(:, ii) is a COLUMN vector!
                        abs_x = abs(out_nlin(:, ii));
                        signs = sign(out_nlin(:, ii));
                        % below ct threshold= abs_x<CtS;
                        % (CtS= ctBM/DRNLa -> abs_x*DRNLa<ctBM)
                        belowThreshold = abs_x<CtS(ii);
                        out_nlin(belowThreshold, ii) = ...
                            pObj.aNonlinPath(ii).*out_nlin(belowThreshold, ii);
                        aboveThreshold = ~belowThreshold;
                        out_nlin(aboveThreshold, ii) = signs(aboveThreshold) *ctBM .* ...
                            exp(pObj.cNonlinPath * log(pObj.aNonlinPath(ii)*abs_x(aboveThreshold)/ctBM) );                                              
                        % nonlinear path GT filtering again afterwards - cascaded
                        % "nCascadeNonlinPathGammatoneFilter2" times
                            out_nlin(:, ii) = ...
                                pObj.GTFilters_nlin2(ii).filter(out_nlin(:, ii));

                    end                 
            end     % end switch                  
                    
            % now add the outputs
            out = out_lin + out_nlin;
        end
        
        function reset(pObj)
            %reset          Order the processor to reset its internal
            %               states, e.g., when some critical parameters in
            %               the processing have been changed
            %USAGE
            %       pObj.reset()
            %       reset(pObj)
            %
            %INPUT ARGUMENT
            %       pObj : Processor object
            
            % number of "CF" channels - check whether this can vary!!!
            nFilter = numel(pObj.cfHz);
            
            % Resetting the internal states of the internal filters
            for ii = 1:nFilter
                pObj.GTFilters_lin(ii).reset();
                pObj.GTFilters_nlin(ii).reset();
                pObj.GTFilters_nlin2(ii).reset();
                if strcmp(pObj.model, 'CASP')
                    pObj.LPFilters_lin(ii).reset();
                    pObj.LPFilters_nlin(ii).reset();
                end
            end
            
        end
        
        function bInBranch = isSuitableForRequest(pObj)
            if strcmp(pObj.parameters.map('fb_type'),'drnl')
                bInBranch = true;
            else
                bInBranch = false;
            end
        end
        
    end
         
    methods (Access=protected)
        
        function verifyParameters(pObj)
            
            % Solve the conflicts between center frequencies, number of channels, and
            % distance between channels
            if ~isempty(pObj.parameters.map('fb_cfHz'))
                % Highest priority case: a vector of channels center 
                %   frequencies is provided
                centerFreq = pObj.parameters.map('fb_cfHz');
                
                pObj.parameters.map('fb_lowFreqHz') = centerFreq(1);
                pObj.parameters.map('fb_highFreqHz') = centerFreq(end);
                pObj.parameters.map('fb_nChannels') = numel(centerFreq);
                pObj.parameters.map('fb_nERBs') = 'n/a';
                
                
            elseif ~isempty(pObj.parameters.map('fb_nChannels'))
                % Medium priority: frequency range and number of channels
                %   are provided
               
                % Build a vector of center ERB frequencies
                ERBS = linspace( freq2erb(pObj.parameters.map('fb_lowFreqHz')), ...
                                freq2erb(pObj.parameters.map('fb_highFreqHz')), ...
                                pObj.parameters.map('fb_nChannels') );  
                centerFreq = erb2freq(ERBS);    % Convert to Hz
                
                pObj.parameters.map('fb_nERBs') = (ERBS(end)-ERBS(1)) ...
                                                / pObj.parameters.map('fb_nChannels');
                pObj.parameters.map('fb_cfHz') = centerFreq;
                
                
            else
                % Lowest (default) priority: frequency range and distance 
                %   between channels is provided (or taken by default)
                
                % Build vector of center ERB frequencies
                ERBS = freq2erb(pObj.parameters.map('fb_lowFreqHz')): ...
                                double(pObj.parameters.map('fb_nERBs')): ...
                                            freq2erb(pObj.parameters.map('fb_highFreqHz'));
                centerFreq = erb2freq(ERBS);    % Convert to Hz
                
                pObj.parameters.map('fb_nChannels') = numel(centerFreq);
                pObj.parameters.map('fb_cfHz') = centerFreq;
                
            end
            
        end
        
    end
    
    methods (Hidden = true)
        
        function prepareForProcessing(pObj)
           
            fs = pObj.FsHzIn;
            
            % Compute internal parameters and instantiate filters
            
            % Switch bUnityComp of preProc to zero to make the middle ear
            % filtering fully effective
            % (By default this is set to 1 for the use with gammatone
            % filterbank, enabling unity gain middle ear filtering)
            % The cell array given by the .LowerDependencies of a processor 
            % contains handle(s) to the processor(s) just below 
            % in the processing tree.
            pObj.LowerDependencies{1}.parameters.map('pp_bUnityComp') = false;
            
            % mocIpsi, mocContra can be a scalar or vector with the same
                % length as cfHz (individual nonlinear gain per channel)
                if isscalar(pObj.mocIpsi)            % single mocIpsi across all channels
                    pObj.parameters.map('fb_mocIpsi') ...
                                    = pObj.mocIpsi*ones(size(pObj.cfHz));
                else                            % mocIpsi given as vector
                    if size(pObj.mocIpsi) ~= size(pObj.cfHz)
                        error('mocIpsi must be a scalar or of the same dimension as cfHz');
                    end
                end
                if isscalar(pObj.mocContra)          % single mocIpsi across all channels
                    pObj.parameters.map('fb_mocContra') ...
                                    = pObj.mocContra*ones(size(pObj.cfHz));
                else                            % mocIpsi given as vector
                    if size(pObj.mocContra) ~= size(pObj.cfHz)
                        error('mocContra must be a scalar or of the same dimension as cfHz');
                    end
                end

                % set the other internal parameters and initialise (depending
                % on model)
                switch pObj.model
                    case 'CASP'
                        % grab default DRNL parameters here
                        % the default parameters follow Jepsen model's definition
                        % linear path
                        pObj.gainLinearPath = 10.^(4.20405 -.47909*log10(pObj.cfHz)); %g
                        pObj.fcLinPathGammatoneFilter = 10.^(-0.06762+1.01679*log10(pObj.cfHz)); % Hz, CF_lin
                        pObj.nCascadeLinPathGammatoneFilter = 2; % number of cascaded gammatone filters
                        pObj.bwLinPathGammatoneFilter = 10.^(.03728+.75*log10(pObj.cfHz)); % Hz, bwLinPathGammatoneFilter
                        pObj.cutoffLinPathLowPassFilter = 10.^(-0.06762+1.01*log10(pObj.cfHz)); % Hz, LP_lin cutoff
                        pObj.nCascadeLinPathLowPassFilter = 4; % no. of cascaded LP filters
                        % nonlinear path
                        pObj.fcNonlinPathGammatoneFilter = 10.^(-0.05252+1.01650*log10(pObj.cfHz)); % Hz, CF_nlin
                        pObj.nCascadeNonlinPathGammatoneFilter = 2; % number of cascaded gammatone filters
                        % RCK 21.10.2014, the 2008 paper uses 0.77 for m
                        % instead of 0.70 below for BW_nlin
                        pObj.bwNonlinPathGammatoneFilter = 10.^(-0.03193+.70*log10(pObj.cfHz)); % Hz, bwNonlinPathGammatoneFilter
                        % Warning: note that cfHz can be a vector now!!
                        for ii=1:length(pObj.cfHz)                  
                            if pObj.cfHz(ii)<=1000
                                % SE 03.02.2011, the 2008 paper states <= 1500 Hz
                                % 06/04/2011 CI: answer from Morten regarding the discontinuity:
                                % This is imprecisely described in the paper. It was simulated as
                                % described with parameter a, having the value for 1500 Hz, for CFs
                                % above 1000 Hz. I do recognize the discontinuity in the derived
                                % parameter, but I think this is not critical
                                pObj.aNonlinPath(ii) = 10.^(1.40298+.81916*log10(pObj.cfHz(ii))); % a, the 1500 assumption is no good for compressionat low freq filters
                                pObj.bNonlinPath(ii) = 10.^(1.61912-.81867*log10(pObj.cfHz(ii))); % b [(m/s)^(1-c)]
                            else
                                pObj.aNonlinPath(ii) = 10.^(1.40298+.81916*log10(1500)); % a, the 1500 assumption is no good for compressionat low freq filters
                                pObj.bNonlinPath(ii) = 10.^(1.61912-.81867*log10(1500)); % b [(m/s)^(1-c)]
                            end
                        end
                        pObj.nCascadeNonlinPathGammatoneFilter2 = 2; % number of cascaded gammatone filters AFTER BROKEN STICK NONLINEARITY STAGE
                        pObj.cNonlinPath = 10^(-.60206); % c, compression coeff
                        pObj.cutoffNonlinPathLowPassFilter = 10.^(-0.05252+1.01*log10(pObj.cfHz)); % LP_nlincutoff
                        pObj.nCascadeNonlinPathLowPassFilter = 1; % no. of cascaded LP filters in nlin path

                        % CASP2008 uses LPF cutoff frequencies for the GTF
                        % centre frequencies (first parameter in function)
                        % cutoff frequency and bandwidth are in Hz (note the
                        % difference from gammatoneProc.m where the bandwidth
                        % is given in ERBs)
                        pObj.GTFilters_lin = pObj.populateGTFilters(pObj.cutoffLinPathLowPassFilter, fs,...
                            pObj.bwLinPathGammatoneFilter, pObj.nCascadeLinPathGammatoneFilter);
                        pObj.GTFilters_nlin = pObj.populateGTFilters(pObj.cutoffNonlinPathLowPassFilter, fs,...
                            pObj.bwNonlinPathGammatoneFilter, pObj.nCascadeNonlinPathGammatoneFilter); 
                        pObj.GTFilters_nlin2 = pObj.populateGTFilters(pObj.cutoffNonlinPathLowPassFilter, fs,...
                            pObj.bwNonlinPathGammatoneFilter, pObj.nCascadeNonlinPathGammatoneFilter2); 
                        % Instantiating the LPFs
                        pObj.LPFilters_lin = pObj.populateLPFilters(pObj.cutoffLinPathLowPassFilter, fs, pObj.nCascadeLinPathLowPassFilter);
                        pObj.LPFilters_nlin = pObj.populateLPFilters(pObj.cutoffNonlinPathLowPassFilter, fs, pObj.nCascadeNonlinPathLowPassFilter);           

                    case 'MAP' 
                        % set parameters based on MAP1_14h implementation
                        % (MAPparamsNormal)

                        % Middle Ear filter (specific to MAP)
                        % pressure to displacement conversion using smoothing filter (50 Hz cutoff)
                        tau = 1/(2*pi*50);
                        dt = 1/fs;
                        a1 = dt/tau-1; a0 = 1;
                        b0 = 1+a1;
                        TMdisp_b = b0; TMdisp_a = [a0 a1]; % filter coeffs
                        pObj.mapTMLowpassFilter = genericFilter(...
                            TMdisp_b, TMdisp_a, fs, 1);
                        % figure(9), freqz(TMdisp_b, TMdisp_a)
    %                     OME_TMdisplacementBndry=[];                 % saved boundary
                        % OME high pass (simulates poor low frequency stapes response)
                        G = 1/(1+tan(pi*pObj.mapStapesHPcutoff*dt));
                        H = (1-tan(pi*pObj.mapStapesHPcutoff*dt))...
                            /(1+tan(pi*pObj.mapStapesHPcutoff*dt));
                        stapesDisp_b=[G -G];
                        stapesDisp_a=[1 -H];
                        pObj.mapMEHighpassFilter = genericFilter(...
                            stapesDisp_b, stapesDisp_a, fs, 1);
    %                     OMEhighPassBndry=[];                        % saved boundary
                        % figure(10), freqz(stapesDisp_b, stapesDisp_a)

                        % linear path parameters
                        pObj.gainLinearPath = 500; % linear path gain g, from MAP1.14h
                        pObj.fcLinPathGammatoneFilter = 0.62*pObj.cfHz + 266; % Hz, CF_lin, from MAP1.14h
                        pObj.nCascadeLinPathGammatoneFilter = 3; % number of cascaded gammatone filters (termed "Order" in MAP? - needs double checking)
                        pObj.bwLinPathGammatoneFilter = 0.2*pObj.cfHz + 235; % Hz, bwLinPathGammatoneFilter, MAP1.14h defines in a new way
    %                     % the following two parameters do not appear in MAP1_14h 
    %                     % (LPF parameters) but appear in previous versions
    %                     pObj.cutoffLinPathLowPassFilter = 10^(-0.06762+1.01*log10(cfHz)); % Hz, LP_lin cutoff
    %                     pObj.nCascadeLinPathLowPassFilter = 4; % no. of cascaded LP filters

                        % nonlinear path parameters
                        pObj.fcNonlinPathGammatoneFilter = pObj.cfHz; % Hz, CF_nlin, grabbed from MAP
                        pObj.nCascadeNonlinPathGammatoneFilter = 3; % number of cascaded gammatone filters (termed "Order" in MAP? - needs double checking)
                        pObj.bwNonlinPathGammatoneFilter = 0.14*pObj.cfHz + 180; % Hz, bwNonlinPathGammatoneFilter, MAP defines in a new way
                        % broken stick compression - note that MAP has changed the
                        % formula from CASP2008 version!! 
                        pObj.aNonlinPath = 4e3*ones(size(pObj.cfHz)); % a
                        pObj.bNonlinPath = 25; % Using b for ctBMdB of MAP
                        pObj.cNonlinPath = .25; % c, compression coeff
                        pObj.nCascadeNonlinPathGammatoneFilter2 = 3; % number of cascaded gammatone filters AFTER BROKEN STICK NONLINEARITY STAGE
    %                     % the following two parameters do not appear in MAP1_14h 
    %                     % (LPF parameters) but appear in previous versions
    %                     pObj.cutoffNonlinPathLowPassFilter = 10^(-0.05252+1.01*log10(cfHz)); % LP_nlincutoff
    %                     pObj.nCascadeNonlinPathLowPassFilter = 3; % no. of cascaded LP filters in nlin path 

                        % initialise GTFs (using corresponding centre freqs)
                        pObj.GTFilters_lin = pObj.populateGTFilters(pObj.fcLinPathGammatoneFilter, fs,...
                            pObj.bwLinPathGammatoneFilter, pObj.nCascadeLinPathGammatoneFilter);
                        pObj.GTFilters_nlin = pObj.populateGTFilters(pObj.fcNonlinPathGammatoneFilter, fs,...
                            pObj.bwNonlinPathGammatoneFilter, pObj.nCascadeNonlinPathGammatoneFilter); 
                        pObj.GTFilters_nlin2 = pObj.populateGTFilters(pObj.fcNonlinPathGammatoneFilter, fs,...
                            pObj.bwNonlinPathGammatoneFilter, pObj.nCascadeNonlinPathGammatoneFilter2); 

                    otherwise
                        error('Model not recognised - CASP or MAP supported only');
                end

    %             % Instantiating the LPFs
    %             pObj.LPFilters_lin = pObj.populateLPFilters(pObj.cutoffLinPathLowPassFilter, fs, pObj.nCascadeLinPathLowPassFilter);
    %             pObj.LPFilters_nlin = pObj.populateLPFilters(pObj.cutoffNonlinPathLowPassFilter, fs, pObj.nCascadeNonlinPathLowPassFilter);           

            
        end
        
    end
    
    methods (Access = private)
        function obj = populateGTFilters(pObj,cfHz,fs,bwHz,cascadeOrder)
            % This function is a workaround to assign an array of objects
            % as one of the processor's property, should remain private
              
            nFilter = numel(cfHz);
         
            % Use genericFilter object to exactly copy CASP2008
            % implementation
            % In this case only cfHz, fs, and bw are necessary
            % bw here should indicate the bandwidth in Hz (compare to the
            % use of bw below)
            % Also note that bw is supposed to be a function of cf (could
            % be a vector!)
            
            theta = 2*pi*cfHz(:)/fs;        % convert cfHz to column
            phi   = 2*pi*bwHz(:)/fs;        % bw should be in Hz!!!
            alpha = -exp(-phi).*cos(theta);

            b1 = 2*alpha;
            b2 = exp(-2*phi);
            a0 = abs( (1+b1.*cos(theta)-1i*b1.*sin(theta)+b2.*cos(2*theta)-1i*b2.*sin(2*theta)) ./ (1+alpha.*cos(theta)-1i*alpha.*sin(theta))  );
            a1 = alpha.*a0;

            % adapt to matlab filter terminology
            B=[a0, a1];
            A=[ones(length(theta), 1), b1, b2];
            
            % Preallocate memory by instantiating last filter
            obj(1,nFilter) = genericFilter(B(nFilter,:), A(nFilter, :), fs, cascadeOrder);
            % Instantiating remaining filters
            for ii = 1:nFilter-1
                obj(1,ii) = genericFilter(B(ii,:), A(ii,:), fs, cascadeOrder);
            end                                  
                        
        end
        
        function obj = populateLPFilters(pObj,cfHz,fs,cascadeOrder)
            % This function is a workaround to assign an array of objects
            % as one of the processor's property, should remain private

            nFilter = numel(cfHz);

            % remember! cfHz can be a vector!
            % so convert cfHz to a column vector first before proceeding
            theta = pi*cfHz(:)/fs;
            % now theta is a column vector regardless of what cfHz was

            C = 1./(1+sqrt(2)*cot(theta)+(cot(theta)).^2);
            D = 2*C.*(1-(cot(theta)).^2);
            E = C.*(1-sqrt(2)*cot(theta)+(cot(theta)).^2);

            B = [C, 2*C, C];
            A = [ones(length(theta), 1), D, E];
                                    
            % Preallocate memory by instantiating last filter
            obj(1,nFilter) = genericFilter(B(nFilter,:), A(nFilter, :), fs,cascadeOrder);
            % Instantiating remaining filters
            for ii = 1:nFilter-1
                obj(1,ii) = genericFilter(B(ii,:), A(ii,:), fs,cascadeOrder);
            end                        

        end
    end
        
    % "Getter" methods
    methods 
        
        function cfHz = get.cfHz(pObj)
            cfHz = pObj.parameters.map('fb_cfHz');
        end
        
        function mocIpsi = get.mocIpsi(pObj)
            mocIpsi = pObj.parameters.map('fb_mocIpsi');
        end
        
        function mocContra = get.mocContra(pObj)
            mocContra = pObj.parameters.map('fb_mocContra');
        end
        
        function model = get.model(pObj)
            model = pObj.parameters.map('fb_model');
        end
        
    end
    
    methods (Static)
        
        function dep = getDependency()
            dep = 'time';
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
            
            
            names = {'fb_type',...
                    'fb_lowFreqHz',...
                    'fb_highFreqHz',...
                    'fb_nERBs',...
                    'fb_nChannels',...
                    'fb_cfHz',...
                    'fb_mocIpsi',...
                    'fb_mocContra',...
                    'fb_model'};
            
            descriptions = {'Filterbank type (''gammatone'' or ''drnl'')',...
                    'Lowest center frequency (Hz)',...
                    'Highest center frequency (Hz)',...
                    'Distance between neighbor filters in ERBs',...
                    'Number of channels',...
                    'Channels characteristic frequencies (Hz)',...
                    'Ipsilateral MOC feedback factor as DRNL nonlinear path gain',...
                    'Contralateral MOC feedback factor as DRNL nonlinear path gain',...
                    'DRNL implementation model (''CASP'' or ''MAP'')'};
            
            defaultValues = {'gammatone',...
                            80,...
                            8000,...
                            1,...
                            [],...
                            [],...
                            1,...
                            1,...
                            'CASP'};
                
        end
        
        function pInfo = getProcessorInfo
            
            pInfo = struct;
            
            pInfo.name = 'Dual resonance non-linear filterbank';
            pInfo.label = 'DRNL filterbank';
            pInfo.requestName = 'filterbank';
            pInfo.requestLabel = 'DRNL output';
            pInfo.outputType = 'TimeFrequencySignal';
            pInfo.isBinaural = false;
            
        end
        
    end
    
end