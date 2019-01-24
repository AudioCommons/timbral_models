classdef adaptationProc < Processor
%ADAPTATIONPROC Adaptation loop processor.
%   The Adaptation loop models corresponds to the adaptive response of 
%   the auditory nerve fibers, in which abrupt changes in the input 
%   result in emphasised overshoots followed by gradual decay to 
%   compressed steady-state level. This is a chain of feedback
%   loops in series, each of which consists of a low-pass filter and a
%   division operator [1,2]. The input to the processor is a time-frequency signal
%   from the inner hair cell, and the output is a time-frequency signal.
%
%   ADAPTATIONPROC properties:
%       overshootLimit      - a limit to the overshoot caused by the loops
%       minLeveldB          - the lowest audible threshhold of the signal (dB)
%       tau                 - Time constants of the loops
%       model (optional)    - implementation model as in various related
%                               studies. If this is given, the first three
%                               parameters are overwritten according to the
%                               description below:
%             'adt_dau'        Choose the parameters as in the Dau 1996 and 1997
%                              models. This consists of 5 adaptation loops with
%                              an overshoot limit of 10 and a minimum level of
%                              1e-5. This is a correction in regard to the model
%                              described in Dau et al. (1996a), which did not use 
%                              overshoot limiting. The adaptation loops have an 
%                              exponential spacing. The default values of the
%                              overshootLimit, minLeveldB, and tau correspond
%                              to this model.
%             'adt_puschel'    Choose the parameters as in the original Puschel 1988
%                              model. This consists of 5 adaptation loops without
%                              overshoot limiting. The adapation loops have a linear spacing.
%             'adt_breebaart'  As 'puschel', but with overshoot limiting.
%             'adt_vandorpschuitman'
%                              5 adaptation loops, same tau spacing as
%                              'adt_dau', no overshoot limiting,
%                              FREQUENCY-DEPENDENT Absloute Threshold of Hearing (ATH)
%                              according to (Terhardt 1979). See [3].

%
%   See also: Processor, ihcProc
%
%   Reference:
%   [1] Puschel, D. (1988). Prinzipien der zeitlichen Analyse beim H?ren. University of G?ttingen.
%   [2] Dau, T., Puschel, D., & Kohlrausch, A. (1996). 
%       A quantitative model of the "effective" signal processing 
%       in the auditory system. I. Model structure. 
%       The Journal of the Acoustical Society of America, 99(6), 3615-3622. 
%   [3] van Dorp Schuitman, J., de Vries, D., & Lindau, A. (2013). 
%       Deriving content-specific measures of room acoustic perception using 
%       a binaural, nonlinear auditory model. 
%       The Journal of the Acoustical Society of America, 133(March), 1572–1585. 

    properties (Dependent = true)
     overshootLimit      % limit to the overshoot of the output
     minLeveldB          % the lowest audible threshhold of the signal 
     tau                 % time constants involved in the adaptation loops 
                         % the number of adaptation loops is determined
                         % by the length of tau
     model               % implementation model - if not empty this setting
                         % overrides the other three parameters
    end

    properties (GetAccess = private)
     stateStore         % cell to store previous output
                        % each element has the same length as tau
                        % # of elements depends on the freq. channels
     minLevel           % minLevel converted from dB to signal value
     cfHz                % CF values to be copied from filterbank
    end
     
    methods
        function pObj = adaptationProc(fs,parObj)
        %adaptationProc   Construct an adaptation loop processor
        %
        % USAGE:
        %   pObj = adaptationProc(fs, parObj)
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
        
        if nargin<2||isempty(parObj); parObj = Parameters; end
        if nargin<1; fs = []; end

        % Call super-constructor
        pObj = pObj@Processor(fs,fs,'adaptationProc',parObj);

        % Prepare to convert minLeveldB following level convention
        % convention: 100 dB SPL corresponds to rms 1
        % calibration factor (see Jepsen et al. 2008)
        dBSPLCal = 100;         % signal amplitude 1 should correspond to max SPL 100 dB
        ampCal = 1;             % signal amplitude to correspond to dBSPLRef

        pObj.minLevel = ampCal*10.^((pObj.minLeveldB-dBSPLCal)/20);

        % initialise the stateStore cell
        % the sizes are unknown at this point - determined by the
        % length of cf (given from the input time-frequency signal)
        pObj.stateStore = {};


        end
         
        function out = processChunk(pObj,in)
            % On-line chunk-based processing is considered
            % in: time-frequency signal (time (row) x frequency (column))
            % The input level is assumed to follow the "convention" that
            % 100 dB SPL corresponds to signal amplitude of 1
            % TODO: To allow for flexibility, configuration should be possible
            % at the beginning...

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for DRNL - CASP2008: needs some additional steps between IHC and adaptation when
            % DRNL is used (not needed when gammatone filterbank is used)
            % Check whether drnlProc is in the dependency list
            % This may need to change if there are more than one stages
            % between filterbank and adaptation loop (in any future)!
            if strcmp(pObj.LowerDependencies{1}.LowerDependencies{1}.Type, ...
                    'drnl filterbank')
                % linear gain to fit ADloop operating point
                in = in*10^(50/20);
                % expansion
                in = in.^2;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            sigLen=size(in, 1);         % signal length = column length
            nChan = size(in, 2);        % # of frequency channels
            nLoops = length(pObj.tau);  % number of loops 
            
            % Extract CF values from filterbank
            % Note that this may need to change if there are more than one
            % dependent processing stages between this and filterbank
            pObj.cfHz = pObj.LowerDependencies{1}.LowerDependencies{1}.cfHz;

            % If stateStore is not defined yet (for the very
            % first chunk of signal)
            % then preallocate storage space
            if isempty(pObj.stateStore)
                pObj.stateStore = cell(1, nChan);
            end

            % b0 from RC-lowpass recursion relation y(n)=b0*x(n)+a1*y(n-1)
            % a1 coefficient of the upper IIR-filter
            % These have same dimension as tau
            b0 = 1./(pObj.tau*pObj.FsHzIn);
            a1 = exp(-b0);
            b0 = 1-a1;
            
            % Check if van Dorp Schuitman model was used
            if strcmp(pObj.model, 'adt_vandorpschuitman')
                % van Dorp Schuitman model uses 
                % frequency-dependent ATH from Terhardt 1979
                % ath_terhardt function is inside /Tools folder
                athdB = ath_terhardt(pObj.cfHz);      % dB SPL (1 x Frequency)
                
                % The following athMU is calculated such that the ATH is
                % returned in MUs instead of dB, as exactly implemented in 
                % the original van Dorp Schuitman model (uses a separate function 
                % adpt_thresholdmu.m), but due to a discrepancy in the
                % final adaptation loop output, it is not used, but left
                % here for future reference
                % athMU = adpt_thresholdmu(pObj.cfHz);    % threshold in MU to be used for the final scaling
                
                dBSPLCal = 100;         % signal amplitude 1 should correspond to max SPL 100 dB
                ampCal = 1;             % signal amplitude to correspond to dBSPLRef
                pObj.minLevel = ampCal*10.^((athdB-dBSPLCal)/20);

                % Repmat extends the dimensions of pObj.minLevel to be
                % (Time x Frequency) for max comparison
                out = max(in, repmat(pObj.minLevel, sigLen, 1));
                
                % Scaling values to get a range from 0 to 100 model units
                % corr and mult are fixed throughout the loops
                % dimensions to be matched to output signal (Time x Freq)
                % Steady state output ~= I.^(1/2^nLoops)
                
                % freq. dependent threshold application (conventional)
                corr = repmat((pObj.minLevel).^(1/(2^nLoops)), sigLen, 1);
                mult = 100./(1 - corr);
                
                % original vd Schuitman version, when athMU is calculated above
                % corr = repmat(athMU, sigLen, 1);
                % mult = 100./...
                %    ((10^5)^(1/(2^nLoops)) - repmat(athMU, sigLen, 1));                
                
                % Determine steady-state levels 
                % Dimension: (nLoops x Frequency)
                state = repmat(pObj.minLevel, nLoops, 1).^...
                    repmat((1./(2.^((1:nLoops).'))), 1, nChan);    
                
            else
                % pObj.minLevel is a scalar here
                % scaling values to get a range from 0 to 100 model units
                % corr and mult are fixed throughout the loops
                % in this case corr and mult are scalars
                corr = pObj.minLevel^(1/(2^nLoops));		
                mult = 100/(1-corr);

                % Apply minimum level to the input
                out = max(in, pObj.minLevel);       % dimension same as in

                % Determine steady-state levels
                % Initial dimension: (nLoops x 1)
                % use repmat to extend to (nLoops x Frequency)
                state = repmat( (pObj.minLevel.^(1./(2.^((1:nLoops).')))), ...
                    1, nChan );    
               
            end


            % Back up the value, because state is overwritten
            stateinit=state;

            if pObj.overshootLimit <=1 
            % No overshoot limitation
                for ch=1:nChan
                    state(:, ch) = stateinit(:, ch);
                    % If there are state values stored from previous call
                    % of the function, overwrite the starting values with
                    % the stored values
                    if ~isempty(pObj.stateStore{ch})
                        state(:, ch) = pObj.stateStore{ch};
                    end
                    for ii=1:sigLen
                        tmp1=out(ii,ch);
                        % Compute the adaptation loops
                        for jj=1:nLoops
                            tmp1=tmp1/state(jj, ch);
                            state(jj, ch) = a1(jj)*state(jj, ch) + b0(jj)*tmp1;         
                        end   
                        % store the result
                        out(ii,ch)=tmp1;
                    end
                    % Now back up the last state (per freq channel)
                    pObj.stateStore{ch} = state(:, ch);
                end

            else 
            % Overshoot Limitation

                % Max. possible output value, dimension: (nLoops x Frequency)
                maxvalue = (1 - state.^2) * pObj.overshootLimit - 1;
                % Factor in formula to speed it up 
                factor = maxvalue * 2; 			
                % Exponential factor in output limiting function
                expfac = -2./maxvalue;
                offset = maxvalue - 1;

                for ch=1:nChan
                    state(:, ch) = stateinit(:, ch);
                    % If there are state values stored from previous call
                    % of the function, overwrite the starting values with
                    % the stored values
                    if ~isempty(pObj.stateStore{ch})
                        state(:, ch) = pObj.stateStore{ch};
                    end
                    for ii=1:sigLen
                        tmp1=out(ii,ch);
                        for jj=1:nLoops
                            tmp1=tmp1/state(jj, ch);
                            if ( tmp1 > 1 )
                                tmp1 = factor(jj, ch)/(1+exp(expfac(jj, ch)*(tmp1-1)))-offset(jj, ch);
                            end
                            state(jj, ch) = a1(jj)*state(jj, ch) + b0(jj)*tmp1;
                        end
                    % store the result
                    out(ii,ch)=tmp1;    
                    end
                    % Now back up the last state (per freq channel)
                    pObj.stateStore{ch} = state(:, ch);                    
                end
            end
            
            % Scale to model units
                out = (out-corr).*mult;
        end
         
        function reset(pObj)
             %reset     Resets the internal states 
             %
             %USAGE
             %      pObj.reset
             %
             %INPUT ARGUMENTS
             %  pObj : adaptation processor instance
             
             % empty the stateStore cell
            if(~isempty(pObj.stateStore))
                pObj.stateStore = {};
            end
        end
         
    end
     
    methods (Access=protected)
        function verifyParameters(pObj)
            % Solve any potential conflicts and/or control the validity of parameters
            % here.
            if ~isempty(pObj.parameters.map('adpt_model'))
                % the adpt_model parameter takes priority if given with
                % other parameters
                modelType = pObj.parameters.map('adpt_model');
                switch(modelType)
                            case 'adt_dau'
                                % This means adpt_model was either not
                                % given or given with any other parameter
                                pObj.parameters.map('adpt_lim') = 10;
                                pObj.parameters.map('adpt_mindB') = 0;
                                pObj.parameters.map('adpt_tau') = ...
                                    [0.005 0.050 0.129 0.253 0.500];                                
                            case 'adt_puschel'
                                % 5 adaptation loops, no overshoot limiting, 
                                % linear tau spacing
                                pObj.parameters.map('adpt_lim') = 0;
                                pObj.parameters.map('adpt_mindB') = 0;
                                pObj.parameters.map('adpt_tau') = ...
                                    linspace(0.005,0.5,5);
                            case 'adt_breebaart'
                                % 5 adaptation loops, with [default] overshoot limiting,
                                % linear tau spacing
                                pObj.parameters.map('adpt_lim') = 10;
                                pObj.parameters.map('adpt_mindB') = 0;
                                pObj.parameters.map('adpt_tau') = ...
                                    linspace(0.005,0.5,5);
                            case 'adt_vandorpschuitman'
                                % 5 adaptation loops, no overshoot limiting,
                                % FREQUENCY-DEPENDENT ATH (Terhardt 1979),
                                % tau spacing same as 'adt_dau',
                                pObj.parameters.map('adpt_lim') = 0;
                                pObj.parameters.map('adpt_mindB') = 0;
                                % THIS WILL CHANGE AT PROCESSCHUNK, because
                                % with this model mindB will be decided at
                                % different cfHz channels!
                                pObj.parameters.map('adpt_tau') = ...
                                    [0.005 0.050 0.129 0.253 0.500]; 

                            otherwise
                                % not supported
                                error('%s: Unsupported adaptation loop model',upper(mfilename));
                end
                fprintf('Note: adaptation loop model name %s was given. The model will override other parameters.', modelType);
            end
            
        end
        
    end
     
     % "Getter" methods
     methods
         function overshootLimit = get.overshootLimit(pObj)
             overshootLimit = pObj.parameters.map('adpt_lim');
         end
         
         function minLeveldB = get.minLeveldB(pObj)
             minLeveldB = pObj.parameters.map('adpt_mindB');
         end
         
         function tau = get.tau(pObj)
             tau = pObj.parameters.map('adpt_tau');
         end
         
         function model = get.model(pObj)
             model = pObj.parameters.map('adpt_model');
         end
     end
     
     methods (Static)
         function dep = getDependency()
             dep = 'innerhaircell';
         end
         
         function [names, defaultValues, descriptions] = getParameterInfo()
            %getParameterInfo   Returns the parameter names, default values
            %                   and descriptions for that processor
            %
            %USAGE:
            %  [names, defaultValues, description] =  ihcProc.getParameterInfo;
            %
            %OUTPUT ARGUMENTS:
            %         names : Parameter names
            % defaultValues : Parameter default values
            %  descriptions : Parameter descriptions
            
            
            names = {'adpt_lim',...
                     'adpt_mindB',...
                     'adpt_tau', ...
                     'adpt_model'};
            
            descriptions = {'Adaptation loop overshoot limit',...
                            'Adaptation loop lowest signal level (dB)',...
                            'Adaptation loop time constants (s)',...
                            'Adaptation loop optional implementation model'};
            
            defaultValues = {10,...
                             0,...
                             [0.005 0.05 0.129 0.253 0.5],...
                             ''};
                
          end
         
          function pInfo = getProcessorInfo
            
             pInfo = struct;
             
             pInfo.name = 'Adaptation loop';
             pInfo.label = 'Neural adaptation model';
             pInfo.requestName = 'adaptation';
             pInfo.requestLabel = 'Adaptation loop output';
             pInfo.outputType = 'TimeFrequencySignal';
             pInfo.isBinaural = 0;
             
         end
         
     end
     
 end