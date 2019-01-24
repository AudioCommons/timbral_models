classdef gammatoneFilter < filterObj
    
    properties
        CenterFrequency     % Center frequency for the filter (Hz)
        FilterOrder         % Gammatone slope order
        IRduration          % Duration of the impulse response for truncation (s)
        delay               % Delay in samples for time alignment
    end   
    
    methods
        function obj = gammatoneFilter(cf,fs,n,bwERB,do_align,cascade)
            %gammatoneFilter    Construct a gammatone filter object
            %
            %USAGE
            %           F = gammatoneFilter(fc,fs)
            %           F = gammatoneFilter(fc,fs,type,n,bw,do_align,durSec,cascade)
            % 
            %INPUT ARGUMENTS
            %          cf : center frequency of the filter (Hz)
            %          fs : sampling frequency (Hz)
            %        type : 'fir' for finite impulse response or 'iir' for
            %               infinite (default: type = 'fir')
            %           n : Gammatone rising slope order (default, n=4)
            %          bw : Bandwidth of the filter in ERBS 
            %               (default: bw = 1.08 ERBS) 
            %    do_align : If true, applies phase compensation and compute
            %               delay for time alignment (default : false)
            %      durSec : Duration of the impulse response in seconds 
            %               (default: durSec = 0.128)
            %     cascade : Cascading order of the filter (default : 1)
            %
            %OUTPUT
            %           F : Gammatone filter object
            
            % TO DO : Instantiating an IIR gammatone filter should not
            % populate the IRDuration property
            % TO DO : Add more code commenting/references for 'iir' case
            % TODO: Implement alignment?
            
            if nargin > 0   % Prevent error when constructors is called 
                            %   without arguments
                % Check for input arguments
                if nargin < 2 || nargin > 8
                    error('Wrong number of input arguments')
                end
                
                % Set default parameters
                if nargin < 6 || isempty(cascade)
                    cascade = 1;
                end
                
                if nargin < 5 || isempty(do_align)
                    do_align = false;
                end
                if nargin < 4 || isempty(bwERB)
                    bwERB = 1.018;
                end
                if nargin < 3 || isempty(n)
                    n = 4;
                end
                
                % One ERB value in Hz at this center frequency
                ERBHz = 24.7 + 0.108 * cf;

                % Bandwidth of the filter in Hertz
                bwHz = bwERB * ERBHz;

                % This is when the function peaks
                delay = 3./(2*pi*bwHz);
                
                % Generate an IIR Gammatone filter
                if do_align
                    % Compute the position of the pole
                    atilde = exp(-2*pi*bwHz/fs - 1i*2*pi*cf/fs);
                    
                    % Repeat the pole n times, and expand the polynomial
                    a = poly(atilde*ones(1,n));
                    
                    btmp = 1-exp(-2*pi*bwHz/fs);
                    b = btmp.^n;
                    b = b*exp(2*pi*1i*cf*delay);
                else
                    btmp=1-exp(-2*pi*bwHz/fs);
                    atmp=[1, -exp(-(2*pi*bwHz + 1i*2*pi*cf)/fs)];
                    
                    b=1;
                    a=1;
                    
                    for jj=1:n
                        b=conv(btmp,b);
                        a=conv(atmp,a);
                    end
                end
                
                % The transfer function is complex-valued
                realTF = false;
                
                % Populate filter Object properties
                %   1- Global properties
                obj = populateProperties(obj,'Type','Gammatone Filter','FsHz',fs,...
                    'b',b,'a',a,'RealTF',realTF);
                %   2- Specific properties
                obj.CenterFrequency = cf;
                obj.FilterOrder = n;
                obj.delay = delay;
                obj.CascadeOrder = cascade;
            end
        end
        
    end
end
    