function p = genParStruct(varargin)
%genParStruct       Generates a valid parameter structure for Two!Ears
%                   Auditory Front-End processor requests.
%
%USAGE:
%     p = genParStruct(par_name,par_value)
%     p = genParStruct(par_name1,par_value,...,par_nameN,par_valueN)
%
%INPUT PARAMETERS:
%  name : Name tag of parameter to be changed
% value : Value for changed parameter
%
%OUTPUT PARAMETER:
%     p : Instance of parameter object class containing non-default parameters
%
%EXAMPLE:
%     p = genParStruct('nERBs',1/2,'IHCMethod','hilbert') will generate a
%     parameter structure corresponding to a request where the gammatone
%     filterbank has a resolution of 0.5 ERBs per channel and where the
%     inner hair-cell envelope is simply the Hilbert envelope.
%
%SEE ALSO: parameterHelper.m

% First check on the inputs
if mod(size(varargin,2),2)==1
    % Then there is an odd number of input parameters, something is wrong
    warning(['Incorrect number of input arguments. Arguments need to be '...
        'provided in pairs of names and values. The whole parameter '...
        'request is disregarded.'])
    
    return
else
    n_par = size(varargin,2)/2;
end

keys = cell(1,n_par);
values = cell(1,n_par);

% Loop on the number of parameters
for ii = 0:n_par-1
   
    % Incorrect parameter name will be picked-up when instantiating the parameter object
    keys{ii+1} = varargin{2*ii+1};
    values{ii+1} = varargin{2*ii+2}; 
    
end

p = Parameters(keys,values);

end