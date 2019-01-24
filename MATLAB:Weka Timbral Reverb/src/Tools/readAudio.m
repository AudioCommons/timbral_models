function signals = readAudio(audioFiles,fsRef)

% Check for proper input arguments
if nargin ~= 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Number of audio files
nFiles = numel(audioFiles);

% Allocate memory
nSamples = zeros(nFiles,2);
fsOrig   = zeros(nFiles,1);

% Loop over number of audio files
for ii = 1 : nFiles
    [nSamples(ii,:),fsOrig(ii)] = wavread(audioFiles{ii},'size'); %#ok
end

% Overall duration
nSamplesMax = round(max(nSamples(:,1)) * (fsRef/max(fsOrig)));

% Allocate memory for signals
signals = zeros(nSamplesMax,nFiles);

% Loop over number of audio files
for ii = 1 : nFiles
    % Read ii-th signal
    [currSig,fs] = wavread(audioFiles{ii}); %#ok
    
    % Resampling, if required
    currSig = resample(currSig,fsRef,fs);

    % Replicate signal to match longest signal
    currSig = repmat(currSig,[ceil(nSamplesMax/size(currSig,1)) 1]);
    
    % Trim edges
    signals(:,ii) = currSig(1:nSamplesMax);
    
    % Normalize ii-th signal by its RMS value
    signals(:,ii) = signals(:,ii) / rms(signals(:,ii));
end

