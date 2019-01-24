function [ y, h_rir ] = generate_rev_speech( Tr, fs, speech )

% [ y, h_rir ] = generate_rev_speech( Tr, fs, speech )
%  generates reverberant speech by means of a statistical room impulse
%  response (RIR)
%
% INPUT
% Tr     : reverberation time in seconds
% fs     : sampling frequency of input speech signal
% speech : non-reverberant input speech
%
% OUTPUT
% y      : reverberant speech
% h_rir  : statistical room impules response
%
% author: Heinrich Loellmann, IND, RWTH Aachen
% 
% created: August 2011


Td = 0.1;        % delay time in sec
Tall = Tr+Td;    % duration of the IR in mssec.
L = Tall*fs;     % length of the RIR
No = round(Td*fs)+1;  % index for first non-zero sample

sigma = 0.04;             % power
n0 = sigma*randn( L,1 );  % white Gaussian noise sequence for model
h_rir = zeros(L,1);
ep = zeros(L,1);

tau = Tr / ( 3*log(10) );                  % decay rate
ep(No:L) = exp( -(0:L-No)' / ( tau *fs) ); % exp. core function

h_rir(No:end) = n0(No:end).*ep(No:end);    % actual RIR

[ n,m ] = size(speech);
if n>m
    speech = speech.';        % ensure right dimension
end

y = conv( speech, h_rir );  % create rev. speech