function amplitude=db2amp(dBlevel);

% function amplitude=db2amp(dBlevel);
%
% transform relative soundpressure
% level into amplitude
%
% Jonas Braasch
% Institut fuer Kommunikationsakustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% e-mail: braasch@ika.ruhr-uni-bochum.de
% 
% 18.09.99


amplitude=exp(dBlevel/20*log(10));