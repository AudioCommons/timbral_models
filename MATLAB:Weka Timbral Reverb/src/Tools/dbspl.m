function y = dbspl(insig,varargin)
%DBSPL RMS value of signal (in dB)
%   Usage: y = dbspl(insig);
%          y = dbspl(insig,'ac');
%
%   `dbspl(insig)` computes the average SPL (sound pressure level) of the
%   input signal measured in dB, using the convention that a pure tone at
%   100 dB SPL has an RMS value of 1.
%
%   `dbspl(insig,'dboffset',dboffset)` specifies a reference level to
%   convert between RMS and dB SPL. Default value is 100. Some commonly
%   used values are:
%
%   * $dboffset = 0$. This convention was used for the development of the
%     DRNL and Breebaart models.
%
%   * $dboffset = -20*log10(20e-6) = 93.98$. This corresponds to the common
%     convention of the reference being 20 micro Pa. Using this
%     convention, the numerical values of signals is the sound pressure
%     measured in Pascal.
%
%   * $dboffset = 100$. This convention was used for the development of the
%     Dau and Jepsen models.
%
%   Globally changing the reference level for all functions in the
%   toolbox can be done by the following code::
%
%     ltfatsetdefaults('dbspl','dboffset',desired_dboffset);
%
%   and the currently used reference level is obtained by::
%
%     current_dboffset = dbspl(1);
%  
%   DBSPL takes the following flags at the end of the line of input
%   parameters:
%
%     'ac'      Consider only the AC component of the signal (i.e. the mean is
%               removed).
%
%     'dim',d   Work along specified dimension. The default value of []
%               means to work along the first non-singleton one.
%
%     'dboffset',dboffset
%               Specify offset in dB. Default value is 100.
%
%   See also: setdbspl
%
%   References: moore2003introduction

%   AUTHOR : Hagen Wierstorf

definput.keyvals.dim=[];
definput.flags.mean={'noac','ac'};
definput.keyvals.dboffset=100;
[flags,kv]=ltfatarghelper({'dim','dboffset'},definput,varargin);

  
% ------ Computation --------------------------

% The level of a signal in dB SPL is given by the following formula:
% level = 20*log10(p/p_0)
% To get to the standard used in the toolbox.
y = 20*log10( rms(insig,flags.mean) )+kv.dboffset;

