% Sample coutinuous-time noise variance
function   v_fSigW = v_fSampleVar(v_fT, s_fDC)
% Output:
%   v_fSigW -  discrete-time variance
% input arguments: 
%   v_fT - time indexes
%   s_fDC - DC cycle

% This function currently implementy periodic continuous pulse. 
% In the future it can be modified to account for other profiles

s_fLow = 0.2;
s_fHigh = 5;
s_fRiseTime = 0.01;

% Get modulo values
v_fTMod = mod(v_fT,1);

% Set output:
% Low interval
v_fSigW = s_fLow*ones(size(v_fTMod));

% Rising interval
v_fSigW(find((v_fTMod <= s_fRiseTime))) = s_fLow + (s_fHigh - s_fLow) *...
                                    v_fTMod(find((v_fTMod <= s_fRiseTime)))/s_fRiseTime;

% High interval
v_fSigW(find((v_fTMod > s_fRiseTime).*(v_fTMod < s_fDC + s_fRiseTime))) = s_fHigh;


% Falling interval
v_fSigW(find((v_fTMod <= s_fDC + 2*s_fRiseTime).*(v_fTMod >= s_fDC + s_fRiseTime))) = ...
                     s_fHigh - (s_fHigh - s_fLow) *...
                      (v_fTMod(find((v_fTMod <= s_fDC + 2*s_fRiseTime) .* (v_fTMod >= s_fDC + s_fRiseTime))) - (s_fDC + s_fRiseTime)) ...
                            /(s_fRiseTime);

