% evaluate rate distortion function of cyclostationary process
function   s_fR = s_fGetrateDD(v_fSigSn, Dist)
%s_fR is the RDF at a specified Dist
%v_fSigSn is the vector of DT cyclostationary processes

s_fPeriod = length(v_fSigSn); %extracts the cyclostationary period

%%
%%%%%%%%%%%%%%Reverse waterfilling%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%NB the distortion and the RDF are normalized
% sigma_sort=[sort(v_fSigSn),inf]; 
% for k=1:s_fPeriod
%    s_fDelta=(s_fPeriod*Dist-(sum(sigma_sort(1:k))))/(s_fPeriod-k);
%    if s_fDelta<sigma_sort(k) 
%        break
%    else
%        continue 
%    end
% end
% 
v_Dk=Dist;


% Evaluate RDF
s_fR = (1/(2*s_fPeriod))*sum(max((log2(v_fSigSn./v_Dk)),0));
