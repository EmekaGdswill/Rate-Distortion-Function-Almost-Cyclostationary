close all
clear all
clc
%% paramaeter definition
s_fTpw = 5e-6;     %original signal duration in continuous time
max_sig=1;         % sets maximum of CT variance sigma_sc 
%v_Dist=[0.125:0.0125:1]*4.6161; %set of distortion values

% v_Dist=3 ;
v_fDC = [0.2, 0.45, 0.75, 0.9];   % DC cycle values
v_fSampOffset = [0, 0.0625];           % Sampling offsets
v_fDistrtnDB = -1:0.0125:6;
% default values
s_fDC = 0.75;           % DC cylce
s_fSampOffset = 0;      % Sampling offset
s_fDist=0.18;

s_nXaxis = 2;       % 1 - Index n
                    % 2 - Sampling interval to signal period ratio
                    % 3 - Signal power
                    switch   s_nXaxis
                        case 1
                            v_fX = 1:500;
                            %v_fY = v_Dist;
                            v_fY=v_fDC;
                            v_fZ = v_fSampOffset;
                            s_nP = 2;               % Integer sampling mismatch i.e. P
                            s_fEps = pi/7;          % Async. sampling mismatch
                            v_fPn = v_fX*s_nP + floor(v_fX*s_fEps); %creates vector of periods Pn for each index n
                            v_fTs = s_fTpw*v_fX ./ v_fPn;   % Sampling period
                            stXaxis = 'Index $n$';
                           
                        case 2
                            v_fX = 2:0.01:4; % P+\eps
                            %v_fY = s_fDist;
                            v_fY = v_fDC;
                            v_fZ = v_fSampOffset;
                            [~, v_fN] = rat(v_fX);   % get index n for each ratio i.e. denominator
                            v_fPn = v_fX.*v_fN;
                            v_fTs = s_fTpw./v_fX;    % Sampling duration
                            stXaxis = 'Sampling interval to signal period ratio $\frac{T_{\rm pw}}{T_{\rm s}}$';
                            
                        case 3
%                             v_fY = 2 + [0, pi/1000, 0.1001, 0.2];
                            v_fY = 2 + [0.5, 5*pi/32, 0.6];
%                             v_fX = v_fDistrtnDB;
%                             v_fX = 10.^(0.1*v_fDistrtnDB);
                            v_fX = 0.05:0.001:0.19;
                            v_fZ = v_fSampOffset;
                            [~, v_fN] = rat(v_fY);
                            v_fPn = v_fY.*v_fN;
                            v_fTs = s_fTpw./v_fY;    % Sampling duration
                            stXaxis = 'Distortion Constraint $D$';
                    end
                    
%% Simulations section
R_D=zeros(length(v_fX),length(v_fY), length(v_fZ));

%loop over sampling offset
for aa=1:length(v_fZ)
    tic
    %loop over the DC values
    for ll=1:length(v_fY)
        %loop over the distortion values
        %for kk=1:length(v_fY)
            %loop over the sampling intervals
            for ii=1:length(v_fX)
                switch s_nXaxis
                    case 1
                         s_fTs=v_fTs(ii);
                         s_fPn=v_fPn(ii);
                         %s_fDist=v_fY(kk);
                         s_fDC=v_fY(ll);
                         s_fsampling_offset=v_fZ(aa);
                    case 2
                         s_fTs=v_fTs(ii);
                         s_fPn=v_fPn(ii);
                         %s_fDist=v_fY(kk);
                         s_fDC=v_fY(ll);
                         s_fsampling_offset=v_fZ(aa);
                    case 3
                         s_fTs=v_fTs(ll);
                         s_fPn=v_fPn(ll);
                         s_fDist=v_fX(ii);
                         s_fsampling_offset=v_fZ(aa);
                         
                        
                        
                end 
                %Evaluate discrete time variance Sn
                 %v_fsigma_Sn=sgnshift(s_fPn,s_fQ,s_fTs,s_fZ); 
                 v_fsigma_Sn=v_fSampleVar((1:s_fPn)*(s_fTs/s_fTpw) - s_fsampling_offset, s_fDC);
                 
                 %Evaluate RDF
                 R_D(ii,ll,aa)= s_fGetrateDD(v_fsigma_Sn, s_fDist);            
                 
            end %loop over the index n
            
        %end %loop over the distortion values
        
    end %loop over the DC values
    toc
end%loop over sampling offset
%% Printing results
%% Print results
v_stPlotType = strvcat( '-k' , '-m',  '-b', '-g', '--r', '-r',...
                        '-r^' , '--g*',  '--cx', '-g*', '-rs', '-c<',...
                        '-k*' , '-m+',  '-bs', '-gv', '-ro', '-c^');
                    
for aa = 1:length(v_fZ)
     fig1 = figure;
    set(fig1, 'WindowStyle', 'docked');
    v_stLegend = [];
    for ll = 1:length(v_fY)
           switch s_nXaxis
               case 1
                   stCurLegend = ['t_d_c = ', num2str(v_fDC(ll)),  ', \phi = ', num2str(v_fSampOffset(aa))];
                   stTitle=['$R_{n}(D)$ Versus Index $n$ for $\phi=$', num2str(v_fSampOffset(aa))];
               case 2
                    stCurLegend = ['t_d_c = ', num2str(v_fDC(ll)),  ', \phi = ', num2str(v_fSampOffset(aa))];
                    stTitle='$R_{n}(D)$ Versus Normalized Sampling Rate $\frac{T_{\rm ps}}{T_{\rm s}}$';
               case 3
                    stCurLegend = ['T_s = ', num2str(v_fTs(ll),'%10.3e\n'),  ', \phi = ', num2str(v_fSampOffset(aa))];
                    stTitle='$R_{n}(D)$ Versus $D$';
           end
          v_stLegend = strvcat(v_stLegend, stCurLegend);
          v_fTemp = R_D(:,ll,aa);
          plot(v_fX, v_fTemp, v_stPlotType(ll,:),'LineWidth',1,'MarkerSize',10);
          hold on;
    end
    xlabel(stXaxis, 'Interpreter', 'latex');
    ylabel('RDF $R_n(D)$', 'Interpreter', 'latex');
    title(stTitle,'Interpreter', 'latex');
    grid on;
    legend(v_stLegend,'Location','NorthWest');
    hold off;    
end

