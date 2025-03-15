clear all;close all;clc
% % tv OIR - Estimation procedure over number of realizations R=[10, 50,
% 100, 200] and over different values of forgetting factor c=[0, 0.5, 0.7, 0.9, 0.94, 0.96, 0.97, 0.98, 0.99, 0.999]
%%% parameters
%%% parameters
c1max      =  0.3;    % max amplitude
fs=100;    % sampling frequency
nfft=500; % number of points on frequency axis (total)
f_osc     =  0.25;    % frequency of oscillation 
nobs = 1000;                        % 10 s of observation (data samples)
k = 1:nobs;                         % vector of time steps
t = (k)/fs;                         % vector of times
DC=50;                              % duty-cycle
c1 = -c1max*square(2*pi*f_osc*t,DC);  % causal coefficient varying as square periodic waveform
c1(c1==-c1max)=0;
c2=1-c1;
M=4; %%% number of processes
Mv=[1,1,1,1]; % structure of blocks
%%% set poles and self-oscillations at each time step
par.poles{2}=([0.85 0.1]);% 10 Hz
par.poles{3}=([0.85 0.1]);
par.poles{4}=([0.85 0.35]);% 35 Hz

for n=1:length(t)
    %%% set the interactions
    C_31=c1(n);
    C_21=c2(n);
    par.coup=[2 1 1 C_21;3 1 2 C_31;4 3 2 1.5;4 2 1 1.5];
    par.Su=ones(1,M); %variance of innovation processes
    [Am1,Su1]=theoreticalVAR(M,par); %% VAR parameters
    Am1=Am1';
    AmT1(:,:,n)=Am1;
    SuT1(:,:,n)=Su1;
    Am1=reshape(Am1,[M,M,size(Am1,2)/M]);
    AmT(:,:,:,n)=Am1;
    SuT(:,:,n)=Su1;
end

AmT_SQ=AmT;
SuT_SQ=SuT;
%% estimation of time-resolved and time-varying OIR
% generate time series
M=4; %%% number of processes
Mv=[1,1,1,1]; % structure of blocks
M1=length(Mv); % Note: Mv contains the dimension of each block
SIGT=eye(M);
R=50; %number of realizations
% forgetting factor
C=0.98;
% generate non stationary time-series
X_SQ =var_to_tsdata_nonstat(AmT_SQ,SuT_SQ,R);
% time-var identification. X to be N.series x R x time steps
X_SQ=permute(X_SQ,[1,3,2]);
% time-varying identification (RLS) with forgetting factor
% RLS model order estimation
[popt_RLS_SQ] = mos_id_tv_VAR_RLS(X_SQ,10,C);
% RLS identification
[Am_n_SQ,Su_n_SQ]=tvID_VAR_RLS(X_SQ,popt_RLS_SQ,C);
% instantaneous PSD
[S_n_SQ,H,f] = tv_fdVAR(Am_n_SQ,Su_n_SQ,nfft,fs);
% time varying OIR - RLS
out_SQ = tv_oir(S_n_SQ,Mv,fs);

%% plot of the results 
load('colmap.mat');
LAB={'\Omega_{Y_1;Y_2;Y_3}','\Omega_{Y_1;Y_2;Y_4}','\Omega_{Y_1;Y_3;Y_4}','\Omega_{Y_2;Y_3;Y_4}','\Omega_{Y_1;Y_2;Y_3;Y_4}'};
LAB_f={'\nu_{Y_1;Y_2;Y_3}','\nu_{Y_1;Y_2;Y_4}','\nu_{Y_1;Y_3;Y_4}','\nu_{Y_2;Y_3;Y_4}','\nu_{Y_1;Y_2;Y_3;Y_4}'};
t_in=popt_RLS_SQ+5;
figure
% theoretical square periodic waveform (TIME-RESOLVED)
for i=1:M
    % time-resolved
    subplot (2,5,i)
    plot(t(popt_RLS_SQ+1:end),squeeze(out_SQ.OIR{3,1}(i,popt_RLS_SQ+1:end)));
    xlabel('time steps [n]')
    title(LAB{1,i})
    % frequency-specific
    subplot(2,5,i+5)
    imagesc(t(popt_RLS_SQ+1:end),f,squeeze(out_SQ.OIRf{3,1}(i,popt_RLS_SQ+1:end,:))');
    title(LAB_f{1,i})
    colormap(VRVmap)
    colorbar
    xlabel('Time [s]');ylabel('f [Hz]');
end
% order 4 time-resolved
subplot (2,5,5)
plot(t(popt_RLS_SQ+1:end),squeeze(out_SQ.OIR{4,1}(:,popt_RLS_SQ+1:end)));
title(LAB{1,5})
xlabel('time steps [n]')
% order 4 time-frequency
subplot (2,5,10)
imagesc(t(popt_RLS_SQ+1:end),f,squeeze(out_SQ.OIRf{4,1}(:,popt_RLS_SQ+1:end,:))');
title(LAB_f{1,5})
colormap(VRVmap)
colorbar
xlabel('Time [s]');ylabel('f [Hz]');
set(gcf, 'Position', get(0, 'Screensize'));