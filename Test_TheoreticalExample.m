%% Theoretical Example (Periodic square waveform)
close all;clear all;clc
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
p_opt=2;
% theoretical time-frequency and time-varying repsentation of VAR
for n=1:length(t)
    %%% set the interactions
    C_31=c1(n);
    C_21=c2(n);
    par.coup=[2 1 1 C_21;3 1 2 C_31;4 3 2 1.5;4 2 1 1.5];
    par.Su=ones(1,M); %variance of innovation processes
    [Am1,Su1]=theoreticalVAR(M,par); %% VAR parameters
    Am1=Am1';
    AmT(:,:,n)=Am1;
    SuT(:,:,n)=Su1;
end
% plot theoretical coupling coefficients
figure
plot(t,squeeze(AmT(1,2,:)));ylim([0.6 1.1]); grid on;
figure
plot(t,squeeze(AmT(1,7,:)));ylim([-0.1 0.4]); grid on
%% Theoretical OIR
% time-frequency representation of TV-VAR model
[S_T,H,f] = tv_fdVAR(AmT,SuT,nfft,fs);
% time-frequency version
out_thSq = tv_oir(S_T,Mv);
% time-resolved version
q=30;
out_thSq_time = tv_lrp_oir(AmT,SuT,q,Mv);
%% plot of the results 
load('colmap.mat');
LAB={'\Omega_{Y_1;Y_2;Y_3}','\Omega_{Y_1;Y_2;Y_4}','\Omega_{Y_1;Y_3;Y_4}','\Omega_{Y_2;Y_3;Y_4}','\Omega_{Y_1;Y_2;Y_3;Y_4}'};
LAB_f={'\nu_{Y_1;Y_2;Y_3}','\nu_{Y_1;Y_2;Y_4}','\nu_{Y_1;Y_3;Y_4}','\nu_{Y_2;Y_3;Y_4}','\nu_{Y_1;Y_2;Y_3;Y_4}'};
t_in=p_opt+5;
figure
% theoretical square periodic waveform (TIME-RESOLVED)
for i=1:M
    % time-resolved
    subplot (2,5,i)
    plot(t(p_opt+1:end),squeeze(out_thSq_time.OIR{3,1}(i,p_opt+1:end)));
    xlabel('time steps [n]')
    title(LAB{1,i})
    % frequency-specific
    subplot(2,5,i+5)
    imagesc(t(t_in:end),f,squeeze(out_thSq.OIRf{3,1}(i,t_in:end,:))');
    title(LAB_f{1,i})
    colormap(VRVmap)
    colorbar
    xlabel('Time [s]');ylabel('f [Hz]');
end
% order 4 time-resolved
subplot (2,5,5)
plot(t(p_opt+1:end),squeeze(out_thSq_time.OIR{4,1}(:,p_opt+1:end)));
title(LAB{1,5})
xlabel('time steps [n]')
% order 4 time-frequency
subplot (2,5,10)
imagesc(t(p_opt+1:end),f,squeeze(out_thSq.OIRf{4,1}(:,p_opt+1:end,:))');
title(LAB_f{1,5})
colormap(VRVmap)
colorbar
xlabel('Time [s]');ylabel('f [Hz]');

set(gcf, 'Position', get(0, 'Screensize'));