function [A,Su,phi_w,Zn]=tvID_VAR_RLS(XM,p,c)
% Recursive-Least-Square identification for Vector AR model
% the function is also implemented for the analysis of multiple realizations 
% input:
% XM --> Time-series [M x R x N] - number of processes x realizations x time steps
% p--> model oreder
% c --> forgetting factor [0.97 - 0.99] defined as 1-c in the main paper
% output:
% A    --> estimate of AR parmeters (M x Mp x N) 
% Su   --> Residual variance (M x M x N)
% phi_w --> time-varying estimate of the correlation matrix between past
% states of Y stored in W (Mp x Mp x N)
%% References : 
% [1] - Y. Antonacci et al. (2025) - Time-varying information measures: 
% an adaptive estimation of information storage with application to brain-heart interactions
% Frontiers in Network Physiology 
%
% [2] - E. Moller et al. (2001) - Instantaneous multivariate EEG coherence analysis by means of adaptive
% high-dimensional autoregressive models.
%
% [3] - Grieszbach et al. (1994) - Dynamic description of stochastic signal
% by adaptive momentary power and momentary frequency estimation and its application in analysis of
% biological signals. Medical and Biological Engineering and Computing 32 

c=1-c; % example 1-c=0.95--> c=1-0.95=0.05 - adaption factor as defined in [2]

[M,R,N]=size(XM);

Ao=zeros(M,M*p); %initial condition for AR parameters
if c>0
    A=zeros(M,M*p,N);
end % init sequence parameter estimation
Z=zeros(M,R);  % initial instantaneous predicition error 
Su_n=zeros(M,M); % init innovation variance 
if c>0
    Su=zeros(M,M,N);
end   

W=zeros(M*p,R); % init past states of Y - Wn main document
phi_wn=eye(M*p); 	% init correlation matrix of lagged terms
for n=1:N    	% begin of the estimation loop
    if n>p
        phi_wn=(1-c)*phi_wn+W*W'; % correlation matrix of the past of Y - Eq.(7)
        K=phi_wn\W;          % gain vector as defined in Eq.(10)
        Z=squeeze(XM(:,:,n))-Ao*W;  % a-priori estimation error Eq.(14c)
        Zn(n)=mean(mean(Z.^2,1),2); % average residual sum of squares
        Ao=Ao+Z*K';  % recursive estimation of AR parameters Eq.(13)
        phi_w(:,:,n)=phi_wn; % correlation of W over time
        if c>0
            A(:,:,n)=Ao;  % update sequence of parameter estimation
            Su_n=Su_n+c*(Z*Z'/R-Su_n);  %  Innovation variance estimation Eq.(15)
            % based on the results of [3], see [Eq. (6)] - [2], see Eq.(14) 
            Su(:,:,n)=Su_n;  % save innovation variance
        else  Su_n=Su_n+(Z*Z'/R-Su_n)/n; % equation (15) Moeller et al.
            Su(:,:,n)=Su_n; A(:,:,n)=Ao;
        end	% update the mean squared prediction error
        
    end
    if p>1, W(M+1:M*p,:)=W(1:(p-1)*M,:); end
    if p>0, W(1:M,:)=XM(:,:,n); end  %update observation
end				     % end of the estimation loop
