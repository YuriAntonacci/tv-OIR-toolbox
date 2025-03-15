
function [S,H,f] = tv_fdVAR(Am_n,Su_n,nfft,fs)
% time varying Frequency domain MVAR Analysis
% inputs:
% Am_n=[A_n(1)...A_n(p)]: M*pM x N matrix of the MVAR model coefficients
% (strictly causal model) at ech time step
% Su_n: M*M x N covariance matrix of the input noises at each time step
% nfft= number of points for calculation of the spectral functions (nfft)
% Fs= sampling frequency
% outputs:
% H= Tranfer Function Matrix at each time step
% S= Spectral Matrix at each time step
%% References:
% [1] - Y. Antonacci et al. (2025) - TBME
if iscell(Am_n)~=1
    N=size(Am_n,3); % time steps
    M= size(Am_n,1); % Am has dim M*pM x N
    p = size(Am_n,2)/M; % p is the order of the MVAR model
    
    if all(size(nfft)==1)	 %if nfft is scalar
        f = (0:nfft-1)*(fs/(2*nfft)); % frequency axis
    else            % if nfft is a vector, we assume that it is the vector of the frequencies
        f = nfft; nfft = length(nfft);
    end
    
    z = 1i*2*pi/fs;
    % Initializations: spectral matrices have M rows, M columns and are calculated at each of the nfft frequencies
    H=zeros(M,M,nfft,N); % Transfer Matrix
    S=zeros(M,M,nfft,N); % Spectral Matrix
    
    for n=p+1:N
        A = [eye(M) -Am_n(:,:,n)]; % matrix from which M*M blocks are selected to calculate spectral functions
        
        %% computation of spectral functions
        for ff=1:nfft % at each frequency
            
            %%% Coefficient matrix in the frequency domain
            As = zeros(M,M); % matrix As(z)=I-sum(A(k))
            for k = 1:p+1
                As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(ff));
                % indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix A (A(1) is in the second block, and so on)
            end
            
            %%% Transfer matrix
            H(:,:,ff,n)  = inv(As);
            
            %%% Spectral matrix
            S(:,:,ff,n)  = H(:,:,ff,n)*Su_n(:,:,n)*H(:,:,ff,n)'; % ' stands for Hermitian transpose
            
        end
    end
else
    N=size(Am_n,2); % time steps
    M= size(Am_n{end},1); % Am has dim M*pM x N
    
    if all(size(nfft)==1)	 %if nfft is scalar
        f = (0:nfft-1)*(fs/(2*nfft)); % frequency axis
    else            % if nfft is a vector, we assume that it is the vector of the frequencies
        f = nfft; nfft = length(nfft);
    end
    
    z = 1i*2*pi/fs;
    % Initializations: spectral matrices have M rows, M columns and are calculated at each of the nfft frequencies
    H=zeros(M,M,nfft,N); % Transfer Matrix
    S=zeros(M,M,nfft,N); % Spectral Matrix
    
    for n=1:N
        p = size(Am_n{n},2)/M; % p is the order of the MVAR model
        if p==0
            n=n+1;
            p = size(Am_n{n},2)/M; % p is the order of the MVAR model
        end
        A = [eye(M) -Am_n{n}]; % matrix from which M*M blocks are selected to calculate spectral functions
        
        %% computation of spectral functions
        for ff=1:nfft % at each frequency
            
            %%% Coefficient matrix in the frequency domain
            As = zeros(M,M); % matrix As(z)=I-sum(A(k))
            for k = 1:p+1
                As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(ff));
                % indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix A (A(1) is in the second block, and so on)
            end
            
            %%% Transfer matrix
            H(:,:,ff,n)  = inv(As);
            
            %%% Spectral matrix
            S(:,:,ff,n)  = H(:,:,ff,n)*Su_n{n}*H(:,:,ff,n)'; % ' stands for Hermitian transpose
            
        end
        clear A As k p
    end
end
end

