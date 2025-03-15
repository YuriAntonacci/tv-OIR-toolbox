%% Computation of the gradient of OIR and OIR for blocks of processes given the spectral matrix 

%%% INPUT
% S: spectral matrix at each time step QxQx nfft xN
% Mv: vector specifying how many processes for each block
% nrange: matrix with rows=number of ranges, columns=2 (limits of the bands measured in Hz)

%%% OUTPUT
% dO, dOf: time domain and spectral delta OIR at each time step
% OIR, OIRf: time domain and spectral OIR at each time step
% OIR_band: integrated OIR values in the given spectral ranges at each time
% step

function out = tv_oir(S_n,Mv,fs,range_i)

% check inputs
if nargin < 2
    error('Not enough input arguments')
end

M = length(Mv); % number of the blocks
nfft = size(S_n,3);

allM=cell(M,1); % number of combination for all orders
dO=cell(M,1); % deltaOIR
dOf=cell(M,1); % spectral deltaOIR 

for N=3:M
    comb=nchoosek(1:M,N); % all multiplets of size N
    allM{N}=comb;
    for n=1:size(comb,1) % cycle on all the combinations of that order
        ij=comb(n,:); % select the combination
        for k=1:size(comb,2) % vary target inside the multiplet
            for t=1:size(S_n,4)
            out1=sir_deltaO(S_n(:,:,:,t),Mv,ij,ij(k)); % this is delta OIR
            dO12(n,k,t)=out1.dO;
            dO12f(n,k,t,:)=out1.dOf; 
            end
        end
    end
    dO{N-1}=dO12;
    dOf{N-1}=dO12f;
end

OIR=cell(M,1); % OIR
OIRf=cell(M,1); % spectral OIR
OIR{3}=squeeze(dO{2}(:,1,:)); % other columns are the same
OIRf{3}=squeeze(dOf{2}(:,1,:,:));

for N=4:M
    t1=1; % index of Xj (can be any number btw 1 and N)
    t2=setdiff(1:N,t1); % index of X-j
    for cnt=1:size(allM{N},1) % cycle on all combinations of order N
        tmp=allM{N}(cnt,:);
        iN1=find(sum(tmp(t2)==allM{N-1},2)==N-1); % position of X-j in the O-info one step back
        OIR{N}(cnt,:) = squeeze(OIR{N-1}(iN1,:)) + squeeze(dO{N-1}(cnt,t1,:))';
        OIRf{N}(cnt,:,:) = OIRf{N-1}(iN1,:,:) + shiftdim(dOf{N-1}(cnt,t1,:,:),1);
    end
end

%%% integral inside spectral bands
if nargin == 4
    nfft=size(S_n,3);
    for nn=1:size(range_i,1)
        nrange(nn,:)=round((nfft*2/fs)*range_i(nn,:));
    end
    nr = size(nrange,1);
    for ir=1:nr
        for N=3:M
            ncomb = size(OIRf{N},1);
            for cnt=1:ncomb
                OIR_band{N,ir}(cnt,:) = sum(OIRf{N}(cnt,:,nrange(ir,1):nrange(ir,2)),3)/nfft;
            end
        end
    end
    out.OIR_band=OIR_band; % integrated values
end

%% OUTPUT
out.dO=dO; out.dOf=dOf; % delta OIR
out.OIR=OIR; out.OIRf=OIRf; % OIR
out.allM=allM;

end

