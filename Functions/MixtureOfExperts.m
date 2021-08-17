%%%%%%%% MoE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MoE data division: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TotalOutput,TotalVar] = MixtureOfExperts(X,Y,Xt,Models,NumberExperts,IDX1)

p1 = cell(1,NumberExperts);
for n=1:NumberExperts
    p1{1,n} = IDX1{1,n};
end

X1      = cell(1,NumberExperts);
Y1      = cell(1,NumberExperts);
yTot    = cell(1,NumberExperts);
sig2Tot = cell(1,NumberExperts);
netTot  = cell(1,NumberExperts);
for n=1:NumberExperts
    X1{1,n}=X(p1{1,n},:);
    Y1{1,n}=Y(p1{1,n},:);
    Model = Models{1,n};
    
    [Ypred0,Ypred_std0,M] = UFPmodelling(X1{1,n},Y1{1,n},X1{1,n},Model);
       
    yTot{1,n}=Ypred0;
    sig2Tot{1,n}=Ypred_std0; 
    netTot{1,n}=M;
end


MaxIterations=100;
[Prior, Mu, Sigma,Likelihood] =...
    TrainMoE2(X1,Y1,yTot,sig2Tot,...
    NumberExperts,MaxIterations);

for n=1:3
    Prior0=Prior;
    Mu0=Mu;
    Sigma0=Sigma;
    Likelihood0=Likelihood;
    
    if  Likelihood0(end)>Likelihood(end)
        Prior=Prior0;
        Mu=Mu0;
        Sigma=Sigma0;
        Likelihood=Likelihood0;
    end
end

[TotalOutput,TotalVar,ExpertOutput,GateProbability] = ...
    OutputMoE3(Prior,Mu,Sigma,X1,Y1,Xt,netTot,Models);
