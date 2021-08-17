function [TotalOutput,TotalVar,ExpertOutput,GateProbability] = ...
    OutputMoE3(Prior,Mu,Sigma,Xtr,Ttr,Xte,netTot,Models)

% This function provides the output and output components of the mixture of 
% experts code for linear regression using a Gaussian mixture model 
% as the gating function 
%
% Inputs -----------------------------------------------------------------%
% PHI: N x M array representing N data points of M dimensions that 
% consist of some nonlinear transformation of the array X  
% X: N x D array representing N datapoints of D dimensions   
% Gate Parameters 
% Priors:  1 x NumberExperts array representing the prior probabilities of the K GMM 
% components.
% Mu:      D x NumberExperts array representing the centers of the K GMM components.
% Sigma:   D x D x NumberExperts array representing the covariance matrices of the 
% Expert Parameters 
% West: MxnumberExperts array Representing the number of parameters, each column represents the parameters for a 
% different expert, analogues to the parameters in linear regression   

% Ouput-------------------------------------------------------------------%
% TotalOutput: Column vector representing the output of the network, 
% where each row represents the output of a different observation 
% ExpertOutput:(Number of Observations)x(Number of sample) matrix where 
% each column represents the output of each expert and 
% the row corresponds to an observation
% GateProbability:(Number of Observations)x(Number of sample) matrix where 
% each column represents expert and 
% the row corresponds to the probability observation   
% Author:Martha Arbayani Zaidan (19/03/2016) 
% modified the code of Joseph Santarcangelo, 2014 
% Contributions :Sylvain Calinon, 2
 
M=length(Xte(:,1));
NumberExperts=length(Prior);
ExpertOutput=zeros(M,NumberExperts);
ExpertOutputsig2=zeros(M,NumberExperts);

GateProbability=zeros(M,NumberExperts);
TotalOutput=zeros(M,1);
TotalVar=zeros(M,1);

Normalization=zeros(M,1);
%Output of each expert
%ExpertOutput= TargetEst; % ExpertOutput=PHI*West; % perhaps, this should be prob

for n=1:NumberExperts
    %%%[ExpertOutput(:,n), ExpertOutputsig2(:,n)] = ...
    %%%    netevfwd(mlppak(net{1,n}), net{1,n}, Xtr{1,n}, Ttr{1,n}, Xte);
    % [Ypred0,Ypred_std0] = predict(mdl,X1{1,n}); % [ypred,yci] = predict(mdl,Xnew)
    
    Model = Models{1,n};
    
    switch Model
        case {'LM1'}
           [Y0,Ystd0] = predict(netTot{1,n},Xte);
           ExpertOutput(:,n)    = Y0;
           ExpertOutputsig2(:,n)= Ystd0(:,2);
        case {'ANN1'}
            net = netTot{1,n};
            outputs = net(Xte');
            ExpertOutput(:,n)     = outputs';
            ExpertOutputsig2(:,n) = ones(size(outputs'));
        otherwise
            warning('Unexpected plot type. No plot created.')
    end
    
    %Probability of each gate for every sample given the Gaussian model
    GateProbability(:,n)=gaussPDFOut(Xte, Mu(:,n), Sigma(:,:,n))*Prior(n);
    
    % Total output of each expert waited by 
    % the probability of that observations:
    TotalOutput=TotalOutput+ExpertOutput(:,n).*GateProbability(:,n);
    TotalVar=TotalVar+ExpertOutputsig2(:,n).*GateProbability(:,n);
    
    %Normalization values
    Normalization= Normalization+GateProbability(:,n);
end

% Gate  normalized gate values for each expert and each samples
GateProbability=GateProbability./repmat(Normalization,1,NumberExperts);
% Normalize output
TotalOutput=TotalOutput./Normalization;
TotalVar=TotalVar./Normalization;
end

function prob = gaussPDFOut(Data, Mu, Sigma)
%
% This function computes the Probability Density Function (PDF) of a
% multivariate Gaussian represented by means and covariance matrix.
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% Inputs -----------------------------------------------------------------
%   o Data:  N x D array representing N datapoints of D dimensions.
%   o Mu:    D x 1 array representing the centers of the K GMM components.
%   o Sigma: D x D x K array representing the covariance matrices of the
%            K GMM components.
% Outputs ----------------------------------------------------------------
%   o prob:  1 x N array representing the probabilities for the
%            N datapoints.

[nbData,nbVar] = size(Data);

Data = Data - repmat(Mu',nbData,1);
prob = sum((Data*pinv(Sigma)).*Data, 2);
prob = exp(-0.5*prob) / sqrt((2*pi)^nbVar * (abs(det(Sigma))+realmin));
end

