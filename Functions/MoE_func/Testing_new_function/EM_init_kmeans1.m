function [Priors,Mu,Sigma,Data_id] = EM_init_kmeans1(Data, nbStates)
%
% This function initializes the parameters of a Gaussian Mixture Model
% (GMM) by using k-means clustering algorithm.
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% Inputs -----------------------------------------------------------------
%   o Data:     D x N array representing N datapoints of D dimensions.
%   o nbStates: Number K of GMM components.
% Outputs ----------------------------------------------------------------
%   o Priors:   1 x K array representing the prior probabilities of the
%               K GMM components.
%   o Mu:       D x K array representing the centers of the K GMM components.
%   o Sigma:    D x D x K array representing the covariance matrices of the
%               K GMM components.
% Comments ---------------------------------------------------------------
%   o This function uses the 'kmeans' function from the MATLAB Statistics
%     toolbox. If you are using a version of the 'netlab' toolbox that also
%     uses a function named 'kmeans', please rename the netlab function to
%     'kmeans_netlab.m' to avoid conflicts.

[nbVar, nbData] = size(Data);

% Below we use 'kmeans' function from the MATLAB Statistics toolbox:
[Data_id, Centers] = kmeans(Data', nbStates);
Mu = Centers';
for i=1:nbStates
    idtmp = find(Data_id==i);
    %Priors(i) = length(idtmp);
    Priors(i) = length(idtmp);
    Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
    %Sigma{1,1,i} = cov([Data(:,idtmp) Data(:,idtmp)]');
    % Add a tiny variance to avoid numerical instability
    Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
    %Sigma{1,1,i} = Sigma{1,1,i} + 1E-5.*diag(ones(nbVar,1));
end
Priors = Priors ./ sum(Priors);
end


% function prob = gaussPDF1(Data, Mu, Sigma)
% %
% % This function computes the Probability Density Function (PDF) of a
% % multivariate Gaussian represented by means and covariance matrix.
% %
% % Author:	Sylvain Calinon, 2009
% %			http://programming-by-demonstration.org
% %
% % Inputs -----------------------------------------------------------------
% %   o Data:  N x D array representing N datapoints of D dimensions.
% %   o Mu:    D x 1 array representing the centers of the K GMM components.
% %   o Sigma: D x D x K array representing the covariance matrices of the 
% %            K GMM components.
% % Outputs ----------------------------------------------------------------
% %   o prob:  1 x N array representing the probabilities for the 
% %            N datapoints.     
% 
% [nbData,nbVar] = size(Data);
% 
% Data = Data - repmat(Mu',nbData,1);
% prob = sum((Data*pinv(Sigma)).*Data, 2);
% prob = exp(-0.5*prob) / sqrt((2*pi)^nbVar * (abs(det(Sigma))+realmin));
% end
% 
% function [ Out ] = ConditionalPDFExpert1 (RealValue,EstimatedValue,Var)
% % This function produces a conditional mixture for each mixture of experts
% %
% % Input
% % RealValue     : Column vector of Target value
% % EstimatedValue: Column vector of Estimated value
% % Var           : Variance of data
% %
% % Output
% % Out:Probability of each expert
% 
% Out = (RealValue-EstimatedValue).^2;
% Out=Out/(Var);
% Out = exp(-0.5*Out) / (sqrt(2*pi*Var)+realmin);
% end
% 
% 
% 
% 
