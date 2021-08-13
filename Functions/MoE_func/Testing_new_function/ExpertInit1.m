function [Prior, Mu, Sigma,loglik] = ...
    ExpertInit1(NumberExperts,X,T,TargetEst,VarTargetEst)
%
% This function initiate expert parameters first by
% Initializes the mixture of experts by using k-means to average (means) to
% initialize the gate , then using the hard labels of sample to train expert
%
% Inputs
% X            : Feature space for expert
% PHI          : Feature space for regression
% T            : Target for expert
% NumberExperts: Number eof experts
%
% Outputs
% Prior: Parameter for gaining function, analogous to posterior in Bayesian likelihood
% Mu: mean of liclyhood
% Sigma: Covariance matrix for each gate
% West: Linear repression parameter of each expert
% Var:

% K-means for gating function
%%[Prior, Mu, Sigma,DataIndex] = EM_init_kmeans(X', NumberExperts);

for expert=1:NumberExperts
    %[Prior, Mu(:,expert), Sigma(:,:,expert),DataIndex] = ...
    %    EM_init_kmeans((X{1,expert})', NumberExperts);
    [Prior, Mu{1,expert}, Sigma{1,expert}(:,:,expert),DataIndex] = ...
        EM_init_kmeans((X{1,expert})', NumberExperts);
    % Compute probability p(x|z)
    %%Pxz(:,expert) = gaussPDF(X{1,expert}, Mu(:,expert), Sigma(:,:,expert));
    Pxz{1,expert} = gaussPDF(X{1,expert}, Mu{1,expert}, Sigma{1,expert}(:,:,expert));    
    % Compute probability p(y|yhat,x)
    % Pyzx(:,expert)= ConditionalPDFExpert (T,PHI*West(:,expert),Var(expert));
    %%Pyzx(:,expert)= ConditionalPDFExpert (T{1,expert},TargetEst(:,expert),VarTargetEst(expert));
    Pyzx{1,expert}= ConditionalPDFExpert (T{1,expert},TargetEst{1,expert},VarTargetEst{1,expert});
    
    ProductLikelihood(1,expert) =Pyzx{1,expert}*Pxz{1,expert}*Prior';

end

% Compute the log likelihood
%%ProductLikelihood =Pyzx.*Pxz*Prior';
loglik = mean(log(ProductLikelihood));

% % % Error of  expert
% % Var=rand(1,NumberExperts);
% % % Parameters for expert
% % West=zeros(length(PHI(1,:)),NumberExperts);
% %
% % NumberDataPoint=length(X(:,1));

% % -> All below parts are changed to be 'TargetEst'
% % Matrix used for hard thresholds
% % IndictorMatrix=zeros(NumberDataPoint,NumberDataPoint);
% % for expert=1:NumberExperts
% %     IndictorMatrix(:,:)=diag(DataIndex==expert);
% %     West(:,expert)=(((PHI'*IndictorMatrix*PHI))^-1)*PHI'*IndictorMatrix*T;
% %     Var(1,expert)=sum((T(DataIndex==expert)...
% %         -PHI(DataIndex==expert,:)*...
% %         West(:,expert)).^2)/sum(DataIndex==expert);
% % end

end

function [Priors,Mu,Sigma,Data_id] = EM_init_kmeans(Data, nbStates)
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
    %Sigma(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
    Sigma{1,i}(:,:,i) = cov([Data(:,idtmp) Data(:,idtmp)]');
    % Add a tiny variance to avoid numerical instability
    %Sigma(:,:,i) = Sigma(:,:,i) + 1E-5.*diag(ones(nbVar,1));
    Sigma{1,i}(:,:,i) = Sigma{1,i}(:,:,i) + 1E-5.*diag(ones(nbVar,1));
end
Priors = Priors ./ sum(Priors);
end


function prob = gaussPDF(Data, Mu, Sigma)
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

function [ Out ] = ConditionalPDFExpert (RealValue,EstimatedValue,Var)
% This function produces a conditional mixture for each mixture of experts
%
% Input
% RealValue     : Column vector of Target value
% EstimatedValue: Column vector of Estimated value
% Var           : Variance of data
%
% Output
% Out:Probability of each expert

Out = (RealValue-EstimatedValue).^2;
Out=Out/(Var);
Out = exp(-0.5*Out) / (sqrt(2*pi*Var)+realmin);
end




