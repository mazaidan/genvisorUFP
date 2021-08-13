function prob = gaussPDF1(Data, Mu, Sigma)
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
%Data = Data - repmat(Mu,nbData,1);
prob = sum((Data*pinv(Sigma)).*Data, 2);
prob = exp(-0.5*prob) / sqrt((2*pi)^nbVar * (abs(det(Sigma))+realmin));
end