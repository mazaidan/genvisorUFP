function [Prior, Mu, Sigma,Likelihood] =...
    TrainMoE2(Xtr1,Ttr1,yTot,sig2Tot,...
    NumberExperts,MaxIterations)

% This code implements an alternative model for mixtures of experts using
% a Gaussian for the gating function.
% The special form is used so that the maximization with respect to
% the parameters  of the gating network can be handled analytically.
% Thus, a single-loop EM can be used, and no learning stepsize is required
% to guarantee convergence.
% The method is based on:
% "An Alternative Model for Mixtures of Experts"
% by Lei Xu, Michael!. Jordan, Geoffrey E. Hinton
% The code also uses Sylvain Calinon Expectation-Maximization (EM) algorithm
% On Learning, Representing and Generalizing a Task in a Humanoid Robot,
% by S. Calinon and F. Guenter and A. Billard,

% Inputs -----------------------------------------------------------------
% T: N x 1 array representing N data points of targets that are meant
% to be Estimated.
% PHI: N x M array representing N data points of M dimensions that
% consist of some nonlinear transformation of the array X
% X: N x D array representing N datapoints of D dimensions
% NumberExperts: Integer indicating number of experts
% RegulationTermDesign: Value for diagonal regularisation term for the design matrix
% RegulationTermCovariance: Regularization term for design matrix
% MaxIterations: maxmum number of maximum number of iterations: Default 100
% Outputs ----------------------------------------------------------------
% Gate Parameters
% Priors:  1 x NumberExperts array representing the prior probabilities of
%          the K GMM components.
% Mu    :  D x NumberExperts array representing the centers of
%          the K GMM components.
% Sigma :  D x D x NumberExperts array representing the covariance matrices of
%          the K GMM components.

% Author: Martha Arbayani Zaidan (19/03/1986)
% modified from Joseph Santarcangelo (29/10/2014)
% Version one
% Also see: ExpertInit,gaussPDF,

if (nargin<=6)
    MaxIterations=100;
end

RegulationTermCovariance=0.01;

% The number of data samples and dimension may vary for each data set
NumberDataPoint=cell(1,NumberExperts);
DimensionsofData=cell(1,NumberExperts);
for n=1:NumberExperts
    NumberDataPoint{1,n}=length(Xtr1{1,n}(:,1));
    DimensionsofData{1,n}=length(Xtr1{1,n}(1,:));
end
%% Criterion to stop the EM iterative update
loglik_threshold = 1e-10;

% Initialization of Gate Parameters --------------------------------------%
[Prior, Mu, Sigma,loglik] = ...
    ExpertInit2(NumberExperts,Xtr1,Ttr1,yTot,sig2Tot);

loglik_old = loglik ;
nStep=1;

% Initialization for E-step
% Likelihood of X given of expert each row represents 
% the probability of a data point

% Probabilities for each gate for every sample xn, each row represents a
% sample, each column represents a sample, for column k.
% [p(x1|z=k),p(x2|z=k),..,p(xN|z=k)]:
Likelihood=[];
while 1
    %-------------- E-step ---------------%
    Pxz=cell(1,NumberExperts);
    Pyzx=cell(1,NumberExperts);
    Pyx=cell(1,NumberExperts);
    for expert=1:NumberExperts
        % Compute probability p(x|z) - gating net num eq(10)
        Pxz{1,expert} = gaussPDF1(Xtr1{1,expert}, Mu(:,expert), Sigma(:,:,expert));
        
        % Compute probability p(y|yhat,x) expert parameter - BNN num eq(10)
        Pyzx{1,expert}= ConditionalPDFExpert1(Ttr1{1,expert},...
            yTot{1,expert},sig2Tot{1,expert}(expert));
        
        % Compute posterior probability p(z|x,y) or:
        % h(y|x)=p(z|x,y) -> eq(10):
        % num eq(10):
        Pyx0=Pyzx{1,expert}.*Pxz{1,expert}.*repmat(Prior(1,expert),NumberDataPoint{1,expert},1);
        Pyx{1,expert}=Pyx0./(sum(Pyx0')'+realmin);
        Pyx1=Pyx{1,expert};
        E(expert) = sum(Pyx1); % -> eq(11)
    end
    
    
    %--- M-step ---%
    for expert=1:NumberExperts
        % Update the priors -> Priors = alpha -> eq(13)
        Priors(expert) = E(expert) / NumberDataPoint{1,expert};
        % Update the centers -> eq(14a)
        Mu(:,expert) = Pyx{1,expert}'*Xtr1{1,expert} /(E(expert)+realmin);
        
        % Update the covariance matrices - eq(14b)
        Xtemp = Xtr1{1,expert} - repmat(Mu(:,expert)',NumberDataPoint{1,expert},1);
        Sigma(:,:,expert) = repmat(Pyx{1,expert}',DimensionsofData{1,expert},1).* Xtemp'* Xtemp ...
            / (E(expert)+realmin);
        % Add a tiny variance to avoid numerical instability
        IregCov=RegulationTermCovariance*eye(length(Xtr1{1,expert}(1,:)));
        Sigma(:,:,expert) = Sigma(:,:,expert) +IregCov;
    end
    
    ProductLikelihood2=cell(1,NumberExperts);
    ProductLikelihood3=[];
    for expert=1:NumberExperts
        % Compute probability p(x|z) - gating net num eq(10)
        Pxz{1,expert} = gaussPDF1(Xtr1{1,expert}, Mu(:,expert), Sigma(:,:,expert));
        
        % Compute probability p(y|yhat,x) - BNN num eq(10)
        Pyzx{1,expert}= ConditionalPDFExpert1(Ttr1{1,expert},yTot{1,expert},sig2Tot{1,expert}(expert));        
        
        % Compute the log likelihood (the matrix operation takes care of the sum)
        ProductLikelihood2{1,expert} = Pyzx{1,expert}.*Pxz{1,expert}*Prior(1,expert) ; % eq(9)
        
        % replace smallest number by realmin (to prevent resulting zero):
        ProductLikelihood2{1,expert}(find(ProductLikelihood2{1,expert}<realmin)) = realmin;
        ProductLikelihood2{1,expert}=log(ProductLikelihood2{1,expert});
        ProductLikelihood3=[ProductLikelihood3;ProductLikelihood2{1,expert}];
    end
    
    loglik=mean(ProductLikelihood3); %between eq(9) and eq(10)
    Likelihood(nStep)=loglik;
    
    % Stop the process depending on the increase of the log likelihood
    if abs((loglik/loglik_old)-1) < loglik_threshold || nStep==MaxIterations
        break;
    end
    loglik_old = loglik;
    nStep = nStep+1;
end

end

%%
function [Prior, Mu, Sigma,loglik] = ...
    ExpertInit2(NumberExperts,Xtr1,Ttr1,yTot,sig2Tot)

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

% We tabulate all data to be clustered by K-means clustering
Xtr1Tot=[];
for expert=1:NumberExperts
    Xtr1Tot=[Xtr1Tot;Xtr1{1,expert}];
end

% K-means for gating function
[Prior, Mu, Sigma,DataIndex] = ...
    EM_init_kmeans1(Xtr1Tot', NumberExperts);

ProductLikelihood1=[];
Pxz=cell(1,NumberExperts);
Pyzx=cell(1,NumberExperts);
for expert=1:NumberExperts
    
    % Compute probability p(x|z)
    Pxz{1,expert} = gaussPDF1(Xtr1{1,expert}, Mu(:,expert), Sigma(:,:,expert));
    
    % Compute probability p(y|yhat,x):
    Pyzx{1,expert}= ConditionalPDFExpert1 (Ttr1{1,expert},yTot{1,expert},sig2Tot{1,expert}(expert));
    
    % Compute the log likelihood:
    ProductLikelihood =log(((Pyzx{1,expert}).*Pxz{1,expert})*Prior(1,expert));
    ProductLikelihood1 =[ProductLikelihood1;ProductLikelihood];
    
end

% Compute the log likelihood (final):
loglik = mean(ProductLikelihood1);

end





