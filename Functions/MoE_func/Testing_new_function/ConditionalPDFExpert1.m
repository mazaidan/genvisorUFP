function [ Out ] = ConditionalPDFExpert1 (RealValue,EstimatedValue,Var)
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