function [R,L,pvalue] = nanxcorr_ms(s1,s2,Lag,Corr_Type)
% function [L, R,pvalue] = nanxcorr (s1, s2, Lag);
% Function that allows obtaining the cross-correlation
% of a pair of time series containing gaps
%
% Input:
% s1 time series [vector]
% s2 time series [vector]
% lag number of lags to correlate (ex. 20)
%
% output
% L lag
% R correlation coefficient
% pvalue
%
% sam 04/16/2013
% Marco Sandoval Belmar 4/1/2018
[r,p]=corr(s1',s2','Type',Corr_Type,'rows','complete');  % correlation a lag == 0
% Performs the correlation for the different lags
L = 0; R = r; pvalue=p;
for i1 =1:1:Lag
    s11 = s1(1:end-i1);
    s21 = s2(i1+1:end);
    [c,pp] = corr(s11',s21','Type',Corr_Type,'rows','complete');
    R = [c;R];
    pvalue = [pp;pvalue];
    L = [-i1;L];
    
    clear s11 s21 c pp
    
    s21 = s2(1:end-i1);
    s11 = s1(i1+1:end);
    [c,pp] = corr(s11',s21','Type',Corr_Type,'rows','complete');
    R = [R;c];
    pvalue = [pvalue;pp];
    L = [L;i1];
    
    clear s21 s11 c
end
end