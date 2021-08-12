function [Xnorm] = normalize_UFPsensors(X,vars_type)

switch vars_type
    case {'Temp'}
        disp('Normalizing Temp');inmin = 10; inmax = 40;
    case {'RH'}
        disp('Normalizing RH');inmin = 10; inmax = 50;
    case {'P'}
        disp('Normalizing P');inmin = 890; inmax = 910;
    case {'PM25'}
        disp('Normalizing PM2.5');
    otherwise
        warning('Unexpected vars_type.')
end

switch vars_type
    case {'Temp','RH','P'}
        l = -1; u = 1; % For normalization
        Xnorm = l + [(X-inmin)./(inmax-inmin)].*(u-l);
    case {'PM25'}
        Xnorm = log10(X);
    otherwise
        warning('Unexpected vars_type.')
end