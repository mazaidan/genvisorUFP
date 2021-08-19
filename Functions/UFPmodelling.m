function [Ypred,Ypred_std,M] = UFPmodelling(X,Y,Xt,Model)

switch Model
    case {'LM1'}
        disp('Standard Linear Model')
        mdl = fitlm(X,Y);
        %mdl = fitlm(X{1,n},Y{1,n});
        %Ypred = predict(mdl,Xt);
        [Ypred,Ypred_std0] = predict(mdl,Xt);
        Ypred_std = Ypred_std0(:,2); % because Ypred_std0 generates two cols
        M = mdl;
    case {'ANN1'}
        disp('Shallow Neural Network')
        % Create a Fitting Network
        hiddenLayerSize = [10 5];% 20;50;100;25;20;%15;
        net = fitnet(hiddenLayerSize);
        net.trainfcn = 'trainlm';%'trainbr'
        % Set up Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = 95/100;%70/100;
        net.divideParam.valRatio = 5/100;
        net.divideParam.testRatio = 0/100;
        % Train the Network
        [net,tr] = train(net,X',Y');
        % Test the Network
        outputs = net(Xt');
        Ypred   = outputs';
        Ypred_std = ones(size(Ypred));
        M = net;
    case {'ANN2'}
        disp('Generalized Regression Neural Networks')
        % https://www.mathworks.com/help/deeplearning/ug/generalized-regression-neural-networks.html
        net = newgrnn(X',Y');
        %net = train(net,X',Y');
        outputs = sim(net,Xt');
        Ypred   = outputs';
        Ypred_std = ones(size(Ypred));
        M = net;
    case {'BNN1'}
        disp('Bayesian Neural Networks')
        % Set up network parameters.
        nin = size(X,2);		% Number of inputs.
        nhidden = 3;		% Number of hidden units.
        nout = 1;		% Number of outputs.
        alpha = 0.01;		% Initial prior hyperparameter.
        beta_init = 50.0;	% Initial noise hyperparameter.
        
        % Create and initialize network weight vector.
        net = mlp(nin, nhidden, nout, 'linear', alpha, beta_init);
        
        % Set up vector of options for the optimiser.
        nouter = 3;			% Number of outer loops.
        ninner = 1;			% Number of innter loops.
        options = zeros(1,18);		% Default options vector.
        options(1) = 1;			% This provides display of error values.
        options(2) = 1.0e-7;		% Absolute precision for weights.
        options(3) = 1.0e-7;		% Precision for objective function.
        options(14) = 500;		% Number of training cycles in inner loop.
        
        % Train using scaled conjugate gradients, re-estimating alpha and beta.
        for k = 1:nouter
            net = netopt(net, options, X, Y, 'scg');
            [net, gamma] = evidence(net, X, Y, ninner);
            fprintf(1, '\nRe-estimation cycle %d:\n', k);
            fprintf(1, '  alpha =  %8.5f\n', net.alpha);
            fprintf(1, '  beta  =  %8.5f\n', net.beta);
            fprintf(1, '  gamma =  %8.5f\n\n', gamma);
            disp(' ')
        end
        
        %fprintf(1, 'true beta: %f\n', 1/(noise*noise));
        
        % Evaluate error bars.
        %[y, sig2] = netevfwd(mlppak(net), net, x, t, plotvals);\
        [y, sig2] = netevfwd(mlppak(net), net, X, Y, Xt);
        sig = sqrt(sig2);
        
        Ypred   = y;
        Ypred_std = sig;
        M = net;
        
    case {'BLM'}
        disp('Bayesian Modelling')
        % https://www.mathworks.com/help/econ/bayesian-linear-regression-models.html
        p = size(X,2);
        %PriorMdl = bayeslm(p,'ModelType','lasso');
        %PriorMdl.Lambda = 10.*ones(p+1,1);%[10; 1e5; 10];
        %%%PriorMdl = bayeslm(p,'ModelType','semiconjugate'); PriorMdl.A = 5; PriorMdl.A = 0.1;
        %%%PriorMdl = bayeslm(p,'ModelType','semiconjugate','Intercept',false); PriorMdl.A = 5; PriorMdl.A = 0.1;
        %PriorMdl = bayeslm(p,'ModelType','mixconjugate');
        %PriorMdl = bayeslm(p,'ModelType','mixsemiconjugate');
        LS = 0;
        if LS ==1
            if p ==4
                PriorMdl = bayeslm(p,'ModelType','lasso','Intercept',false,'VarNames',["PM25" "T" "RH" "P"]);
                PriorMdl.Lambda = [0.1 1e2 1e2 1e3];
            else
                PriorMdl = bayeslm(p,'ModelType','lasso','Intercept',false);
            end
        else
            %PriorMdl = bayeslm(p,'ModelType','diffuse','Intercept',false);
            %PriorMdl = bayeslm(p,'ModelType','diffuse','Intercept',true);
            %PriorMdl = bayeslm(p,'ModelType','lasso','Intercept',true);
            %PriorMdl = bayeslm(p,'ModelType','mixsemiconjugate');
            %PriorMdl = bayeslm(p,'ModelType','mixconjugate');
            PriorMdl = bayeslm(p,'ModelType','conjugateblm','Intercept',true);
            %PriorMdl = bayeslm(p,'ModelType','semiconjugateblm','Intercept',true);
        end
        
        rng(1);
        PosteriorMdl = estimate(PriorMdl,X,Y);
        [yfit,ycov] = forecast(PosteriorMdl ,Xt);
        
        Ypred   = yfit;
        Ypred_std = sqrt(ycov);
        M = PosteriorMdl;
        
    otherwise
        warning('Please choose a Model')
end