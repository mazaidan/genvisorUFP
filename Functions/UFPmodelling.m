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
        hiddenLayerSize = [15 5];% 20;50;100;25;20;%15;
        net = fitnet(hiddenLayerSize);
        net.trainfcn = 'trainbr';'trainlm';%'trainbr'
        % Set up Division of Data for Training, Validation, Testing
        net.divideParam.trainRatio = 100/100;%70/100;
        net.divideParam.valRatio = 0/100;
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
    otherwise
        warning('Please choose a Model')
end