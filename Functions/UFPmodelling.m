function [Ypred] = UFPmodelling(X,Y,Xt,Model)

switch Model
    case {'LM1'}
        disp('Standard Linear Model')
        mdl = fitlm(X,Y);
        Ypred = predict(mdl,Xt);
    case {'ANN1'}
        disp('Shallow Neural Network')
        % Create a Fitting Network
        hiddenLayerSize = 20;50;100;25;20;%15;
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
    otherwise
        warning('Please choose a Model')
end