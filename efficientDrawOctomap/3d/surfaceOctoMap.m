function [ surfacePoints ] = surfaceOctoMap( prior, dataPoints, iterations)
%surfaceOctoMap Summary of this function goes here
%   Detailed explanation goes here

    X = dataPoints(1:3, :);
    m = length(X); 
    
    data = zeros(1,length(dataPoints));
    data = [data ; dataPoints(4:6,:)]; 
    f = reshape(data, [1,4*length(X)])';

    sigma = prior.sigma;
    gamma = prior.gamma;
    noiseVal = prior.noiseVal;
    noiseGrad = prior.noiseGrad;
    mean = prior.mean;
    
    K = ComputeFullKder(sigma, gamma, X, noiseVal, noiseGrad);
    for i = 1:m
        mu(((i-1)*4+1):((i-1)*4 +4)) = mean(X(:,i));
    end

    %% Efficiend draw
    xLimits = [min(X(1,:)), max(X(1,:))]*1.2;
    yLimits = [min(X(2,:)), max(X(2,:))]*1.2;
    zLimits = [min(X(3,:)), max(X(3,:))]*1.2;
    centroid = [sum(xLimits)/2, sum(yLimits)/2,sum(zLimits)/2];

    % 2 valid, 1 new, 0 invalid.
    root = {2,  {}, xLimits, yLimits, centroid, -1, zLimits};
    Qmat = K\(f - mu');

    evalFun = @(x) mean(x) + ComputeKderX1X2(sigma, gamma, x, X)*Qmat;

    root = expandCell(root,evalFun, 2);
    for iter=1:iterations
        % Expand tree
        root = expandCell(root, evalFun, 1); 
        % Validate branches
        root = validatePoints(root, root, evalFun);
    end
    root = validatePoints(root, root, evalFun);
  
    surfacePoints = getPointsTree(points, cols, root, true);
end

