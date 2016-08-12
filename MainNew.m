%% DATA
close all
addpath('efficientDrawOctomap/3D')
addpath('surfacePlottingGradient/covMat')
addpath('grasping')

load 'data/crawlerHandlerNoFloorNoWheels_cpp'

data = [PartMeans;SurfNormals];

%% GPIS parameters
prior.sigma = Prior.Sigma;
prior.gamma = 5;%Prior.Gamma;
prior.noiseVal = Prior.noiseVals;
prior.noiseGrad = Prior.noiseGrad;

% R = 0.5;
% cen = (sum(PartMeans')/length(PartMeans))';
% prior.mean = @(x) [   1/2/R*((x-cen)'*(x-cen) - R^2);...
%                 1/R*((x(1)-cen(1)));...
%                 1/R*((x(2)-cen(2)));...
%                 1/R*((x(3)-cen(3)))];

Prior.pos = [0.1650, 0.0766, -3.8024]';
Prior.param = [0.6,0.4,0.8];
Prior.rot = [0,0,0];
cs2 = @(r) (cos(norm(r)/2)^2);
si2 = @(r) (sin(norm(r)/2)^2);
cs = @(r) (cos(norm(r)/2));
si = @(r) (sin(norm(r)/2));
R = @(r) ([     (r(1)^2 - r(2)^2 - r(3)^2) * si2(r)  + norm(r)^2 * cs2(r)    , 2 * si(r)* (r(1)*r(2)*si(r) + norm(r) * r(3) * cs(r)) , 2 * si(r)* (r(1)*r(3)*si(r) - norm(r) * r(2) * cs(r)) ;
             2 * si(r)* (r(1)*r(2)*si(r) - norm(r) * r(3) * cs(r)) ,  (r(2)^2 - r(3)^2 - r(1)^2) * si2(r)  + norm(r)^2 * cs2(r)    ,  2 * si(r)* (r(2)*r(3)*si(r) + norm(r) * r(1) * cs(r));
             , 2 * si(r)* (r(1)*r(3)*si(r) + norm(r) * r(2) * cs(r)) , 2 * si(r)* (r(2)*r(3)*si(r) - norm(r) * r(1) * cs(r)),   (r(3)^2 - r(1)^2 - r(2)^2) * si2(r)  + norm(r)^2 * cs2(r)  ]/norm(r)^2);
priorRotation = R(Prior.rot);

A = diag([1/Prior.param(1)^2, 1/Prior.param(2)^2, 1/Prior.param(3)^2]);
prior.mean = @(x) [Prior.param(1)/2 * ((x-Prior.pos)'* priorRotation' * A * priorRotation * (x-Prior.pos) - 1);...
                Prior.param(1) *  priorRotation' * A * priorRotation * (x-Prior.pos)];
            
%% Surface plotting
surface = surfaceOctoMap(prior, data, 3, 2, true);
%%
figure();
hold on;
points = [];
cols = [];
plot3(surface(:,1), surface(:,2),surface(:,3), 'go');
plot3(data(1,:), data(2,:), data(3,:), 'r.', 'MarkerSize',10);
quiver3(data(1,:), data(2,:), data(3,:), data(4,:), data(5,:),data(6,:));
grid;
axis equal

%% Sample grasping points and evaluate points

% Get initial point.
% candidatePoint = [0.0874, 0.0218, -3.277]';
candidatePoint = PartMeans(:,10);
% Converge to surface.
X = data(1:3, :);    
m = length(X);     
f = zeros(1,length(data));
f = [f ; data(4:6,:)]; 
f = reshape(f, [1,4*length(X)])';
sigma = prior.sigma;
gamma = prior.gamma;
noiseVal = prior.noiseVal;
noiseGrad = prior.noiseGrad;
mean = prior.mean;
K = ComputeFullKder(sigma, gamma, X, noiseVal, noiseGrad);
for i = 1:m
    mu(((i-1)*4+1):((i-1)*4 +4)) = mean(X(:,i));
end
Qmat = K\(f - mu');
evalFun = @(x) mean(x) + ComputeKderX1X2(sigma, gamma, x, X)*Qmat;
% 
% fGradCurrent = evalFun(x);
% x = x - fGradCurrent(1)/norm(fGradCurrent(2:end)) * fGradCurrent(2:end)*0.1;
% plot3(x(1), x(2), x(3), 'x','MarkerSize',20)

% Compute contact points. By now because the movement of the gripper is
% unkown, take a circle of radius Rg and discretice it to get the point in
% the surface.

nPointsInCircle = 20;
points = pointsInCircle(0.1, candidatePoint, -pi/4,pi/2,nPointsInCircle);
plot3(points(1,:), points(2,:), points(3,:));
hold on;

valuesCircle = [];
for i =1:nPointsInCircle
    newEval = evalFun(points(:,i));
    valuesCircle  = [valuesCircle , newEval];
    if newEval(1) < 0
        plot3(points(1,i), points(2,i), points(3,i), 'gx','MarkerSize',20);
    else
        plot3(points(1,i), points(2,i), points(3,i), 'rx','MarkerSize',20);
    end
end

% Evaluate grasping point. Evaluate grasping point using the Grasp Wrench
% Space.


