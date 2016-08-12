%% DATA
close all
addpath('efficientDrawOctomap/3D')

load 'data/crawlerHandlerNoFloorNoWheels_cpp'

data = [PartMeans;SurfNormals];

prior.sigma = Prior.Sigma;
prior.gamma = 10;%Prior.Gamma;
prior.noiseVal = Prior.noiseVals;
prior.noiseGrad = Prior.noiseGrad;

R = 0.5;
cen = (sum(PartMeans')/length(PartMeans))';
prior.mean = @(x) [   1/2/R*((x-cen)'*(x-cen) - R^2);...
                1/R*((x(1)-cen(1)));...
                1/R*((x(2)-cen(2)));...
                1/R*((x(3)-cen(3)))];
surface = surfaceOctoMap(prior, data, 3, 3, true);

figure();
hold on;
points = [];
cols = [];
plot3(surface(:,1), surface(:,2),surface(:,3), 'go');
plot3(data(1,:), data(2,:), data(3,:), 'r.', 'MarkerSize',40);
quiver3(data(1,:), data(2,:), data(3,:), data(4,:), data(5,:),data(6,:));
grid;
axis equal
