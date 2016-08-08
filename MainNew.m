%% DATA
close all
addpath('covMat')
addpath('grasping')
addpath('efficientDrawOctomap/3D')

load 'crawlerHandlerNoFloorNoWheels_cpp'

data = [PartMeans;SurfNormals];

prior.sigma = Prior.Sigma;
prior.gamma = 10;%Prior.Gamma;
prior.noiseVal = Prior.noiseVals;
prior.noiseGrad = Prior.noiseGrad;

R = Prior.param(1);
cen = (sum(PartMeans')/length(PartMeans))';
prior.mean = @(x) [   1/2/R*((x-cen)'*(x-cen) - R^2);...
                1/R*((x(1)-cen(1)));...
                1/R*((x(2)-cen(2)));...
                1/R*((x(3)-cen(3)))];
% prior.mean = @(x) [0.5,0,0,0]';
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

[x,y,z] = sphere(20);
x = x * R;
y = y * R;
z = z * R;

for i = 1:21
    for j=1:21
        x(i,j) = x(i,j)+ cen(1);
    end
end
for i = 1:21
    for j=1:21
        y(i,j) = y(i,j)+ cen(2);
    end
end
for i = 1:21
    for j=1:21
        z(i,j) = z(i,j)+ cen(3);
    end
end
surf(x,y,z);

% % For bmw_11
% cx = -6.1655;
% cy = -0.0472;
% cz = -3.6693;
% 
% a = 1.0763;
% b = 2.3691;
% c = 1.3254;
% r = [2 * pi , 0, 0];
% Prior.type = 'N';
% Prior.param(1) = 0.1;
% Prior.Gamma = 2;
% Prior.Sigma = 0.05;
% Prior.noiseGrad = 0.05;
% Prior.noiseVals = 0.002;
% % 
% % Prior = struct('pos',[cx cy cz],'type',prior_type, 'param', [a b c], 'rot', r);
% %% SURFACE
% [meanValue, meanGrad] = computePriorFunctions(Prior)
% 
% dist = 0.005;
% initPoints = locations;%r * [-4;0;-2];

% 
% [faces, vertices] = computeSurface(locations, surfNormals, ...
%     Prior, ...
%     meanValue, meanGrad, initPoints(:,1:5:end), dist, false);
% 
% figure
% hold on
% set(gca,'dataaspectratio',[1 1 1])
% 
% plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',5);
% quiver3(locations(1,:),locations(2,:),locations(3,:),...
%     surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',1,'color','r');
% quiver3(locations(1,:),locations(2,:),locations(3,:),...
%     surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',1,'color','r');
% 
% patch('faces',faces,'vertices',vertices,...
%     'facecolor',[0.5 0.5 0.5], ...
%     'edgecolor', 'none', ...
%     'facelighting','gouraud');
% camlight('headlight')
% set(gca,'view',[46.8000   18.8000]);
% light('Position',[-1 -1 0])

% %% Grasping
% sigma = Prior.Sigma;
% gamma = Prior.Gamma;
% noiseVals = Prior.noiseVals;
% noiseGrad = Prior.noiseGrad;
% 
% covMatData = ComputeCovMatFull(sigma,gamma,locations,noiseVals,noiseGrad);
% invCovMatData = inv(covMatData);
% covPoint = @(x) ComputeCovMatFull(sigma, gamma, x, noiseVals, noiseGrad);
% uncertaintyGp = @(x) covPoint(x) - CovMatStar(sigma, gamma, x, locations)*invCovMatData*CovMatStar(sigma, gamma, x, locations)';
% 
% fPlusData = ComputeFplus(locations, surfNormals, meanValue, meanGrad);
% RVector = covMatData\fPlusData;
% fPlusGP = @(x)(CovMatStar(sigma, gamma, x, locations) * RVector);
% fPlus = @(x)([meanValue(x);meanGrad(x)] + fPlusGP(x));
% 
% 
% 
% uncertainties = [];
% values = [];
% for i = 1:length(vertices)
%     sig = uncertaintyGp(vertices(i,:)');
%     uncertainties = [uncertainties; sqrt(sig(1,1))];
%     values = [values ; fPlus(vertices(i,:)')'];
% end
% 
% figure
% hold on
% set(gca,'dataaspectratio',[1 1 1])
% patch('faces',faces,'vertices',vertices,...
%     'facecolor',[0.5 0.5 0.5], ...
%     'edgecolor', 'none', ...
%     'facelighting','gouraud');
% camlight('headlight')
% set(gca,'view',[46.8000   18.8000]);
% 
% 
% minUncer = min(uncertainties);
% maxUncer = max(uncertainties-minUncer);
% 
% color = (uncertainties-minUncer)/maxUncer;
% for i = 1:length(vertices)
%     plot3(vertices(i,1),vertices(i,2),vertices(i,3), 'Color',[color(i),0,0] ,'Marker','+','MarkerSize', 10)
% end
% 
% quiver3(vertices(:,1),vertices(:,2),vertices(:,3),...
%     values(:,2),values(:,3),values(:,4),'linewidth',1,'color','r');
% 

