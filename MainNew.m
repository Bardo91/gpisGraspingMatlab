%% DATA
close all
addpath('covMat')
addpath('grasping')

noiseVals = 0.00000;
noiseGrad = 0.037;
sigma = 0.1243;%0.2844; %0.1
L0 = .2823; 
gamma = 2.5977;%0.6525;%1/L0^2;

% % load '/Users/Yannick/Coding/SurfaceConstruction/SurfaceConstruction/bmw_total'
% load  '/Users/Yannick/Coding/SurfaceConstruction/SurfaceConstruction/bmw_11'
load 'crawlerHandlerNoFloor'
% load  'bmw_11'
locations = PartMeans;
surfNormals = SurfNormals;

% % For bmw_11
% cx = -6.1655;
% cy = -0.0472;
% cz = -3.6693;
% 
% a = 1.0763;
% b = 2.3691;
% c = 1.3254;
% r = [2 * pi , 0, 0];
Prior.type = 'N';
Prior.param(1) = 0.5;
Prior.Gamma = 7;
Prior.Sigma = 0.5;
Prior.noiseGrad = 0.01;
Prior.noiseVals = 0.001;
% 
% Prior = struct('pos',[cx cy cz],'type',prior_type, 'param', [a b c], 'rot', r);
%% SURFACE
[meanValue, meanGrad] = computePriorFunctions(Prior)

dist = 0.01;
initPoints = locations;%r * [-4;0;-2];


[faces, vertices] = computeSurface(locations, surfNormals, ...
    Prior, ...
    meanValue, meanGrad, initPoints(:,1:5), dist, false);

figure
hold on
set(gca,'dataaspectratio',[1 1 1])

plot3(locations(1,:),locations(2,:),locations(3,:),'r.','markersize',5);
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',1,'color','r');
quiver3(locations(1,:),locations(2,:),locations(3,:),...
    surfNormals(1,:),surfNormals(2,:),surfNormals(3,:),'linewidth',1,'color','r');

patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','gouraud');
camlight('headlight')
set(gca,'view',[46.8000   18.8000]);
% light('Position',[-1 -1 0])

%% Grasping
sigma = Prior.Sigma;
gamma = Prior.Gamma;
noiseVals = Prior.noiseVals;
noiseGrad = Prior.noiseGrad;

covMatData = ComputeCovMatFull(sigma,gamma,locations,noiseVals,noiseGrad);
invCovMatData = inv(covMatData);
covPoint = @(x) ComputeCovMatFull(sigma, gamma, x, noiseVals, noiseGrad);
uncertaintyGp = @(x) covPoint(x) - CovMatStar(sigma, gamma, x, locations)*invCovMatData*CovMatStar(sigma, gamma, x, locations)';

fPlusData = ComputeFplus(locations, surfNormals, meanValue, meanGrad);
RVector = covMatData\fPlusData;
fPlusGP = @(x)(CovMatStar(sigma, gamma, x, locations) * RVector);
fPlus = @(x)([meanValue(x);meanGrad(x)] + fPlusGP(x));



uncertainties = [];
values = [];
for i = 1:length(vertices)
    sig = uncertaintyGp(vertices(i,:)');
    uncertainties = [uncertainties; sqrt(sig(1,1))];
    values = [values ; fPlus(vertices(i,:)')'];
end

figure
hold on
set(gca,'dataaspectratio',[1 1 1])
patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','gouraud');
camlight('headlight')
set(gca,'view',[46.8000   18.8000]);


minUncer = min(uncertainties);
maxUncer = max(uncertainties-minUncer);

color = (uncertainties-minUncer)/maxUncer;
for i = 1:length(vertices)
    plot3(vertices(i,1),vertices(i,2),vertices(i,3), 'Color',[color(i),0,0] ,'Marker','+','MarkerSize', 10)
end

quiver3(vertices(:,1),vertices(:,2),vertices(:,3),...
    values(:,2),values(:,3),values(:,4),'linewidth',1,'color','r');


