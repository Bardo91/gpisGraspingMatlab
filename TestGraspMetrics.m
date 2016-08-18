%%
%   Contact point 
%       |--> pos: Position (vec3)
%       |--> normal: Normal (Vec3)
%       |--> nu: Friccion coefficient (double)
%       |--> frame: Coordinate Frame (rotation matrix)

%clc; close all; clear all;

addpath('grasping');

% ContactMatrix = [   0;...
%                     0;...
%                     1 ;...
%                     0;...
%                     0;...
%                     0];

ContactMatrix = [   1 0 0 0;...
                    0 1 0 0;...
                    0 0 1 0;...
                    0 0 0 0;...
                    0 0 0 0;...
                    0 0 0 1];
                    

% ContactMatrix = [   1 0 0 0 0 0;...
%                     0 1 0 0 0 0;...
%                     0 0 1 0 0 0;...
%                     0 0 0 0 0 0;...
%                     0 0 0 0 0 0;...
%                     0 0 0 0 0 1];

contactPoints = {};
%Generate random points on sphere.
% for i = 1:1000
%     x = random('uni', -1,1);
%     limY = sqrt(1-x*x);
%     y = random('uni', -limY,limY);
%     z = sign(random('uni', -1,1))*sqrt((1-x^2-y^2));
%     contactPoints{i}.pos = [x;y;z];
%     contactPoints{i}.normal = [x;y;z];
%     contactPoints{i}.nu = 0.1;
%     contactPoints{i}.contactFrame = contactFrame(contactPoints{i}.normal);
% end

contactPoints{1}.pos = [1,0.5,0]';
contactPoints{1}.normal = [-1,0,0]';
contactPoints{1}.nu = 0.1;
contactPoints{1}.contactFrame = contactFrame(contactPoints{1}.normal);
contactPoints{1}.selectionMatrix = ContactMatrix;

contactPoints{2}.pos = [-1,-0.5,0]';
contactPoints{2}.normal = [1,0,0]';
contactPoints{2}.nu = 0.1;
contactPoints{2}.contactFrame = contactFrame(contactPoints{2}.normal);
contactPoints{2}.selectionMatrix = ContactMatrix;

% contactPoints{3}.pos = [-0.5,-1,0]';
% contactPoints{3}.normal = [0,1,0]';
% contactPoints{3}.nu = 0.1;
% contactPoints{3}.contactFrame = contactFrame(contactPoints{3}.normal);
% contactPoints{3}.selectionMatrix = ContactMatrix;
% 
% 
% contactPoints{4}.pos = [0.5,1,0]';
% contactPoints{4}.normal = [0,-1,0]';
% contactPoints{4}.nu = 0.1;
% contactPoints{4}.contactFrame = contactFrame(contactPoints{4}.normal);
% contactPoints{4}.selectionMatrix = ContactMatrix;

figure();
N = length(contactPoints);
for i = 1:N
   showAxis(contactPoints{i}.pos, contactPoints{i}.contactFrame); 
end
grid;
axis equal;

objCenter = [0,0,0]';

plot3(objCenter(1), objCenter (2), objCenter(3), 'r.', 'MarkerSize', 10);

G = graspMatrix(contactPoints, objCenter);
% eigVals = eig(G*G');
% %Get only N positive eigen values
% nPosEigVals = eigVals(end-5:end);
% singularValues = sqrt(nPosEigVals);
singularValues = svd(G)'

Qmin = min(singularValues);
Qmax = max(singularValues);
Qelli = prod(singularValues);
Qiso = min(singularValues)/max(singularValues);
title(['Qmin:', num2str(Qmin), ', ','Qmax:', num2str(Qmax), ', ', 'Qelli:', num2str(Qelli), ', ', 'Qiso:', num2str(Qiso)]);
