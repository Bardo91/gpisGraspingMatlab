%%
%   Contact point 
%       |--> pos: Position (vec3)
%       |--> normal: Normal (Vec3)
%       |--> nu: Friccion coefficient (double)
%       |--> frame: Coordinate Frame (rotation matrix)


addpath('grasping');

contactPoints = {};
for i = 1:6
    contactPoints{i}.pos = random('uni',-1,1, [3,1]);
    contactPoints{i}.normal = random('uni',-1,1, [3,1]);
    contactPoints{i}.normal = contactPoints{i}.normal/norm(contactPoints{i}.normal);
    contactPoints{i}.nu = 0.1;
    contactPoints{i}.contactFrame = contactFrame(contactPoints{i}.normal);
end

N = length(contactPoints);
for i = 1:N
   showAxis(contactPoints{i}.pos, contactPoints{i}.contactFrame); 
end
grid;
axis equal;

G = graspMatrix(contactPoints);
% eigVals = eig(G*G');
% %Get only N positive eigen values
% nPosEigVals = eigVals(end-5:end);
% singularValues = sqrt(nPosEigVals);
singularValues = svd(G);

Q = min(singularValues)/max(singularValues)
