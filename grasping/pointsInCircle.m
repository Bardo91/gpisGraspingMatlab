function [ points ] = pointsInCircle( radius,center, azimutAngle, zenitAngle, nPoints )
%POINTSINCIRCLE Summary of this function goes here
%   Detailed explanation goes here

n = [cos(azimutAngle)*sin(zenitAngle);sin(zenitAngle)*sin(azimutAngle); cos(azimutAngle)];
u = [-sin(azimutAngle); cos(azimutAngle); 0];

nxu = [cos(zenitAngle)*cos(azimutAngle); cos(zenitAngle)*sin(azimutAngle); -sin(zenitAngle)];

pCircle = @(t) radius*cos(t)*u + radius*sin(t)*nxu + center;

points = [];
for angle=0:(2*pi/nPoints):2*pi
    points = [ points , pCircle(angle)];
end
points = points(:,1:end-1);
end

