% GPIS learning
% 2D example, Sphere mean and normals

close all; clear all; clc;

% X = [0,-0.5;
%     -0.3,-0.1;
%     0.5,-1.0;
%      -0.5,0.5;
%      0.5,0.5;
%     1,0]';
% 
% f = [   0,-cos(20/180*pi),-sin(20/180*pi),...
%         0,-cos(45/180*pi),-sin(45/180*pi),...
%         0,-cos(120/180*pi),-sin(120/180*pi),...
%         0, -cos(45/180*pi), sin(45/180*pi),...
%         0, cos(45/180*pi), sin(45/180*pi),...
%         0,1,0]';

X = [0,-1]';

f = [   0,0,-1]';
data = f ;
m  = 1; 

lim = 1.2;
step = 1.2/15;

[Xg,Yg] = meshgrid(-lim:step:lim,-lim:step:lim);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
n = length(Xs);

sigma = 1;
gamma = 1;

display('Computing covariance matrix K');
K = ComputeFullKder(sigma, gamma, X, 0.0, 0.0);

display('Computing covariance matrix Ks');
Ks = ComputeKderX1X2(sigma, gamma, Xs, X);

display('Computing covariance matrix Kss');
Kss = ComputeFullKder(sigma, gamma, Xs, 0.0, 0.0);

% figure();
% imagesc(K);
% colorbar;
% 
% figure();
% imagesc(Ks);
% colorbar;
% 
% figure();
% imagesc(Kss);
% colorbar;


display('Computing means');

R = 1;
cen = [0.0, 0.0]';
mean = @(x) 1/2/R*((x-cen)'*(x-cen) - R^2);
meandx = @(x) 1/R*((x(1)-cen(1)));
meandy = @(x) 1/R*((x(2)-cen(2)));

for i = 1:m
    mu((i-1)*3 +1) = mean(X(:,i));
    mu((i-1)*3 +2) = meandx(X(:,i));
    mu((i-1)*3 +3) = meandy(X(:,i));
end

for i = 1:n
    mus((i-1)*3 +1) = mean(Xs(:,i));
    mus((i-1)*3 +2) = meandx(Xs(:,i));
    mus((i-1)*3 +3) = meandy(Xs(:,i));
end

mu = mu';
mus = mus';


display('Computing Regression');
kinv = inv(K);
fs = mus + Ks*kinv*(f - mu);
sig = Kss - Ks*kinv*Ks';

% figure();
% imagesc(sig);
% colorbar;

d = diag(sig);
devFsP = fs + d;
devFsP = reshape(devFsP(1:3:end),d1,d2);

devFsN = fs - d;
devFsN = reshape(devFsN(1:3:end),d1,d2);

Fs = reshape(fs(1:3:end),d1,d2);

% figure();
% hold on;
% plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
% contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
% quiver(X(1,:), X(2,:), data(2,:), data(3,:));
% axis equal;

% figure();
% hold on;
% plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
% surf(Xg,Yg,Fs);
% contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
% axis equal;

% figure();
% hold on;
% plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
% surf(Xg,Yg,Fs);
% contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
% surf(Xg,Yg,devFsN);
% surf(Xg,Yg,devFsP);
% axis equal;

% Compute the inside probability.
cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

D = reshape(d(1:3:end), d1,d1);
prob = cdfBell((0-Fs)./D);

% figure();
% hold on;
% surf(Xg,Yg,prob);
% axis equal;

figure();
hold on;
contour(Xg,Yg,prob);
contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
axis equal;

% figure;
% hold on;
% contourf(Xg,Yg,D);     
% contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
% colorbar;
% axis equal;

%% Plot cost function

% Create mask with circle to compute cost
a = 90;
b = 180;
efficiency = @(x) exp(-(x-a)./(b-x));
angles=a:1:b;
nu = efficiency(angles);

R = @(x) 0.5*exp((x-a)./(b-x));

EA = zeros(length(nu),1)';
index = 1;
for angle=angles
    mask = zeros(n,1);
    r = R(angle);
    c = [0, -1 + r];
    for i=1:n
        x = Xs(1,i);
        y = Xs(2,i);
        val = (x-c(1))^2 + (y-c(2))^2;
        if(val < r^2)
            mask(i) = 1;
        end
    end
    mask = reshape(mask, d1,d1);
    normFactor = sum(sum(mask));
    if normFactor ~= 0
        EA(index) = sum(sum(prob.*mask))/normFactor;
    else
        EA(index) = sum(sum(prob.*mask));
    end
    index = index +1;
end

% mask = reshape(mask, d1,d1);
% figure();
% hold on;
% contour(Xg,Yg,mask);
% colorbar;
% axis equal;

figure();
hold on;
plot(angles, nu);
plot(angles, EA);
plot(angles, nu.*EA);
legend('nu', 'EA', 'U');
