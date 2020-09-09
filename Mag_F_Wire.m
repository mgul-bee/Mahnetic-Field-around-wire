%% Initialize
clc
clear all 
close all 
%% parameters
L=-10:1:10;
dl = [zeros(length(L'),1),zeros(length(L'),1),L']; % x,y,z 

I = 1; % filament current [A]
dG = 1e9; % filament max discretization step
n = length(dl(:,1)); % number of filaments on wire, 21

% Field points (where we want to calculate the field)
x_M = linspace(-10,10,21); % x [m]
y_M = linspace(-10,10,21); % y [m]
z_M = linspace(-10,10,21); % z [m]
[X,Y,Z]=meshgrid(x_M,y_M,z_M);

HX = zeros(size(X,1),size(X,2),size(X,3));
HY = zeros(size(X,1),size(X,2),size(X,3));
HZ = zeros(size(X,1),size(X,2),size(X,3));
%% Loop on each filament
for nF = 1:n     
    % Discretization of dl
    x_P = []; y_P = []; z_P = [];
    N = size(dl,1)-1; % Number of points defining dl
    for i = 1:N % Loop on the segments defining dl
        L_dl_i = norm(dl(i,:)-dl(i+1,:));
        NP = ceil(L_dl_i/dG); % Number of points required to have a discretization step smaller than dG
        x_P = [x_P,linspace(dl(i,1), dl(i+1,1), NP)]; % discretization of dl for x component
        y_P = [y_P,linspace(dl(i,2), dl(i+1,2), NP)]; % discretization of dl for y component
        z_P = [z_P,linspace(dl(i,3), dl(i+1,3), NP)]; % discretization of dl for z component
    end
    % Add contribution of each source point P on each field point M (where we want to calculate the field)
    for m = 1:size(X,1)
        for n = 1:size(X,2)
            for p = 1:size(X,3)
            % M is the field point
            x_M = X(m,n,p);
            y_M = Y(m,n,p);
            z_M = Z(m,n,p);
            % Loop on each discretized segment of dl PkPk+1
            for k = 1:length(x_P)-1
                PkM3 = (sqrt((x_M-x_P(k))^2 + (y_M-y_P(k))^2 + (z_M-z_P(k))^2))^3;
                DHx(k) = ((y_P(k+1)-y_P(k))*(z_M-z_P(k))-(z_P(k+1)-z_P(k))*(y_M-y_P(k)))/PkM3;
                DHy(k) = ((z_P(k+1)-z_P(k))*(x_M-x_P(k))-(x_P(k+1)-x_P(k))*(z_M-z_P(k)))/PkM3;
                DHz(k) = ((x_P(k+1)-x_P(k))*(y_M-y_P(k))-(y_P(k+1)-y_P(k))*(x_M-x_P(k)))/PkM3;
            end
            % Sum
            HX(m,n,p) = HX(m,n,p) + I/4/pi*sum(DHx);
            HY(m,n,p) = HY(m,n,p) + I/4/pi*sum(DHy);
            HZ(m,n,p) = HZ(m,n,p) + I/4/pi*sum(DHz);
            end
        end
    end
end

%% Plots
figure 
plot3(dl(:,1),dl(:,2),dl(:,3),'.-r', 'LineWidth',2) % plot the wire along z-axis in red
axis tight
axis square
hold on
title('3D Plot of magnetic field of current in wires')
xlabel('X')
ylabel('Y')
zlabel('Z')
normH=sqrt(HX.^2+HY.^2+HZ.^2);
quiver3(X,Y,Z,HX./normH,HY./normH,HZ./normH,1)