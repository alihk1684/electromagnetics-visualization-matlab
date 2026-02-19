clear all; close all; clc;

% Parameters
mu0 = 4*pi*1e-7;
I = 4;
R = 0.01;

% Biot-Savart Function
function [Bx, By, Bz] = biot_savart_loop(X, Y, Z, R, I, mu0)
    N = 100;
    theta = linspace(0, 2*pi, N+1);
    theta = theta(1:end-1);
    
    dl_points = [R*cos(theta); R*sin(theta); zeros(1, N)]';
    dl_vectors = zeros(N, 3);
    dl_vectors(:,1) = -R * sin(theta)' * (2*pi/N);  % dl_x
    dl_vectors(:,2) =  R * cos(theta)' * (2*pi/N);  % dl_y
    
    % Biot-Savart calculation
    Bx = zeros(size(X));
    By = zeros(size(Y));
    Bz = zeros(size(Z));
    
    for i = 1:numel(X)
        r_obs = [X(i), Y(i), Z(i)];   % obsercation point
        B_total = [0, 0, 0];
    
        for s = 1:N
            r_prime = dl_points(s, :);  % source point
            dl = dl_vectors(s, :);
    
            R_vec = r_obs - r_prime;
            R_mag = norm(R_vec);
    
            if R_mag < 1e-10
                continue;
            end
    
            dB = (mu0/(4*pi)) * I * cross(dl, R_vec) / (R_mag^3);
            B_total = B_total + dB;
        end
    
        Bx(i) = B_total(1);
        By(i) = B_total(2);
        Bz(i) = B_total(3);
    end
end

%% 3D plot of loop and B Lines
figure;
set(gcf, 'Name', '3D plot of loop and B Lines', 'NumberTitle', 'off');
hold on; grid on; axis([-0.0175,0.0175,-0.0175,0.0175,-0.0175,0.0175])
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Magnetic Field Lines of a Circular Current Loop (Biot–Savart)');

% Create Cartesian grid
x_range = linspace(-0.0175, 0.0175, 31);
y_range = linspace(-0.0175, 0.0175, 31);
z_range = linspace(-0.0175, 0.0175, 31);

[X_cart, Y_cart, Z_cart] = meshgrid(x_range, y_range, z_range);

[Bx_cart, By_cart, Bz_cart] = biot_savart_loop(X_cart, Y_cart, Z_cart, R, I, mu0);

% Build unit direction field (so magnitude doesn't mess up the look)
Bmag_cart = sqrt(Bx_cart.^2 + By_cart.^2 + Bz_cart.^2);
Bmag_cart(Bmag_cart == 0) = 1;
Bxn_cart = Bx_cart ./ Bmag_cart;
Byn_cart = By_cart ./ Bmag_cart;
Bzn_cart = Bz_cart ./ Bmag_cart;

% Plot the current loop
theta = linspace(0, 2*pi, 100);
plot3(R*cos(theta), R*sin(theta), zeros(size(theta)), 'r', 'LineWidth', 3);

% Tangent direction:
theta0 = 17*pi/10;
p0 = [R*cos(theta0), R*sin(theta0), 0];
tCW = [ sin(theta0), -cos(theta0), 0];
tCW = tCW / norm(tCW);
arrowLenLoop = 0.0025;
quiver3(p0(1), p0(2), p0(3), ...
        arrowLenLoop*tCW(1), arrowLenLoop*tCW(2), arrowLenLoop*tCW(3), ...
        0, 'g', 'LineWidth', 3, 'MaxHeadSize', 25, 'AutoScale', 'off');

% Plot the field lines
seedL = 0.01;
nSeed = 3;
u = linspace(-seedL, seedL, nSeed);

% x = 0 plane
[YY, ZZ] = meshgrid(u, u);
sx0 = zeros(numel(YY),1);
sy0 = YY(:);
sz0 = ZZ(:);

% y = 0
[XX, ZZ] = meshgrid(u, u);
sx1 = XX(:);
sy1 = zeros(numel(XX),1);
sz1 = ZZ(:);

% x = y
[SS, ZZ] = meshgrid(u, u);
sx2 = SS(:);
sy2 = SS(:);
sz2 = ZZ(:);

% x = -y
sx3 = SS(:);
sy3 = -SS(:);
sz3 = ZZ(:);

% Combine all seeds
sx = [sx0; sx1; sx2; sx3];
sy = [sy0; sy1; sy2; sy3];
sz = [sz0; sz1; sz2; sz3];

% Plot streamlines (field lines)
h = streamline(X_cart, Y_cart, Z_cart, Bxn_cart, Byn_cart, Bzn_cart, sx, sy, sz);
set(h, 'Color', 'b', 'LineWidth', 1.0);

view(40,25);

%% 3D plot of loop and B vectors
figure;
set(gcf, 'Name', '3D plot of loop and B vectors', 'NumberTitle', 'off');
hold on; grid on; axis([-0.0175,0.0175,-0.0175,0.0175,-0.0175,0.0175]);
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Magnetic Field Vectors (quiver3)');

% Create cylindrical grid
rho_range = linspace(0.00001, 0.0175, 25);
phi_range = linspace(0, 2*pi - 2*pi/41, 40);
z_range = linspace(-0.0175, 0.0175, 25);

[Rho, Phi, Z_cyl] = meshgrid(rho_range, phi_range, z_range);
X_cyl = Rho .* cos(Phi);
Y_cyl = Rho .* sin(Phi);

% Calculate B-field on cylindrical grid
[Bx_cyl, By_cyl, Bz_cyl] = biot_savart_loop(X_cyl, Y_cyl, Z_cyl, R, I, mu0);

% Plot the current loop
plot3(R*cos(theta), R*sin(theta), zeros(size(theta)), 'r', 'LineWidth', 3);

% Tangent direction:
theta0 = 17*pi/10;
p0 = [R*cos(theta0), R*sin(theta0), 0];
tCW = [ sin(theta0), -cos(theta0), 0];
tCW = tCW / norm(tCW);
arrowLenLoop = 0.0025;
quiver3(p0(1), p0(2), p0(3), ...
        arrowLenLoop*tCW(1), arrowLenLoop*tCW(2), arrowLenLoop*tCW(3), ...
        0, 'g', 'LineWidth', 3, 'MaxHeadSize', 23, 'AutoScale', 'off');

% Downsample the grid for plotting arrows
step = 4;
Xs = X_cyl(1:step:end, 1:step:end, 1:step:end);
Ys = Y_cyl(1:step:end, 1:step:end, 1:step:end);
Zs = Z_cyl(1:step:end, 1:step:end, 1:step:end);

Bxs = Bx_cyl(1:step:end, 1:step:end, 1:step:end);
Bys = By_cyl(1:step:end, 1:step:end, 1:step:end);
Bzs = Bz_cyl(1:step:end, 1:step:end, 1:step:end);

% Mask points too close to the wire
rho_s = sqrt(Xs.^2 + Ys.^2);
delta = 1.0e-3;
mask = ~(abs(rho_s - R) < delta & abs(Zs) < delta);

Xs = Xs(mask); Ys = Ys(mask); Zs = Zs(mask);
Bxs = Bxs(mask); Bys = Bys(mask); Bzs = Bzs(mask);

% Normalize vectors
Bms = sqrt(Bxs.^2 + Bys.^2 + Bzs.^2);
Bms(Bms==0) = 1;
Bxs_n = Bxs ./ Bms;
Bys_n = Bys ./ Bms;
Bzs_n = Bzs ./ Bms;

% Log-scaled magnitude for coloring
logB = log10(Bms + eps);
% Normalize to [0,1] for colormap indexing
logB_norm = (logB - min(logB(:))) ./ (max(logB(:)) - min(logB(:)));

cmap = jet(256);
color_idx = round(logB_norm * (size(cmap,1)-1)) + 1;
colors = cmap(color_idx, :);

% Plot arrows
arrowLen = 0.001;
for k = 1:numel(Xs)
    quiver3(Xs(k), Ys(k), Zs(k), ...
        arrowLen*Bxs_n(k), arrowLen*Bys_n(k), arrowLen*Bzs_n(k), ...
        0, 'Color', colors(k,:), ...
        'LineWidth', 0.9, 'MaxHeadSize', 3);
end

% Colorbar
colormap(jet);
cb = colorbar;
cb.Label.String = 'log_{10}(|B|)  [T]';
cb.Label.FontSize = 11;

view(40,25);

%% ANALYTICAL COMPARISON
% Analytical solution
z_axis = linspace(-0.05, 0.05, 200);
Bz_analytic = (mu0 * I * R^2) ./ (2 * (R^2 + z_axis.^2).^(3/2));

% Numerical from our calculation
z_test = [-0.0175, -0.015, -0.0125, -0.01, -0.0075, -0.005, -0.0025, 0, 0.0025, 0.005, 0.0075, 0.01, 0.012, 0.015, 0.0175];
Bz_numeric = zeros(size(z_test));

for i = 1:length(z_test)
    % Find closest point in cylindrical grid
    [~, idx_z] = min(abs(z_range - z_test(i)));
    [~, idx_rho] = min(abs(rho_range - 0));
    
    Bz_numeric(i) = Bz_cyl(1, idx_rho, idx_z);
end

figure;
set(gcf, 'Name', 'Magnetic Field on Z-axis', 'NumberTitle', 'off');
plot(z_axis*100, Bz_analytic*1e6, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(z_test*100, Bz_numeric*1e6, 'ro', 'MarkerSize', 4, 'LineWidth', 2, 'DisplayName', 'Numerical');
grid on;
xlabel('z (cm)');
ylabel('B_z (\muT)');
title('Magnetic Field on Z-axis: Analytical vs Numerical');
legend('Location', 'best');
