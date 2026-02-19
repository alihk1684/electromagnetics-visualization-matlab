clear all; close all; clc;

% Dipole parameters
Q = 4e-9;
d = 0.01;
epsilon0 = 8.854e-12;

% Charge positions: dipole along z-axis
z_plus = d/2;    % +Q at (0,0,d/2)
z_minus = -d/2;  % -Q at (0,0,-d/2)

% Grid for XZ plane (y=0)
[x, z] = meshgrid(-0.01:0.0002:0.01, -0.01:0.0002:0.01);
y = zeros(size(x));

% Calculate distances
r_plus = sqrt(x.^2 + y.^2 + (z - z_plus).^2);
r_minus = sqrt(x.^2 + y.^2 + (z - z_minus).^2);

% Calculate potential
V = (Q/(4*pi*epsilon0)) * (1./r_plus - 1./r_minus);

% Calculate electric field
[Ex, Ez] = gradient(-V, 0.0002, 0.0002);

%% Part A: Plot in XZ plane only
figure(1)
set(gcf, 'Name', 'Part A: Potential Contours & Electric field', 'NumberTitle', 'off');

% XZ Plane - Potential contours
subplot(1,2,1)
contour(x, z, V, 100, 'LineWidth', 1);
hold on;
% charge markers
plot(0, z_plus, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
plot(0, z_minus, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
axis([-0.01 0.01 -0.01 0.01]);
xlabel('X (m)'); ylabel('Z (m)');
title('Potential Contours in XZ Plane');
grid on; hold off;

% XZ Plane - Electric field vectors
subplot(1,2,2)
% Create mask to exclude region very near charges
mask = sqrt(x.^2 + (z - z_plus).^2) > 0.001 & sqrt(x.^2 + (z - z_minus).^2) > 0.001;

% Apply mask to arrow data
x_masked = x;
z_masked = z;
Ex_masked = Ex;
Ez_masked = Ez;

x_masked(~mask) = NaN;
z_masked(~mask) = NaN;
Ex_masked(~mask) = 0;
Ez_masked(~mask) = 0;

% Plot with masked data (cleaner near charges)
quiver(x_masked(1:2:end, 1:2:end), z_masked(1:2:end, 1:2:end), ...
       Ex_masked(1:2:end, 1:2:end), Ez_masked(1:2:end, 1:2:end), 'b', 'LineWidth', 1);
hold on;

% charge markers
plot(0, z_plus, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
plot(0, z_minus, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
axis([-0.01 0.01 -0.01 0.01]);
xlabel('X (m)'); ylabel('Z (m)');
title('Electric Field Vectors in XZ Plane');
grid on; hold off;

%% Part B: Calculate |E| on X and Y axes correctly
figure(2)
set(gcf, 'Name', 'Part B: |E| on X & Y axes', 'NumberTitle', 'off');

% X-axis
x_points = linspace(-0.01, 0.01, 1000);
y_axis = zeros(size(x_points));
z_axis = zeros(size(x_points));

% Calculate field COMPONENTS along X-axis
% Distance vectors
r_plus_x = sqrt(x_points.^2 + y_axis.^2 + (z_axis - z_plus).^2);
r_minus_x = sqrt(x_points.^2 + y_axis.^2 + (z_axis - z_minus).^2);

% Field from +Q
Ex_plus = (Q/(4*pi*epsilon0)) * (x_points - 0) ./ (r_plus_x.^3);
Ey_plus = (Q/(4*pi*epsilon0)) * (y_axis - 0) ./ (r_plus_x.^3);
Ez_plus = (Q/(4*pi*epsilon0)) * (z_axis - z_plus) ./ (r_plus_x.^3);

% Field from -Q
Ex_minus = (-Q/(4*pi*epsilon0)) * (x_points - 0) ./ (r_minus_x.^3);
Ey_minus = (-Q/(4*pi*epsilon0)) * (y_axis - 0) ./ (r_minus_x.^3);
Ez_minus = (-Q/(4*pi*epsilon0)) * (z_axis - z_minus) ./ (r_minus_x.^3);

% Total field components
Ex_total = Ex_plus + Ex_minus;
Ey_total = Ey_plus + Ey_minus;
Ez_total = Ez_plus + Ez_minus;

% Magnitude along X-axis
E_mag_x = sqrt(Ex_total.^2 + Ey_total.^2 + Ez_total.^2);

% Y-axis
y_axis = linspace(-0.01, 0.01, 1000);
x_axis = zeros(size(y_axis));
z_axis = zeros(size(y_axis));

% Calculate field COMPONENTS along Y-axis
r_plus_y = sqrt(x_axis.^2 + y_axis.^2 + (z_axis - z_plus).^2);
r_minus_y = sqrt(x_axis.^2 + y_axis.^2 + (z_axis - z_minus).^2);

% Field from +Q
Ex_plus = (Q/(4*pi*epsilon0)) * (x_axis - 0) ./ (r_plus_y.^3);
Ey_plus = (Q/(4*pi*epsilon0)) * (y_axis - 0) ./ (r_plus_y.^3);
Ez_plus = (Q/(4*pi*epsilon0)) * (z_axis - z_plus) ./ (r_plus_y.^3);

% Field from -Q
Ex_minus = (-Q/(4*pi*epsilon0)) * (x_axis - 0) ./ (r_minus_y.^3);
Ey_minus = (-Q/(4*pi*epsilon0)) * (y_axis - 0) ./ (r_minus_y.^3);
Ez_minus = (-Q/(4*pi*epsilon0)) * (z_axis - z_minus) ./ (r_minus_y.^3);

% Total field components
Ex_total = Ex_plus + Ex_minus;
Ey_total = Ey_plus + Ey_minus;
Ez_total = Ez_plus + Ez_minus;

% Magnitude along Y-axis
E_mag_y = sqrt(Ex_total.^2 + Ey_total.^2 + Ez_total.^2);

% Plot |E| on X-axis
subplot(2,1,1)
plot(x_points, E_mag_x, 'b-', 'LineWidth', 2);
xlabel('X position (m)'); ylabel('|E| (V/m)');
title('Electric Field Magnitude along X-axis (y=0, z=0)');
grid on;

% Plot |E| on Y-axis
subplot(2,1,2)
plot(y_axis, E_mag_y, 'r-', 'LineWidth', 2);
xlabel('Y position (m)'); ylabel('|E| (V/m)');
title('Electric Field Magnitude along Y-axis (x=0, z=0)');
grid on;