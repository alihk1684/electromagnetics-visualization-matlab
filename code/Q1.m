clear all; close all; clc;

% Parameters
R = 3;
rho_s = 3;
epsilon_0 = 8.854e-12;

% 3D grid
[x, y, z] = meshgrid(-4:0.5:4, -4:0.5:4, -4:0.5:4);
r = sqrt(x.^2 + y.^2 + z.^2);

%% Part A
Ex = zeros(size(x));
Ey = zeros(size(y));
Ez = zeros(size(z));

for i = 1:size(x,1)
    for j = 1:size(x,2)
        for k = 1:size(x,3)
            if r(i,j,k) <= R
                magnitude = (rho_s * r(i,j,k)) / (3*epsilon_0);
            else
                magnitude = (rho_s * R^3) / (3*epsilon_0 * r(i,j,k)^2);
            end

            if r(i,j,k) ~= 0
                Ex(i,j,k) = magnitude * x(i,j,k) / r(i,j,k);
                Ey(i,j,k) = magnitude * y(i,j,k) / r(i,j,k);
                Ez(i,j,k) = magnitude * z(i,j,k) / r(i,j,k);
            end
        end
    end
end

% plot sphere and electric field vectors inside only
figure(1)
set(gcf, 'Name', 'Part A: Sphere electric field', 'NumberTitle', 'off');

% transparent sphere
sphere_points = 50;
[Xs, Ys, Zs] = sphere(sphere_points);
surf(R*Xs, R*Ys, R*Zs, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
axis equal; xlabel('X'); ylabel('Y'); zlabel('Z');
title('solid sphere (Radius = 3m) with Electric Field');
grid on; hold on;

% Electric field vectors for indside points
inside = r <= R-0.2;

quiver3(x(inside), y(inside), z(inside), ...
    Ex(inside), Ey(inside), Ez(inside), 'r');
hold off;

%% part B
figure(2)
set(gcf, 'Name', 'Part B: Plots', 'NumberTitle', 'off');

r_vals = linspace(0, 10, 1000);
Q_enc = zeros(size(r_vals));
for i = 1:length(r_vals)
    if r_vals(i) <= R
        Q_enc(i) = (4/3) * pi * r_vals(i)^3 * rho_s;
    else
        Q_enc(i) = (4/3) * pi * R^3 * rho_s;
    end
end

E_mag = zeros(size(r_vals));
for i = 1:length(r_vals)
    if r_vals(i) <= R
        E_mag(i) = (rho_s * r_vals(i)) / (3*epsilon_0);
    else
        E_mag(i) = (rho_s * R^3) / (3*epsilon_0 * r_vals(i)^2);
    end
end

flux = Q_enc / epsilon_0;

subplot(3,1,1)
plot(r_vals, Q_enc, 'b-', 'LineWidth', 2);
hold on;
plot([R, R], [0, max(Q_enc)], 'r--', 'LineWidth', 1.5);
xlabel('Radius r (m)'); ylabel('Enclosed Charge Q_{enc} (C)');
title('Enclosed Charge vs Radius');
legend('Q_{enc}', 'Sphere Boundary', 'Location', 'northwest');
grid on;

subplot(3,1,2)
plot(r_vals, E_mag, 'g-', 'LineWidth', 2);
hold on;
plot([R, R], [0, max(E_mag)], 'r--', 'LineWidth', 1.5);
xlabel('Radius r (m)'); ylabel('Electric Field |E| (V/m)');
title('Electric Field Magnitude vs Radius');
legend('|E|', 'Sphere boundary', 'Location', 'northwest');
grid on;

subplot(3,1,3)
plot(r_vals, flux, 'm-', 'LineWidth', 2);
hold on;
plot([R, R], [0, max(flux)], 'r--', 'LineWidth', 1.5);
xlabel('Radius r (m)'); ylabel('Flux Φ_E (V·m)');
title('Electric Flux vs Radius');
legend('Φ_E', 'Sphere Boundary', 'Location', 'northwest');
grid on;

%% Part C
figure(3)
set(gcf, 'Name', 'Part C: Sphere shell', 'NumberTitle', 'off');
Q_total = (4/3) * pi * R^3 * rho_s;
sigma = Q_total / (4 * pi * R^2);

r_vals = linspace(0, 10, 1000);

E_mag_shell = zeros(size(r_vals));
Q_enc_shell = zeros(size(r_vals));

for i = 1:length(r_vals)
    if r_vals(i) < R
        E_mag_shell(i) = 0;
        Q_enc_shell(i) = 0;
    elseif r_vals(i) == R
        E_mag_shell(i) = sigma / (2*epsilon_0);
        Q_enc_shell(i) = Q_total;
    else
        E_mag_shell(i) = Q_total / (4*pi*epsilon_0 * r_vals(i)^2);
        Q_enc_shell(i) = Q_total;
    end
end

flux_shell = Q_enc_shell / epsilon_0;

subplot(3,1,1)
plot(r_vals, Q_enc, 'b-', 'LineWidth', 2); hold on;
plot(r_vals, Q_enc_shell, 'r--', 'LineWidth', 2);
plot([R, R], [0, max(Q_enc)], 'k:', 'LineWidth', 1);
xlabel('Radius r (m)'); ylabel('Enclosed Charge Q_{enc} (C)');
title('Q_{enc}: Solid Sphere vs Spherical Shell');
legend('Solid Sphere', 'Spherical Shell', 'Boundary', 'Location', 'northwest');
grid on;

subplot(3,1,2)
plot(r_vals, E_mag, 'b-', 'LineWidth', 2); hold on;
plot(r_vals, E_mag_shell, 'r--', 'LineWidth', 2);
plot([R, R], [0, max(E_mag)], 'k:', 'LineWidth', 1);
xlabel('Radius r (m)'); ylabel('Electric Field |E| (V/m)');
title('|E|: Solid Sphere vs Spherical Shell');
legend('Solid Sphere', 'Spherical Shell', 'Boundary', 'Location', 'northeast');
grid on;

subplot(3,1,3)
plot(r_vals, flux, 'b-', 'LineWidth', 2); hold on;
plot(r_vals, flux_shell, 'r--', 'LineWidth', 2);
plot([R, R], [0, max(flux)], 'k:', 'LineWidth', 1);
xlabel('Radius r (m)'); ylabel('Flux Φ_E (V·m)');
title('Flux: Solid Sphere vs Spherical Shell');
legend('Solid Sphere', 'Spherical Shell', 'Boundary', 'Location', 'northwest');
grid on;