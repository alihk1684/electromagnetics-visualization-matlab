clear all; close all; clc;

% parameters
teta0_deg = 15;
teta1_deg = 85;
V0 = -2;
V1 = 10;
% Radians
teta0 = deg2rad(teta0_deg);
teta1 = deg2rad(teta1_deg);

% Space mesh: teta direction & Time mesh
x = linspace(teta0, teta1, 100);
t = linspace(0, 1, 100);

% PDE
m = 0;
sol = pdepe(m, @pdefun, @icfun, @bcfun, x, t);
V = sol;

%% Part A: Time-dependent solution
% 3D Surface plot of V
[T, X] = meshgrid(t, rad2deg(x));
figure(1);
set(gcf, 'Name', 'Part A: 3D Surface Plot', 'NumberTitle', 'off');
surf(T, X, V', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
colormap(jet);
colorbar;
xlabel('Time (s)');
ylabel('\theta (degrees)');
zlabel('V(\theta, t) (V)');
title('3D Plot of Potential');
grid on;
view(45, 30);

figure(2);
set(gcf, 'Name', 'Part A: Contour Plot', 'NumberTitle', 'off');
contourf(T, X, V', 30, 'LineColor', 'k', 'LineWidth', 0.01);
colorbar;
xlabel('Time (s)');
ylabel('\theta (degrees)');
title('Contour Plot of V(\theta, t)');

%% Part B: Steady-state comparison
V_numerical_steady = V(end, :);
theta_deg = rad2deg(x);

syms A_sym B_sym

eq1 = A_sym * log(tan(teta0/2)) + B_sym == V0;  % V(teta0) = V0
eq2 = A_sym * log(tan(teta1/2)) + B_sym == V1;  % V(teta1) = V1

sol_AB = solve([eq1, eq2], [A_sym, B_sym]);
A_val = double(sol_AB.A_sym);
B_val = double(sol_AB.B_sym);

% Analytical function
V_analytical = @(theta) A_val * log(tan(theta/2)) + B_val;

% Calculate analytical values at our grid points
V_analytical_steady = arrayfun(V_analytical, x);

% Plot comparison
figure(3);
set(gcf, 'Name', 'Part B: Steady-State', 'NumberTitle', 'off');
plot(theta_deg, V_numerical_steady, 'b-', 'LineWidth', 2, 'DisplayName', 'Numerical (t = 1 s)');
hold on;
plot(theta_deg, V_analytical_steady, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold off;

xlabel('\theta (degrees)');
ylabel('V(\theta) (V)');
title('Part B: Steady-State Potential Comparison');
legend('Location', 'best');
grid on;

% Text box
error = max(abs(V_numerical_steady - V_analytical_steady));
Text = sprintf('Analytical coefficients:\nA = %.4f\nB = %.4f\nMax error: %.4f V', A_val, B_val, error);
annotation('textbox', [0.15, 0.7, 0.2, 0.15], 'String', Text, ...
           'BackgroundColor', 'none', 'FontSize', 10, 'EdgeColor', 'none');

% Function Definitions

% PDE function
function [c, f, s] = pdefun(x, ~, ~, dudx)
    c = sin(x);
    f = sin(x) * dudx;
    s = 0;
end

% initial condition function
function u0 = icfun(~)
    u0 = 0;
end

% Boundary condition function
function [pl, ql, pr, qr] = bcfun(~, ul, ~, ur, ~)
    %Left
    pl = ul - (-2);
    ql = 0;
    %Right
    pr = ur - 10;
    qr = 0;
end

%% Part C: Ploting electric field between the shells
% Use the A value we already calculated from Part B
r = linspace(0.5, 2, 8);
theta = linspace(teta0, teta1, 12);
theta = theta(2:end-1);
phi = linspace(0, 2*pi, 16);

X = []; Y = []; Z = [];
Ex = []; Ey = []; Ez = [];
E_mag = [];  % Store actual magnitudes

for i = 1:length(r)
    for j = 1:length(theta)
        for k = 1:length(phi)
            r_curr = r(i);
            theta_curr = theta(j);
            phi_curr = phi(k);
            
            % Position
            x = r_curr * sin(theta_curr) * cos(phi_curr);
            y = r_curr * sin(theta_curr) * sin(phi_curr);
            z = r_curr * cos(theta_curr);
            
            % Electric field
            E_theta = -A_val / (r_curr * sin(theta_curr));
            
            % Convert to Cartesian
            Ex_curr = E_theta * cos(theta_curr) * cos(phi_curr);
            Ey_curr = E_theta * cos(theta_curr) * sin(phi_curr);
            Ez_curr = -E_theta * sin(theta_curr);
            
            % Magnitude
            E_mag_curr = abs(E_theta);
            
            % Store
            X = [X; x]; Y = [Y; y]; Z = [Z; z];
            Ex = [Ex; Ex_curr]; Ey = [Ey; Ey_curr]; Ez = [Ez; Ez_curr];
            E_mag = [E_mag; E_mag_curr];
        end
    end
end

figure(4);
set(gcf, 'Name', 'Part C: Electric Field - Normalized Arrows with Color-Coded Magnitude', ...
         'NumberTitle', 'off');
hold on;

% normalize arrows
vector_magnitudes = sqrt(Ex.^2 + Ey.^2 + Ez.^2);
vector_magnitudes(vector_magnitudes == 0) = 1;  % Avoid division by zero

Ex_norm = Ex ./ vector_magnitudes;
Ey_norm = Ey ./ vector_magnitudes;
Ez_norm = Ez ./ vector_magnitudes;

scale_factor = 0.2;  % Fixed arrow length

% COLOR by log-scaled magnitude
% Calculate log-scaled magnitude for coloring
log_E_mag = log10(E_mag + 1);  % +1 to avoid log(0)
min_log = min(log_E_mag);
max_log = max(log_E_mag);

% Normalize to [0, 1] for colormap
color_values = (log_E_mag - min_log) / (max_log - min_log);

% Colormap
cmap = jet(256);
color_indices = round(color_values * (size(cmap, 1) - 1)) + 1;
colors = cmap(color_indices, :);

% Ploting arrows one by one to assign different colors
for i = 1:length(X)
    quiver3(X(i), Y(i), Z(i), ... 
            scale_factor * Ex_norm(i), ... 
            scale_factor * Ey_norm(i), ...
            scale_factor * Ez_norm(i), ...
            0, 'Color', colors(i,:), 'LineWidth', 0.6, ...
            'MaxHeadSize', 2, 'HandleVisibility', 'off');
end

% Create a dummy plot for legend
plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 2, 'DisplayName', 'Electric Field Vectors');

% Add transparent conical shells
r_shell = linspace(0, 2.2, 40);
phi_shell = linspace(0, 2*pi, 60);
[R_shell, PHI_shell] = meshgrid(r_shell, phi_shell);

% Inner cone
X_inner = R_shell * sin(teta0) .* cos(PHI_shell);
Y_inner = R_shell * sin(teta0) .* sin(PHI_shell);
Z_inner = R_shell * cos(teta0);
surf(X_inner, Y_inner, Z_inner, ...
     'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.05, ...
     'DisplayName', sprintf('\\theta = %.0f° (V = %.0f V)', teta0_deg, V0));

% Outer cone
X_outer = R_shell * sin(teta1) .* cos(PHI_shell);
Y_outer = R_shell * sin(teta1) .* sin(PHI_shell);
Z_outer = R_shell * cos(teta1);
surf(X_outer, Y_outer, Z_outer, ...
     'FaceColor', 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'g', 'EdgeAlpha', 0.05, ...
     'DisplayName', sprintf('\\theta = %.0f° (V = %.0f V)', teta1_deg, V1));

% Create colorbar
colormap(jet);
c = colorbar;
c.Label.String = 'Log_{10}(|E| + 1)  [V/m]';
c.Label.FontSize = 10;

% Add vertex
plot3(0, 0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k', ...
      'DisplayName', 'Vertex');

% Add axes
axis_limits = [-2.2, 2.2];
plot3(axis_limits, [0, 0], [0, 0], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
plot3([0, 0], axis_limits, [0, 0], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
plot3([0, 0], [0, 0], axis_limits, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');

text(2.2, 0, 0, 'x', 'FontSize', 12, 'FontWeight', 'bold');
text(0, 2.2, 0, 'y', 'FontSize', 12, 'FontWeight', 'bold');
text(0, 0, 2.2, 'z', 'FontSize', 12, 'FontWeight', 'bold');

% Set plot properties
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(sprintf('Part C: Electric Field - Normalized Arrows, Color-Coded Magnitude\n\\theta: %.0f° to %.0f°, V: %.0f V to %.0f V', ...
      teta0_deg, teta1_deg, V0, V1));
legend('Location', 'best', 'FontSize', 9);
grid on; axis equal; view(135, 25);
xlim([-2.2, 2.2]); ylim([-2.2, 2.2]); zlim([-2.2, 2.2]);

camlight left; lighting gouraud;
hold off;