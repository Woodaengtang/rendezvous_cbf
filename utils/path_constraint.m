close all; clear; clc;

offset = 1;
x_min = 1;
x_max = 10 + offset;
x = linspace(x_min, x_max, 100);
theta = linspace(0, 2*pi, 50);

a_h = 0.1;
delta_h = 1;
[X, Theta] = meshgrid(x, theta);
R = sqrt(a_h*(X - delta_h).^3);
Y = R .* cos(Theta);
Z = R .* sin(Theta);

rtPlot = figure();
rtPlot.Theme = 'light';
hold on; grid on;
surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor', 'r');

meshData = readSurfaceMesh('SmallSat.glb');

V = double(meshData.Vertices);
F = double(meshData.Faces);
V = (eul2rotm([-deg2rad(90), 0 ,0])*V')';

scale_factor = 0.8; 
V = V * scale_factor;

trisurf(F, V(:,1), V(:,2), V(:,3), ...
      'FaceColor', [0.9, 0.9, 0.9], ...
      'EdgeColor', 'none', ...
      'FaceLighting', 'gouraud');

camlight('headlight');

view(3);
xlabel('x_t [m]'); ylabel('y_t [m]'); zlabel('z_t [m]');
axis equal;
xlim([-5, 10]); ylim([-10, 10]); 
saveas(gcf, 'assets/rt_plot.png');