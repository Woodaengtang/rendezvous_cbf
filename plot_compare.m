close all; clear all; clc;

addpath(genpath(pwd));

fig_size = [600, 450];
line_width = 1.5;

load('assets\nominal_clfqp\simLogger_20260326_085923.mat');
nominal = struct('sig', simLogger.loggerRelative.state.log(1:3, :),...
                 'rho', simLogger.loggerRelative.state.log(7:9, :),...
                 'h', simLogger.loggerBarrier.h.log(1,:));


load('assets\ccbf\simLogger_20260326_090134.mat');
ccbf = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
              'force', simLogger.loggerControl.force.log,...
              'sig', simLogger.loggerRelative.state.log(1:3, :),...
              'rho', simLogger.loggerRelative.state.log(7:9, :),...
              'h', simLogger.loggerBarrier.h.log(1,:));

load('assets\hocbf\simLogger_20260326_090025.mat');
hocbf = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
               'force', simLogger.loggerControl.force.log,...
               'sig', simLogger.loggerRelative.state.log(1:3, :),...
               'rho', simLogger.loggerRelative.state.log(7:9, :),...
               'h', simLogger.loggerBarrier.h.log(1,:));

simCfg = SimCfg();

lyapunovFig = figure();
lyapunovFig.Theme = 'light';
lyapunovFig.Position(3:4) = fig_size;
subplot(2, 1, 1);
hold on; grid on;
CcbfLyp = plot(simCfg.sim_time, ccbf.lyapunov(1,:), 'LineWidth', line_width);
HocbfLyp = plot(simCfg.sim_time, hocbf.lyapunov(1,:), 'LineWidth', line_width);
ylabel('V_\rho'); legend([CcbfLyp, HocbfLyp], {'Cascaded CBF', 'HOCBF'}, 'Location', 'northeast');

subplot(2, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, ccbf.lyapunov(2,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.lyapunov(2,:), 'LineWidth', line_width);
ylabel('V_v');

inputFig = figure();
inputFig.Theme = 'light';
inputFig.Position(3:4) = fig_size;
subplot(3, 1, 1);
hold on; grid on;
plot(simCfg.sim_time, ccbf.force(1,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.force(1,:), 'LineWidth', line_width);
ylabel('F_x');
legend({'Cascaded CBF', 'HOCBF'}, 'Location', 'northeast');

subplot(3, 1, 2);
hold on; grid on;
plot(simCfg.sim_time, ccbf.force(2,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.force(2,:), 'LineWidth', line_width);
ylabel('F_y');

subplot(3, 1, 3);
hold on; grid on;
plot(simCfg.sim_time, ccbf.force(3,:), 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.force(3,:), 'LineWidth', line_width);
ylabel('F_z');

nominal_traj = NaN([3, simCfg.sim_len]);
ccbf_traj = NaN([3, simCfg.sim_len]);
hocbf_traj = NaN([3, simCfg.sim_len]);
for i = 1:simCfg.sim_len
    nominal_sigma = nominal.sig(:, i);
    nominal_rho = nominal.rho(:, i);
    nominal_traj(:, i) = (get_R_tc(nominal_sigma))'*nominal_rho;

    ccbf_sigma = ccbf.sig(:, i);
    ccbf_rho = ccbf.rho(:, i);
    ccbf_traj(:, i) = (get_R_tc(ccbf_sigma))'*ccbf_rho;

    hocbf_sigma = hocbf.sig(:, i);
    hocbf_rho = hocbf.rho(:, i);
    hocbf_traj(:, i) = (get_R_tc(hocbf_sigma))'*hocbf_rho;
end

meshData = readSurfaceMesh('SmallSat.glb');
V = double(meshData.Vertices);
F = double(meshData.Faces);
V = (eul2rotm([-deg2rad(90), 0 ,0])*V')';
scale_factor = 0.8; 
V = V * scale_factor;

rtPlot = figure();
rtPlot.Theme = 'light';
rtPlot.Position(3:4) = [1000, 800]; 

t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

view_angles = {[-46.5, 25], [0, 0], [0, 90], [90, 0]};

for i = 1:4
    nexttile;
    hold on; grid on;
    
    plot3(ccbf_traj(1, :), ccbf_traj(2, :), ccbf_traj(3, :), 'LineWidth', line_width);
    plot3(hocbf_traj(1, :), hocbf_traj(2, :), hocbf_traj(3, :), 'LineWidth', line_width);
    plot3(nominal_traj(1, :), nominal_traj(2, :), nominal_traj(3, :), 'LineWidth', line_width);
    
    trisurf(F, V(:,1), V(:,2), V(:,3), ...
          'FaceColor', [0.7 0.7 0.7], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
    view(view_angles{i});
    camlight('headlight');
    
    xlabel('x_t [m]'); ylabel('y_t [m]'); zlabel('z_t [m]');
    if i ~= 4
        axis equal;
    end
    
    if i == 1
        legend({'Cascaded CBF', 'HOCBF', 'Nominal Controller'}, 'Location', 'northwest');
    end
end

hPlot = figure();
hPlot.Theme = 'light';
hPlot.Position(3:4) = fig_size;
hold on; grid on;
plot(simCfg.sim_time, ccbf.h, 'LineWidth', line_width);
plot(simCfg.sim_time, hocbf.h, 'LineWidth', line_width);
plot(simCfg.sim_time, nominal.h, 'LineWidth', line_width);
xlabel('Time (s)'); ylabel('h');
ylim([-0.5, 5]);
legend({'Cascaded CBF', 'HOCBF', 'Nominal Controller'}, 'Location', 'northeast');

%%
load("assets\s_PCCBF\simLogger_20260326_092938.mat");
PCCBF = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
              'force', simLogger.loggerControl.force.log,...
              'sig', simLogger.loggerRelative.state.log(1:3, :),...
              'rho', simLogger.loggerRelative.state.log(7:9, :),...
              'h', simLogger.loggerBarrier.h.log(1,:));

load("assets\s_CCBF\simLogger_20260326_092805.mat");
CCBF = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
              'force', simLogger.loggerControl.force.log,...
              'sig', simLogger.loggerRelative.state.log(1:3, :),...
              'rho', simLogger.loggerRelative.state.log(7:9, :),...
              'h', simLogger.loggerBarrier.h.log(1,:));

phPlot = figure();
phPlot.Theme = 'light';
phPlot.Position(3:4) = [fig_size(1), fig_size(2)*1.5];
t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile([2, 1]); 
hold(ax1, 'on'); grid(ax1, 'on');
p1 = plot(ax1, simCfg.sim_time, PCCBF.h, 'LineWidth', line_width);
p2 = plot(ax1, simCfg.sim_time, CCBF.h, 'LineWidth', line_width);
margin1 = plot(ax1, simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
uistack(margin1, 'bottom');
ylabel(ax1, 'h');
legend(ax1, [p1, p2], {'PCCBF', 'CCBF'}, 'Location', 'northeast');
ax2 = nexttile([1, 1]);
hold(ax2, 'on'); grid(ax2, 'on');
p1 = plot(ax2, simCfg.sim_time, PCCBF.h, 'LineWidth', line_width);
p2 = plot(ax2, simCfg.sim_time, CCBF.h, 'LineWidth', line_width);
margin1 = plot(ax2, simCfg.sim_time, zeros(simCfg.sim_len, 1), 'r--', 'LineWidth', line_width);
ylim([-0.2, 1]);
ylabel(ax2, 'h');
xlabel(ax2, 'Time (s)');

pccbf_traj = NaN([3, simCfg.sim_len]);
ccbf_traj = NaN([3, simCfg.sim_len]);
for i = 1:simCfg.sim_len
    pccbf_sigma = PCCBF.sig(:, i);
    pccbf_rho = PCCBF.rho(:, i);
    pccbf_traj(:, i) = (get_R_tc(pccbf_sigma))'*pccbf_rho;

    ccbf_sigma = CCBF.sig(:, i);
    ccbf_rho = CCBF.rho(:, i);
    ccbf_traj(:, i) = (get_R_tc(ccbf_sigma))'*ccbf_rho;
end

rtSPlot = figure();
rtSPlot.Theme = 'light';
rtSPlot.Position(3:4) = [1000, 800]; 

t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

view_angles = {[-46.5, 25], [0, 0], [0, 90], [90, 0]};

for i = 1:4
    nexttile;
    hold on; grid on;
    
    plot3(pccbf_traj(1, :), pccbf_traj(2, :), pccbf_traj(3, :), 'LineWidth', line_width);
    plot3(ccbf_traj(1, :), ccbf_traj(2, :), ccbf_traj(3, :), 'LineWidth', line_width);
    
    trisurf(F, V(:,1), V(:,2), V(:,3), ...
          'FaceColor', [0.7 0.7 0.7], ...
          'EdgeColor', 'none', ...
          'FaceLighting', 'gouraud');
    view(view_angles{i});
    camlight('headlight');
    
    xlabel('x_t [m]'); ylabel('y_t [m]'); zlabel('z_t [m]');
    if i ~= 4
        axis equal;
    end
    
    if i == 1
        legend({'PCCBF', 'CCBF'}, 'Location', 'northwest');
    end
end

%%

function R = get_R_tc(mrp)
    mrp_sq = mrp' * mrp;
    skew = [      0, -mrp(3),  mrp(2);...
             mrp(3),       0, -mrp(1);...
            -mrp(2),  mrp(1),       0];
    denom = (1 + mrp_sq)^2;
    R = eye(3) - (4*(1-mrp_sq)/denom)*skew + (8/denom)*(skew*skew);
end
