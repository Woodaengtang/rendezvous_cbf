close all; clear; clc;
addpath(genpath(pwd));

fig_size = [600, 675]; 
line_width = 1.5;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 12); 

color_pccbf = '#0072BD';  style_pccbf = '-';   

load('assets\s_PCCBF\simLogger_20260326_092938.mat');
PCCBF = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
               'force', simLogger.loggerControl.force.log,...
               'moment', simLogger.loggerControl.moment.log,...
               'sig', simLogger.loggerRelative.state.log(1:3, :),...
               'omg', simLogger.loggerRelative.state.log(4:6, :),...
               'rho', simLogger.loggerRelative.state.log(7:9, :),...
               'vel', simLogger.loggerRelative.state.log(10:12, :),...
               'h', simLogger.loggerBarrier.h.log(1,:));
              
simCfg = SimCfg();
f_bound = 20; % Control saturation bound (N)
m_bound = 5;  % Control Moment saturation bound (Nm)

pccbf_traj = NaN([3, simCfg.sim_len]);
for i = 1:simCfg.sim_len
    pccbf_traj(:, i) = (get_R_tc(PCCBF.sig(:, i)))' * PCCBF.rho(:, i);
end

meshData = readSurfaceMesh('SmallSat.glb');
V_mesh = double(meshData.Vertices);
F_mesh = double(meshData.Faces);
V_mesh = (eul2rotm([-deg2rad(90), 0 ,0])*V_mesh')';
V_mesh = V_mesh * 0.8; 

%%
rtPlot = figure('Theme', 'light', 'Position', [50, 50, 1200, 450]); 
t_rt = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
view_angles = {[-35.5, 25], [0, 0]};
for i = 1:2
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot3(ax, pccbf_traj(1, :), pccbf_traj(2, :), pccbf_traj(3, :), style_pccbf, 'Color', color_pccbf, 'LineWidth', line_width);
    trisurf(F_mesh, V_mesh(:,1), V_mesh(:,2), V_mesh(:,3), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'Parent', ax);
    view(ax, view_angles{i}); 
    xlabel(ax, '$x_t$ (m)'); ylabel(ax, '$y_t$ (m)'); zlabel(ax, '$z_t$ (m)');
    if i ~= 2, axis(ax, 'equal'); end
    if i == 1, legend(ax, {'PCCBF'}, 'Location', 'northwest'); end
end

rho_targets = [1, 0, 0];
sig_targets = [0, 0, 0];
vel_targets = [0, 0, 0];
omg_targets = [0, 0, 0];

rhoFig = plot_3x1_alone(simCfg.sim_time, PCCBF.rho, {'$\rho_x$ (m)', '$\rho_y$ (m)', '$\rho_z$ (m)'}, 'PCCBF', color_pccbf, line_width, fig_size, rho_targets);
velFig = plot_3x1_alone(simCfg.sim_time, PCCBF.vel, {'$v_x$ (m/s)', '$v_y$ (m/s)', '$v_z$ (m/s)'}, 'PCCBF', color_pccbf, line_width, fig_size, vel_targets);
sigFig = plot_3x1_alone(simCfg.sim_time, PCCBF.sig, {'$\sigma_1$', '$\sigma_2$', '$\sigma_3$'}, 'PCCBF', color_pccbf, line_width, fig_size, sig_targets);
omgFig = plot_3x1_alone(simCfg.sim_time, PCCBF.omg, {'$\omega_x$ (rad/s)', '$\omega_y$ (rad/s)', '$\omega_z$ (rad/s)'}, 'PCCBF', color_pccbf, line_width, fig_size, omg_targets);

forceFig  = plot_3x1_alone_sat(simCfg.sim_time, PCCBF.force, {'$F_x$ (N)', '$F_y$ (N)', '$F_z$ (N)'}, 'PCCBF', color_pccbf, line_width, fig_size, f_bound);
momentFig = plot_3x1_alone_sat(simCfg.sim_time, PCCBF.moment, {'$M_x$ (Nm)', '$M_y$ (Nm)', '$M_z$ (Nm)'}, 'PCCBF', color_pccbf, line_width, fig_size, m_bound);

%%
lypFig = figure('Theme', 'light', 'Position', [400, 400, fig_size]);
t_lyp = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_v1 = nexttile; hold(ax_v1, 'on'); grid(ax_v1, 'on');
plot(ax_v1, simCfg.sim_time, PCCBF.lyapunov(1,:), style_pccbf, 'Color', color_pccbf, 'LineWidth', line_width);
ylabel(ax_v1, '$V_\rho$'); 
legend(ax_v1, {'PCCBF'}, 'Location', 'northeast');

ax_v2 = nexttile; hold(ax_v2, 'on'); grid(ax_v2, 'on');
plot(ax_v2, simCfg.sim_time, PCCBF.lyapunov(2,:), style_pccbf, 'Color', color_pccbf, 'LineWidth', line_width);
ylabel(ax_v2, '$V_v$');
xlabel(ax_v2, 'Time (s)');


%%
hFig = figure('Theme', 'light', 'Position', [450, 450, fig_size]);
t_h = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top plot (2/3 영역 차지)
ax_h1 = nexttile([2, 1]); hold(ax_h1, 'on'); grid(ax_h1, 'on');
plot(ax_h1, simCfg.sim_time, PCCBF.h, style_pccbf, 'Color', color_pccbf, 'LineWidth', line_width);
yline(ax_h1, 0, 'r-', 'LineWidth', 1.5, 'Alpha', 0.5); 
ylabel(ax_h1, '$h(x)$');
ylim(ax_h1, [-1000, 20000]);
ax_h1.YAxis.Exponent = 3;
legend(ax_h1, {'PCCBF', 'Safety Bound'}, 'Location', 'northeast');

% Bottom plot (1/3 영역 차지, Zoom)
ax_h2 = nexttile; hold(ax_h2, 'on'); grid(ax_h2, 'on');
plot(ax_h2, simCfg.sim_time, PCCBF.h, style_pccbf, 'Color', color_pccbf, 'LineWidth', line_width);
yline(ax_h2, 0, 'r-', 'LineWidth', 1.5, 'Alpha', 0.5); 
ylim(ax_h2, [-1, 5]); 
xlabel(ax_h2, 'Time (s)');
ylabel(ax_h2, '$h(x)$ (Zoom)');

% =========================================================================
% 10. Export All Figures to Assets
% =========================================================================
save_dir = 'assets';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

figs_to_save = { ...
    rtPlot,    's01_PCCBF_3D_Trajectory'; ...
    rhoFig,    's02_PCCBF_Relative_Position'; ...
    velFig,    's03_PCCBF_Relative_Velocity'; ...
    sigFig,    's04_PCCBF_MRP'; ...
    omgFig,    's05_PCCBF_Relative_Angular_Velocity'; ...
    forceFig,  's06_PCCBF_Control_Force'; ...
    momentFig, 's07_PCCBF_Control_Moment'; ...
    lypFig,    's08_PCCBF_Lyapunov_Function'; ...
    hFig,      's09_PCCBF_Barrier_Function' ...
};

fprintf('Exporting Stochastic PCCBF Figures...\n');
for idx = 1:size(figs_to_save, 1)
    f_handle = figs_to_save{idx, 1};
    f_name   = figs_to_save{idx, 2};
    
    % PNG (3D 플롯 포함 모두 기본 저장)
    exportgraphics(f_handle, fullfile(save_dir, [f_name, '.png']), 'Resolution', 300);
    
    % % PDF (3D 플롯을 제외한 2D 그래프들만 저장)
    % if ~contains(f_name, '_3D')
    %     exportgraphics(f_handle, fullfile(save_dir, [f_name, '.pdf']), 'ContentType', 'vector');
    % end
    fprintf('Saved: %s\n', f_name);
end
fprintf('All stochastic standalone figures exported successfully!\n\n');

%% Helper Functions
function R = get_R_tc(mrp)
    mrp_sq = mrp' * mrp;
    skew = [0, -mrp(3), mrp(2); mrp(3), 0, -mrp(1); -mrp(2), mrp(1), 0];
    denom = (1 + mrp_sq)^2;
    R = eye(3) - (4*(1-mrp_sq)/denom)*skew + (8/denom)*(skew*skew);
end

function fig = plot_3x1_alone(t, data, y_labels, leg_str, color, lw, f_size, targets)
    % targets 인자가 들어오지 않았을 경우 빈 배열로 처리
    if nargin < 8
        targets = [];
    end
    
    fig = figure('Theme', 'light', 'Position', [randi([200,400]), randi([200,400]), f_size]);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for i = 1:3
        ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
        plot(ax, t, data(i,:), '-', 'Color', color, 'LineWidth', lw);
        
        % 타겟 값이 존재하면 빨간 점선 추가
        if ~isempty(targets)
            yline(ax, targets(i), 'r--', 'LineWidth', 1.5, 'Alpha', 0.8);
        end
        
        ylabel(ax, y_labels{i});
        
        % 첫 번째 서브플롯에만 범례 추가
        if i == 1
            if isempty(targets)
                legend(ax, {leg_str}, 'Location', 'northeast');
            else
                legend(ax, {leg_str, 'Target'}, 'Location', 'northeast', 'NumColumns', 2);
            end
        end
        
        if i == 3
            xlabel(ax, 'Time (s)'); 
        end
    end
end

function fig = plot_3x1_alone_sat(t, data, y_labels, leg_str, color, lw, f_size, bound)
    fig = figure('Theme', 'light', 'Position', [randi([200,400]), randi([200,400]), f_size]);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:3
        ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
        plot(ax, t, data(i,:), '-', 'Color', color, 'LineWidth', lw);
        yline(ax, bound, 'r:', 'LineWidth', 1.5, 'Alpha', 0.8);
        yline(ax, -bound, 'r:', 'LineWidth', 1.5, 'Alpha', 0.8);
        ylim(ax, [-bound*1.2, bound*1.2]); ylabel(ax, y_labels{i});
        if i == 1, legend(ax, {leg_str, 'Saturation'}, 'Location', 'northeast', 'NumColumns', 2); end
        if i == 3, xlabel(ax, 'Time (s)'); end
    end
end