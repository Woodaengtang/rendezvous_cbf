close all; clear all; clc;
addpath(genpath(pwd));
fig_size = [600, 675]; % 3x1 서브플롯에 최적화된 세로로 긴 사이즈
line_width = 1.5;

% --- 학술 논문용 전역 그래픽 설정 ---
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesFontSize', 12); 

% 색상 및 선 스타일 지정
color_ccbf = '#0072BD';  style_ccbf = '-';   
color_hocbf = '#D95319'; style_hocbf = '--'; 
color_nom = '#7E2F8E';   style_nom = '-.';   

% 데이터 로드
load('assets\nominal_clfqp\simLogger_20260326_085923.mat');
nominal = struct('sig', simLogger.loggerRelative.state.log(1:3, :),...
                 'omg', simLogger.loggerRelative.state.log(4:6, :),...
                 'rho', simLogger.loggerRelative.state.log(7:9, :),...
                 'vel', simLogger.loggerRelative.state.log(10:12, :),...
                 'h', simLogger.loggerBarrier.h.log(1,:));
                 
load('assets\ccbf\simLogger_20260326_090134.mat');
ccbf = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
              'force', simLogger.loggerControl.force.log,...
              'moment', simLogger.loggerControl.moment.log,...
              'sig', simLogger.loggerRelative.state.log(1:3, :),...
              'omg', simLogger.loggerRelative.state.log(4:6, :),...
              'rho', simLogger.loggerRelative.state.log(7:9, :),...
              'vel', simLogger.loggerRelative.state.log(10:12, :),...
              'h', simLogger.loggerBarrier.h.log(1,:));
              
load('assets\hocbf\simLogger_20260326_090025.mat');
hocbf = struct('lyapunov', simLogger.loggerLyapunov.V.log,...
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

% --- Trajectory Data Computation ---
nominal_traj = NaN([3, simCfg.sim_len]);
ccbf_traj = NaN([3, simCfg.sim_len]);
hocbf_traj = NaN([3, simCfg.sim_len]);
for i = 1:simCfg.sim_len
    nominal_traj(:, i) = (get_R_tc(nominal.sig(:, i)))' * nominal.rho(:, i);
    ccbf_traj(:, i) = (get_R_tc(ccbf.sig(:, i)))' * ccbf.rho(:, i);
    hocbf_traj(:, i) = (get_R_tc(hocbf.sig(:, i)))' * hocbf.rho(:, i);
end

% =========================================================================
% 1. 3D Trajectory (1x2)
% =========================================================================
meshData = readSurfaceMesh('SmallSat.glb');
V = double(meshData.Vertices);
F = double(meshData.Faces);
V = (eul2rotm([-deg2rad(90), 0 ,0])*V')';
V = V * 0.8; 

rtPlot = figure('Theme', 'light', 'Position', [50, 50, 1200, 450]); 
t_rt = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
view_angles = {[-35.5, 25], [0, 0]};
for i = 1:2
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    p_c = plot3(ax, ccbf_traj(1, :), ccbf_traj(2, :), ccbf_traj(3, :), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
    p_h = plot3(ax, hocbf_traj(1, :), hocbf_traj(2, :), hocbf_traj(3, :), '-', 'Color', color_hocbf, 'LineWidth', line_width);
    p_n = plot3(ax, nominal_traj(1, :), nominal_traj(2, :), nominal_traj(3, :), style_nom, 'Color', color_nom, 'LineWidth', line_width);
    trisurf(F, V(:,1), V(:,2), V(:,3), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'Parent', ax);
    view(ax, view_angles{i}); 
    xlabel(ax, '$x_t$ (m)'); ylabel(ax, '$y_t$ (m)'); zlabel(ax, '$z_t$ (m)');
    if i ~= 4, axis(ax, 'equal'); end
    if i == 1, legend(ax, [p_c, p_h, p_n], {'Cascaded CBF', 'HOCBF', 'Nominal'}, 'Location', 'northwest'); end
end

% =========================================================================
% 2. Relative Position (Rho) (3x1)
% =========================================================================
rhoFig = figure('Theme', 'light', 'Position', [100, 100, fig_size]);
t_rho = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylabels_rho = {'$\rho_x$ (m)', '$\rho_y$ (m)', '$\rho_z$ (m)'};
rho_targets = [1, 0, 0]; % 목표 랑데부 위치 설정 (\rho_x = 1m offset)

for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, simCfg.sim_time, ccbf.rho(ax_idx,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, hocbf.rho(ax_idx,:), '-', 'Color', color_hocbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, nominal.rho(ax_idx,:), style_nom, 'Color', color_nom, 'LineWidth', line_width);
    
    % 타겟 기준선 추가
    yline(ax, rho_targets(ax_idx), 'r--', 'LineWidth', 1.5, 'Alpha', 0.8);
    
    ylabel(ax, ylabels_rho{ax_idx});
    if ax_idx == 1 
        % 항목이 4개가 되었으므로 가독성을 위해 NumColumns를 2로 설정 (2x2 배열)
        legend(ax, {'Cascaded CBF', 'HOCBF', 'Nominal', 'Target'}, 'Location', 'northeast', 'NumColumns', 2); 
    end
    if ax_idx == 3, xlabel(ax, 'Time (s)'); end
end

% =========================================================================
% 3. Relative Velocity (Vel) (3x1)
% =========================================================================
velFig = figure('Theme', 'light', 'Position', [150, 150, fig_size]);
t_vel = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylabels_vel = {'$v_x$ (m/s)', '$v_y$ (m/s)', '$v_z$ (m/s)'};
vel_targets = [0, 0, 0]; % 목표 상대 속도 설정 (정지 상태)

for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, simCfg.sim_time, ccbf.vel(ax_idx,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
    
    % 사용자가 이전에 '-'로 설정했던 hocbf 선 스타일을 반영
    plot(ax, simCfg.sim_time, hocbf.vel(ax_idx,:), '-', 'Color', color_hocbf, 'LineWidth', line_width); 
    
    plot(ax, simCfg.sim_time, nominal.vel(ax_idx,:), style_nom, 'Color', color_nom, 'LineWidth', line_width);
    
    % 타겟 기준선 추가
    yline(ax, vel_targets(ax_idx), 'r--', 'LineWidth', 1.5, 'Alpha', 0.8);
    
    ylabel(ax, ylabels_vel{ax_idx});
    
    % 첫 번째 서브플롯에만 범례 추가 (Target 포함 2x2 배열)
    if ax_idx == 1 
        legend(ax, {'Cascaded CBF', 'HOCBF', 'Nominal', 'Target'}, 'Location', 'southeast', 'NumColumns', 2); 
    end
    if ax_idx == 3, xlabel(ax, 'Time (s)'); end
end

% =========================================================================
% 4. MRP (Sigma) (3x1)
% =========================================================================
sigFig = figure('Theme', 'light', 'Position', [200, 200, fig_size]);
t_sig = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylabels_sig = {'$\sigma_x$', '$\sigma_y$', '$\sigma_z$'};
sig_targets = [0, 0, 0]; % 목표 자세 설정 (MRP = 0)

for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, simCfg.sim_time, ccbf.sig(ax_idx,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, hocbf.sig(ax_idx,:), style_hocbf, 'Color', color_hocbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, nominal.sig(ax_idx,:), style_nom, 'Color', color_nom, 'LineWidth', line_width);
    
    % 타겟 기준선 추가
    yline(ax, sig_targets(ax_idx), 'r--', 'LineWidth', 1.5, 'Alpha', 0.8);
    
    ylabel(ax, ylabels_sig{ax_idx});
    if ax_idx == 1 
        legend(ax, {'Cascaded CBF', 'HOCBF', 'Nominal', 'Target'}, 'Location', 'northeast', 'NumColumns', 2); 
    end
    if ax_idx == 3, xlabel(ax, 'Time (s)'); end
end

% =========================================================================
% 5. Relative Angular Velocity (Omega) (3x1)
% =========================================================================
omgFig = figure('Theme', 'light', 'Position', [250, 250, fig_size]);
t_omg = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylabels_omg = {'$\omega_x$ (rad/s)', '$\omega_y$ (rad/s)', '$\omega_z$ (rad/s)'};
omg_targets = [0, 0, 0]; % 목표 상대 각속도 설정 (동기화 상태)

for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, simCfg.sim_time, ccbf.omg(ax_idx,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, hocbf.omg(ax_idx,:), style_hocbf, 'Color', color_hocbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, nominal.omg(ax_idx,:), style_nom, 'Color', color_nom, 'LineWidth', line_width);
    
    % 타겟 기준선 추가
    yline(ax, omg_targets(ax_idx), 'r--', 'LineWidth', 1.5, 'Alpha', 0.8);
    
    ylabel(ax, ylabels_omg{ax_idx});
    
    % 첫 번째 서브플롯에만 범례 추가 (Target 포함 2x2 배열)
    if ax_idx == 1 
        legend(ax, {'Cascaded CBF', 'HOCBF', 'Nominal', 'Target'}, 'Location', 'northeast', 'NumColumns', 2); 
    end
    if ax_idx == 3, xlabel(ax, 'Time (s)'); end
end

% =========================================================================
% 6. Control Force (3x1)
% =========================================================================
forceFig = figure('Theme', 'light', 'Position', [300, 300, fig_size]);
t_f = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylabels_f = {'$F_x$ (N)', '$F_y$ (N)', '$F_z$ (N)'};
for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, simCfg.sim_time, ccbf.force(ax_idx,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, hocbf.force(ax_idx,:), '-', 'Color', color_hocbf, 'LineWidth', line_width);
    yline(ax, f_bound, 'r:', 'LineWidth', 1.5, 'Alpha', 0.8);
    yline(ax, -f_bound, 'r:', 'LineWidth', 1.5, 'Alpha', 0.8);
    ylim(ax, [-f_bound*1.2, f_bound*1.2]);
    ylabel(ax, ylabels_f{ax_idx});
    if ax_idx == 1, legend(ax, {'Cascaded CBF', 'HOCBF', 'Saturation'}, 'Location', 'northeast', 'NumColumns', 3); end
    if ax_idx == 3, xlabel(ax, 'Time (s)'); end
end

% =========================================================================
% 7. Control Moment (3x1)
% =========================================================================
momentFig = figure('Theme', 'light', 'Position', [350, 350, fig_size]);
t_m = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ylabels_m = {'$M_x$ (Nm)', '$M_y$ (Nm)', '$M_z$ (Nm)'};
for ax_idx = 1:3
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, simCfg.sim_time, ccbf.moment(ax_idx,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
    plot(ax, simCfg.sim_time, hocbf.moment(ax_idx,:), style_hocbf, 'Color', color_hocbf, 'LineWidth', line_width);
    yline(ax, m_bound, 'r:', 'LineWidth', 1.5, 'Alpha', 0.8);
    yline(ax, -m_bound, 'r:', 'LineWidth', 1.5, 'Alpha', 0.8);
    ylim(ax, [-m_bound*1.2, m_bound*1.2]);
    ylabel(ax, ylabels_m{ax_idx});
    if ax_idx == 3, xlabel(ax, 'Time (s)'); end
end

% =========================================================================
% 8. Lyapunov Function Value (2x1)
% =========================================================================
% 기존 [600, 450] 고정 크기에서 다른 피겨들과 동일한 fig_size로 변경
lypFig = figure('Theme', 'light', 'Position', [400, 400, fig_size]);
t_lyp = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_v1 = nexttile; hold(ax_v1, 'on'); grid(ax_v1, 'on');
plot(ax_v1, simCfg.sim_time, ccbf.lyapunov(1,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
plot(ax_v1, simCfg.sim_time, hocbf.lyapunov(1,:), style_hocbf, 'Color', color_hocbf, 'LineWidth', line_width);
ylabel(ax_v1, '$V_\rho$'); 
legend(ax_v1, {'Cascaded CBF', 'HOCBF'}, 'Location', 'northeast');

ax_v2 = nexttile; hold(ax_v2, 'on'); grid(ax_v2, 'on');
plot(ax_v2, simCfg.sim_time, ccbf.lyapunov(2,:), style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
plot(ax_v2, simCfg.sim_time, hocbf.lyapunov(2,:), style_hocbf, 'Color', color_hocbf, 'LineWidth', line_width);
ylabel(ax_v2, '$V_v$');
xlabel(ax_v2, 'Time (s)');

% =========================================================================
% 9. Control Barrier Function Value (3x1: Full + Zoom)
% =========================================================================
hFig = figure('Theme', 'light', 'Position', [450, 450, fig_size]);
t_h = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Top plot (2/3 영역 차지)
ax_h1 = nexttile([2, 1]); hold(ax_h1, 'on'); grid(ax_h1, 'on');
plot(ax_h1, simCfg.sim_time, ccbf.h, style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
plot(ax_h1, simCfg.sim_time, hocbf.h, style_hocbf, 'Color', color_hocbf, 'LineWidth', line_width);
plot(ax_h1, simCfg.sim_time, nominal.h, style_nom, 'Color', color_nom, 'LineWidth', line_width);
yline(ax_h1, 0, 'r-', 'LineWidth', 1.5, 'Alpha', 0.5); 
ylabel(ax_h1, '$h(x)$');
ylim(ax_h1, [-1000, 20000]);
ax_h1.YAxis.Exponent = 3;
legend(ax_h1, {'Cascaded CBF', 'HOCBF', 'Nominal'}, 'Location', 'northeast');

% Bottom plot (1/3 영역 차지)
ax_h2 = nexttile; hold(ax_h2, 'on'); grid(ax_h2, 'on');
plot(ax_h2, simCfg.sim_time, ccbf.h, style_ccbf, 'Color', color_ccbf, 'LineWidth', line_width);
plot(ax_h2, simCfg.sim_time, hocbf.h, style_hocbf, 'Color', color_hocbf, 'LineWidth', line_width);
plot(ax_h2, simCfg.sim_time, nominal.h, style_nom, 'Color', color_nom, 'LineWidth', line_width);
yline(ax_h2, 0, 'r-', 'LineWidth', 1.5, 'Alpha', 0.5); 
ylim(ax_h2, [-1, 5]); 
xlabel(ax_h2, 'Time (s)');
ylabel(ax_h2, '$h(x)$ (Zoom)');

% =========================================================================
% 10. Export Figures (PDF & PNG)
% =========================================================================
save_dir = 'assets';
if ~exist(save_dir, 'dir')
    mkdir(save_dir); % 혹시 폴더가 없다면 생성
end

% 저장할 피겨 핸들과 파일명 지정
figs_to_save = { ...
    rtPlot,    '01_3D_Trajectory'; ...
    rhoFig,    '02_Relative_Position'; ...
    velFig,    '03_Relative_Velocity'; ...
    sigFig,    '04_MRP'; ...
    omgFig,    '05_Relative_Angular_Velocity'; ...
    forceFig,  '06_Control_Force'; ...
    momentFig, '07_Control_Moment'; ...
    lypFig,    '08_Lyapunov_Function'; ...
    hFig,      '09_Barrier_Function' ...
};

fprintf('Saving figures to %s folder...\n', save_dir);
for idx = 1:size(figs_to_save, 1)
    fig_handle = figs_to_save{idx, 1};
    file_name  = figs_to_save{idx, 2};

    % pdf_path = fullfile(save_dir, [file_name, '.pdf']);
    png_path = fullfile(save_dir, [file_name, '.png']);

    % % PDF 저장 (벡터 그래픽 유지, LaTeX 삽입 최적화)
    % exportgraphics(fig_handle, pdf_path, 'ContentType', 'vector');

    % PNG 저장 (300 DPI 고해상도)
    exportgraphics(fig_handle, png_path, 'Resolution', 300);

    fprintf('Saved: %s\n', file_name);
end
fprintf('All figures exported successfully!\n\n');

%% Helper Function
function R = get_R_tc(mrp)
    mrp_sq = mrp' * mrp;
    skew = [      0, -mrp(3),  mrp(2);...
             mrp(3),       0, -mrp(1);...
            -mrp(2),  mrp(1),       0];
    denom = (1 + mrp_sq)^2;
    R = eye(3) - (4*(1-mrp_sq)/denom)*skew + (8/denom)*(skew*skew);
end