classdef NominalRho < handle
    % NominalRho Implements CLF-CBF-QP based nominal controller
    % Incorporates coordinate transformation for CBF defined in Target Frame
    
    properties
        gamma   % CLF decay rate
        alpha   % CBF decay rate (Class K function)
        p_relax % High cost for slack variable
        
        qp_options % Quadprog options
        
        V       % Value of CLF
        LfV     % Lie derivative of V along f
        LgV     % Lie derivative of V along g
        h       % Value of CBF
        Lfh     % Lie derivative of h along f
        Lgh     % Lie derivative of h along g
        
        Alpha   % CBF Parameter: Cone shape factor
        Delta   % CBF Parameter: Shift along x-axis (Target frame)

        RelativeChaser  RelativeDynamics
    end
    
    methods
        function obj = NominalRho(control_setting, relativeDynamics)
            obj.gamma = control_setting.gamma;
            obj.alpha = control_setting.alpha;
            obj.p_relax = control_setting.p_relax;
            
            obj.qp_options = control_setting.qp_options;
            
            % Initialize properties
            obj.V   = 0;
            obj.LfV = 0;
            obj.LgV = zeros([1, 3]);
            obj.h   = 0;
            obj.Lfh = 0;
            obj.Lgh = zeros([1, 3]);
            
            % Geometric Parameters for the Barrier Function
            obj.Alpha = 0.02; 
            obj.Delta = 0.1;

            obj.RelativeChaser = relativeDynamics;
        end

        function [Vd, slack, feas] = command(obj)
            obj.update_clf();
            obj.update_cbf();
            
            quadprog_H = blkdiag(eye(3), obj.p_relax);
            quadprog_f = zeros([4, 1]);
            quadprog_f(4) = obj.p_relax;
            quadprog_f(1:3) = 0;

            A = [obj.LgV, -1;...
                 -obj.Lgh, 0];
            b = [-obj.gamma * obj.V - obj.LfV;...
                 obj.alpha * obj.h - obj.Lfh];
            [output, ~, exitflag, ~] = quadprog(quadprog_H, quadprog_f, A, b, [], [], [], [], [], obj.qp_options);
            Vd = output(1:3);
            slack = output(end);
            if exitflag == -2
                feas = 0;
                disp("Infeasible QP. CBF constraint is conflicting with input constraints.");
            else
                feas = 1;
            end
        end
        
        % ========================================================
        %%%%%%%%%% CLF Functions (V = 0.5 * rho' * rho) %%%%%%%%%%
        % ========================================================
        function update_clf(obj)
            state = obj.RelativeChaser.state;
            
            % Compute clf V
            rho = state(7:9);
            obj.V = 0.5 * (rho' * rho);

            % Compute  LfV
            S_omg_c = obj.RelativeChaser.skew(obj.RelativeChaser.omg_c);
            obj.LfV = -rho' * S_omg_c * rho;
            
            % Compute LgV
            obj.LgV = rho';
        end

        % =========================================================
        % CBF Functions (h in Target Frame)
        % =========================================================

        function update_cbf(obj)
            % Compute cbf h
            state = obj.RelativeChaser.state;
            R_tc = obj.RelativeChaser.get_Rt_c(state(1:3));
            r_t = R_tc' * state(7:9); 
            
            obj.h = obj.Alpha * (r_t(1) + obj.Delta)^3 - r_t(2)^2 - r_t(3)^2;
            
            % Compute Lfh and Lgh
            dh_dxt = 3 * obj.Alpha * (r_t(1) + obj.Delta)^2;
            dh_dyt = -2 * r_t(2);
            dh_dzt = -2 * r_t(3);
            grad_h_rt = [dh_dxt, dh_dyt, dh_dzt];
            
            S_omega = obj.RelativeChaser.skew(state(4:6));
            S_omega_c = obj.RelativeChaser.skew(obj.RelativeChaser.omg_c);
            rho = state(7:9);

            obj.Lfh = grad_h_rt * (-S_omega * R_tc' * rho - R_tc' * S_omega_c * rho);
            obj.Lgh = grad_h_rt * R_tc';
        end
    end
end