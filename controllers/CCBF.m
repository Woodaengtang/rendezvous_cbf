classdef CCBF < ClfQp
    % Backstepping-based Control Lyapunov Function and Control Barrier Function Quadratic Program (CLF-CBF-QP) controller
    properties
        alpha_rho
        alpha_vel
        alpha_sig
        alpha_omg

        h_rho
        h_vel
        h_sig
        h_omg

        vel_nom
        vel_safe
        f_nom
        f_safe

        partial
    end

    methods
        function obj = CCBF(Cfg, relativeDynamics)
            obj@ClfQp(Cfg, relativeDynamics);

            obj.alpha_rho = Cfg.alpha_rho;
            obj.alpha_vel = Cfg.alpha_vel;
            obj.alpha_sig = Cfg.alpha_sig;
            obj.alpha_omg = Cfg.alpha_omg;

            obj.h_rho = NaN;
            obj.h_vel = NaN;
            obj.h_sig = NaN;
            obj.h_omg = NaN;

            obj.vel_nom = NaN;
            obj.f_nom = NaN;

            obj.partial = struct('dh_drho', [], 'd2h_drho2', []);
            obj.RD = relativeDynamics;
        end
        
        function u_ctrl = command(obj)
            % 현재 연결된 동역학 모델(Plant)이 관측 상태(state_obs)를 지원하는지 확인
            has_obs = isprop(obj.RD, 'state_obs') && ~isempty(obj.RD.state_obs);

            if has_obs
                % [State Swapping 트릭] 참값 백업 및 관측 상태 덮어쓰기
                true_state = obj.RD.state;
                true_R_tc = obj.RD.R_tc;

                obj.RD.state = obj.RD.state_obs;
                obj.RD.R_tc = obj.RD.get_R_tc(obj.RD.state_obs(1:3));
            end

            % 편미분 및 제어 명령 계산 (has_obs가 true면 관측값 기반으로 계산됨!)
            obj.partial_rho_t();
            obj.ref_vel_cal();
            obj.ref_omg_cal();
            u_ctrl = struct('f', obj.command_force(),...
                'tau', obj.command_torque());

            if has_obs
                % 계산 종료 후 시뮬레이터 물리 엔진을 위해 참값 원상복구
                obj.RD.state = true_state;
                obj.RD.R_tc = true_R_tc;
            end
        end

        function ref_vel_cal(obj)
            ref_vel_cal@ClfQp(obj);
            obj.vel_nom = obj.ref_vel;
            obj.h_rho = obj.barrier_value();

            Lgh = obj.partial.dh_drho * obj.RD.R_tc';
            w_t = obj.RD.Target.stateECI(10:12);
            rho = obj.RD.state(7:9);
            
            Lfh = -obj.partial.dh_drho * obj.RD.R_tc' * obj.RD.skew(obj.RD.R_tc*w_t) * rho;
            
            H = eye(3);
            f = -obj.vel_nom; 
            
            A = -Lgh;
            b = Lfh + obj.alpha_rho * obj.h_rho;
            % Solve
            [v_safe, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], obj.qp_option);
            if exitflag ~= 1
                warning('Velocity CBF-QP did not converge, using nominal control');
                v_safe = obj.vel_nom;
            end
            obj.ref_vel = v_safe;
        end

        function input_F = command_force(obj)
            obj.f_nom = command_force@ClfQp(obj);
            w_t = obj.RD.Target.stateECI(10:12);
            rho = obj.RD.state(7:9);
            vel = obj.RD.state(10:12);
            
            h = obj.barrier_value();
            dot_h = obj.partial.dh_drho * obj.RD.R_tc' * (vel - obj.RD.skew(obj.RD.R_tc*w_t)*rho);

            rho_t = obj.RD.R_tc' * rho;
            D2 = obj.RD.gravitational_force() - obj.RD.R_tc * obj.RD.Target.gravitational_force();

            dot_rho_t = obj.RD.R_tc' * (vel - obj.RD.skew(obj.RD.R_tc*w_t)*rho);
            
            Lfh = dot_rho_t' * obj.partial.d2h_drho2 * dot_rho_t...
                + obj.partial.dh_drho * (-2 * obj.RD.skew(w_t) * dot_rho_t - obj.RD.skew(w_t)^2 * rho_t + obj.RD.R_tc' * D2)...
                + obj.alpha_rho * dot_h;
            Lgh = obj.partial.dh_drho * obj.RD.R_tc' ./ obj.RD.m_c;
            h_2 = dot_h + obj.alpha_rho * h;
            
            H = eye(3);
            f = -obj.f_nom;

            A = -Lgh;
            b = Lfh + obj.alpha_vel * h_2;
            
            % Solve
            [force_safe, ~, exitflag] = quadprog(H, f, A, b, [], [], obj.force_lb(1:3), obj.force_ub(1:3), [], obj.qp_option);
        
            if exitflag ~= 1
                warning('Force CBF-QP did not converge, using nominal control');
                force_safe = obj.f_nom;
            end
            input_F = force_safe;
            obj.f_safe = input_F;
        end

        function partial_rho_t(obj)
            rho = obj.RD.state(7:9);
            rho_t = obj.RD.R_tc' * rho;
            
            grad_rho_t = [3*obj.a_h*(rho_t(1) - obj.delta_h)^2,...
                       -2*rho_t(2),...
                       -2*rho_t(3)];
            hessian_rho_t = [6*obj.a_h*(rho_t(1) - obj.delta_h), 0, 0;...
                         0, -2, 0;...
                         0, 0, -2];
            obj.partial.dh_drho = grad_rho_t;
            obj.partial.d2h_drho2 = hessian_rho_t;
        end
    end
end