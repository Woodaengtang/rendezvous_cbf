classdef SCCBF < CCBF
    properties
        sig_rho
        sig_vel
        SigmaRho
        SigmaVel
        mu_rho
        mu_vel
        eta_rho
        eta_vel
    end
    methods
        function obj = SCCBF(Cfg, stochasticDynamics)
            obj@CCBF(Cfg, stochasticDynamics);
            obj.sig_rho = stochasticDynamics.sig_rho;
            obj.sig_vel = stochasticDynamics.sig_vel;
            obj.SigmaRho = diag(Cfg.sigma.rho.^2);
            obj.SigmaVel = diag(Cfg.sigma.vel.^2);
            obj.mu_rho = Cfg.mu.rho;
            obj.mu_vel = Cfg.mu.vel;
            
            obj.eta_rho = Cfg.eta_rho; 
            obj.eta_vel = Cfg.eta_vel;
        end
        
        function u_ctrl = command(obj)
            % =========================================================
            % [상태 스왑 (State Swapping) 트릭]
            % 명목 제어기와 편미분 함수들이 관측 상태(Observed state)를 
            % 기반으로 계산되도록, 잠시 참값을 관측값으로 덮어씁니다.
            % =========================================================
            true_state = obj.RD.state;
            true_R_tc = obj.RD.R_tc;
            
            % 제어기가 인식하는 관측 상태 대입
            obj.RD.state = obj.RD.state_obs;
            obj.RD.R_tc = obj.RD.get_R_tc(obj.RD.state_obs(1:3));
            
            % 모든 계산은 이제 관측값(Noisy data) 기반으로 이루어짐
            obj.partial_rho_t();
            obj.ref_vel_cal();
            obj.ref_omg_cal();
            u_ctrl = struct('f', obj.command_force(),...
                            'tau', obj.command_torque());
                            
            % 계산 종료 후, 시뮬레이터 로깅을 위해 다시 참값으로 원상복구
            obj.RD.state = true_state;
            obj.RD.R_tc = true_R_tc;
        end
        
        function ref_vel_cal(obj)
            ref_vel_cal@ClfQp(obj);
            obj.vel_nom = obj.ref_vel;
            obj.h_rho = obj.barrier_value();
            Lgh = obj.partial.dh_drho * obj.RD.R_tc';
            w_t = obj.RD.Target.stateECI(10:12);
            rho = obj.RD.state(7:9); % 여기서의 state는 이미 Swap된 관측 상태!
            
            % [1st Layer] Cantelli's Margin Calculation
            c = obj.partial.dh_drho * obj.RD.R_tc';
            mu_w1 = c * obj.mu_rho;
            sigma_w1 = sqrt(c * obj.SigmaRho * c');
            
            k_rho = sqrt(obj.eta_rho / (1 - obj.eta_rho));
            cantelli_margin_1 = mu_w1 - k_rho * sigma_w1;
            
            Lfh = -obj.partial.dh_drho * obj.RD.R_tc' * obj.RD.skew(obj.RD.R_tc*w_t) * rho;
            
            H = eye(3); f = -obj.vel_nom; A = -Lgh;
            b = Lfh + obj.alpha_rho * obj.h_rho + cantelli_margin_1; 
            
            [v_safe, ~, exitflag] = quadprog(H, f, A, b, [], [], [], [], [], obj.qp_option);
            if exitflag ~= 1, v_safe = obj.vel_nom; end
            obj.ref_vel = v_safe;
        end
        
        function input_F = command_force(obj)
            obj.f_nom = command_force@ClfQp(obj);
            w_t = obj.RD.Target.stateECI(10:12);
            rho = obj.RD.state(7:9);
            vel = obj.RD.state(10:12);
            omg = obj.RD.state(4:6);

            h = obj.barrier_value();
            dot_rho_t = obj.RD.R_tc' * (vel - obj.RD.skew(obj.RD.R_tc*w_t)*rho);
            dot_h = obj.partial.dh_drho * dot_rho_t;
            rho_t = obj.RD.R_tc' * rho;
            D2 = obj.RD.gravitational_force() - obj.RD.R_tc * obj.RD.Target.gravitational_force();

            % [2nd Layer] Cantelli's Margin Calculation (Isserlis' Theorem)
            Q_mat = obj.RD.R_tc * obj.partial.d2h_drho2 * obj.RD.R_tc';
            M_mat = obj.RD.R_tc' * (obj.RD.skew(omg) - obj.RD.skew(obj.RD.R_tc*w_t));

            L_lin = 2 * dot_rho_t' * obj.partial.d2h_drho2 * obj.RD.R_tc' ...
                + obj.partial.dh_drho * M_mat ...
                + obj.alpha_rho * obj.partial.dh_drho * obj.RD.R_tc';
            K_lin = obj.partial.dh_drho * obj.RD.R_tc';

            % Exact Moments of the Polynomial Gaussian
            mu_quad = trace(Q_mat * obj.SigmaRho) + obj.mu_rho' * Q_mat * obj.mu_rho;
            mu_lin = L_lin * obj.mu_rho + K_lin * obj.mu_vel;
            mu_w2 = mu_quad + mu_lin;

            var_quad = 2 * trace((Q_mat * obj.SigmaRho)^2) + 4 * obj.mu_rho' * Q_mat * obj.SigmaRho * Q_mat * obj.mu_rho;
            var_lin = L_lin * obj.SigmaRho * L_lin' + K_lin * obj.SigmaVel * K_lin';
            cov_lin_quad = 2 * L_lin * obj.SigmaRho * Q_mat * obj.mu_rho;
            sigma_w2 = sqrt(var_quad + var_lin + cov_lin_quad);

            k_vel = sqrt(obj.eta_vel / (1 - obj.eta_vel));
            cantelli_margin_2 = mu_w2 - k_vel * sigma_w2;

            % Nominal terms
            Lfh = dot_rho_t' * obj.partial.d2h_drho2 * dot_rho_t...
                + obj.partial.dh_drho * (-2 * obj.RD.skew(w_t) * dot_rho_t - obj.RD.skew(w_t)^2 * rho_t + obj.RD.R_tc' * D2)...
                + obj.alpha_rho * dot_h;
            Lgh = obj.partial.dh_drho * obj.RD.R_tc' ./ obj.RD.m_c;
            h_2 = dot_h + obj.alpha_rho * h;

            H = eye(3);
            f = -obj.f_nom;
            A = -Lgh;
            b = Lfh + obj.alpha_vel * h_2 + cantelli_margin_2;

            [force_safe, ~, exitflag] = quadprog(H, f, A, b, [], [], obj.force_lb(1:3), obj.force_ub(1:3), [], obj.qp_option);
            if exitflag ~= 1, force_safe = obj.f_nom; end
            input_F = force_safe;
            obj.f_safe = input_F;
        end
    end
end