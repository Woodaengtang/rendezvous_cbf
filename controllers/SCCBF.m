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
        function obj = SCCBF(ControlCfg, stochasticDynamics)
            obj@CCBF(ControlCfg, stochasticDynamics);

            obj.sig_rho = stochasticDynamics.sig_rho;
            obj.sig_vel = stochasticDynamics.sig_vel;
            obj.SigmaRho = diag(stochasticDynamics.sig_rho.^2);
            obj.SigmaVel = diag(stochasticDynamics.sig_vel.^2);

            obj.mu_rho = stochasticDynamics.mu_rho;
            obj.mu_vel = stochasticDynamics.mu_vel;
            obj.rho_eta = ControlCfg.eta_rho;
            obj.vel_eta = ControlCfg.eta_vel;
        end

        
        function ref_vel_cal(obj)
            ref_vel_cal@ClfQp(obj);
            obj.vel_nom = obj.ref_vel;
            obj.h_rho = obj.barrier_value();

            Lgh = obj.partial.dh_drho * obj.RD.R_tc';
            w_t = obj.RD.Target.stateECI(10:12);
            rho = obj.RD.state(7:9);

            c = (obj.partial.dh_drho * obj.RD.R_tc')';
            Phi_eta = norminv(obj.rho_eta);
            a_bar = obj.mu_rho;
            sCSCT = sqrt(c' * obj.SigmaRho * c);
            
            Lfh = -obj.partial.dh_drho * obj.RD.R_tc' * obj.RD.skew(obj.RD.R_tc*w_t) * rho - c' * a_bar - Phi_eta * sCSCT;
            
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
    end
end