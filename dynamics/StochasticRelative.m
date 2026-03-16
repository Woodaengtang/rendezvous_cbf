classdef StochasticRelative < RelativeDynamics

    properties
        % sigma denotes standarad deviation of Normal Distribution
        sig_rho
        sig_vel
        sig_mrp
        sig_omg

        mu_rho
        mu_vel
        mu_mrp
        mu_omg
    end

    methods
        function obj = StochasticRelative(sim_cfg, targetSatellite)
            obj@RelativeDynamics(sim_cfg, targetSatellite);
            
            obj.sig_rho = sim_cfg.sigma.rho;
            obj.sig_vel = sim_cfg.sigma.vel;
            obj.sig_mrp = sim_cfg.sigma.mrp;
            obj.sig_omg = sim_cfg.sigma.omg;

            obj.mu_rho = sim_cfg.mu.rho;
            obj.mu_vel = sim_cfg.mu.vel;
            obj.mu_mrp = sim_cfg.mu.mrp;
            obj.mu_omg = sim_cfg.mu.omg;
        end

        function dstate = dynamics(obj, x, u_ctrl, u_dist)
            dstate_ = dynamics@RelativeDynamics(obj, x, u_ctrl, u_dist);
            dstate = dstate_ + [normrnd(obj.mu_mrp, obj.sig_mrp);...
                                normrnd(obj.mu_omg, obj.sig_omg);...
                                normrnd(obj.mu_rho, obj.sig_rho);...
                                normrnd(obj.mu_vel, obj.sig_vel)];
        end

        function step(obj, u_ctrl, u_dist)
            % STEP Performs RK4 integration and updates internal state
            %
            % Usage: obj.step(dt, u_ctrl, u_dist, target_state)
            
            x_curr = obj.state;
            
            % RK4 Integration
            k1 = obj.dynamics(              x_curr, u_ctrl, u_dist);
            k2 = obj.dynamics(x_curr + obj.dt/2*k1, u_ctrl, u_dist);
            k3 = obj.dynamics(x_curr + obj.dt/2*k2, u_ctrl, u_dist);
            k4 = obj.dynamics(  x_curr + obj.dt*k3, u_ctrl, u_dist);
            
            % Update State
            obj.state = x_curr + (obj.dt/6) * (k1 + 2*k2 + 2*k3 + k4);

            obj.get_chaser_pos();
            obj.get_chaser_euler();
            obj.R_tc = obj.get_R_tc(obj.state(1:3));
        end
    end
end