classdef StochasticRelative < RelativeDynamics
    properties
        sig_rho
        sig_vel
        sig_mrp
        sig_omg
        mu_rho
        mu_vel
        mu_mrp
        mu_omg
        
        % LPF 객체 및 백색 노이즈 스펙
        lpf_rho
        lpf_vel
        lpf_mrp
        lpf_omg
        w_sig_rho
        w_sig_vel
        w_sig_mrp
        w_sig_omg
        
        % 제어기가 보게 될 관측 상태
        state_obs 
    end
    methods
        function obj = StochasticRelative(Cfg, targetSatellite)
            obj@RelativeDynamics(Cfg, targetSatellite);
            
            obj.sig_rho = Cfg.sigma.rho;
            obj.sig_vel = Cfg.sigma.vel;
            obj.sig_mrp = Cfg.sigma.mrp;
            obj.sig_omg = Cfg.sigma.omg;
            obj.mu_rho = Cfg.mu.rho;
            obj.mu_vel = Cfg.mu.vel;
            obj.mu_mrp = Cfg.mu.mrp;
            obj.mu_omg = Cfg.mu.omg;
            
            % 1차 LPF 초기화 (Cut-off Freq: 10 rad/s 가정)
            omg_c = 10;
            obj.lpf_rho = LPF(omg_c, Cfg);
            obj.lpf_vel = LPF(omg_c, Cfg);
            obj.lpf_mrp = LPF(omg_c, Cfg);
            obj.lpf_omg = LPF(omg_c, Cfg);
            
            % LPF의 분산 감쇠를 역산하여 입력 백색 노이즈의 분산 증폭
            alpha = obj.lpf_rho.alpha; 
            scale_factor = sqrt((1 + alpha) / (1 - alpha));
            obj.w_sig_rho = obj.sig_rho * scale_factor;
            obj.w_sig_vel = obj.sig_vel * scale_factor;
            obj.w_sig_mrp = obj.sig_mrp * scale_factor;
            obj.w_sig_omg = obj.sig_omg * scale_factor;
            
            obj.update_observation();
        end
        
        function update_observation(obj)
            % 1. 증폭된 분산을 가진 순수 백색 노이즈 생성
            w_rho = normrnd(obj.mu_rho, obj.w_sig_rho, [3, 1]);
            w_vel = normrnd(obj.mu_vel, obj.w_sig_vel, [3, 1]);
            w_mrp = normrnd(obj.mu_rho, obj.w_sig_rho, [3, 1]);
            w_omg = normrnd(obj.mu_vel, obj.w_sig_vel, [3, 1]);

            
            % 2. LPF 통과하여 논문과 동일한 Colored Noise 생성
            c_rho = obj.lpf_rho.forward(w_rho);
            c_vel = obj.lpf_vel.forward(w_vel);
            c_mrp = obj.lpf_mrp.forward(w_mrp);
            c_omg = obj.lpf_omg.forward(w_omg);
            
            % 3. 관측 상태 업데이트 (참값 + 유색 노이즈)
            obj.state_obs = obj.state; 
            % obj.state_obs(1:3) = obj.state_obs(1:3) + c_mrp;
            % obj.state_obs(4:6) = obj.state_obs(4:6) + c_omg;
            obj.state_obs(7:9) = obj.state_obs(7:9) + c_rho;
            obj.state_obs(10:12) = obj.state_obs(10:12) + c_vel;
        end
        
        function step(obj, u_ctrl, u_dist)
            % [핵심] 실제 우주선 물리 동역학은 참값(True)으로 업데이트 (노이즈 X)
            step@RelativeDynamics(obj, u_ctrl, u_dist);
            
            % 물리 업데이트가 끝난 참값에 센서 노이즈를 씌워 관측 상태 생성
            obj.update_observation();
        end
    end
end