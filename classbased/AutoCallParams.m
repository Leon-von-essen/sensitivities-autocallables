classdef AutoCallParams
    properties
        s_0
        risk_free_rate
        sigma
        barrier_level
        corr_chol
        t_steps
        rebates
        survival_payoff
        dim
        t_len
        t_0
        obs_times
    end
    
    
    methods
        function obj = AutoCallParams(s_0, r, sigma, B, correlation_matrix, t_steps, rebates, survival_payoff)
            obj.s_0 = s_0;
            obj.risk_free_rate = r;
            obj.sigma = sigma;
            obj.barrier_level = B;
            obj.corr_chol = chol(correlation_matrix)';
            obj.t_steps = t_steps;
            obj.t_0 = 0;
            obj.rebates = rebates;
            obj.survival_payoff = survival_payoff;
            obj.dim = length(s_0);
            obj.t_len = length(t_steps);
            obj.obs_times = t_steps * triu(ones(obj.t_len));
        end
        
        function ret = disp(o)
            s0 = o.s_0
            risk_free = o.risk_free_rate
            sigma = o.sigma
            barrier = o.barrier_level
            correlation_matrix = o.corr_chol * o.corr_chol'
            rebates = o.rebates
            payoff = o.survival_payoff
            windows = o.obs_times
        end
        
        function obj = adj_t_init(obj)
            obj.t_steps(1) = obj.t_steps(1) - obj.t_0; 
        end
        
        function obj = shift_s(obj, dim, h)
            scaled_unit_vec = h * eye(obj.dim, dim);
            obj.s_0 = obj.s_0 + scaled_unit_vec;
        end
        
        function obj = shift_sigma(obj, dim, h)
            scaled_unit_vec = h * eye(obj.dim, dim);
            obj.sigma = obj.sigma + scaled_unit_vec;
        end
        
        function obj = shift_t(obj, h)
            obj.t_0 = obj.t_0 + h;
            obj = obj.adj_t_init();
            obj.obs_times(2:end) =  obj.obs_times(2:end) - h;
        end
        
        function obj = shift_r(obj, h)
            obj.risk_free_rate = obj.risk_free_rate + h;
        end
        
        function alt = shift(obj, shift_var, h, dim)
            switch shift_var
                case 's_0'
                    alt = obj.shift_s(dim, h);
                case 'sigma'
                    alt = obj.shift_sigma(dim, h);
                case 't_0'
                    alt = obj.shift_t(h);
                case 'risk_free_rate'
                    alt = obj.shift_r(h);
            end
        end
        
        
        function test = test_shift(obj, shift_var)    
            try
                copy_obj = obj.shift(shift_var, 0.1, 1);
                assert(isa(copy_obj, 'AutoCallParams'), ...
                        'TestAutoCallParams:CopyError',...
                        'Copy was initialized in the wrong class')
                assert(abs(copy_obj.(shift_var)(1) - 0.1 - obj.(shift_var)(1)) < 1e-16, ...
                    'TestAutoCallParams:ShiftError' ,...
                    'Shift distance was wrong')


            catch ME
                switch ME.identifier
                    case 'TestAutoCallParams:CopyError'
                        class(obj)
                        class(copy_obj)
                        rethrow(ME)
                    case 'TestAutoCallParams:ShiftError'
                        obj
                        copy_obj
                        obj.(shift_var)
                        copy_obj.(shift_var)
                        rethrow(ME)
                    otherwise
                        rethrow(ME)
                end
            end
            test = 1;
        end
    end
    
    
    methods(Static)
        function obj = test_instance_dim_3()
            payoff = @(s) (sum(s - 1));
            obj = AutoCallParams([1, 1, 1]' , 0.05, [0.1, 0.2, 0.05]', 1.05, ...
                                        [1, 0.9, 0.5; 0.9, 1, 0.4; 0.5, 0.4, 1], [1],... 
                                        [0.05], payoff);
        end
        
        function obj = test_instance_dim_2()
            payoff = @(s) (sum(s - 1));
            obj = AutoCallParams([1, 1]' , 0.05, [0.1, 0.2]', 1.05, ...
                                        [1, 0.9; 0.9, 1], [1, 1],... 
                                        [0.05, 0.1], payoff);
        end
        
        function test = run_class_tests()
            obj = AutoCallParams.test_instance();
            props = ["s_0", "sigma", "risk_free_rate", "t_0"];
            for shift = 1 : 4 
                shift_var = props(shift);
                test = obj.test_shift(shift_var);
            end
        end
    end
end
            
        