classdef SmoothPathGenerator < pathGenerator
    properties
        R
        Qs
        T
        proj_Rs
        proj_bs
        surv_mat
        non_surv
        sample_surv_payoff
    end
    methods
        function obj = SmoothPathGenerator(params)
            obj@pathGenerator(params);
            [obj.T, obj.R, obj.proj_Rs, obj.proj_bs, obj.Qs] = loadProjectedMatrices(params.corr_chol);
            [obj.surv_mat, obj.non_surv] = probEvalMatrix(params.dim, params.t_len + 1);
            obj.sample_surv_payoff = @(M) params.survival_payoff(M'); 
        end
        
        function [payoff, sum_probs] = incremental_payoff(obj, m, s)
            p = obj.params;
            r = p.risk_free_rate;
            K = p.dim;
            M = p.t_len;
            payoff = 0;
            U = rand(K, K);
            B = (log(p.barrier_level ./ s') - ...
                 (p.risk_free_rate - (p.sigma' .^ 2) / 2) * p.t_steps(m)) ...
                 ./ (p.sigma' * sqrt(p.t_steps(m)));
            [samples, probs] = sampleConeComplement(obj.proj_Rs, obj.proj_bs, obj.Qs, ...
                                                       B, U);                                  
            sum_probs = 0;
            for k = 1:K
                final_payoff = 0;
                next_payoff = 0;
                inc_prob = 0;
                end_probs = 0;
                sample = samples(k, :);
                prob = probs(k, :);
                s_new = s .* exp((p.risk_free_rate - (p.sigma.^2)  / 2) * ...
                                 p.t_steps(m) + sqrt(p.t_steps(m)) * ... 
                                 p.sigma .* obj.params.corr_chol * obj.T * sample');
                early_probs = prod(prob(1:end-1)) * (1 - prob(end));
                early_payoff = exp(-r * p.obs_times(m)) * early_probs * p.rebates(m);
                if m == M
                    end_probs =  prod(prob);
                    final_payoff = exp(- r * p.obs_times(m))* end_probs * p.survival_payoff(s_new);
                else
                    [inc_p, inc_prob] = obj.incremental_payoff(m + 1, s_new);
                    next_payoff = prod(prob) * inc_p;
                end
                payoff = payoff + early_payoff + next_payoff + final_payoff;
                sum_probs = sum_probs + prod(prob) * inc_prob + early_probs + end_probs;
            end
            
        end
        
        function [payoff, new_j] = pr_payoff(obj, m, s, V, j)
            %pathrecycled version of the incremental payoff
            p = obj.params;
            r = p.risk_free_rate;
            K = p.dim;
            M = p.t_len;
            payoff = 0;
            U = V(:, (j - 1) * K  + 1:j * K)';
            B = (log(p.barrier_level ./ s') - ...
                 (p.risk_free_rate - (p.sigma' .^ 2) / 2) * p.t_steps(m)) ...
                 ./ (p.sigma' * sqrt(p.t_steps(m)));
            [samples, probs] = sampleConeComplement(obj.proj_Rs, obj.proj_bs, obj.Qs, ...
                                                       B, U);                                  
            new_j = j + 1;
            for k = 1:K
                final_payoff = 0;
                next_payoff = 0;
                sample = samples(k, :);
                prob = probs(k, :);
                s_new = s .* exp((p.risk_free_rate - (p.sigma.^2)  / 2) * ...
                                 p.t_steps(m) + sqrt(p.t_steps(m)) * ... 
                                 p.sigma .* obj.params.corr_chol * obj.T * sample');
                early_probs = prod(prob(1:end-1)) * (1 - prob(end));
                early_payoff = exp(-r * p.obs_times(m)) * early_probs * p.rebates(m);
                if m == M
                    end_probs =  prod(prob);
                    final_payoff = exp(- r * p.obs_times(m))* end_probs * p.survival_payoff(s_new);
                else
                    [inc_p, new_j] = obj.incremental_payoff(m + 1, s_new, V, new_j);
                    next_payoff = prod(prob) * inc_p;
                end
                payoff = payoff + early_payoff + next_payoff + final_payoff;
            end
            
        end
        
        
        
        function payoff = samplePathPayoff(obj)
            %[paths, probs] = obj.generatePath();
            %payoff = obj.evaluatePaths(paths, probs);
            payoff = obj.incremental_payoff(1, obj.params.s_0);
        end
        
        function payoff = pr_samplePathPayoff(obj, U)
            %Generates payoff with U as random input variable.
            payoff = obj.pr_payoff(1, obj.params.s_0, U, 1);
        end
        
        function delta = pr_Delta(obj, dim, h, n_sims)
            s_pathGen = obj.shift("s_0", h, dim);
            deltas = zeros(n_sims,1);
            K = obj.params.dim;
            t_steps = sum(obj.params.dim .^ (0:obj.params.t_len));
            parfor n = 1:n_sims
                u = rand(K, t_steps);
                deltas(n) = s_pathGen.pr_samplePathPayoff(u) - obj.pr_samplePathPayoff(u);
            end
            delta = mean(deltas) / h;
        end
        
        function gamma = pr_Gamma(obj, dim, h, n_sims)
            forward = obj.shift("s_0", h, dim);
            backward = obj.shift("s_0", - h, dim);
            gammas = zeros(n_sims,1);
            K = obj.params.dim;
            t_steps = sum(obj.params.dim .^ (0:obj.params.t_len));
            parfor n = 1:n_sims
                u = rand(K, t_steps);
                gammas(n) = forward.pr_samplePathPayoff(u) - 2*obj.pr_samplePathPayoff(u) + backward.pr_samplePathPayoff(u);
            end
            gamma = mean(gammas) / h^2;
        end
        
    end
    
    methods(Static)
        function obj = test_instance()
            params = AutoCallParams.test_instance();
            obj = SmoothPathGenerator(params);
        end
        
        function [test, obj] = run_tests()
            try
                obj = SmoothPathGenerator.test_instance();
                p = obj.params;
                [paths, probs] = obj.generatePath();
                for t = 1:p.t_len
                    loc = (t - 1) * p.dim + 1: t * p.dim;
                    paths
                    assert(all(max(min(paths(:, loc), 2)) < p.barrier_level), ...
                        "TestSmoothPathGenerator:PathBarrierError", ...
                        "Minimum of underlying was not in path")
                    assert(all(0 <= min(probs(:, loc)), 'all') &&...
                            all(1 >= max(probs(:, loc)), 'all'), ...
                            "TestSmoothPathGenerator:ProbsError", ...
                            "Probabilites not bounded by 0, 1")
                end
                payoff = obj.evaluatePaths(paths, probs);
                assert(length(payoff) == 1, ...
                    "TestSmoothPathGenerator:PayoffDimensionError", ...
                    "Dimension of Payoff was greater than 1")
                assert(isnumeric(payoff), ...
                    "TestSmoothPathGenerator:PayoffTypeError", ...
                    "Type of Payoff is Non-numeric")
                test = 1;
            catch ME
                switch ME
                    case "TestSmoothPathGenerator:PathBarrierError"
                        paths
                    case "TestSmoothPathGenerator:ProbsError"
                        probs
                    case "TestSmoothPathGenerator:PayoffDimensionError"
                        payoff
                    case "TestSmoothPathGenerator:PayoffTypeError"
                        payoff
                end
                rethrow(ME)
            end
        end
    end
    
   
    
end
            
        