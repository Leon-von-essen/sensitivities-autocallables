classdef NaivePathGenerator < pathGenerator
    methods
        function payoff = samplePathPayoff(obj)
            p = obj.params;
            payoff = 0;
            N = length(p.t_steps);
            K = length(p.corr_chol);
            s = p.s_0;
            breaker = 0;
            for n = 1:N
                Z = sqrt(p.t_steps(n)) .* (p.sigma .* (p.corr_chol * randn(K,1)));
                s = s .* exp(Z + p.t_steps(n) .* (p.risk_free_rate - (p.sigma .^ 2) / 2));
                if min(s) > p.barrier_level
                    payoff = exp(- p.risk_free_rate * p.obs_times(n)) * p.rebates(n);
                    breaker = 1;
                    break
                end
            end
            if breaker < 1
                payoff = exp(- p.risk_free_rate * p.obs_times(end)) * p.survival_payoff(s);
            end
        end
        
        function payoff = pr_samplePathPayoff(obj, U)
            p = obj.params;
            payoff = 0;
            N = length(p.t_steps);
            K = length(p.corr_chol);
            s = p.s_0;
            breaker = 0;
            for n = 1:N
                Z = sqrt(p.t_steps(n)) .* (p.sigma .* (p.corr_chol * U(:, n)));
                s = s .* exp(Z + p.t_steps(n) .* (p.risk_free_rate - (p.sigma .^ 2) / 2));
                if min(s) > p.barrier_level
                    payoff = exp(- p.risk_free_rate * p.obs_times(n)) * p.rebates(n);
                    breaker = 1;
                    break
                end
            end
            if breaker < 1
                payoff = exp(- p.risk_free_rate * p.obs_times(end)) * p.survival_payoff(s);
            end
        end
        
    end
end
            
        