classdef pathGenerator
    properties
        params
    end
    methods
        function obj = pathGenerator(params)
            obj.params = params;
        end
        
        function payoff = samplePathPayoff(obj)
            payoff = 1;
        end
        
        function time = timeSamplePath(obj)
            run = @() obj.samplePathPayoff();
            time = timeit(run);
        end 
        
        function [pv, stdev] = MCEstimate(obj, n_sims)
            p = zeros(n_sims, 1);
            parfor n = 1:n_sims
                p(n) = obj.samplePathPayoff();
            end
            pv = mean(p);
            stdev = std(p);
        end
        
        function [avg, stand, time] = statsMC(obj, n_sims, n_runs)
            res = zeros(n_runs,1);
            for run = 1:n_runs
                res(run) = obj.MCEstimate(n_sims);
            end
            avg = mean(res);
            stand = std(res);
            time = obj.timeSamplePath;
        end
        
        function delta = Delta(obj, dim, h, n_sims)
            s_pathGen = obj.shift("s_0", h, dim);
            delta = (s_pathGen.MCEstimate(n_sims) - obj.MCEstimate(n_sims)) / h;
        end
        
        function delta = pr_Delta(obj, dim, h, n_sims)
            s_pathGen = obj.shift("s_0", h, dim);
            deltas = zeros(n_sims,1);
            parfor n = 1:n_sims
                u = rand(obj.params.dim, obj.params.t_len);
                deltas(n) = s_pathGen.pr_samplePathPayoff(u) - obj.pr_samplePathPayoff(u);
            end
            delta = mean(deltas) / h;
        end
        
        function gamma = Gamma(obj, dim, h, n_sims)
            forward = obj.shift("s_0", h, dim);
            backward = obj.shift("s_0", - h, dim);
            gamma = (forward.MCEstimate(n_sims) - 2 * obj.MCEstimate(n_sims) + ...
                    backward.MCEstimate(n_sims)) / h^ 2;
        end
        
        function gamma = pr_Gamma(obj, dim, h, n_sims)
            forward = obj.shift("s_0", h, dim);
            backward = obj.shift("s_0", - h, dim);
            gammas = zeros(n_sims,1);
            parfor n = 1:n_sims
                u = rand(obj.params.dim, obj.params.t_len);
                gammas(n) = forward.pr_samplePathPayoff(u) - 2*obj.pr_samplePathPayoff(u) + backward.pr_samplePathPayoff(u);
            end
            gamma = mean(gammas) / h^2;
        end
        
        function alt = shift(obj, shift_var, h, dim)
            p = obj.params.shift(shift_var, h, dim);
            alt = obj;
            alt.params = p;
        end
        
    end
    
end
            
        