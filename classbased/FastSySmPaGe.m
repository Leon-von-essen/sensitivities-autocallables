classdef FastSySmPaGe < pathGenerator
    %Symbolic path generation class
    properties
        R
        T
        Qs
        proj_Rs
        proj_bs
        orig_term
        sym_payoff
        o_delta
        s_delta
        s_gamma
        payoff
        delta
        gamma
        lens
        fdpayoff
        fddelta
    end
    methods
        function obj = FastSySmPaGe(params)
            obj@pathGenerator(params);
            [obj.T, obj.R, obj.proj_Rs, obj.proj_bs, obj.Qs] = loadProjectedMatrices(params.corr_chol);          
            obj.orig_term = obj.built_integrals();
            syms s [params.dim, 1]
            x = symvar(subs(symvar(obj.orig_term), s, ones(params.dim, 1)));
            obj.sym_payoff = obj.orig_term.sub_pdf(x);
            obj.o_delta = obj.orig_term.diff(s(1)).clean();
            obj.s_delta = obj.o_delta.sub_pdf(x);
            obj.s_gamma = obj.o_delta.diff(s(1)).sub_pdf(x);
            u = symvar(obj.sym_payoff);
            obj.payoff = matlabFunction(subs(obj.sym_payoff, s, params.s_0), 'Vars', {u});
            obj.delta = matlabFunction(subs(obj.s_delta, s, params.s_0), 'Vars', {u});
            obj.gamma = matlabFunction(subs(obj.s_gamma, s, params.s_0), 'Vars', {u});
            obj.lens = length(u);
            obj.fdpayoff = matlabFunction(obj.sym_payoff, 'Vars', {u});
            obj.fddelta =  matlabFunction(obj.s_delta, 'Vars', {u});
        end
        
        function term = built_integrals(obj)
            K = obj.params.dim;
            syms s [K, 1]
            term = obj.get_int(K, 1, s, 0, 0);
        end
        
        function [term, next_ind] = get_int(o, k, subcone, s, m, ind,...
                                            ubounds, lbounds, ineqs, y)
            %GET_INT builts the integral structure recursively. Depending
            %on the subcone, the integral boundaries are set. Then a
            %SymTerm object is returned as term, with conditions based on
            %timestep and dimension
            K = o.params.dim;
            r = o.params.risk_free_rate;
            p = o.params;
            sigma = o.params.sigma;
            next_ind = ind + 1;           
            terms = {0};
            %set rebate payments
            if k == K && m == o.params.t_len
                terms(1) = {exp(- r * p.obs_times(p.t_len)) * p.survival_payoff(s)}; 
                term = SymTerm(terms);
            else
                if k == K
                    for d = 1:K
                        %locate transformation 
                        loc = (d - 1) * K + 1:K * d;
                        Q = o.Qs(loc, loc);
                        %get relevant section of integration variables
                        syms x [next_ind + K 1]
                        y = x(next_ind :next_ind + K - 1);
                        %get bounds
                        [ubounds, lbounds, ineqs] = o.get_bounds(d, s, m + 1, y);
                        %sample next s and iterate
                        Z = sqrt(p.t_steps(m + 1)) * p.sigma .* ...
                            (p.corr_chol * o.T * Q' * y);
                        s_new = s .* exp(Z + p.t_steps(m + 1) * (r - sigma .^ 2 / 2)); 
                        [newterm,  provis] = o.get_int(1, d, s_new, m + 1, next_ind, ubounds, lbounds, ineqs, y);
                        terms = [terms, {GaussIntegral(ubounds(1), lbounds(1), x(next_ind), ...
                                           newterm)}]; 
                        next_ind = provis;
                    end
                    term = SymTerm(terms);
                else
                    if k == K - 1
                        terms(1) = {(ineqs(K) * (1 - normcdf(ubounds(K))) + ...
                                     (1-ineqs(K)) *(1 - normcdf(lbounds(K)))) * ...
                                    exp(- r * p.obs_times(m)) * p.rebates(m)};
                    end
                    [newterm,  next_ind] = o.get_int(k + 1, subcone, s, m, next_ind, ubounds, lbounds, ineqs, y);
                    int_var = y(k + 1);
                    term = SymTerm([terms, {GaussIntegral(ubounds(k + 1), lbounds(k + 1), int_var, newterm)}]);
                end
            end
        end
    
        
        function [ubounds, lbounds, ineqs] = get_bounds(o, subcone, s, m, x)
           %GET_BOUNDS pulls the boundary conditions for the integrals
           p = o.params;
           r = p.risk_free_rate;
           sigma = p.sigma;
           b = p.barrier_level;
           [~, pR, pB, T] = o.loc_vars(subcone);
           B = (log(b) - log(s) - (r - sigma .^2 / 2) * p.t_steps(m)) ./ sigma / sqrt(p.t_steps(m)) ;
           ubounds = (diag(pR) < 0) .* 1 ./ diag(pR) .* (pB * B - T * x);
           lbounds = (diag(pR) > 0) .* 1 ./ diag(pR) .* (pB * B - T * x);
           dia = diag(pR);
            for i = 1:length(B)
                if dia(i) > 0
                    ubounds(i) = inf;
                else 
                    lbounds(i) = -inf;
                end
            end
           ineqs = diag(pR) < 0;
        end
        
        function [loc, pR, pB, T] = loc_vars(o, subcone)
            %LOC_VARS returns the respective parts of the equations.
            K = o.params.dim;
            loc = (subcone-1) * K + 1:K*subcone;
            pR = o.proj_Rs(loc, loc);
            T = tril(pR, -1);
            pB = o.proj_bs(loc, loc);
        end
        
        function payoff = samplePathPayoff(obj)
            u = rand(obj.lens, 1);
            payoff = obj.payoff(u');            
        end
        
        function [delta, stdev] = Delta(obj, dim,  h, n_sims)
            d = zeros(n_sims, 1);
            parfor sim = 1:n_sims
                d(sim) = obj.delta(rand(obj.lens, 1)');
            end
            delta = mean(d);
            stdev = std(d);
        end
        
        function [xGamma, stdev] = Gamma(obj, dim, h, n_sims)
            Gamma = zeros(n_sims, 1);
            parfor sim = 1:n_sims
                Gamma(sim) = obj.gamma(rand(obj.lens, 1)');
            end
            xGamma = mean(Gamma);
            stdev = std(Gamma);
        end
        
        function test = testGreeks(obj)
            test = 0;
            payoff = obj.samplePathPayoff()
            delta = obj.Delta(1,0.1, 10)
            %vega = obj.Vega(1,0.1, 10);
            %rho = obj.Rho(0.1, 10);
            %theta = obj.Theta(0.1, 10);
            %Gamma = obj.Gamma(1,0.1,10);
            %xGamma = obj.CrossGamma(1, 2, 0.1, 10);
            test = 1;
        end
        
        function [m, s] = fd_Delta(o, dim, h, nsims)
            s = o.params.s_0;
            K = o.params.dim;
            len =o.lens;
            s1 = s + h * [zeros(dim - 1, 1);1;zeros(K - dim, 1)];
            vals = zeros(1, nsims);
            parfor i = 1:nsims
                u = [rand(1, len - K), s1'];
                v = [rand(1, len - K), s1'];
                vals(i) = o.fdpayoff(u) - o.payoff(v);
            end
            m = mean(vals) / h;
            s = std(vals) / sqrt(h);
        end
        
        function [m, s] = fd_Gamma(o, dim, h, nsims)
            s = o.params.s_0;
            K = o.params.dim;
            len = o.lens;
            s1 = s + h * [zeros(dim - 1, 1);1;zeros(K - dim, 1)];
            vals = zeros(1, nsims);
            parfor i = 1:nsims
                u = [rand(1, len - K), s1'];
                v = [rand(1, len - K), s1'];
                vals(i) = o.fddelta(u) - o.delta(v);
            end
            m = mean(vals) / h;
            s = std(vals) /  sqrt(h);
        end
        function gamma = pr_Gamma(o, dim, h, n_sims)
            s = o.params.s_0;
            K = o.params.dim;
            len = o.lens;
            s1 = s + h * [zeros(dim - 1, 1);1;zeros(K - dim, 1)];
            vals = zeros(1, n_sims);
            parfor i = 1:n_sims
                u = [rand(1, len - K), s1'];
                vals(i) = o.fddelta(u) - o.delta(u);
            end
            gamma = mean(vals) / h;
        end
    end
    
    methods(Static) 
        function [obj, test] = test_instance()
            params = AutoCallParams.test_instance();
            obj = FastSySmPaGe(params);
            test = obj.testGreeks();
        end
        
    end
end

