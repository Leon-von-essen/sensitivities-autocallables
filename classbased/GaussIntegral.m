classdef GaussIntegral
    %GaussIntegral Integral class with implied normal density
    %props:
    %ubound: upper bound, sym or numeric
    %lbound: lower bound, sym or numeric
    %int_var: integration variable, sym
    %term: SymTerm object
    %depth: integer
    
    properties
        ubound
        lbound
        int_var
        term
        
    end
    
    methods
        function obj = GaussIntegral(ubound, lbound, int_var, term)
            %SimpleIntegral Construct an instance of this class
            %   Detailed explanation goes here
            obj.ubound = ubound;
            obj.lbound = lbound;
            obj.int_var = int_var;
            obj.term = term;
        end
        
        function alt = scale(obj, scale)
            %Scales the terms within the integral
            alt = obj;
            alt.term = obj.term.scale(scale);
        end
        
        function alt = sub(obj, var, sub_var)
            %Substitutes the variable into all components of the object;
            %var
            %obj.int_var
            alt = obj;
            %handle edgecase of infinite subs (matlab unable to manage
            %subs if inf is involved)
            if sub_var == inf || sub_var == -inf
                if obj.ubound ~= inf && obj.ubound ~= -inf
                    uvars = symvar(obj.ubound);
                    vars = uvars(uvars~=var);
                    temp = limit(subs(obj.ubound, vars, 0.5 * ones(1, length(vars))), var, sub_var);
                    if temp == inf || temp == -inf
                        alt.ubound = temp;
                    else 
                        alt.ubound = subs(obj.ubound, var, sub_var);
                    end
                end
                if obj.lbound ~= inf && obj.lbound ~= -inf
                    lvars = symvar(obj.lbound);
                    vars = lvars(lvars~=var);
                    temp = limit(subs(obj.lbound, vars, 0.5 * ones(1, length(vars))), var, sub_var);
                    if temp == inf || temp == -inf
                        alt.lbound = temp;
                    else 
                        alt.lbound = subs(obj.lbound, var, sub_var);
                    end
                end
            else
                %standard case
                alt.ubound = subs(obj.ubound, var, sub_var);
                alt.lbound = subs(obj.lbound, var, sub_var);
            end
            alt.term = obj.term.sub(var, sub_var);
            if isnan(alt.ubound) || isnan(alt.lbound)
                alt.ubound
                alt.lbound
            end
        end
        
        function term = sub_pdf(obj, x)
            %pdf_var = obj.int_var
            up = normcdf(obj.ubound);
            low = normcdf(obj.lbound);
            p = up - low;
            i = find(x == obj.int_var);
            syms U [i, 1]
            subbed_term = obj.term.sub(obj.int_var, UNSAFE_norminv(low + p * U(end)));
            scaled_term = subbed_term.scale(p);
            term = scaled_term.sub_pdf(x);
            if isnan(term) 
                term
            end
        end
        
        function boolean = isempty(obj)
            boolean = 1;
            if  obj.term.term{1} == 0 && length(obj.term) == 1
                boolean = 1;
            else    
                for i = 2:length(obj.term)
                    if isa(obj.term.term{i}, 'GaussIntegral')
                        boolean = boolean * obj.term.term{i}.isempty();
                    elseif obj.term.term{i} ~= 0 
                        boolean = 0;
                    end
                end
            end                        
        end
        
        function alt = clean(obj)
            %Cleans up the term. 
            alt = obj;
            alt.term = obj.term.clean();    
        end
                
        function bool = isequiv(obj, other)
            %Checks if other has same integral bounds as obj
            if isa(other, 'GaussIntegral') && (obj.ubound == other.ubound) && (obj.lbound == other.lbound)
                bool = 1;
            else
                bool = 0;
            end
        end
        
         function bool = eq(o, alt)
             %overrides the == operator
            bool = o.isequiv(alt) * (o.term == alt.term);
         end
         
        function terms_o = leibnitz(obj, param)
            %Calculates the symbolic differentiation by employing the
            %leibnitz rule
            terms = {};
            if  ~isAlways(isinf(obj.ubound))
                ubound1 = diff(obj.ubound, param);
                terms = [terms, obj.term.sub(obj.int_var, obj.ubound).scale(normpdf(obj.ubound) * ...
                         ubound1).term];
            end
            if ~isAlways(isinf(obj.lbound))
                lbound1 = diff(obj.lbound, param);
                terms = [terms, obj.term.sub(obj.int_var, obj.lbound).scale(- normpdf(obj.lbound) * ...
                        lbound1).term];
            end
            alt = obj;
            alt.term = obj.term.diff(param);
            terms(end + 1) = {alt};
            terms_o = SymTerm(terms).clean();
        end
        
        function arr = symvar(o)
            %get the symbolic variables for the array.
            arr = symvar([o.int_var, symvar(o.term)]);
        end
        
        function obj = subsindex(o, i)
            obj = o.term(i);
        end
        
        function ret = disp(o)
            upper_bound = o.ubound
            lower_bound = o.lbound
            integration_variable = o.int_var
            term = o.term
        end
    end
                  
        
    methods(Static)
        function t_instance = t_instance()
            syms x y z
            ubound = 0;
            lbound = - inf;
            int_var = y;
            inner_term = {y};
            ub2 = inf;
            lb2 = 0;
            iv2 = z;
            it2 = SymTerm({x * z});
            inner_term = SymTerm([inner_term, {GaussIntegral(ub2, lb2, iv2, it2)}]);
            t_instance = GaussIntegral(ubound, lbound, int_var, inner_term);
        end
        
        function [t, t_subbed, test] = t_subs()
            syms x z
            ubound = x;
            lbound = -inf;
            int_var = z;
            term = SymTerm({x * z});
            t = GaussIntegral(ubound, lbound, int_var, term);
            t_subbed = t.sub(x, 0);
            alt = GaussIntegral(0, -inf, z, 0);
            assert(alt.isequiv(t_subbed), 'SubError:integralbounds not correctly subbed');
            assert(alt.term == t_subbed.term, 'SubError: Integrals not corretly subbed');
        end
        
        function [ts, terms, tests] = t_sub_pdfs()
            syms x z
            ubound = 0;
            lbound = -inf;
            int_var = z;
            term = SymTerm({1});
            t = GaussIntegral(ubound, lbound, int_var, term);
            t_sp = t.sub_pdf(z);
            assert(abs(t_sp(1) - 0.5) < 0, 'Suberror:test_value_fail')
        end
            
            
        
        function [tins, terms, test] = t_leibnitz()
            syms x
            tins = GaussIntegral.t_instance();
            terms = tins.leibnitz(x);
            try
            assert(isa(terms,'SymTerm'), 'TypeError:Postleibnitz')
            assert(isa(terms.term{1},'sym'), 'TypeError:firstObj')
            assert(isa(terms.term{2},'GaussIntegral'), 'TypeError:Deriv2')
            test = 1;
            catch ME
                terms
                tins
                rethrow(ME)
            end
        end
        
        function [tins, tsubs, test] = t_subs_1()
            tins = GaussIntegral.t_instance();
            syms U [1 2]
            tsubs = tins.sub_pdf(U);
            assert(isa(tsubs,'sym'), 'TypeError:subbed')
            test = 1;
        end
        
    end
end

