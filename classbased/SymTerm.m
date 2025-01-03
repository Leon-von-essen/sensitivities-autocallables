classdef SymTerm
    %SymTerm Class to combine GaussIntegrals and syms/numerics
    properties
        term
    end
    
    methods
        function obj = SymTerm(terms)
            %Constructor, takes term array
            if isa(terms, 'SymTerm')
                obj = terms;
            else
                obj.term = terms;
            end
        end
        
        %function ret = disp(o)
         %   for i = 1:length(o.term)
          %      sub_term = o(i)
           % end
        %end
        
        function obj = subsindex(o, i)
            obj = o.term{i};
        end
        
         function bool = eq(o, alt)
             %overrides the == operator
            bool = 1;
            if isa(alt, 'SymTerm')
                for i = 1:length(o)
                    bool = bool * (o(i) == alt(i));
                end
            elseif isa(alt, 'sym') || isnumeric(alt)
                if o(1) ~= alt || o.length() ~= 1
                    bool = 0;
                end
            else 
                bool = 0;
            end
        end
        
        function alt_term = sub(o, var, sub_var)
            alt_terms = {[0]};
            for i = 1:length(o.term)
                if isa(o.term{i}, 'GaussIntegral')
                   alt_terms{end+1} = o.term{i}.sub(var, sub_var);
                else
                    alt_terms{1} = alt_terms{1} + subs(o.term{i}, var, sub_var);
                end
            end
            alt_term = SymTerm(alt_terms);
        end
        
        function val = sub_pdf(o, x)
            val = 0;
            for i = 1:length(o.term)
                xpr = o.term;
                if isa(xpr{i}, 'GaussIntegral')
                    add = xpr{i}.sub_pdf(x);
                    %syms U [k + 1 1]
                    %v = U;
                else 
                    add = xpr{i};
                end
                val = val + add;
            end
        end
        
        function alt = scale(o, factor)
            alt_terms = {0};
            for i = 1:length(o.term)
                if isa(o.term{i}, 'GaussIntegral')
                   alt_terms{i} = o.term{i}.scale(factor);
                else
                    alt_terms{i} = factor * o.term{i};
                end
            end
            alt = SymTerm(alt_terms);
        end
        
        function alt = diff(o, param)
            alt = o;
            for i = 1:length(o.term)
                if isa(o.term{i}, 'GaussIntegral')
                   alt.term{i} = o.term{i}.leibnitz(param);
                else
                   alt.term{i} = diff(o.term{i}, param);
                end
            end
            alt = alt.clean();
        end
        
        function alt = plus(o, other)
            newterm = [o.term, other.term];
            alt = SymArray(newterm).clean();
        end
        
        function len = length(o)
            len = length(o.term);
        end
        
        function arr = symvar(o)
            arr = [];
            for i = 1:o.length()
                arr = [arr, symvar(o.term{i})];
            end
            arr = symvar(arr);
        end
        
        function term = clean(o)
            alt_term = {[0]};
            try
            for t = 1:length(o.term)
                xpr = o.term{t};
                if isa(xpr, 'sym') || isnumeric(xpr) 
                    alt_term{1} = alt_term{1} + xpr;
                elseif isa(xpr, 'SymTerm')
                    new_term = xpr.clean().term;
                    alt_term = [alt_term, new_term];
                else
                    new_term = xpr.clean();
                    alt_term{end + 1} = new_term;
                end
            end
            term = SymTerm(alt_term);
            catch ME
            o
            rethrow(ME);
            end
        end
        
    end
end
                
            
        