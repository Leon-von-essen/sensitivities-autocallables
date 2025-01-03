function val = UNSAFE_norminv(p)
val = - sqrt(2) * erfcinv(2 * p);
%val = subs(val);
end
