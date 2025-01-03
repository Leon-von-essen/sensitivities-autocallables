function prob = probNormSlice(lbound, ubound)
%%Input:
%lbound: MxN Array
%ubound: MxN Array
%output: prob: MxN Array
%Gives the comulative probabiltiy of a slice of the Normal distribution
%depending on the upper and lower bounds of the truncation
try
    assert(all(ubound >= lbound, 'all'), 'InputError:probNormSlice:boundaries',...
       'boundaries do not respect >= inequality')
   assert(~isnan(ubound) && ~isnan(lbound), 'InputError:probNormSlice:nanboundaries',...
       'boundaries are nan')
catch ME
    if ME == 'InputError:probNormSlice:boundaries'
            ubound
            lbound
    elseif ME == 'InputError:probNormSlice:nanboundaries'
        ubound
        lbound
    end
    rethrow(ME)
end
    
norm_ubound = normcdf(ubound);
norm_lbound = normcdf(lbound);
prob = norm_ubound - norm_lbound;

end