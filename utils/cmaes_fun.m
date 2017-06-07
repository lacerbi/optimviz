function y = cmaes_fun(x)
%CMAESFUN Auxiliary optimization function for CMA-ES.

% We need this since CMA-ES has an annoying input format.

persistent fun;

if isa(x,'function_handle')
    fun = x;
    y = [];
else
    y = fun(x');
end

end