function result = function_exact_solution_u1(x,y)
%% 
global solution_case
if solution_case == 1
    result = exp(x*y);
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    result = sin(y^2+6*x);
elseif solution_case == 4
    result = x+x^2-2*x*y+x^3-3*x*y^2+x^2*y;
elseif solution_case == 5
    result = x+x^2-2*x*y+x^3-3*x*y^2+x^2*y;
end

end