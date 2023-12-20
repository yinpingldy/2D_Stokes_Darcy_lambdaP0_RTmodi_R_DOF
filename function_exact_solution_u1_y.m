function result = function_exact_solution_u1_y(x,y)
%% 
global solution_case
if solution_case == 1
    result = x*exp(x*y);
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    result = 2*y*cos(y^2 + 6*x);
elseif solution_case == 4
    result = x^2 - 6*x*y - 2*x;
elseif solution_case == 5
    result = x^2 - 6*x*y - 2*x;
end

end