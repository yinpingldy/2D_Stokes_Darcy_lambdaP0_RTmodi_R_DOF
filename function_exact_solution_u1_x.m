function result = function_exact_solution_u1_x(x,y)
%% 
global solution_case
if solution_case == 1
    result = y*exp(x*y);
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    result = 6*cos(y^2 + 6*x);
elseif solution_case == 4
    result = 3*x^2 + 2*x*y + 2*x - 3*y^2 - 2*y + 1;
elseif solution_case == 5
    result = 3*x^2 + 2*x*y + 2*x - 3*y^2 - 2*y + 1;
end

end