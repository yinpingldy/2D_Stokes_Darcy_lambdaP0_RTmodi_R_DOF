function result = function_exact_solution_u2_y(x,y)
%% 
global solution_case
if solution_case == 1
    if x<=1/2
        result = 2*exp(x+2*y);
    elseif x>1/2
        result = 2*y;
    end
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    if x<=1/2
        result = -4*x^2*sin(4*x^2*y);
    elseif x>1/2
        result = -3*sin(2*x)*sin(3*y);
    end
elseif solution_case == 4
    if x<=1/2
        result = - 3*x^2 - 2*x*y - 2*x + 3*y^2 + 2*y - 1;
    elseif x>1/2
        result = 2*y;
    end
elseif solution_case == 5
    result = - 3*x^2 - 2*x*y - 2*x + 3*y^2 + 2*y - 1;
end
end