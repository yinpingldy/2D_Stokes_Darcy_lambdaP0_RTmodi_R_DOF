function result = function_exact_solution_u2_x(x,y)
%% 
global solution_case
if solution_case == 1
    if x<=1/2
        result = exp(x + 2*y);
    elseif x>1/2
        result = 3;
    end
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    if x<=1/2
        result = -8*x*y*sin(4*x^2*y);
    elseif x>1/2
        result = 2*cos(2*x)*cos(3*y);
    end
elseif solution_case == 4
    if x<=1/2
        result = - 2*y - 6*x*y - y^2;
    elseif x>1/2
        result = 3;
    end
elseif solution_case == 5
    result = - 2*y - 6*x*y - y^2;
end
end