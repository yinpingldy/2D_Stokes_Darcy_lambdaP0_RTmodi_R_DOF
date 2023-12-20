function result = function_exact_solution_u2(x,y)
%% 
global solution_case
if solution_case == 1
    if x<=1/2
        result = exp(x+2*y);
    elseif x>1/2
        result = 3*x + y^2;
    end
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    if x<=1/2
        result = cos(4*x^2*y);
    elseif x>1/2
        result = sin(2*x)*cos(3*y);
    end
elseif solution_case == 4
    if x<=1/2
        result = -y-2*x*y+y^2-3*x^2*y+y^3-x*y^2;
    elseif x>1/2
        result = 3*x + y^2;
    end
elseif solution_case == 5
    if x<=1/2
        result = -y-2*x*y+y^2-3*x^2*y+y^3-x*y^2;
    elseif x>1/2
        result = -y-2*x*y+y^2-3*x^2*y+y^3-x*y^2 + 1;
    end
end
end