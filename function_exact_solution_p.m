function result = function_exact_solution_p(x,y)
%% 
global solution_case Ra
if solution_case == 1
    if x<=1/2
        result = Ra*12*x^2*exp(y);
    elseif x>1/2
        result = Ra*(16*x*y^3-exp(1)-2);
    end
elseif solution_case == 2
    result = (x*y)^3 - 5/80;
elseif solution_case == 3
    if x<=1/2
        result = 2*(y-1)*cos(x)^2;
    elseif x>1/2
        result = y*cos(y^2)+4*x-5/2;
    end
elseif solution_case == 4
    if x<=1/2
        result = 12*x^2*exp(y);
    elseif x>1/2
        result = 16*x*y^3-exp(1)-2;
    end
elseif solution_case == 5
    if x<=1/2
        result = 12*x^2*exp(y);
    elseif x>1/2
        result = 16*x*y^3-exp(1)-2;
    end
end


end