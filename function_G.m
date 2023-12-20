function result = function_G(x,y)
%% the boundary item
global solution_case
if solution_case == 1
    if x<=1/2
        result = y*exp(x*y) + 2*exp(x+2*y);
    else
        result = y*exp(x*y) + 2*y;
    end
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    if x<=1/2
        result = 6*cos(y^2 + 6*x) - 4*x^2*sin(4*x^2*y);
    else
        result = 6*cos(y^2 + 6*x) - 3*sin(2*x)*sin(3*y);
    end
elseif solution_case == 4
    if x<=1/2
        result = 0;
    else
        result = 3*x^2 + 2*x*y + 2*x - 3*y^2 - 2*y + 1 + 2*y;
    end
elseif solution_case == 5
    result = 0;
end

end