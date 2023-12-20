function [result] = function_F2(x, y)
%% the right hand item
global mu solution_case Ra
if solution_case == 1
    if x<=1/2
        result = -9*mu*exp(x+2*y) - mu*(exp(x*y)+x*y*exp(x*y)) + Ra*12*x^2*exp(y);
    else
        result = mu*(3*x+y^2)+Ra*16*3*x*y^2;
    end
elseif solution_case == 2
    result = 3*y^2*x^3;
elseif solution_case == 3
    if x<=1/2
        result = 2*mu*16*x^4*cos(4*x^2*y) + mu*(8*y*sin(4*x^2*y) + 64*x^2*y^2*cos(4*x^2*y)) +mu*12*y*sin(y^2 + 6*x) + 2*cos(x)^2;
    else
        result = mu*sin(2*x)*cos(3*y) + cos(y^2) - 2*y^2*sin(y^2);
    end
elseif solution_case == 4
    if x<=1/2
        result = -2*mu*(6*y - 2*x + 2) -mu*(-6*y) -mu*(2*x - 6*y - 2) + 12*x^2*exp(y);
    else
        result = mu*(3*x+y^2)+16*3*x*y^2;
    end
elseif solution_case == 5
    if x<=1/2
        result = -2*mu*(6*y - 2*x + 2) -mu*(-6*y) -mu*(2*x - 6*y - 2) + 12*x^2*exp(y);
    else
        result = mu*(-y-2*x*y+y^2-3*x^2*y+y^3-x*y^2 + 1)+16*3*x*y^2;
    end
end

end