function [result] = function_F1(x, y)
%% the right hand item
global mu solution_case Ra
if solution_case == 1
    if x<=1/2
        result = -2*mu*y^2*exp(x*y) - mu*x^2*exp(x*y) - 2*mu*exp(x+2*y) + Ra*24*x*exp(y);
    elseif x>1/2
        result = mu*exp(x*y)+Ra*16*y^3;
    end
elseif solution_case == 2
    result = 3*x^2*y^3;
elseif solution_case == 3
    if x<=1/2
        result = 2*mu*36*sin(y^2 + 6*x) - mu*(2*cos(y^2 + 6*x) - 4*y^2*sin(y^2 + 6*x)) + mu*(8*x*sin(4*x^2*y) + 32*x^3*y*cos(4*x^2*y)) -2*cos(x)*sin(x)*(2*y - 2);
    elseif x>1/2
        result = mu*(sin(y^2+6*x)) + 4;
    end
elseif solution_case == 4
    if x<=1/2
        result = -2*mu*( 6*x + 2*y + 2)- mu*(-6*x) -mu*(- 6*x - 2*y - 2) + 24*x*exp(y);
    elseif x>1/2
        result = mu*(x+x^2-2*x*y+x^3-3*x*y^2+x^2*y)+16*y^3;
    end
elseif solution_case == 5
    if x<=1/2
        result = -2*mu*( 6*x + 2*y + 2)- mu*(-6*x) -mu*(- 6*x - 2*y - 2) + 24*x*exp(y);
    elseif x>1/2
        result = mu*(x+x^2-2*x*y+x^3-3*x*y^2+x^2*y)+16*y^3;
    end
end

end