function result = function_neg_theta1(x,y)
%%
global solution_case mu Ra
if solution_case == 1
    result = -Ra*12*x^2*exp(y) + Ra*(16*x*y^(3)-exp(1)-2) + 2*mu*y*exp(x*y);
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    result = -(2*(y-1)*cos(x)^2 - (y*cos(y^2)+4*x-5/2) - 2*mu*(6*cos(y^2 + 6*x)));
elseif solution_case == 4
    result = -12*x^2*exp(y) + 16*x*y^(3)-exp(1)-2 + 2*mu*(3*x^2 + 2*x*y + 2*x - 3*y^2 - 2*y + 1);
elseif solution_case == 5
    result = -12*x^2*exp(y) + 16*x*y^(3)-exp(1)-2 + 2*mu*(3*x^2 + 2*x*y + 2*x - 3*y^2 - 2*y + 1);
end
end

