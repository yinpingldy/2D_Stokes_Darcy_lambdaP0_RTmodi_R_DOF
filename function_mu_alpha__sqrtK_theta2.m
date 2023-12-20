function result = function_mu_alpha__sqrtK_theta2(x, y)

global solution_case mu alpha
k_j = 1;
if solution_case == 1
    result =  mu*alpha*sqrt(k_j)*(exp(x+2*y) + x*exp(x*y) + exp(x+2*y));
elseif solution_case == 2
    result = 0;
elseif solution_case == 3
    result = mu*alpha*sqrt(k_j)*(cos(4*x^2*y) + 2*y*cos(y^2+6*x)-8*x*y*sin(4*x^2*y));
elseif solution_case == 4
    result =  mu*alpha*sqrt(k_j)*( -y-2*x*y+y^2-3*x^2*y+y^3-x*y^2   +   x^2-6*x*y-2*x    - 2*y - 6*x*y - y^2 );
elseif solution_case == 5
    result =  mu*alpha*sqrt(k_j)*( -y-2*x*y+y^2-3*x^2*y+y^3-x*y^2   +   x^2-6*x*y-2*x    - 2*y - 6*x*y - y^2 );
end
end

