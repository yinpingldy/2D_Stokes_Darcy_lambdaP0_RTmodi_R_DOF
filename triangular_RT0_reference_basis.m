function [q_RT] = triangular_RT0_reference_basis(x,y,basis_index,derivative_degree)
%% the reference FE basis function on triangle ABC where A=(0,0), B=(1,0) and C=(0,1).
% x,y: the coordinates of the point where we want to evaluate the reference FE basis function.
% basis_type: the type of the FE.
%basis_type=0: 2D constant FE.
% basis_type=1: 2D linear FE.
% basis_type=2: 2D Lagrange quadratic FE.
% basis_type=10: 2D Crouzeix-Raviart FE.
% basis_index: the index of FE basis function to specify which FE basis function we want to use.
% derivative_degree_x: the derivative degree of the FE basis function with respect to x.
% derivative_degree_y: the derivative degree of the FE basis function with respect to y.
%%
% A = [0 -1 0;sqrt(2)/2 sqrt(2)/2 sqrt(2)/2;-1 0 0];
% q1 = A\[1;0;0];
% q2 = A\[0;1;0];
% q3 = A\[0;0;1];
q_RT = zeros(1,2);
if derivative_degree == 0
    if basis_index == 1
        q_RT(1) = sqrt(2)*x;
        q_RT(2) = sqrt(2)*y;
    elseif basis_index == 2
        q_RT(1) = -1 + x;
        q_RT(2) = y;
    elseif basis_index == 3
        q_RT(1) = x;
        q_RT(2) = -1 + y;
    end
elseif derivative_degree == 1 % 用不到，直接定义了local basis
    if basis_index == 1
        q_RT(1) = sqrt(2)*2;
    elseif basis_index == 2
        q_RT(1) = 2;
    elseif basis_index == 3
        q_RT(1) = 2;
    end
end