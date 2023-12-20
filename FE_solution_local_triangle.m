function [result] = FE_solution_local_triangle(x, y, uh_local, vertices, basis_type, derivative_degree_x, derivative_degree_y)
%% Evaluate a finite element solution at (x,y) which is in a triangular element T.(解函数、插值函数、投影函数)
% uh_local: the values of numerical solution at all the nodes of FE basis functions in the triangular element T.
% vertices: the coordinates of all vertices of the triangular element T.
% basis_type: the type of the FE.
% basis_type=1: 2D linear FE.
% basis_type=2: 2D Lagrange quadratic FE.
% basis_index: the index of basis function to specify which basis function we want to use.
% derivative_degree_x: the derivative degree of the FE basis function with respect to x.
% derivative_degree_y: the derivative degree of the FE basis function with respect to y.
%% 
result = 0;
number_of_local_basis = length(uh_local);
for i = 1:number_of_local_basis
    result = result + uh_local(i) * triangular_local_basis(x, y, vertices, basis_type, i ,derivative_degree_x, derivative_degree_y);
end

end