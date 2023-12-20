function result = Gauss_quadrature_for_surface_integral_test_triangle_RT0(coefficient_function1, coefficient_function2, Gauss_weights_local, Gauss_nodes_local, vertices, test_basis_index, test_derivative_degree)
%% Use Gauss quadrature to compute a volume integral on a local triangular element T for a vector.
%We will use "FE" to replace "finite element" in the comments.
%The integrand of the volume integral must be in the following format:
%a coefficient function * a test FE basis function (or its derivatives).
%coefficient_function_name: the coefficient function of the integrand.
%Gauss_coefficient_local,Gauss_point_local:the Gauss coefficients and Gauss points on the triangular element T.
%vertices: the coordinates of all vertices of the triangular element T.
%test_basis_type:the type of the test FE basis function.
%test_basis_type=1:2D linear FE.  
%test_basis_type=2:2D Lagrange quadratic FE.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.
%%
Gpn = length(Gauss_weights_local); % Gpn: the number of the Gauss points of the Gauss quadrature we are using.
result = 0;
for i = 1:Gpn
      result = result + Gauss_weights_local(i) * dot( [feval(coefficient_function1, Gauss_nodes_local(i,1), Gauss_nodes_local(i,2)) ,feval(coefficient_function2, Gauss_nodes_local(i,1), Gauss_nodes_local(i,2))] , ...
          triangular_RT0_local_basis(Gauss_nodes_local(i,1), Gauss_nodes_local(i,2), vertices, test_basis_index, test_derivative_degree) );
end