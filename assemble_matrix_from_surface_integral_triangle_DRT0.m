function result = assemble_matrix_from_surface_integral_triangle_DRT0(coefficient_function, P_partition, T_partition, N_basis_trial, N_basis_test, Darcy_t, triangle_edge,...
    number_of_trial_local_basis, number_of_test_local_basis, number_of_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    derivative_degree_trial, derivative_degree_test)
%% Assemble a matrix from a surface inegral on the whole domain with a triangular mesh.
%The integrand of the volume integral must be in the following format:
%a coefficient function * a trial FE basis function(or its derivatives) * a test FE basis function (or its derivatives).
%coefficient_function_name: the coefficient function of the integrand.
%M_partition: store the coordinates of all the grid points of the partition,not FE.
%T_partition: store the global indices of the grid points of every element for the partition,not FE.
%T_basis: store the global indices of the nodes of every element for FE,not the partition.
%T_basis_trial: T_basis for the trial basis function.
%T_basis_test: T_basis for the test basis function.
%The explanation for M_partition,T_partition,T_basis is in generate_M_T_triangular.m.
%number_of_trial_local_basis: the number of local FE basis functions for the trial function in a local element.
%number_of_test_local_basis: the number of local FE basis functions for the test function in a local element.
%number_of_element: the number of the local triangluar elements of the partition.
%N1_basis,N2_basis:The N1 and N2 for the basis,not the partition.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangular element.

%trial_basis_type:the type of the trial FE basis function.
%basis_type=0:2D constant FE.
%basis_type=1:2D Lagrange linear FE.
%basis_type=2:2D Lagrange quadratic FE.
%basis_type=10:2D Crouzeix-Raviart FE.
%trial_basis_index: the index of trial FE basis function to specify which trial FE basis function we want to use.
%trial_derivative_degree_x:the derivative degree of the trial FE basis function with respect to x.
%trial_derivative_degree_y:the derivative degree of the trial FE basis function with respect to y.
%test_basis_type:the type of the test FE basis function.
%basis_type=0:2D constant FE.
%basis_type=1:2D Lagrange linear FE.
%basis_type=2:2D Lagrange quadratic FE.
%basis_type=10:2D Crouzeix-Raviart FE.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.

%vertices: the coordinates of all vertices of a triangular element.
%Gauss_coefficient_local_triangle,Gauss_point_local_triangle: the Gauss coefficients and Gauss points on the local triangular element.
%%
result = sparse(N_basis_test,N_basis_trial);

% Go through all elements.
% On each element, compute the volume integrals for all possible combinations of trial and test FE basis functions.
% Assemble the values of those volume integrals into the matrix.
for n=1:number_of_elements
    vertices = P_partition(:,T_partition(:,n));
    [Gauss_weights_local_triangle,Gauss_nodes_local_triangle] = generate_Gauss_local_triangle(Gauss_weights_reference_triangle,Gauss_nodes_reference_triangle,vertices);
    element_now = Darcy_t(n);
    
    for alpha = 1:number_of_trial_local_basis
       for beta = 1:number_of_test_local_basis      
          int_value = Gauss_quadrature_for_surface_integral_trial_test_triangle_RT0(coefficient_function,Gauss_weights_local_triangle, Gauss_nodes_local_triangle, vertices,...
              alpha, derivative_degree_trial, beta, derivative_degree_test);
          i = triangle_edge(beta,element_now);
          j = triangle_edge(alpha,element_now);
          result(i,j) = result(i,j) + int_value;
       end
    end

end