function result = assemble_vector_from_surface_integral_triangle_p(coefficient_function, P_partition, T_partition, N_basis_test, T_basis_p,...
    number_of_test_local_basis, number_of_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_test, derivative_degree_x_test, derivative_degree_y_test)
%% Assemble a vector from a volume inegral on the whole domain with a triangular mesh.
%The integrand of the volume integral must be in the following format:
%a coefficient function * a test FE basis function (or its derivatives).
%coefficient_function_name: the coefficient function of the integrand.
%M_partition: store the coordinates of all the grid points of the partition,not FE.
%T_partition: store the global indices of the grid points of every element for the partition,not FE.
%T_basis: store the global indices of the nodes of every element for FE,not the partition.
%T_basis_test: T_basis for the test basis function.
%The explanation for M_partition,T_partition,T_basis is in generate_M_T_triangular.m
%h_partition: the step size of the partition
%number_of_test_local_basis: the number of local FE basis functions for the test function in a local element.
%number_of_element: the number of the local triangluar elements of the partition.
%N1_basis,N2_basis:The N1 and N2 for the basis,not the partition.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangular element.
%test_basis_type:the type of the test FE basis function.
%test_basis_type=1:2D linear FE.  
%test_basis_type=2:2D Lagrange quadratic FE.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.

%vertices: the coordinates of all vertices of a triangular element.
%Gauss_coefficient_local_triangle,Gauss_point_local_triangle: the Gauss coefficients and Gauss points on the local triangular element.
%%
result = zeros(N_basis_test,1);

%Go through all elements.
%On each element, compute the volume integrals for all test FE basis functions.
%Assemble the values of those volume integrals into the vector.
for n = 1:number_of_elements

    vertices = P_partition(:,T_partition(:,n));
    [Gauss_coefficient_local_triangle,Gauss_point_local_triangle] = generate_Gauss_local_triangle(Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, vertices);

    for beta = 1:number_of_test_local_basis     
       int_value = Gauss_quadrature_for_surface_integral_test_triangle(coefficient_function, Gauss_coefficient_local_triangle, Gauss_point_local_triangle, vertices,...
           basis_type_test, beta, derivative_degree_x_test, derivative_degree_y_test);
       i = T_basis_p(n);
       result(i,1) = result(i,1) + int_value;
    end

end