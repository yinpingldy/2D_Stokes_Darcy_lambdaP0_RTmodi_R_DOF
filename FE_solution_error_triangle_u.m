function [result] = FE_solution_error_triangle_u(exact_solution, uh, index_u, P_partition, T_partition, number_of_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, stokes_t, nodes_stokes, Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle,...
    derivative_degree, basis_type, derivative_degree_x, derivative_degree_y)
%% Numerically compute a norm error of FE solution on the whole domain [left,right]*[bottom,top].
%uh: the values of the FE solution at all nodes of FE in the whole domain. These values are in 1D index of nodes of FE.
%accurate_function: the accurate function in the error.
%When we take the L2 norm,accurate_function is the exact solution.
%When we take the H1 seminorm, accurate_function is the first derivative of the exact solution.
%h_partition: the stepsize of the partition.
%basis_type: the type of the FE.
%basis_type=1:2D linear FE.
%basis_type=2:2D Lagrange quadratic FE.
%derivative_degree_x:the derivative degree of the FE solution with respect to x.
%derivative_degree_y:the derivative degree of the FE solution with respect to y.
%Gauss_point_number: the number of the Gauss points of the Gauss quadrature we want to use.

%N1_partition,N2_partition:The N1 and N2 for the partition,not the FE basis functions.
%N1 is the number of the sub-intervals of the partition in x-direction.
%N2 is the number of the sub-intervals of the partition in y-direction.
%vertices: the coordinates of the vertices of a triangular element.
%uh_local: the values of the FE solution at the nodes of FE in a triangular element.
%% 
result=0;
%Go through all elements and accumulate the error on them.
for n = 1:number_of_elements
    vertices = P_partition(:,T_partition(:,n)); % T:triangle
    
    result = result + Gauss_quadrature_for_surface_integral_u_error_triangle(n, uh, index_u, exact_solution, vertices,...
        T_partition, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, stokes_t, nodes_stokes, Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle,...
        derivative_degree, basis_type, derivative_degree_x, derivative_degree_y);

end
result = sqrt(result);