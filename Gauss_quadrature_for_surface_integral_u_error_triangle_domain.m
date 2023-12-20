function [result] = Gauss_quadrature_for_surface_integral_u_error_triangle_domain(element, uh, index_u, domain, domain_t, exact_function, vertices,...
    triangle, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    derivative_degree, basis_type, derivative_degree_x, derivative_degree_y)
%% Use Gauss quadrature to numerically compute a norm error of FE solution on a local triangular element T.
%accurate_function: the accurate function in the error.
%When we take the L2 norm,accurate_function is the exact solution.
%When we take the H1 seminorm, accurate_function is the first derivative of the exact solution.
%vertices: the coordinates of the vertices of the triangular element T.
%uh_local: the values of the FE solution at the nodes of FE in the triangular element T.
%Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle: the Gauss coefficients and Gauss points on the reference triangle.
%derivative_degree_x:the derivative degree of the FE solution with respect to x.
%derivative_degree_y:the derivative degree of the FE solution with respect to y.
%%
Gpn = length(Gauss_weights_reference_triangle); % Gpn: the Gauss point number.

result = 0;
[Gauss_weights_local_triangle,Gauss_nodes_local_triangle] = generate_Gauss_local_triangle(Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, vertices);

for i = 1:Gpn
    result = result + Gauss_weights_local_triangle(i) * (feval(exact_function, Gauss_nodes_local_triangle(i,1), Gauss_nodes_local_triangle(i,2))...
        - FE_solution_local_triangle_u_domain(element, Gauss_nodes_local_triangle(i,1), Gauss_nodes_local_triangle(i,2),...
        uh, index_u, domain, domain_t, vertices, triangle, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, derivative_degree, basis_type, derivative_degree_x, derivative_degree_y))^2;
end

end