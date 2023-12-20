function [result] = Gauss_quadrature_for_surface_integral_u_mass_triangle(element, uh1, uh2,  coefficient_function, triangle_area, vertices,...
    triangle, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, stokes_t, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
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

Ph_fun = Gauss_quadrature_for_surface_integral(coefficient_function, Gauss_weights_local_triangle, Gauss_nodes_local_triangle)/triangle_area(element);
for i = 1:Gpn
    [~, div_project] = triangular_div_RT0_project(element, Gauss_nodes_local_triangle(i,1), Gauss_nodes_local_triangle(i,2), uh1, uh2, vertices, triangle, triangle_edge, Gamma_inter_num, number_of_edges,...
        nodes_stokes_num, stokes_t, nodes_stokes, derivative_degree, basis_type, derivative_degree_x, derivative_degree_y);
    
    result = result + Gauss_weights_local_triangle(i) * ( div_project - Ph_fun)^2;
end

end