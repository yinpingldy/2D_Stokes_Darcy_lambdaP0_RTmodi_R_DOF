function [result] = Delta_mass_conservation(uh1, uh2, coefficient_function, triangle_area, P_partition, T_partition, number_of_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, stokes_t, nodes_stokes, Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle,...
    derivative_degree, basis_type, derivative_degree_x, derivative_degree_y)
%%
result = 0;
%Go through all elements and accumulate the error on them.
for n = 1:number_of_elements
    vertices = P_partition(:,T_partition(:,n)); % T:triangle
    
    result = result + Gauss_quadrature_for_surface_integral_u_mass_triangle(n, uh1, uh2, coefficient_function, triangle_area, vertices,...
    T_partition, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, stokes_t, nodes_stokes, Gauss_coefficient_reference_triangle, Gauss_point_reference_triangle,...
    derivative_degree, basis_type, derivative_degree_x, derivative_degree_y);

end
result = sqrt(result);
end

