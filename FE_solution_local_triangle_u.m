function [result] = FE_solution_local_triangle_u(element, x, y, uh, index_u, vertices, triangle, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, stokes_t, nodes_stokes, derivative_degree, basis_type, derivative_degree_x, derivative_degree_y)
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
number_of_local_basis = 3;
edge_element = triangle_edge(:,element);

if ismember(element, stokes_t) == 1
    for alpha = 1:number_of_local_basis
        node = triangle(alpha,element);
        i_node = find(nodes_stokes == node) + number_of_edges;
        uh_L = uh(i_node);
        i = edge_element(alpha);
        if any(Gamma_inter_num == i)
            i = find(Gamma_inter_num == i) + number_of_edges + nodes_stokes_num;
        end
        Q_basis = triangular_RT0_local_basis(x, y, vertices, alpha ,derivative_degree);
        result = result + uh(i) * Q_basis(index_u) + uh_L * triangular_local_basis(x, y, vertices, basis_type, alpha ,derivative_degree_x, derivative_degree_y);
    end

else
    for alpha = 1:number_of_local_basis
        Q_basis = triangular_RT0_local_basis(x, y, vertices, alpha ,derivative_degree);
        result = result + uh(edge_element(alpha)) * Q_basis(index_u);
    end
end

end