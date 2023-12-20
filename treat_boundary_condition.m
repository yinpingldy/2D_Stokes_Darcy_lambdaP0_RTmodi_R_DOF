function [A, b] = treat_boundary_condition(exact_solution_u1, exact_solution_u2, A, b, P_partition, N_basis_u, N_basis_p, Gamma_stokes, Gamma_stokes_num, number_of_edges, Gamma_darcy_num, Gamma_darcy, edge_condition, triangle, triangle_area, nodes_stokes, Gauss_weights_reference_1D, Gauss_nodes_reference_1D)
%% Deal with boundary condition.
% boundary_nodes(1,k): specifiy the type of the kth boundary node.
% boundary_nodes(1,k)=-1: Dirichlet boundary node;
% boundary_nodes(1,k)=-2: Neumann boundary node;
% boundary_nodes(1,k)=-3: Robin boundary node. 
% boundary_nodes(2,k): global index of the kth boundary node among all nodes of FE. 
% P_basis(trial): store the coordinates of all the nodes for the FE,not the partition.
% Dirichelet_boundary_fucntion: the name of the Dirichelet boundary function.
%% Darcy boundary treatment
N_Gamma_Darcy = size(Gamma_darcy_num,1);
for k = 1:N_Gamma_Darcy
    Darcy_nodes = Gamma_darcy(k,:);
    vertices = P_partition(:,Darcy_nodes);
%     vertices_center = [ (vertices(1,1) + vertices(1,2))/2,(vertices(2,1) + vertices(2,2))/2];
%     
%     element_connected = edge_condition(1,Gamma_darcy_num(k));
%     element_nodes = triangle(:,element_connected);
%     inside_node = setdiff(element_nodes,Darcy_nodes);
%     vertices_inside_node = P_partition(:,inside_node);
    vector_x = vertices(1,1)-vertices(1,2);
    vector_y = vertices(2,1)-vertices(2,2);
    if vector_x == 0
        normal_vector = [1,0];
        h_e = abs(vector_y);
    elseif vector_y == 0
        normal_vector = [0,1];
        h_e = abs(vector_x);
    else
        normal_vector = [1/vector_x; -1/vector_y]; % The type does not exist in the rule area.
    end
    
    i = Gamma_darcy_num(k);
    A(i,:) = 0;
    A(i,i) = 1;
    exact = Gauss_quadrature_for_line_integral_2D(exact_solution_u1, exact_solution_u2, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, vertices(:,1), vertices(:,2), normal_vector);
    b(i,1) = exact/h_e;
    % RT0 part,Darcy u2 not true, so the following processing is unnecessary.
%     A(i+N_basis_u,:) = 0;
%     A(i+N_basis_u,i+N_basis_u) = 1;
%     b(i+N_basis_u,1) = dot( [feval(function_exact_solution_u1, vertices_center(1), vertices_center(2)), feval(function_exact_solution_u2, vertices_center(1), vertices_center(2))], normal_verctor);
end

%% Interface boundary treatment
% N_Gamma_Interface = size(Gamma_inter_num,1);
% for k = 1:N_Gamma_Interface
%     Interface_nodes = Gamma_inter(k,:);
%     vertices = P_partition(:,Interface_nodes);
%     vertices_center = [ (vertices(1,1) + vertices(1,2))/2,(vertices(2,1) + vertices(2,2))/2];
%     
%     element_connected = edge_condition(1,Gamma_darcy_num(k));
%     element_nodes = triangle(element_connected);
%     inside_node = setdiff(element_nodes,Interface_nodes);
%     vertices_inside_node = P_partition(inside_node);
%     verctor_x = vertices(1,1)-vertices(1,2);
%     verctor_y = vertices(2,1)-vertices(2,2);
%     if verctor_x == 0
%         if vertices_inside_node(1) < vertices(1,1)
%             normal_verctor = [1,0];
%         else
%             normal_verctor = [-1,0];
%         end
%     elseif verctor_y == 0
%         if vertices_inside_node(2) < vertices(2,1)
%             normal_verctor = [0,1];
%         else
%             normal_verctor = [0,-1];
%         end
%     else
%         normal_verctor = [1/verctor_x; -1/verctor_y]; % The type does not exist in the rule area.
%     end
%     
%     i = Gamma_darcy_num(k);
%     A(i,:) = 0;
%     A(i,i) = 1;
%     b(i,1) = dot( [feval(function_exact_solution_u1, vertices_center(1), vertices_center(2)),...
%         feval(function_exact_solution_u2, vertices_center(1), vertices_center(2))], normal_verctor); % =0
%     % RT0 part,Darcy u2 not true, so the following processing is unnecessary.
% %     A(i+N_basis_u,:) = 0;
% %     A(i+N_basis_u,i+N_basis_u) = 1;
% %      b(i+N_basis_u,1) = dot( [feval(function_exact_solution_u1, vertices_center(1), vertices_center(2)),...
% %          feval(function_exact_solution_u2, vertices_center(1), vertices_center(2))], normal_verctor);
% end

%% Stokes boundary treatment: Dirichlet
N_Gamma_Stokes = size(Gamma_stokes_num,1);
for k = 1:N_Gamma_Stokes
    Stokes_nodes = Gamma_stokes(k,:);
    vertices = P_partition(:,Stokes_nodes);
%     vertices_center = [ (vertices(1,1) + vertices(1,2))/2;(vertices(2,1) + vertices(2,2))/2];

    vector_x = vertices(1,1)-vertices(1,2);
    vector_y = vertices(2,1)-vertices(2,2);
    if vector_x == 0
        normal_vector = [1,0];
        h_e = abs(vector_y);
    elseif vector_y == 0
        normal_vector = [0,1];
        h_e = abs(vector_x);
    else
        normal_vector = [1/vector_x; -1/vector_y]; % The type does not exist in the rule area.
    end

%     b(i,1) = feval(exact_solution_u1, vertices_center(1), vertices_center(2));
    % RT0 part u2 not true, so the following processing is unnecessary.
%     A(i+N_basis_u,:) = 0;
%     A(i+N_basis_u,i+N_basis_u) = 1;
%     b(i+N_basis_u,1) = feval(function_exact_solution_u2, vertices_center(1), vertices_center(2));

    j_1 = find(nodes_stokes == Stokes_nodes(1)) + number_of_edges;
    j_2 = find(nodes_stokes == Stokes_nodes(2)) + number_of_edges;
    A(j_1,:) = 0;
    A(j_1,j_1) = 1;
    b(j_1,1) = feval(exact_solution_u1, vertices(1,1), vertices(2,1));
    A(j_2,:) = 0;
    A(j_2,j_2) = 1;
    b(j_2,1) = feval(exact_solution_u1, vertices(1,2), vertices(2,2));
    A(j_1+N_basis_u,:) = 0;
    A(j_1+N_basis_u,j_1+N_basis_u) = 1;
    b(j_1+N_basis_u,1) = feval(exact_solution_u2, vertices(1,1), vertices(2,1));
    A(j_2+N_basis_u,:) = 0;
    A(j_2+N_basis_u,j_2+N_basis_u) = 1;
    b(j_2+N_basis_u,1) = feval(exact_solution_u2, vertices(1,2), vertices(2,2));
    
    i = Gamma_stokes_num(k);
    A(i,:) = 0;
    A(i,i) = 1;
    if vector_x == 0
        A(i,j_1) = 1/2;
        A(i,j_2) = 1/2;
    elseif vector_y == 0
        A(i,j_1+N_basis_u) = 1/2;
        A(i,j_2+N_basis_u) = 1/2;
    end
    exact = Gauss_quadrature_for_line_integral_2D(exact_solution_u1, exact_solution_u2, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, vertices(:,1), vertices(:,2), normal_vector);
    b(i,1) = exact/h_e;
end
%% Determining pressure solution
N_PU = N_basis_u*2 + N_basis_p;
A(N_PU,:)=0;
for i = N_basis_u + N_basis_u + 1 : N_PU
    num = i-(N_basis_u + N_basis_u);
    A(N_PU,i) = triangle_area(num);
end
b(N_PU) = 0;

end