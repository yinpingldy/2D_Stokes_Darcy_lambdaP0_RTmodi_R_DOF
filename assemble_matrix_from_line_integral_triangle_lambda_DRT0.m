function result = assemble_matrix_from_line_integral_triangle_lambda_DRT0(coefficient_function, edges, Gamma_inter_num, Gamma_inter_nodes, edge_condition, Stokes_or_Darcy_t, P_partition, T_partition, N_basis_trial,  N_basis_test, triangle_edge,...
    number_of_trial_local_basis, number_of_test_local_basis, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, trial_basis_type_1D, trial_derivative_degree_1D, test_derivative_degree)
%% Assemble a matrix from a surface inegral on the whole domain with a triangular mesh.
%The integrand of the volume integral must be in the following format:
%a coefficient function * a trial FE basis function(or its derivatives) * a test FE basis function (or its derivatives).
%coefficient_function_name: the coefficient function of the integrand.
%P_partition: store the coordinates of all the grid points of the partition,not FE.
%T_partition: store the global indices of the grid points of every element for the partition,not FE.
%T_basis: store the global indices of the nodes of every element for FE,not the partition.
%T_basis_trial: T_basis for the trial basis function.
%T_basis_test: T_basis for the test basis function.
%The explanation for P_partition,T_partition,T_basis is in generate_P_T_triangular.m.
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

element_use = edge_condition(:,Gamma_inter_num);
% Go through all elements.
% On each element, compute the volume integrals for all possible combinations of trial and test FE basis functions.
% Assemble the values of those volume integrals into the matrix.
for n=1:number_of_Gamma_inter
    element_Gamma = element_use(:,n);
    [element_S_D,~,~]=intersect(element_Gamma,Stokes_or_Darcy_t);
    vertices = P_partition(:,T_partition(:,element_S_D));
    
    j_nodes_h = edges(Gamma_inter_num(n),:);
    endpoint = P_partition(:,j_nodes_h);
    vertices_inside_node = setdiff(vertices',endpoint','rows');
    vector_x = endpoint(1,1)-endpoint(1,2);
    vector_y = endpoint(2,1)-endpoint(2,2);
    if vector_x == 0
        vertices_1D = [endpoint(2,1),endpoint(2,2)];
        if vertices_inside_node(1) < endpoint(1,1)
            normal_vector = [1,0];
        else
            normal_vector = [-1,0];
        end
    elseif vector_y == 0
        vertices_1D = [endpoint(1,1),endpoint(1,2)];
        if vertices_inside_node(2) < endpoint(2,1)
            normal_vector = [0,1];
        else
            normal_vector = [0,-1];
        end
    else
        normal_vector = [1/vector_x; -1/vector_y]; % The type does not exist in the rule area.
    end
    
    j_1 = find(Gamma_inter_nodes == j_nodes_h(1));
    j_2 = find(Gamma_inter_nodes == j_nodes_h(2));
    jj1 = min(j_1,j_2);
    jj2 = max(j_1,j_2);
    for alpha = 1:number_of_trial_local_basis
       for beta = 1:number_of_test_local_basis     
           i = triangle_edge(beta,element_S_D);
           
           if alpha == 1
               j = jj1;
           else
               j = jj2;
           end
          int_value = Gauss_quadrature_for_line_integral_trial_test_lambda_RT0(coefficient_function, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, endpoint(:,1), endpoint(:,2), normal_vector,...
              trial_basis_type_1D, vertices_1D, alpha, trial_derivative_degree_1D, vertices, beta, test_derivative_degree);
          
          result(i,j) = result(i,j) + int_value;
       end
    end

end