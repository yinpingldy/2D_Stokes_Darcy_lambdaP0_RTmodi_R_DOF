function [ uh_1, uh_2, ph, error_u_L2, error_S_u_H1_semi, error_S_u_L2, error_D_u_L2, error_p_L2, lambda, delta_mass] = FE_solver_2D_Stokes_Darcy_triangle(P_partition, T_partition, boundary_edges, basis_type_u, basis_type_p, basis_type_lambda)
%% Finite element solver for Stokes-Darcy equation on a triangular mesh.
% The trial FE functions and test FE functions need to be the same.
% The problem domain is [left,right]*[bottom,top].
% h_partition is the step size of the partition.
% basis_type: the type of the FE basis function.
% basis_type=1: 2D linear FE.  
% basis_type=2: 2D Lagrange quadratic FE.
% basis_type=10: 2D Crouzeix-Raviart FE.

% N1_basis,N2_basis:The N1 and N2 for the FE basis functions,not the partition.
% N1_partition,N2_partition:The N1 and N2 for the partition,not the FE basis functions.
% N1 is the number of the sub-intervals of the partition in x-direction.
% N2 is the number of the sub-intervals of the partition in y-direction.
% h_basis is the step size for the finite element nodes, not the partition.
% P_partition,T_partition, P_basis,T_basis: see the note in "generate_P_T_triangular.m".
% "u" is for the velocity vector function u.
% function_Lamda, function_Mu: the two coefficient functions in Elasticity equation.
% funciton_F1,funciton_F2: the right hand side functions of the Stokes equation.
% function_G1,function_G2: the Dirichelet boundary functions for (u1,u2).
%% Preliminary
% Gamma_stokes, Gamma_darcy, Gamma_inter: 边界边--（stokes,darcy,交界面）两顶点
% Gamma_stokes_num, Gamma_darcy_num, Gamma_inter_num: 各个区域边界边编号（在全部边中的编号）
% stokes_t, darcy_t: 两个区域的单元编号，哪个单元为stokes，哪个单元为darcy
% triangle_stokes, triangle_darcy: 每个区域上的三角形单元顶点，与上述编号对应
% Boundary_num: 边界边在全部边里的编号
% triangle, triangle_area, triangle_edge: 全部三角形单元顶点；每个单元的面积；每个单元的边在全部边中的编号
% point: 全部节点坐标
% edge, Boundary_edge: 全部边集合，边界边集合

number_of_elements = size(T_partition,2);
[triangle, triangle_Stokes, triangle_Darcy, Stokes_t, Darcy_t, number_of_Stokes_elements, number_of_Darcy_elements] = split_T(T_partition);

[edges, triangle_edge] = generate_edges(number_of_elements,T_partition); % triangle_edge：每个单元的边在全部边中的编号
number_of_edges = size(edges,1);
edge_condition = edge_condition_treat(edges,triangle_edge);
triangle_area = cal_triangle_area(number_of_elements,P_partition,T_partition);

[Gamma_stokes, Gamma_darcy, Gamma_inter, Gamma_stokes_num, Gamma_darcy_num, Gamma_inter_num, Gamma_inter_nodes, Boundary_num] = split_boundary_edges(boundary_edges, edges, P_partition);
Gamma_stokes_num_full = [Gamma_stokes_num;Gamma_inter_num];

number_of_Gamma_inter = size(Gamma_inter_num,1);

nodes_stokes = unique(triangle_Stokes);
nodes_stokes_num = size(nodes_stokes,1);
% Mesh information for partition and finite element basis functions.
N_basis_u = number_of_edges + nodes_stokes_num + number_of_Gamma_inter;
T_basis_p = zeros(1,number_of_elements);
for n = 1:number_of_elements
    T_basis_p(1,n) = n;
end
N_basis_p = number_of_elements;

% N_element_lambda = size(Gamma_inter_num,1);
% T_basis_lambda = zeros(2,N_element_lambda);
% for n = 1:N_element_lambda
%    T_basis_lambda(1,n) = n;    
%    T_basis_lambda(2,n) = n+1;
% end
N_basis_lambda = number_of_Gamma_inter;
% N_basis_lambda = (size(Gamma_inter_nodes,1)+1)/2;

number_of_local_basis_u = 3; % P1 or RT0
number_of_local_basis_p = 1; % P0
number_of_local_basis_lambda = 1; % P0

% Guass quadrature's points and coefficient on the refenrece triangle and reference interval.
Gauss_nodes_number = 9;
[Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle] = generate_Gauss_reference_triangle(Gauss_nodes_number);
% The following line is necessary if we need to deal with Neumann or Robin boundary condition. Otherwise, we can comment this line to save time and memory.
[Gauss_weights_reference_1D,Gauss_nodes_reference_1D] = generate_Gauss_reference_1D(4);

%% Assemble the stiffness matrix. 
AA1 = assemble_matrix_from_surface_integral_triangle('function_2mu', P_partition, triangle_Stokes, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u, 1, 0, basis_type_u, 1, 0);
AA2 = assemble_matrix_from_surface_integral_triangle('function_mu', P_partition, triangle_Stokes, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u, 0, 1, basis_type_u, 0, 1);
% AA3 = assemble_matrix_from_surface_integral_triangle_SRT0_0('function_2mu_alphaT', P_partition, triangle_Stokes, N_basis_u, N_basis_u, Stokes_t, triangle_edge, Gamma_stokes_num_full,...
%     number_of_local_basis_u, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, 0, 0); % 需处理 h
AA3 = assemble_matrix_from_surface_integral_triangle_SRT0_div('function_2mu_alphaT', P_partition, triangle_Stokes, N_basis_u, N_basis_u, Stokes_t, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num,...
    number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, 1, 1); %未加边界

AA4 = assemble_matrix_from_surface_integral_triangle_DRT0('function_mu__K', P_partition, triangle_Darcy, N_basis_u, N_basis_u, Darcy_t, triangle_edge,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Darcy_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, 0, 0);
AA5 = assemble_matrix_from_line_integral_triangle('function_mu_alpha__sqrtK', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges, ...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 1, 1, basis_type_u, 0, 0, basis_type_u, 0, 0);

AA6 = assemble_matrix_from_line_integral_triangle_L_RT0('function_2mu', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges, nodes_stokes_num, triangle_edge,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 1, basis_type_u, 1, 0, 0);
AA7 = assemble_matrix_from_line_integral_triangle_RT0_L('function_neg_2mu', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges, nodes_stokes_num, triangle_edge,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 1, 0, basis_type_u, 1, 0);
% AA8 = assemble_matrix_from_line_integral_triangle_RT0('function_mu_alpha__sqrtK', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, triangle_edge, Gamma_stokes_num_full,...
%     number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 0, 0); % 交界面有几项，包括RT0的吗

AB1 = assemble_matrix_from_surface_integral_triangle('function_mu', P_partition, triangle_Stokes, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u, 1, 0, basis_type_u, 0, 1);
AB2 = assemble_matrix_from_line_integral_triangle('function_mu_alpha__sqrtK', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 2, 1, basis_type_u, 0, 0, basis_type_u, 0, 0);
% AB3 = assemble_matrix_from_line_integral_triangle_L_RT0('function_mu_alpha__sqrtK', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges, triangle_edge, Gamma_stokes_num_full,...
%     number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 2, basis_type_u, 0, 0, 0);

AC1 = assemble_matrix_from_surface_integral_triangle_p_u('function_negative_one', P_partition, triangle_Stokes, N_basis_p, N_basis_u, T_basis_p, Stokes_t, nodes_stokes, number_of_edges,...
    number_of_local_basis_p, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_p, 0, 0, basis_type_u, 1, 0); 
AC2 = assemble_matrix_from_surface_integral_triangle_P_SRT0('function_negative_one', P_partition, triangle_Stokes, N_basis_p,  N_basis_u, T_basis_p, Stokes_t, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num,...
    number_of_local_basis_p, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_p, 0, 0, 1);
AC3 = assemble_matrix_from_surface_integral_triangle_P_DRT0('function_negative_one', P_partition, triangle_Darcy, N_basis_p, N_basis_u, T_basis_p, Darcy_t, triangle_edge,...
    number_of_local_basis_p, number_of_local_basis_u, number_of_Darcy_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_p, 0, 0, 1);

AD1 = assemble_matrix_from_line_integral_triangle_lambda('function_one', edges, Gamma_inter_num, Gamma_inter_nodes, edge_condition, Stokes_t, P_partition, triangle, N_basis_lambda, N_basis_u, nodes_stokes, number_of_edges, ...
    number_of_local_basis_lambda, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 1, basis_type_lambda, 0, basis_type_u, 0, 0);
AD2 = assemble_matrix_from_line_integral_triangle_lambda_SRT0('function_one', edges, Gamma_inter_num, Gamma_inter_nodes, edge_condition, Stokes_t, P_partition, triangle, N_basis_lambda,  N_basis_u, triangle_edge, number_of_edges, nodes_stokes_num,...
    number_of_local_basis_lambda, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, basis_type_lambda, 0, 0); % SRT0加吗
AD3 = assemble_matrix_from_line_integral_triangle_lambda_DRT0('function_one', edges, Gamma_inter_num, Gamma_inter_nodes, edge_condition, Darcy_t, P_partition, triangle, N_basis_lambda,  N_basis_u, triangle_edge,...
    number_of_local_basis_lambda, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, basis_type_lambda, 0, 0);

BA1 = assemble_matrix_from_surface_integral_triangle('function_mu', P_partition, triangle_Stokes, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u, 0, 1, basis_type_u, 1, 0); % =AB1'
BA2 = assemble_matrix_from_line_integral_triangle('function_mu_alpha__sqrtK', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 1, 2, basis_type_u, 0, 0, basis_type_u, 0, 0); % =AB2'
% BA3 = assemble_matrix_from_line_integral_triangle_RT0_L('function_mu_alpha__sqrtK', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges, triangle_edge, Gamma_stokes_num_full,...
%     number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 2, 0, basis_type_u, 0, 0);

BB1 = assemble_matrix_from_surface_integral_triangle('function_2mu', P_partition, triangle_Stokes, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u, 0, 1, basis_type_u, 0, 1);
BB2 = assemble_matrix_from_surface_integral_triangle('function_mu', P_partition, triangle_Stokes, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u, 1, 0, basis_type_u, 1, 0);
BB3 = assemble_matrix_from_line_integral_triangle('function_mu_alpha__sqrtK', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 2, 2, basis_type_u, 0, 0, basis_type_u, 0, 0);

BC = assemble_matrix_from_surface_integral_triangle_p_u('function_negative_one', P_partition, triangle_Stokes, N_basis_p, N_basis_u, T_basis_p, Stokes_t, nodes_stokes, number_of_edges,...
    number_of_local_basis_p, number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_p, 0, 0, basis_type_u, 0, 1);

BD = assemble_matrix_from_line_integral_triangle_lambda('function_one', edges, Gamma_inter_num, Gamma_inter_nodes, edge_condition, Stokes_t, P_partition, triangle, N_basis_lambda, N_basis_u, nodes_stokes, number_of_edges, ...
    number_of_local_basis_lambda, number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 2, basis_type_lambda, 0, basis_type_u, 0, 0);

O_33 = sparse(N_basis_p, N_basis_p);
O_34 = sparse(N_basis_p, N_basis_lambda);
O_44 = sparse(N_basis_lambda, N_basis_lambda);
A = [AA1+AA2+AA3+AA4+AA5+AA6+AA7  AB1+AB2  AC1+AC2+AC3 AD1+AD2+AD3;
            BA1+BA2               BB1+BB2+BB3       BC          BD; 
          AC1'+AC2'+AC3'                  BC'          O_33        O_34;
          AD1'+AD2'+AD3'                  BD'          O_34'       O_44];
clear AA1 AA2 AA3 AA4 AA5 AA6 AA7 AB1 AB2 AC1 AC2 AC3 AD1 AD2 AD3 BA1 BA2 BB1 BB2 BB3 BC BD O_33 O_34 O_44 
%% Assemble the load vector.
b11 = assemble_vector_from_surface_integral_triangle_Proj_S('function_F1','function_F2', P_partition, triangle_Stokes, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u,0,0, 1);
b12 = assemble_vector_from_surface_integral_triangle_SRT0('function_F1','function_F2', P_partition, triangle_Stokes, N_basis_u, Stokes_t, triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num,...
    number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, 0);
b13 = assemble_vector_from_surface_integral_triangle_DRT0('function_F1','function_F2', P_partition, triangle_Darcy, N_basis_u, Darcy_t, triangle_edge,...
    number_of_local_basis_u, number_of_Darcy_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, 0);
b14 = assemble_vector_from_line_integral_triangle('function_neg_theta1', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle,  N_basis_u, nodes_stokes, number_of_edges,...
     number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 'n', 1, basis_type_u, 0, 0,  1);
b15 = assemble_vector_from_line_integral_triangle_SRT0('function_neg_theta1', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, triangle_edge, number_of_edges, nodes_stokes_num,...
    number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 'n', 0);
b16 = assemble_vector_from_line_integral_triangle('function_mu_alpha__sqrtK_theta2', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle,  N_basis_u, nodes_stokes, number_of_edges,...
     number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 't', 1, basis_type_u, 0, 0, 1); % 0
% b17 = assemble_vector_from_line_integral_triangle_SRT0('function_mu_alpha__sqrtK_theta2', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle, N_basis_u, triangle_edge, Gamma_stokes_num_full,...
%     number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 't', 0);

% b21 = assemble_vector_from_surface_integral_triangle('function_F2',P_partition, triangle_Stokes, N_basis_u, nodes_stokes, number_of_edges,...
%     number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
%     basis_type_u,0,0);
b21 = assemble_vector_from_surface_integral_triangle_Proj_S('function_F1','function_F2', P_partition, triangle_Stokes, N_basis_u, nodes_stokes, number_of_edges,...
    number_of_local_basis_u, number_of_Stokes_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_u,0,0, 2);
b22 = assemble_vector_from_line_integral_triangle('function_neg_theta1', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle,  N_basis_u, nodes_stokes, number_of_edges,...
     number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 'n', 2, basis_type_u, 0, 0, 2); % 0
b23 = assemble_vector_from_line_integral_triangle('function_mu_alpha__sqrtK_theta2', edges, Gamma_inter_num, edge_condition, Stokes_t, P_partition, triangle,  N_basis_u, nodes_stokes, number_of_edges,...
     number_of_local_basis_u, number_of_Gamma_inter, Gauss_weights_reference_1D, Gauss_nodes_reference_1D, 't', 2, basis_type_u, 0, 0, 2);
 
b3 = assemble_vector_from_surface_integral_triangle_p('function_G',P_partition, triangle, N_basis_p, T_basis_p,...
    number_of_local_basis_p, number_of_elements, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    basis_type_p,0,0);
b4 = zeros(N_basis_lambda,1);
b = [b11+b12+b13+b14+b15+b16; b21+b22+b23; -b3; b4];

clear b11 b12 b13 b14 b15 b16 b21 b22 b23 b3 b4
%% Treat boundary condition

[A, b] = treat_boundary_condition('function_exact_solution_u1', 'function_exact_solution_u2', A, b, P_partition, N_basis_u, N_basis_p,...
    Gamma_stokes, Gamma_stokes_num, number_of_edges, Gamma_darcy_num, Gamma_darcy, edge_condition, triangle, triangle_area, nodes_stokes, Gauss_weights_reference_1D, Gauss_nodes_reference_1D);

%% Compute the numerical solution and Get the finite element solution for u1, u2 and p.

tips_A = find(~any(A));
for i = 1:size(tips_A,1)
    A(tips_A(i),tips_A(i)) = 1;
    b(tips_A(i),:) = 0;
end

solution = A\b; %直接求逆，若存在奇异情况，结果可能出错

% A = full(A);
% solution = pinv(A)*b; %求伪逆，结果准确但需要转换为全矩阵，计算成本和内存占用大
% solution = lsqminnorm(A, b, 1e-12);  % 最小二乘求解，速度快，但有可能会出现结果不准确现象

uh_1 = solution(1:N_basis_u);
uh_2 = solution(N_basis_u+1 : 2*N_basis_u);
uh_2(1:number_of_edges) = uh_1(1:number_of_edges);
uh_2(number_of_edges+nodes_stokes_num+1:end) = uh_1(number_of_edges+nodes_stokes_num+1:end);
ph = solution(2*N_basis_u+1 : 2*N_basis_u+N_basis_p);
lambda = solution(2*N_basis_u+N_basis_p+1:end);
%% mass conservation
delta_mass = Delta_mass_conservation(uh_1, uh_2, 'function_G', triangle_area, P_partition, triangle, number_of_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, Stokes_t, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 0);
%% error computation
error_u1_L2 = FE_solution_error_triangle_u('function_exact_solution_u1', uh_1, 1, P_partition, triangle, number_of_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, Stokes_t, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 0);
error_u2_L2 = FE_solution_error_triangle_u('function_exact_solution_u2', uh_2, 2, P_partition, triangle, number_of_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, Stokes_t, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 0);
error_p_L2 = FE_solution_error_triangle('function_exact_solution_p', ph, P_partition, triangle, T_basis_p, number_of_elements,...
    Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle, basis_type_p, 0, 0);
error_u_L2 = sqrt(error_u1_L2^2+error_u2_L2^2);

% Stokes
% L2
error_S_u1_L2 = FE_solution_error_triangle_u_domain('function_exact_solution_u1', uh_1, 1, 'S', Stokes_t, P_partition, triangle_Stokes, number_of_Stokes_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 0);
error_S_u2_L2 = FE_solution_error_triangle_u_domain('function_exact_solution_u2', uh_2, 2, 'S', Stokes_t, P_partition, triangle_Stokes, number_of_Stokes_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 0);
error_S_u_L2 = sqrt(error_S_u1_L2^2+error_S_u2_L2^2);
% H1-semi
error_S_u1_H1_semi_x = FE_solution_error_triangle_u_domain('function_exact_solution_u1_x', uh_1, 1, 'S', Stokes_t, P_partition, triangle_Stokes, number_of_Stokes_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 1, 0);
error_S_u1_H1_semi_y = FE_solution_error_triangle_u_domain('function_exact_solution_u1_y', uh_1, 1, 'S', Stokes_t, P_partition, triangle_Stokes, number_of_Stokes_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 1);
% error_u1_H1_semi = sqrt(error_u1_H1_semi_x^2 + error_u1_H1_semi_y^2);
error_S_u2_H1_semi_x = FE_solution_error_triangle_u_domain('function_exact_solution_u2_x', uh_2, 2, 'S', Stokes_t, P_partition, triangle_Stokes, number_of_Stokes_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 1, 0);
error_S_u2_H1_semi_y = FE_solution_error_triangle_u_domain('function_exact_solution_u2_y', uh_2, 2, 'S', Stokes_t, P_partition, triangle_Stokes, number_of_Stokes_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 1);
error_S_u_H1_semi = sqrt(error_S_u1_H1_semi_x^2 + error_S_u1_H1_semi_y^2 + error_S_u2_H1_semi_x^2 + error_S_u2_H1_semi_y^2);

% Darcy
error_D_u1_L2 = FE_solution_error_triangle_u_domain('function_exact_solution_u1', uh_1, 1, 'D', Darcy_t, P_partition, triangle_Darcy, number_of_Darcy_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 0);
error_D_u2_L2 = FE_solution_error_triangle_u_domain('function_exact_solution_u2', uh_2, 2, 'D', Darcy_t, P_partition, triangle_Darcy, number_of_Darcy_elements,...
    triangle_edge, Gamma_inter_num, number_of_edges, nodes_stokes_num, nodes_stokes, Gauss_weights_reference_triangle, Gauss_nodes_reference_triangle,...
    0, basis_type_u, 0, 0);
error_D_u_L2 = sqrt(error_D_u1_L2^2+error_D_u2_L2^2);
end
