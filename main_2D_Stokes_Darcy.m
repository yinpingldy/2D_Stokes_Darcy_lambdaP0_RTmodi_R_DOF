clc
clear
%%
global mu solution_case alpha Ra

for n = 1:5
filename = [ 'data' num2str(n) '.mat'];
load(filename)
% Check the 2D linear finite element.
basis_type_u = 1;
basis_type_p = 0;
basis_type_lambda = 100; % 1D
% mu = 1;
mu = 10^(-6);
alpha = 1;
solution_case = 1;
Ra = 1;
% The problem domain is [left,right]*[bottom,top].
Gauss_nodes_number = 9;
tic
[ uh_1, uh_2, ph, error_u_L2, error_S_u_H1_semi, error_S_u_L2, error_D_u_L2, error_p_L2, lambda, delta_mass] = FE_solver_2D_Stokes_Darcy_triangle(P_partition, T_partition, boundary_edges, basis_type_u, basis_type_p, basis_type_lambda);
t=toc;
%% error analysis
fid=fopen('result.txt','at');
fprintf(fid, '%d\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t\n', n, error_u_L2, error_S_u_H1_semi, error_S_u_L2, error_D_u_L2, error_p_L2, delta_mass, t);
fclose(fid);
end