function [Gamma_stokes, Gamma_darcy, Gamma_inter, Gamma_stokes_num, Gamma_darcy_num, Gamma_inter_num, Gamma_inter_nodes, Boundary_num] = split_boundary_edges(boundary_edges,edges,P_partition)
% 这里用到了e矩阵的特点，67行代表的是区域这里01是gamma-s 02是gamma-d 12是gamma-I
number_of_boundary_edges = size(boundary_edges,2);
Boundary_edge = zeros(number_of_boundary_edges,2);
boundary_belong = boundary_edges(6,:)+boundary_edges(7,:);
for i=1:number_of_boundary_edges
    Boundary_edge(i,1)=boundary_edges(1,i);
    Boundary_edge(i,2)=boundary_edges(2,i);
end
Boundary_edge=sort(Boundary_edge,2); % 边界边
[~,~,Boundary_num]=intersect(Boundary_edge,edges,'rows','stable');
gamma_s=boundary_belong==1;
gamma_d=boundary_belong==2;
gamma_i=boundary_belong==3;
Gamma_stokes = Boundary_edge(gamma_s,:);
Gamma_darcy = Boundary_edge(gamma_d,:);
Gamma_inter = Boundary_edge(gamma_i,:);
Gamma_stokes_num=Boundary_num(gamma_s,:);
Gamma_darcy_num=Boundary_num(gamma_d,:);
Gamma_inter_num=Boundary_num(gamma_i,:);

Gamma_inter_nodes = unique(Gamma_inter);
N = size(Gamma_inter_nodes,1);
for i = 1:N
    for j = 1:N-i
        if P_partition(2,Gamma_inter_nodes(j)) > P_partition(2,Gamma_inter_nodes(j+1))
            temp = Gamma_inter_nodes(j);
            Gamma_inter_nodes(j) = Gamma_inter_nodes(j+1);
            Gamma_inter_nodes(j+1) = temp;
        end
    end
end

end

