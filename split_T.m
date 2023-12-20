function [triangle, triangle_Stokes, triangle_Darcy, Stokes_t, Darcy_t, number_of_Stokes_elements, number_of_Darcy_elements] = split_T(T_partition)
triangle = T_partition((1:3),:);
triangle_belong = T_partition(4,:);

%t的第四行表示了在哪个区域 1 是stokes 2是darcy
Stokes_t = triangle_belong==1;
Darcy_t = triangle_belong==2;
triangle_Stokes=triangle(:,Stokes_t);
triangle_Darcy=triangle(:,Darcy_t);
NUM = [1:1:size(triangle,2)];
Stokes_t = NUM(:,Stokes_t); % 那个单元为stokes
Darcy_t = NUM(:,Darcy_t);
number_of_Stokes_elements = size(triangle_Stokes,2);
number_of_Darcy_elements = size(triangle_Darcy,2);
end

