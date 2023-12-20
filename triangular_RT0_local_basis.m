function [q_RT] = triangular_RT0_local_basis(x, y, vertices, basis_index, derivative_degree)
%% the local basis functions of triangular FE.
% x,y: the coordinates of the point where we want to evaluate the local FE basis function.
% basis_index: the index of basis function to specify which basis function we want to use.
% derivative_degree_x: the derivative degree of the FE basis function with respect to x.
% derivative_degree_y: the derivative degree of the FE basis function with respect to y.

% J is the Jacobi matrix of the affine mapping, which has four elememts J_11, J_12, J_21, J_22.
% J_det is the determinant of J.
% x_hat,y_hat: the '\hat{x}' and '\hat{y}' in the affine mapping, which are for the coordinates of the reference element.
%% 
vertices_0 = vertices;
vertices_0(:,4) = vertices(:,1);
mid_point = zeros(size(vertices));
for i=1:3
    mid_point(:,i) = 1/2*(vertices_0(:,i+1) - vertices_0(:,i)) + vertices_0(:,i);
end

normal_vector = generate_normal_verctor (vertices); % Not the outer normal,ensure that the normal direction is continuous, and the same edge has the same direction.
A=zeros(3);
for i=1:3
    A(i,1)=normal_vector(1,i);
    A(i,2)=normal_vector(2,i);
    A(i,3)=normal_vector(1,i)*mid_point(1,i)+normal_vector(2,i)*mid_point(2,i);
end
basic_project_martix(:,1)=A\[0;1;0];
basic_project_martix(:,2)=A\[0;0;1];
basic_project_martix(:,3)=A\[1;0;0];
% basic_project_martix(:,1)=A\[1;0;0];
% basic_project_martix(:,2)=A\[0;1;0];
% basic_project_martix(:,3)=A\[0;0;1];

q_RT = zeros(1,2);
if derivative_degree == 0
    if basis_index == 1
        q_RT(1) = basic_project_martix(1,1)+basic_project_martix(3,1)*x;
        q_RT(2) = basic_project_martix(2,1)+basic_project_martix(3,1)*y;
    elseif basis_index == 2
        q_RT(1) = basic_project_martix(1,2)+basic_project_martix(3,2)*x;
        q_RT(2) = basic_project_martix(2,2)+basic_project_martix(3,2)*y;
    elseif basis_index == 3
        q_RT(1) = basic_project_martix(1,3)+basic_project_martix(3,3)*x;
        q_RT(2) = basic_project_martix(2,3)+basic_project_martix(3,3)*y;
    end
elseif derivative_degree == 1
    if basis_index == 1
        q_RT(1) = 2*basic_project_martix(3,1);
    elseif basis_index == 2
        q_RT(1) = 2*basic_project_martix(3,2);
    elseif basis_index == 3
        q_RT(1) = 2*basic_project_martix(3,3);
    end
end

end