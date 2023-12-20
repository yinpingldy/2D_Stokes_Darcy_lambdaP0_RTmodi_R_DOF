function [normal_vector] = generate_normal_verctor (vertices)
%%
% Not the outer normal,ensure that the normal direction is continuous, and the same edge has the same direction.
%%
normal_vector = zeros(2,3);
vertices_0 = vertices;
vertices_0(:,4) = vertices(:,1);

for i = 1:3
    endpoint = vertices_0(:,i:i+1);
    
    vector_x = endpoint(1,1)-endpoint(1,2);
    vector_y = endpoint(2,1)-endpoint(2,2);
    if vector_x == 0
        normal_vector_0 = [1,0];
    elseif vector_y == 0
        normal_vector_0 = [0,1];
    else
        normal_vector_0 = [1/vector_x; -1/vector_y];
    end
    
    % Specify the same orientation
    if normal_vector_0(1)>0 && normal_vector_0(2)>0
        k = 1;
    elseif normal_vector_0(1)>0 && normal_vector_0(2)<0
        k = 1;
    elseif normal_vector_0(1)<0 && normal_vector_0(2)>0
        k = -1;
    elseif normal_vector_0(1)<0 && normal_vector_0(2)<0
        k = -1;
    else
        k = 1;
    end
    
    normal_vector(:,i) = k * normal_vector_0/norm(normal_vector_0,2);
end
end

