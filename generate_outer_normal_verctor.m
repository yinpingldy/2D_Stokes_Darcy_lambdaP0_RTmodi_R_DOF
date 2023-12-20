function [normal_vector] = generate_outer_normal_verctor (vertices)
%%
% Not the outer normal,ensure that the normal direction is continuous, and the same edge has the same direction.
%%
normal_vector = zeros(2,3);
vertices_0 = vertices;
vertices_0(:,4) = vertices(:,1);

for i = 1:3
    endpoint = vertices_0(:,i:i+1);
    
    vertices_inside_node = setdiff(vertices',endpoint','rows');
    vector_x = endpoint(1,1)-endpoint(1,2);
    vector_y = endpoint(2,1)-endpoint(2,2);
    if vector_x == 0
        if vertices_inside_node(1) < endpoint(1,1)
            normal_vector(:,i) = [1,0];
        else
            normal_vector(:,i) = [-1,0];
        end
    elseif vector_y == 0
        if vertices_inside_node(2) < endpoint(2,1)
            normal_vector(:,i) = [0,1];
        else
            normal_vector(:,i) = [0,-1];
        end
    else
        normal_vector_0 = [1/vector_x; -1/vector_y];
        normal_vector(:,i) = normal_vector_0/norm(normal_vector_0,2);
        k = vector_y / vector_x;
        b = endpoint(2,1) - k*endpoint(1,1);
        if vertices_inside_node(2) < k*vertices_inside_node(1)+b
            if normal_vector(2,i) < 0
                normal_vector(:,i) = -normal_vector(:,i);
            end
        elseif vertices_inside_node(2) > k*vertices_inside_node(1)+b
            if normal_vector(2,i) > 0
                normal_vector(:,i) = -normal_vector(:,i);
            end
        end
    end
end
end

