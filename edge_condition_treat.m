function edge_condition=edge_condition_treat(edges,triangle_edge)
N=size(edges,1);
edge_condition=zeros(2,N);
N=size(triangle_edge,2);
for i =1:N
    for j=1:3
        num=triangle_edge(j,i);
        if edge_condition(1,num)~=0
            edge_condition(2,num)=i;
        else
            edge_condition(1,num)=i;
        end
    end
end