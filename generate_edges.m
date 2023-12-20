function [edges, triangle_edge] = generate_edges(number_of_elements,T_partition)
% 首先找到每一个单元的边,并起来去重
E1=zeros(number_of_elements,2);
E2=zeros(number_of_elements,2);
E3=zeros(number_of_elements,2);
for i=1:number_of_elements
    E1(i,1)=T_partition(1,i);
    E1(i,2)=T_partition(2,i);
    E2(i,1)=T_partition(2,i);  
    E2(i,2)=T_partition(3,i);  
    E3(i,1)=T_partition(1,i); 
    E3(i,2)=T_partition(3,i);
end
Sum_E = [E1;E2;E3];
Sum_E = sort(Sum_E,2); % 对每行中的元素进行排序
edges = unique(Sum_E,'rows'); % 单元边的集合

triangle_edge=zeros(3,number_of_elements);
for i=1:number_of_elements
    E11=sort(E1(i,:),2);
    E22=sort(E2(i,:),2);
    E33=sort(E3(i,:),2);
    [~,~,ib1]=intersect(E11,edges,'rows');
    [~,~,ib2]=intersect(E22,edges,'rows');
    [~,~,ib3]=intersect(E33,edges,'rows');   
    triangle_edge(:,i)=[ib2;ib3;ib1];   % 每个单元的边在全部边中的编号
end

end

