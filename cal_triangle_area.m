function triangle_area = cal_triangle_area(number_of_elements,P_partition,T_partition)
%%
triangle_area=zeros(1,number_of_elements);

for i=1:number_of_elements
    tri=T_partition(:,i);
    p1=P_partition(:,tri(1));
    p2=P_partition(:,tri(2));
    p3=P_partition(:,tri(3));
    A=[[p1,p2,p3];[1,1,1]];
    triangle_area(i)=det(A)/2;  % 单元面积信息*2
end
end

