function h_T = diameter(vertices)
A = vertices(:,1);
B = vertices(:,2);
C = vertices(:,3);

h_T1 = norm(A-B);
h_T2 = norm(A-C);
h_T3 = norm(B-C);
h_T = max(h_T1 , max(h_T2,h_T3));
end

