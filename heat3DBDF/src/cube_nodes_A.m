function To=cube_nodes_A(m,n,l)
j1=(1:n)';
j2=repmat(j1,1,m);
j3=repmat(sort(j2(:)),1,8);
j4=repmat(j3,l,1);
% disp(j4)
clear j1 j2 j3
%-------------
I1=[-1 -1 0 0 -1 -1 0 0];
I2=repmat(I1,l*m*n,1);
% disp(I2)
clear I1
%-------------
P1=[-1 -1 -1 -1 0 0 0 0];
P2=repmat(P1,l*m*n,1);
% disp(P2)
clear P1
%-------------
K1=repmat(repmat([1 2],1,4),l*m*n,1);
% disp(K1)
%-------------
h1=(1:l);
h2=repmat(h1,n*m,1);
h3=h2(:);
h4=repmat(h3,1,8);
% disp(h4)
clear h1 h2 h3
%--------------
M1=ones(1,8)*(m+1);
M2=repmat(M1,l*m*n,1);
% disp(M2)
clear M1
%--------------
MN1=repmat(ones(1,8)*(m+1)*(n+1),l*m*n,1);
% disp(MN1)
%---------------
R=repmat((0:1:(m-1))',l*n,8);
% disp(R)
%---------------
To=(j4+I2).*(M2)+K1+MN1.*(h4+P2)+R;
end