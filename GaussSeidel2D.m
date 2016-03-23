function [ ite, r ] = GaussSeidel2D( F,M,Tol,N )
% Gauss Seidel iterative method using Red and Black Ordering

h=1/(N+1);
U1=zeros(N+2,N+2);
U2=zeros(N+2,N+2);
U3=zeros(N+2,N+2);
V=zeros(N+2,N+2);


r=zeros(3,1);
ite=zeros(3,1);


% 1-Norm
for k=1:M
%Red Points
%Found using 2 Loops and counting by 2's. The first loop fills in every
%other row and every other column. The second loop goes back to the skipped
%rows and fills them in for every other column, but it starts at the red
%node in column 3, since column two should have a black node. i's are rows,
%j's are columns
i=2;   
while i<=N+1
    j=2;
    while j<=N+1
        U1(i,j)=[U1(i-1,j)+U1(i+1,j)+U1(i,j-1)+U1(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
i=3;
while i<=N+1
    j=3;
    while j<=N+1
        U1(i,j)=[U1(i-1,j)+U1(i+1,j)+U1(i,j-1)+U1(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
%Black Points
%These are found similarly to the red nodes, except we start at the 3rd
%column since it is a black node, and on the second loop we start at the
%first column since it is a black node as well. 
i=2;   
while i<=N+1
    j=3;
    while j<=N+1
        U1(i,j)=[U1(i-1,j)+U1(i+1,j)+U1(i,j-1)+U1(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
i=3;
while i<=N+1
    j=2;
    while j<=N+1
        U1(i,j)=[U1(i-1,j)+U1(i+1,j)+U1(i,j-1)+U1(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end

%The 1-Norm is used to set tolerance
  r(1)= norm(U1-V);
    if r(1)<Tol
       break
    end
    ite(1)=k;
    V=U1; %I keep a copy of the last value to use in the norm calculations

end


V=zeros(N+2,N+2);

% 2-Norm
for k=1:M
%Red Nodes
i=2;   
while i<=N+1
    j=2;
    while j<=N+1
        U2(i,j)=[U2(i-1,j)+U2(i+1,j)+U2(i,j-1)+U2(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
i=3;
while i<=N+1
    j=3;
    while j<=N+1
        U2(i,j)=[U2(i-1,j)+U2(i+1,j)+U2(i,j-1)+U2(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
%Black Points

i=2;   
while i<=N+1
    j=3;
    while j<=N+1
        U2(i,j)=[U2(i-1,j)+U2(i+1,j)+U2(i,j-1)+U2(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
i=3;
while i<=N+1
    j=2;
    while j<=N+1
        U2(i,j)=[U2(i-1,j)+U2(i+1,j)+U2(i,j-1)+U2(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end

%The 2-Norm is used to set tolerance
  r(2)=norm(U2-V);
    if r(2)<Tol
       break
    end
    ite(2)=k;
    V=U2; %I keep a copy of the last value to use in the norm calculations

end

% Inf-Norm
for k=1:M
%Red Nodes
i=2;   
while i<=N+1
    j=2;
    while j<=N+1
        U3(i,j)=[U3(i-1,j)+U3(i+1,j)+U3(i,j-1)+U3(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
i=3;
while i<=N+1
    j=3;
    while j<=N+1
        U3(i,j)=[U3(i-1,j)+U3(i+1,j)+U3(i,j-1)+U3(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
%Black Points

i=2;   
while i<=N+1
    j=3;
    while j<=N+1
        U3(i,j)=[U3(i-1,j)+U3(i+1,j)+U3(i,j-1)+U3(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end
i=3;
while i<=N+1
    j=2;
    while j<=N+1
        U3(i,j)=[U3(i-1,j)+U3(i+1,j)+U3(i,j-1)+U3(i,j+1)+(h^2)*F(i,j)]/4;
        j=j+2;
    end
    i=i+2;
end

%The 2-Norm is used to set tolerance
 r(3)=norm(U3-V,inf);
    if r(3)<Tol
       break
    end
    ite(3)=k;
    V=U3; %I keep a copy of the last value to use in the norm calculations

end

end

