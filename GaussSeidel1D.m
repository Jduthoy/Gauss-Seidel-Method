function [x] = GaussSeidel1D(A,b,x0,M,Tol)
%1D Gauss Seidel Iterative Method on Ax=b
[n,m]=size(A);
x=x0;
xx=A\b % answer using matlab command


for k=1:M
    r=norm(x0-xx);
    if r<Tol
        break;
    end
        x(1)=[b(1)-A(1,2:m)*x(2:m)]/A(1,1); %Sets value for first approximation
        for i=1:n-1
            x(i)=[b(i)-(A(i,1:i-1)*x(1:i-1)+A(i,i+1:m)*x(i+1:m))]/A(i,i); 
        end
        x(n)=[b(n)-A(n,1:n-1)*x(1:n-1)]/A(n,n); %Sets value for last
   
end
r
k
end

