clear all
 
tic 
xmatlab1 = fsolve(@NonlinearExample,[1;1;1]);
toc 

 options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true);
 tic
 [xmatlab2,F,exitflag,output,JAC] = fsolve(@NonlinearExampleWithJacobian,[1;1;1],options);
 toc
 
 
 x=[1;1;1]
 [f,J]=NonlinearExampleWithJacobian(x);
 while max(abs(f)) > 1.0e-8
    [L,U]=LUDecomposition(J);
    R=LowerTriangularSolver(L,-f);
    diff=UpperTriangularSolver(U,R);
    x=diff+x
    [f,J]=NonlinearExampleWithJacobian(x);
 end
  
 
 
 function [f,J] = NonlinearExampleWithJacobian(x)
f(1) = 3*x(1)+cos(x(2)*x(3))-1/2;
f(2) = x(1)^2-81*(x(2)+0.1)^2 +sin(x(3)) + 1.06;
f(3)=exp(-x(1)*x(2)) +20*x(3)+(10*pi-3)/3;

J(1,1)=3;J(1,2)=-x(3)*sin(x(2)*x(3));J(1,3)=-x(2)*sin(x(2)*x(3));
J(2,1)=2*x(1);J(2,2)=-16.2-162*x(2);J(2,3)=cos(x(3));
J(3,1)=-x(2)*exp(-x(1)*x(2));J(3,2)=-x(1)*exp(-x(1)*x(2));J(3,3)=20;
end
 
function f = NonlinearExample(x)
f(1) = 3*x(1)+cos(x(2)*x(3))-1/2;
f(2) = x(1)^2-81*(x(2)+0.1)^2 +sin(x(3)) + 1.06;
f(3)=exp(-x(1)*x(2)) +20*x(3)+(10*pi-3)/3
end
 
 
function [L,U]=LUDecomposition(A)

[n,n]=size(A);
L=zeros(n,n);
U=zeros(n,n);

for i=1:n
    L(i,1)=A(i,1);
end
for j=2:n
    U(1,j)=A(1,j)/L(1,1);
end

for j=2:n
    for i=j:n
        sum=0;
        for k=1:j-1
            sum=sum+L(i,k)*U(k,j);
        end
        L(i,j)=A(i,j)-sum;
    end
        
    for k=j+1:n
        sum=0;
        for i=1:j-1
            sum=sum+L(j,i)*U(i,k);
        end
        U(j,k)=(A(j,k)-sum)/L(j,j);
    end
end

for i=1:n
   U(i,i)=1.0; 
end
    
end

function x=UpperTriangularSolver(A,C)
[n,n]=size(A); % Finding the size of the problem
x=zeros(n,1); % set up  x vector with elements initially set to zero
x(n)=C(n)/A(n,n);
 for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+A(i,j)*x(j);
    end
    x(i)=(C(i)-sum)/A(i,i);
end
end

function x=LowerTriangularSolver(A,C)
[n,n]=size(A); % Finding the size of the problem
x=zeros(n,1); % set up  x vector with elements initially set to zero
x(1)=C(1)/A(1,1);
for i=2:n
    sum=0;
    for j=1:i-1
        sum=sum+A(i,j)*x(j);
    end
    x(i)=(C(i)-sum)/A(i,i);
end
end
        