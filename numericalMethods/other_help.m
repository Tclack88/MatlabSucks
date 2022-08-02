% check if iterative method converges
A3 = [1 2 3; 3 2 2; 2 3 2];
% triangular matrix gauss-seidel
M3 = tril(A3);
N3 = M3-A3;
abs(eig(inv(M3)*N3));
% diagonal matrix point jacobi
M3 = diag(diag(A3));
N3 = M3-A3;
abs(eig(inv(M3)*N3));


% nonlinear system of equations; jacobian; Lecture07M05.m

 x=[1;1;1]
 [f,J]=NonlinearExampleWithJacobian(x);
 while max(abs(f)) > 1.0e-8
    [f,J]=NonlinearExampleWithJacobian(x);
    [L,U]=LUDecomposition(J);
    R=LowerTriangularSolver(L,-f);
    diff=UpperTriangularSolver(U,R);
    x=diff+x
 end

function f = NonlinearExample(x)
f(1) = 3*x(1)+cos(x(2)*x(3))-1/2;
f(2) = x(1)^2-81*(x(2)+0.1)^2 +sin(x(3)) + 1.06;
f(3)=exp(-x(1)*x(2)) +20*x(3)+(10*pi-3)/3
end

function [f,J] = NonlinearExampleWithJacobian(x)
f(1) = 3*x(1)+cos(x(2)*x(3))-1/2;
f(2) = x(1)^2-81*(x(2)+0.1)^2 +sin(x(3)) + 1.06;
f(3)=exp(-x(1)*x(2)) +20*x(3)+(10*pi-3)/3;

J(1,1)=3;J(1,2)=-x(3)*sin(x(2)*x(3));J(1,3)=-x(2)*sin(x(2)*x(3)); 
J(2,1)=2*x(1);J(2,2)=-16.2-162*x(2);J(2,3)=cos(x(3)); 
J(3,1)=-x(2)*exp(-x(1)*x(2));J(3,2)=-x(1)*exp(-x(1)*x(2));J(3,3)=20; 
end 




% polynomial least squares regression ; least square
% this is a quadratic fit
xdata=[0.075 0.5 1.0 1.2 1.7 2.0 2.3]
ydata=[600 800 1200 1400 2050 2650 3750]
N=length(xdata);
A=[N sum(xdata) sum(xdata.^2) ;
   sum(xdata) sum(xdata.^2) sum(xdata.^3);
   sum(xdata.^2) sum(xdata.^3) sum(xdata.^4)];
C=[sum(ydata);
    sum(xdata.*ydata);
    sum(xdata.^2.*ydata)];
%
%Note that for brevity I have used the MATLAB linsolve() function
%solve the resulting equation.  You can use any of the
%solvers from previous lectures.
%
a=linsolve(A,C);
x=linspace(0,2.5);
y=a(1)+a(2)*x+a(3)*x.^2;
hold off



% y = p^x fit; exponential fit ; lecture08
% Lecture08M04.m
x=[0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0]
y=[6.00 4.83 3.70 3.15 2.41 1.83 1.49 1.21 0.96 0.73 0.64]
xdata=x
ydata=log(y)
A=[sum(xdata.^2) sum(xdata);
    sum(xdata) length(xdata)];
C=[sum(xdata.*ydata);
    sum(ydata)];
vec=inv(A)*C;
a=vec(1); %this should be m
b=vec(2); % this should be ln(p)
a=m;
p=exp(b);
xplot=linspace(0,5);
yplot=p*exp(m*xplot);
hold off
plot(x,y,'ko','MarkerSize',15,'MarkerFaceColor','g')
hold on
plot(xplot,yplot,'k-','LineWidth',2)
xlabel('x');
ylabel('y');


% differentiation ; derivative ; finite difference ; CDS ; FDS ; BDS ; differencing
% lecture 13
MyFunc=@(x) x.^2.*cos(x)
MyDeriv=@(x)x.*(2.*cos(x)-x.*sin(x))

% Define the vector of Delta values
DeltaVec=[0.001 0.002 0.005 0.01 0.1 0.2 0.4 1.0] ;

%Defining data point
xi=2.0;
for k=1:length(DeltaVec)
    Delta=DeltaVec(k)
    xim2=xi-2*Delta
    xim1=xi-Delta;
    xip1=xi+Delta;
    xip2=xi+2*Delta;

    %2nd order centred
    dfdx2ndCentered=(1./(2*Delta))*(-MyFunc(xim1)+MyFunc(xip1));
    error2ndCentered(k)=abs(MyDeriv(xi)-dfdx2ndCentered);

    %4th order centred
    dfdx4thCentered=(1./(12*Delta))*(MyFunc(xim2)-8*MyFunc(xim1)+8*MyFunc(xip1)-MyFunc(xip2));
    error4thCentered(k)=abs(MyDeriv(xi)-dfdx4thCentered);
end

hold off
plot(log10(DeltaVec),log10(error2ndCentered),'ko-')
hold on
plot(log10(DeltaVec),log10(error4thCentered),'go-')
xlabel('log(\Delta)','FontSize',14)
ylabel('log(Error)','FontSize',14)


% derivative unequally spaced data ; lecture13
clear all

xveryfine=linspace(0,5,1000);
fxveryfine=exp(-xveryfine.^2).*sin(10*xveryfine);
dfactualveryfine=exp(-xveryfine.^2).*(-2*xveryfine.*sin(10*xveryfine)+10*cos(10*xveryfine))

x=GenGrid(0,5,30,1.0);
fx=exp(-x.^2).*sin(10*x);
df=DerivCDA(fx,x);

figure(1)
subplot(211)
hold off
plot(x,fx,'bo-',xveryfine,fxveryfine,'k-')
xlabel('x');ylabel('f(x)')

subplot(212)
hold off
plot(x,df,'bo-',xveryfine,dfactualveryfine,'k-')
xlabel('x');ylabel('df(x)/dx')

function df=DerivCDA(f,x)
df=zeros(size(f))
N=length(f)
df(1)=(f(2)-f(1))/(x(2)-x(1))
for i=2:N-1
    df(i)=(f(i+1)-f(i-1))/(x(i+1)-x(i-1));
end
df(N)=(f(N)-f(N-1))/(x(N)-x(N-1))
end

function x=GenGrid(a,b,N,r)
x=linspace(a,b,N);
if abs(r-1) > 1.0e-8
    dx0=(b-a)*(1-r)/(1-r.^(N-1));
    x(1)=a;
    dx=dx0;
    for i=2:N
        x(i)=x(i-1)+dx;
        dx=r*dx;
    end
end
end


% boundary value problem; finite difference; equally spaced data;  
% lecture 15M01.m
% original interval 1-2. BC y(1) = 5 y(2) = 3; dirichlet; 2nd order central difference approximation
clear all
clear all;

Delta=0.2;
x=1:Delta:2;
alpha=5.0;
beta=3.0;

n=length(x)-1;
A=zeros(n-1,n-1);
C=zeros(n-1,1);
ysol=zeros(size(x));

A(1,1)=-(2/Delta.^2)-2./(x(2).^2);
A(1,2)=(1/Delta.^2)+(2./x(2))*(1/(2*Delta));
for i=2:n-2
    A(i,i-1)=(1/Delta.^2)-(2./x(i+1))*(1/(2*Delta));
    A(i,i)=-(2/Delta.^2)-2./(x(i+1).^2);
    A(i,i+1)=(1/Delta.^2)+(2./x(i+1))*(1/(2*Delta));
end
A(n-1,n-2)=(1/Delta.^2)-(2./x(n))*(1/(2*Delta));
A(n-1,n-1)=-(2/Delta.^2)-2./(x(n).^2);

C(1)=-alpha*((1/Delta.^2)-(2./x(2))*(1/(2*Delta)));
C(n-1)=-beta*((1/Delta.^2)+(2./x(n))*(1/(2*Delta)));

%
% This is a very inefficient way of solving this system of equations
% Matrix [A] has got lots of zeros
% Should use Thomas algorithm and not store all the zeros.
%
yinner=A\C;

ysol(1)=alpha;
ysol(2:n)=yinner;
ysol(n+1)=beta;

hold off;
plot(x,ysol,'bo-','Markersize',15,'Markerfacecolor',[0 0 1]);
hold on
plot(x,x+4./x.^2,'k-','linewidth',2);

xlabel('x','FontSize',14)
ylabel('y(x)','FontSize',14)
