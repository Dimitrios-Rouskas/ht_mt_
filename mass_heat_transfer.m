% Askisi 3.1.2 (Convection-Diffusion)
clear; clc;

n = 21;
nel=n-1;

%Define matrices
a = zeros(n,n);
A = zeros (n,n);
b = zeros(n,1);
B = zeros (n,1);
T = zeros(n,1);
C = zeros(n,1);
Tm = zeros(nel,1);
K = zeros(n,1);
z = zeros(nel,2);
ml = zeros(2,2);
ll = zeros(2,2);
cl = zeros(2,2);
%Parameters
Pe=100;
Bi=50;
Da=1000;
bi=10;
T0=1;
C0=1;
qL=0;
CL=0;
L=1;
dz = L/nel;
for i=1:n
    K(i)=(i-1)*dz;
end
%Heat transfer 
%Form connectivity matrix
for ie=1:nel
    z(ie,1)=ie;
    z(ie,2)=ie+1;
end

%Form local element matrices
ll(1,1)=1/dz;
ll(1,2)=-1/dz;
ll(2,1)=-1/dz;
ll(2,2)=1/dz;

ml(1,1)=dz/3;
ml(1,2)=dz/6;
ml(2,1)=dz/6;
ml(2,2)=dz/3;

cl(1,1)=-1;
cl(1,2)=1;
cl(2,1)=-1;
cl(2,2)=1;
%Form linear system
%Local element matrix assembly
for ie=1:nel
    for i=1:2
        i1=z(ie,i);
        for j=1:2
            j1=z(ie,j);
            a(i1,j1)=a(i1,j1)+Pe*cl(i,j)+ll(i,j)+Bi*ml(i,j);
        end
    end
end

%Boundary conditions
%First line (Dirichlet condition)
a(1,:)=0;
a(1,1)=1;
b(1)=T0;

%Last line (Neumann or natural boundary condition)
b(n)=b(n)-qL;

% Solve linear system A.T = b
%metafora mazas
T = a\b;
for j = 1:nel
    Tm(j)=(T(j)+T(j+1))/2;
end
%Form linear system
%Local element matrix assembly
for ie=1:nel
    for i=1:2
        i1=z(ie,i);
        for j=1:2
            j1=z(ie,j);
            A(i1,j1)=A(i1,j1)+Pe*cl(i,j)+ll(i,j)+Da*exp(-bi*Tm(ie))*ml(i,j);
        end
    end
end

%Boundary conditions
%First line (Dirichlet condition)
A(1,:)=0;
A(1,1)=1;
B(1)=C0;

%Last line (Neumann or natural boundary condition)
B(n)=B(n)-CL;
C=A\B;
%
figure(1)
plot(K,T)
xlabel('position (z)');
ylabel('Temprature');
title('Temprature vs Position');
legend('Location', 'best');
grid on;
hold off;
figure(2)
plot(K,C)
xlabel('position (z)');
ylabel('Concentration');
title('Concentration vs Position');
legend('Location', 'best');
grid on;
hold off;
