clc;
clear all;
%参数定义
n=20; M=20;       
b=ones(n,1);
%定义A1
A = zeros(20,20);
for i =1:20
    for j=1:20
        A(i,j)=1/(i+j-1);
    end
end  
%定义A2
%A =diag(2*ones(n,1),0)+diag(-1*ones(n-1,1),-1)+diag(ones(n-1,1),1); 
%定义A3
% Au = diag(2*ones(n/2,1),0)+diag(ones(n/2-1,1),-1)+diag(-1*ones(n/2-1,1),1);
% Ad = diag(-2*ones(n/2,1),0)+diag(-1*ones(n/2-1,1),-1)+diag(ones(n/2-1,1),1);
% A =[Au,zeros(n/2,n/2);zeros(n/2,n/2),Ad];
x0 = [1;zeros(n-1,1)];
r0=b-A*x0;
r = zeros(n,M);
%利用arnoldi过程计算H与V
H=zeros(M+1,M);
V =zeros(n,M); 
V(:,1)=r0/norm(r0);
for j=1:M
    H(1:j,j)=V(:,1:j)'*(A*V(:,j));
    q=zeros(n,1)+ V(:,1:j)*H(1:j,j);
    r(:,j)=A*V(:,j)-q;
    H(j+1,j)=norm(r(:,j));
    V(:,j+1)=r(:,j)/H(j+1,j);    
end
V = V(:,1:M);
%GMRES
u = zeros(M,1);
p = zeros(M,1);
for m = 1:M
    v=V(:,1:m);
    h=H(1:m+1,1:m);
    g=[norm(r0);zeros(m,1)];
    for i=1:m
        W=eye(m+1);
        c=1/sqrt(1+(h(i+1,i)/h(i,i))^2);
        s=h(i+1,i)/h(i,i)*c;
        W(i,i:i+1)=[c,s];
        W(i+1,i:i+1)=[-s,c];
        h=W*h;
        g=W*g;
    end
    y=h(1:m,1:m)\g(1:m,:);
    x=x0 + v*y;
    u(m,1)=norm(b-A*x); %计算误差
    [L,U]=lu(H(1:m,1:m));
    p(m,1)=norm(V(:,1:m)/U);
end
figure;
subplot(1,2,1);semilogy(u,'k.');title('余量模变化');
xlabel('m');ylabel('||r||');
subplot(1,2,2);plot(p,'k.');title('Pm模变化');
xlabel('m');ylabel('||Pm||');
        