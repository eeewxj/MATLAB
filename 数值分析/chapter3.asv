%arnoldi¹ý³Ì
m=4;
n=4;
b=[1;1;1;1];
x0=zeros(4,1);
h=zeros(n,n);
r=zeros(n,n);
A=[1 0 0 0
   1 1 0 0
   1 1 1 0
   1 1 1 1];
r0=b-A*(x0);
v(:,1)=r0/norm(r0);
for j=1:n
    for i=1:j
        h(i,j)=v(:,i)'*A*v(:,j);
    end
    q=zeros(n,1);
    for i=1:j
        q=q+h(i,j)*v(:,i);
    end
    r(:,j)=A*v(:,j)-q;
    h(j+1,j)=norm(r(:,j));
    v(:,j+1)=r(:,j)/h(j+1,j);
end
v = v(:,1:m);
h = h(1:m,:);
