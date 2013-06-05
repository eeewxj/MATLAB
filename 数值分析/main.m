clc;clear all;
M = 500;
N = 400;
NUM = 1;
A = 10*rand(M,N,NUM);

U1=zeros(M,N,NUM); U2=zeros(M,M,NUM); U3=zeros(M,N,NUM);
V1=zeros(N,N,NUM); V2=zeros(N,N,NUM); V3=zeros(N,N,NUM);
S1=zeros(N,N,NUM); S2=zeros(M,N,NUM); S3=zeros(N,N,NUM);

tic;
for i=1:NUM
    [U1(:,:,i),S1(:,:,i),V1(:,:,i)] =lansvd(A(:,:,i),N,'L');
end 
toc;

tic;
for i=1:NUM
    [U2(:,:,i),S2(:,:,i),V2(:,:,i)] = svd(A(:,:,i));
end
toc;

tic;
for i=1:NUM
    [U3(:,:,i),S3(:,:,i),V3(:,:,i)] = jacobi_svd(A(:,:,i));
end
 toc;

 norm(U1(:,:,1)*S1(:,:,1)*V1(:,:,1)'-A(:,:,1))/norm(A(:,:,1))
 max(diag(S1(:,:,1)-S2(1:N,:,1)))

 norm(U3(:,:,1)*S3(:,:,1)*V3(:,:,1)'-A(:,:,1))/norm(A(:,:,1))
 max(diag(S3(:,:,1)-S2(1:N,:,1)))
 