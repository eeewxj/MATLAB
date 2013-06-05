clear all;
clc;
% 常量定义
SIZE = 100;
EPS = 1e-8;
TIMES = 99;
NUM1 = 200;
NUM2 = 100;
% 构造5个100阶对角阵D(:,:,1) -- D(:,:,5)
D = zeros(SIZE,SIZE,5);
D(:,:,1) = diag([100,80*rand(1,SIZE-2)+10,1]);%中间特征值分布范围大
D(:,:,2) = diag([100,42*rand(1,SIZE-2)+10,1]);%中间特征值偏向最小值，分布范围大
D(:,:,3) = diag([100,42*rand(1,SIZE-2)+50,1]);%中间特征值偏向最大值，分布范围大
D(:,:,4) = diag([100,8*rand(1,SIZE-2)+80,1]);%中间特征值偏向最大值，分布范围小
D(:,:,5) = diag([100,8*rand(1,SIZE-2)+10,1]);%中间特征值偏向最小值，分布范围小

% 随机构造10个100阶方阵M(:,:,1) --M(:,:,10)
% 并将M(:,:,i)QR分解为Q(:,:,i)R(:,:,i)
M = zeros(SIZE,SIZE,10);
Q = zeros(SIZE,SIZE,10);
R = zeros(SIZE,SIZE,10);
for i = 1:10
   M(:,:,i) = TIMES*rand(SIZE,SIZE);
   [Q(:,:,i),R(:,:,i)] = qr(M(:,:,i));
end
%计算得到A(:,:,i,j)
A = zeros(SIZE,SIZE,5,10);
for i = 1:5
    for j = 1:10
        A(:,:,i,j) = Q(:,:,j)*D(:,:,i)*(Q(:,:,j).');        
    end
end
%随机生成一个向量b，计算得到精确解 X(:,:,i,j) = inv(A(:,:,i,j))*b
b = TIMES*rand(SIZE,1);
X = zeros(SIZE,1,5,10);
for i = 1:5
    for j = 1:10
        X(:,:,i,j) = Q(:,:,j)*(inv(D(:,:,i)))*(Q(:,:,j).')*b;        
    end
end

%最速下降法解A(:,:,i,j)*x1=b
x1 = zeros(SIZE,1,5,10); 
e1 = zeros(SIZE,NUM1,5,10);
for i = 1:5
    for j = 1:10
        r = b;
        k = 1;
        e1(:,1,i,j) = x1(:,:,i,j) - X(:,:,i,j);
        while norm(r)>EPS
            k = k+1;
            x1(:,:,i,j) = x1(:,:,i,j) + (r'*r)/(r'*(A(:,:,i,j)*r))*r; 
            if k < NUM1+1
                e1(:,k,i,j) = x1(:,:,i,j) - X(:,:,i,j);
            end
            r = b - A(:,:,i,j)*x1(:,:,i,j);    
        end
    end
end
%g共轭梯度法解A(:,:,i,j)*x2=b
x2 = zeros(SIZE,1,5,10); 
e2 = zeros(SIZE,NUM2,5,10);
for i = 1:5
    for j = 1:10
        rtmp = b;    
        p = rtmp;
        e2(:,1,i,j) = x2(:,:,i,j) - X(:,:,i,j);
        x2(:,:,i,j) =  x2(:,:,i,j) + (rtmp'*p)/(p'*A(:,:,i,j)*p)*p;
        e2(:,2,i,j) = x2(:,:,i,j) - X(:,:,i,j);
        r = rtmp - (rtmp'*p)/(p'*A(:,:,i,j)*p)*A(:,:,i,j)*p;
        k = 2;
        while norm(r)>EPS
            k = k+1;
            p = r + p*(r'*r)/(rtmp'*rtmp);
            x2(:,:,i,j) = x2(:,:,i,j) + (rtmp'*p)/(p'*(A(:,:,i,j)*p))*p;
            if k<NUM2+1
                e2(:,k,i,j) = x2(:,:,i,j) - X(:,:,i,j);
            end
            rtmp = r;
            r = rtmp - (rtmp'*p)/(p'*A(:,:,i,j)*p)*A(:,:,i,j)*p;    
        end
    end
end
%计算最速下降法收敛曲线
t1 = zeros(1,NUM1-1,5,10);
for i = 1:5
    for j = 1:10
        for k = 1:NUM1-1
            t1(:,k,i,j) = ((e1(:,k+1,i,j)'*A(:,:,i,j)*e1(:,k+1,i,j))^0.5)/...
            ((e1(:,k,i,j)'*A(:,:,i,j)*e1(:,k,i,j))^0.5);
        end
    end
end
%计算共轭梯度法收敛曲线
t2 = zeros(1,NUM2-1,5,10);
for i = 1:5
    for j = 1:10
        for k = 1:NUM2-1
            if e2(:,k+1,i,j) ~= 0
                t2(:,k,i,j) = ((e2(:,k+1,i,j)'*A(:,:,i,j)*e2(:,k+1,i,j))^0.5)/...
                ((e2(:,k,i,j)'*A(:,:,i,j)*e2(:,k,i,j))^0.5);
            end
        end
    end
end
%收敛曲线图
for i = 1:5
    figure;
    subplot(2,2,1),plot(t1(:,:,i,1));
    subplot(2,2,2),plot(t2(:,:,i,1));
    subplot(2,2,3:4),plot(D(:,:,i),ones(1,100),'*');
end
for i = 1:5
    fprintf('对A[%i][1]，收敛速率均值：%e\n',i,sum(t1(:,:,i,1))/size(t1(:,:,i,1),2));
    fprintf('对A[%i][1]，最速下降法计算%i次的差值的2范数：%e\n',i,NUM1,(e1(:,NUM1,i,1)'*e1(:,NUM1,i,1))^0.5);
end
for i = 1:5
    fprintf('对A[%i][1]，收敛速率均值：%e\n',i,sum(t2(:,:,i,1))/size(t2(:,:,i,1),2));
    fprintf('对A[%i][1]，共轭梯度法计算%i次的差值的2范数：%e\n',i,NUM2,(e2(:,NUM2,i,1)'*e2(:,NUM2,i,1))^0.5);
end
