function [U, S, V] = jacobi_svd(A)
B = A'*A;
[S, V]=jacobi(B);
S = abs(diag(S)).^0.5;
% 奇异值按递减排序，对应的V中的特征向量也要重排序
[~, in] = sort(S,'descend');
S = diag(S(in));
V = V(:,in);
U = A*V*diag(1./(diag(S)))';
end

function [matrix,J] = jacobi(matrix)
    THRESHOLD = 1e-8;
    ITERATION = 1e8; 
    J = eye(size(matrix,1));
    iteration = 0;
    while iteration <ITERATION
        pass = true;
%方法1：遍历非对角元，对每一对非对角元做jacobi旋转变换
%         for i = 1:size(matrix,1)
%             for j = i+1:size(matrix,1)
%                 [matrix,J,pass] = rotate(matrix, i, j, pass, J);
%             end
%         end
%方法2：仅搜索得到绝对值最大的一对非对角元，对其做jacobi旋转变换
        [x,y]=max(abs(matrix-diag(diag(matrix))));
        [~,j]=max(x);
        i = y(j);
        [matrix,J,pass] = rotate(matrix, i, j, pass, J);
         %当非对角元素全部小于阈值时停止迭代
        if pass  
            break;
        end
        %将绝对值小于阈值的非对角元置零
        matrix = matrix - diag(diag(matrix).*(diag(matrix)<THRESHOLD));
        iteration = iteration + 1; 
    end

end

function [a, J, pass] = rotate(a, i, j, pass, J)
%计算Jacobi旋转矩阵J，以及J'*B*J
    THRESHOLD = 1e-8;
    if abs(a(i,j)) < THRESHOLD
        return;
    end
    pass = false;
    c = (a(i,i) - a(j,j)) / (2 * a(i,j));
    t = sign((c>=0)-1/2) / (abs(c) + sqrt(1 + c^2));
    Cos = 1 / sqrt(1 + t^2);
    Sin = Cos * t;
    G = eye(size(a,1));
    G([i j],[i j])=[Cos -Sin;Sin Cos];
    a = G' * a * G;
    J = J * G;
end


  

