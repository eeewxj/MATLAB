clear all;
clc;
cycle = 1e3;
T = zeros(1,cycle);
for i = 1:cycle
    N = randi(100,1,1e6);
    con = zeros(1,100);
    j=0;
    while sum(con)<100
        j=j+1;
        if con(N(j))==1
        else
           con(N(j)) = 1;
         end
    end
    T(i) = j;    
end
disp(sum(T)/cycle);
f = zeros(1,1e6);
for i = 1:1e6
    f(i) = sum(T<i)/cycle;
end
semilogx(f(1:1e6));