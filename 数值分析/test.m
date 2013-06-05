clear all;
clc;
x1=0;
x2=1.5;
iratations = 0;
while  max(abs(x1-1),abs(x2-1))> 1e-4
    tmp = (x1^2+x2^2+8)/10;
    x2 = (x1*x2^2+x1+8)/10;
    x1 = tmp;
    iratations = iratations + 1;
end
disp(iratations);