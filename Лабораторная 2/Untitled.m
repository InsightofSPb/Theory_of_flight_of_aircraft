clc;
clear all
n(1)=0
for i=1:1:6
    x=15;
    n(i+1)=n(i)+1
    Sum(i+1)=(3*40+x*n(i+1))/(3+n(i+1));
end