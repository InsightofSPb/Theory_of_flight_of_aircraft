clc
clear all
n=23;
m0=28000+20*n;
msyx=4000+10*n;
S=2.14+0.001*n;
Sa=2.01+0.001*n;
ms=212+n;
tk1=115;
tk2=(m0-msyx)/ms;
tk=min(tk1,tk2);
t1=5;
t2=80;
Teta0=90;
Tetapr=35+n;
Rpus=500310;
h=0.01;
g=9.81;
%Параметры стандартной атмосферы Земли
step1=0:1:10^5;
H=[0,0.5,1,2,2.5,3,4,5,6,7,8,9,10,11,14,18,24,28,32,36,40,50,60,80,100]*10^3;
T=[288.2,284.9,281.7,275.2,271.9,268.7,262.2,255.7,249.2,242.7,236.2,229.7,223.3,216.8,216.7,216.7,220.6,224.5,228.5,239.3,250.4,270.7,247,198.6,196.6];
P=[101330,95464,89877,79499,74690,70123,61661,54052,47217,41106,35653,30801,26500,22700,14170,7565,2971,1616,889,499,287,80,22,1,3.19*10^(-2)];
f1=interp1(H,T,step1,'spline');
f2=interp1(H,P,step1,'spline');
%1 участок. Вертикальный
i=0;
V(1)=0; x(1)=0; y(1)=0; t(1)=0;
teta(1)=deg2rad(Teta0);
while t(end)<t1
    i=i+1;
    T(i)=f1(find(step1==round(y(i))));
    P(i)=f2(find(step1==round(y(i))));
    a(i)=sqrt(1.4*287*T(i));
    M=V(i)/a(i);
    if M>=0 & M<=0.8
        cx=0.29;
    elseif 0.8<=M & M<=1.068
        cx=M-0.51;
    else
        cx=0.091+0.5*M^(-1);
    end
    ro(i)=1.22*exp(-y(i)/9800);
    R(i)=Rpus-Sa*P(i);
    
    dx(i)=h*V(i)*cos(teta(i));
    x(i+1)=x(i)+dx(i);
    
    dy(i)=h*V(i)*sin(teta(i));
    y(i+1)=y(i)+dy(i);
    
    dV(i)=h*(R(i)/(m0-ms*t(i))-(cx*ro(i)*S*V(i)^2)/(2*(m0-ms*t(i)))-g*sin(teta(i)));
    V(i+1)=V(i)+dV(i);
    
    t(i+1)=t(i)+h;
    teta(i+1)=teta(i);
end
%2 участок. Участок "завала" и "разворота" (криволинейный участок)
A=[1 t1 t1^2 t1^3;
    0 1 2*t1 3*t1^2;
    1 t2 t2^2 t2^3;
    0 1 2*t2 3*t2^2];
B=[deg2rad(Teta0);0;deg2rad(Tetapr);0];
C=A\B;
while t<=t2
    i=i+1;
    teta(i)=C(1)+C(2)*t(i)+C(3)*t(i)^2+C(4)*t(i)^3;
    T(i)=f1(find(step1==round(y(i))));
    P(i)=f2(find(step1==round(y(i))));
    a(i)=sqrt(1.4*287*T(i));
    M=V(i)/a(i);
    if M>=0 & M<=0.8
        cx=0.29;
    elseif 0.8<=M & M<=1.068
        cx=M-0.51;
    else
        cx=0.091+0.5*M^(-1);
    end
    ro(i)=1.22*exp(-y(i)/9800);
    R(i)=Rpus-Sa*P(i);
   
    dx(i)=h*V(i)*cos(teta(i));
    x(i+1)=x(i)+dx(i);
    
    dy(i)=h*V(i)*sin(teta(i));
    y(i+1)=y(i)+dy(i);
    
    dV(i)=h*(R(i)/(m0-ms*t(i))-(cx*ro(i)*S*V(i)^2)/...
            (2*(m0-ms*t(i)))-g*sin(teta(i)));
    V(i+1)=V(i)+dV(i);
    
    t(i+1)=t(i)+h;
end
teta(end+1)=teta(end);
%3 участок. Участок "наведения", прямолинейный участок
while t<tk
    i=i+1;
    T(i)=f1(find(step1==round(y(i))));
    P(i)=f2(find(step1==round(y(i))));
    a(i)=sqrt(1.4*287*T(i));
    M=V(i)/a(i);
    if M>=0 & M<=0.8
        cx=0.29;
    elseif 0.8<=M & M<=1.068
        cx=M-0.51;
    else
        cx=0.091+0.5*M^(-1);
    end
    ro(i)=1.22*exp(-y(i)/9800);
    R(i)=Rpus-Sa*P(i);
    
    dx(i)=h*V(i)*cos(teta(i));
    x(i+1)=x(i)+dx(i);
    
    dy(i)=h*V(i)*sin(teta(i));
    y(i+1)=y(i)+dy(i);
    
    dV(i)=h*(R(i)/(m0-ms*t(i))-(cx*ro(i)*S*V(i)^2)/(2*(m0-ms*t(i)))-g*sin(teta(i)));
    V(i+1)=V(i)+dV(i);
    
    t(i+1)=t(i)+h;
    teta(i+1)=teta(i);
end
%Графики
figure(1)
plot(t,rad2deg(teta))
grid on
axis([0 tk+10 Tetapr-10 Teta0+10])
hold on
line([t1 t1],[Tetapr-10 Teta0+10],'linestyle','--')
hold on
line([t2 t2],[Tetapr-10 Teta0+10],'linestyle','--')
hold on
line([tk tk],[Tetapr-10 Teta0+10],'linestyle','--')
hold off
xlabel('t, с')
ylabel('\theta, градусы')
title('Зависимость угла возвышения от времени')
delete Traectory1.jpg
saveas(figure(1),'Traectory1.jpg');

figure(2)
plot(x,y)
grid on
axis ([0 27000 0 45000]);
xlabel('x, м')
ylabel('y, м')
title({'Траектория'; 'активного участка БР'})
delete Traectory2.jpg
saveas(figure(2),'Traectory2.jpg');

t1=t(1:100:end);
t2=t(end);
t=[t1 t2];
y1=y(1:100:end);
y2=y(end);
y=[y1 y2];
x1=x(1:100:end);
x2=x(end);
x=[x1 x2];
teta1=teta(1:100:end);
teta2=teta(end);
teta=[teta1 teta2];
V1=V(1:100:end);
V2=V(end);
V=[V1 V2];
delete trajectory_RB.xlsx
filename='trajectory_RB.xlsx';
xlswrite(filename,t','sheet1','A2');
xlswrite(filename,y','sheet1','B2');
xlswrite(filename,x','sheet1','C2');
xlswrite(filename,rad2deg(teta)','sheet1','D2');
xlswrite(filename,V','sheet1','E2');