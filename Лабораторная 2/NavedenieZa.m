clc
clear all
close all
n=23;
%%
x0=0;
y0=0;
xc=0;
yc=500+10*n;
Vc=80;
tetac=0;
xkp=0;
ykp=0;
h=0.01;
m0=170+n;
msyx=50+n;
ms=8.3;
S=0.03+0.001*n;
Cy=0.3;
Cx0=0.3;
A=2;
Snapr=3;
ftr=0.15;
Iud=1200;
R=Iud*ms;
g=9.81;

%%
beta(1)=atan(xc/yc);
tetad=pi/2-beta(1);
td=sqrt(2*Snapr/(R/m0-g*(sin(tetad)+ftr*cos(tetad))));
Vd=(R/m0-g*(sin(tetad)+ftr*cos(tetad)))*td;
xd=x0+Snapr*cos(tetad);
yd=y0+Snapr*sin(tetad);
r(1)=sqrt((xc-xd)^2+(yc-yd)^2);
rastold(1)=r;
alpha(1)=0;

%%
x(1)=xd;
y(1)=yd;
V(1)=Vd;
t(1)=td;
teta(1)=tetad;
m=m0-ms*td;
rc(1)=sqrt((xc-xkp)^2+(yc-ykp)^2);
i=0;

%%
while r(end)<=rastold
    i=i+1;
    rastold(i)=r(i);
    t(i+1)=t(i)+h;
    m=m0-ms*t(i+1);
    if m<msyx
        m=msyx;
        R=0;
    else
    end
    xc(i+1)=xc(i)+Vc*h;
    yc(i+1)=yc(i);
    rc(i+1)=sqrt((xc(i+1)-xkp)^2+(yc(i+1)-ykp)^2);
    beta(i+1)=beta(i)+h*(Vc*cos(beta(i)+tetac)/rc(i+1));
    x(i+1)=x(i)+h*(V(i)*cos(teta(i)));
    y(i+1)=y(i)+h*(V(i)*sin(teta(i)));
    rla=sqrt((x(i+1)-xkp)^2+(y(i+1)-ykp)^2);
    ro(i)=1.22*exp(-y(i+1)/10000);
    V(i+1)=V(i)+h*(R/m-(Cx0+A*alpha(i)^2)*S*ro(i)*V(i)^2/(2*m)-g*sin(teta(i)));
    teta(i+1)=teta(i)+h*(R*alpha(i)/(m*V(i))+Cy*alpha(i)*S*ro(i)*V(i)/(2*m)-g*cos(teta(i))/V(i));
    teta1(i)=acos(Vc*cos(beta(i+1)+tetac)/V(i+1)*rla/rc(i+1))-beta(i+1);
    alpha(i+1)=((teta1(i)-teta(i+1))/h+g*cos(teta(i+1))/V(i+1))/(R/(m*V(i+1))+Cy*S*ro(i)*V(i+1)/(2*m));
    r(i+1)=sqrt((xc(i+1)-x(i+1))^2+(yc(i+1)-y(i+1))^2);
end

if abs(r(end))>abs(r(end-1))
    x(end)=[]; y(end)=[]; r(end)=[];
    xc(end)=[]; yc(end)=[];
else
end
while xc(end)>x(end)
    xc(end)=[]; yc(end)=[];
end
r(end)

%%
time1=round(i/4);
time2=round(i/2);
time3=round(3*i/4);
figure(1)
plot(x,y);
hold on
plot(xc,yc);
hold on
l1=line([xkp xc(time1)],[ykp yc(time1)],'linestyle','-.');
l2=line([xkp xc(time2)],[ykp yc(time2)],'linestyle','-.');
l3=line([xkp xc(time3)],[ykp yc(time3)],'linestyle','-.');
grid on
xlabel('x, м')
ylabel('y, м')
title('Траектория ЛА, если цель удаляется')
legend('Траектория ЛА','Траектория Цели','Линия визирования 1','Линия визирования 2','Линия визирования 3','location','southeast')
saveas(figure(1),'Trajectory1.jpg');

tviz=[t(time1) t(time2) t(time3)];
xviz=[x(time1) x(time2) x(time3)];
yviz=[y(time1) y(time2) y(time3)];
xcviz=[xc(time1) xc(time2) xc(time3)];
ycviz=[yc(time1) yc(time2) yc(time3)];
Vviz=[V(time1) V(time2) V(time3)];

V(11:1:end-10)=[];
t(11:1:end-10)=[];
x(11:1:end-10)=[];
y(11:1:end-10)=[];
xc(11:1:end-10)=[];
yc(11:1:end-10)=[];


%% 
delete NavedenieZa.xlsx
filename='NavedenieZa.xlsx';
%% 
xlswrite(filename,t','sheet1','A2');
xlswrite(filename,x','sheet1','B2');
xlswrite(filename,y','sheet1','C2');
xlswrite(filename,xc','sheet1','D2');
xlswrite(filename,yc','sheet1','E2');
xlswrite(filename,V','sheet1','F2');
xlswrite(filename,tviz','sheet2','A2');
xlswrite(filename,xviz','sheet2','B2');
xlswrite(filename,yviz','sheet2','C2');
xlswrite(filename,xcviz','sheet2','D2');
xlswrite(filename,ycviz','sheet2','E2');
xlswrite(filename,Vviz','sheet2','F2');



      