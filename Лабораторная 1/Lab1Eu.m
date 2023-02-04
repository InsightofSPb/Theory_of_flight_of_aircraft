clc
clear all
%Dano
d=0.150; m0=60; mtop=13;ta=0.7;Lh=5.88; tetad=pi/2; ue=2500; ro0=1.23; 
ftr=0.15; g=9.81; ij=1.15;
%Движение по направляющим
Q=mtop/ta; %расход
R=Q*ue; %тяга
td=sqrt(2*Lh/((R/m0)-g*(sin(tetad)+ftr*cos(tetad)))); %время схода с напр.
Vd=(R/m0-g*(sin(tetad)+ftr*cos(tetad)))*td; %скорость схода с напр.
%Расчёт активного участка траектории ЛА методом Эйлера
S=pi/4*d^2; %характерная площадь
h=0:1:1000;
Vet=[100,150,200,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,500,550,600,650,700,750,800,850,900,950,1000];
cxet=[0.255,0.257,0.260,0.268,0.274,0.280,0.295,0.321,0.361,0.411,0.460,0.500,0.542,0.574,0.603,0.628,0.648,0.665,0.680,0.690,0.700,0.709,0.715,0.719,0.735,0.732,0.721,0.704,0.680,0.665,0.645,0.624,0.605,0.583,0.567];
Cx=cxet*ij; % Сx с учётом несоответсвия эталонному ЛА
gr1=pchip(Vet,Cx,h);
%
figure(1)
plot(Vet,Cx,'x',h,gr1)
grid on
xlabel('Эталонная скорость V,м/с');
ylabel('Реальный Сх для Ла');
%
t(1)=td; %начальное время для активного участка
h1=0.001;
i=0;
V(1)=Vd; %начальная скорость для активного участка
teta(1)=tetad; %начальный угол для активного участка
y(1)=Lh*sin(tetad); %начальная высота для активного участка
x(1)=Lh*cos(tetad); %начальная координата х для активного участка
while t<=ta 
    
    i=i+1;
    %
    dx(i)=h1*(V(i)*cos(teta(i)));
    x(i+1)=x(i)+dx(i);
    %
    dy(i)=h1*(V(i)*sin(teta(i)));
    y(i+1)=y(i)+dy(i);
    %
    Cx(i)=gr1(round(V(i)));
    ro(i)=ro0*exp(-y(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V(i)^2;
    %
    dV(i)=h1*(R/(m0-Q*t(i))-g*sin(teta(i))-X(i)/(m0-Q*t(i)));
    V(i+1)=V(i)+dV(i);
    %
    dteta(i)=(h1/V(i))*(-g)*cos(teta(i));
    teta(i+1)=teta(i)+dteta(i);
    %
    t(i+1)=t(i)+h1;
end
%Расчёт пассивного участка траектории ЛА методом Эйлера
m=m0-mtop;
while y>0
    i=i+1;
    Cx(i)=gr1(round(V(i)));
    ro(i)=ro0*exp(-y(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V(i)^2;
    %
    dx(i)=h1*(V(i)*cos(teta(i)));
    x(i+1)=x(i)+dx(i);
    %
    dy(i)=h1*V(i)*sin(teta(i));
    y(i+1)=y(i)+dy(i);
    %
    dV(i)=h1*(-X(i)/m-m*g*sin(teta(i))/m);
    V(i+1)=V(i)+dV(i);
    %
    dteta(i)=h1/V(i)*(-g*cos(teta(i)));
    teta(i+1)=teta(i)+dteta(i);
    %
    t(i+1)=t(i)+h1;
end
figure(2)
h2=0:10:x(end);
gr2=pchip(x,y,h2);
plot(x,y,h2,gr2)
grid on
xlim([0,11500]);
xticks(0:1000:11500);
xlabel('Дальность полёта,м');
ylabel('Высота полёта, м');
%
figure(3)
tetagrad=rad2deg(teta);
h3=0:0.1:t(end);
gr3=spline(t,tetagrad,h3);
plot(t,tetagrad,h3,gr3)
grid on
xlabel('Время полёта,с');
ylabel('Угол полёта,град');
%
figure(4)
gr4=spline(t,V,h3);
plot(t,V,h3,gr4)
grid on
xlabel('Время полёта,с');
ylabel('Скорость полёта,м/с');


