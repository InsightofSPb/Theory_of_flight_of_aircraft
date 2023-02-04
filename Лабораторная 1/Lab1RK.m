clear all
%Dano
d=0.150; m0=60; mtop=13;ta=0.7;Lh=4; tetad=1.55334; ue=2500; ro0=1.23; 
ftr=0.15; g=9.81; ij=1.15; 
%Движение по направляющим
mpas=m0;
Q=mtop/ta; %расход
R=Q*ue; %тяга
td=sqrt(2*Lh/((R/m0)-g*(sin(tetad)+ftr*cos(tetad)))); %время схода с напр.
Vd=(R/m0-g*(sin(tetad)+ftr*cos(tetad)))*td; %скорость схода с напр.
%Расчёт траектории на активном участке методом Рунге-Кутта 4-го порядка
t(1)=td;
V(1)=Vd;
teta(1)=tetad;
y(1)=Lh*sin(tetad);
x(1)=Lh*cos(tetad);
ro0=1.23;
i=0;
S=pi/4*d^2;
h=0:1:1000;
Vet=[100,150,200,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,500,550,600,650,700,750,800,850,900,950,1000];
cxet=[0.255,0.257,0.260,0.268,0.274,0.280,0.295,0.321,0.361,0.411,0.460,0.500,0.542,0.574,0.603,0.628,0.648,0.665,0.680,0.690,0.700,0.709,0.715,0.719,0.735,0.732,0.721,0.704,0.680,0.665,0.645,0.624,0.605,0.583,0.567];
Cx=cxet*ij; % Сx с учётом несоответсвия эталонному ЛА
gr1=pchip(Vet,Cx,h);
h2=0.001;
while t(end)<=ta
    i=i+1;
    dx1(i)=(V(i)*cos(teta(i)));
    x2(i)=x(i)+h2/2*dx1(i);
    %
    dy1(i)=(V(i)*sin(teta(i)));
    y2(i)=y(i)+h2/2*dy1(i);
    %
    dteta1(i)=(1/V(i))*(-g)*cos(teta(i));
    teta2(i)=teta(i)+h2/2*dteta1(i);
    %
    t2(i+1)=t(i)+h2/2;
    %
    Cx(i)=gr1(round(V(i)));
    ro(i)=ro0*exp(-y(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V(i)^2;
    %
    dV1(i)=(R/(m0-Q*t(i))-g*sin(teta(i))-X(i)/(m0-Q*t(i)));
    V2(i)=V(i)+h2/2*dV1(i);
    %%%
    dx2(i)=(V2(i)*cos(teta2(i)));
    x3(i)=x(i)+h2/2*dx2(i);
    %
    dy2(i)=(V2(i)*sin(teta2(i)));
    y3(i)=y(i)+h2/2*dy2(i);
    %
    dteta2(i)=(1/V2(i))*(-g)*cos(teta2(i));
    teta3(i)=teta(i)+h2/2*dteta2(i);
    %
    t3(i+1)=t(i)+h2/2;
    %
    Cx(i)=gr1(round(V2(i)));
    ro(i)=ro0*exp(-y2(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V2(i)^2; 
    %
    dV2(i)=(R/(m0-Q*t(i))-g*sin(teta2(i))-X(i)/(m0-Q*t(i)));
    V3(i)=V(i)+h2/2*dV2(i);
    %%%
    dx3(i)=(V3(i)*cos(teta3(i)));
    x4(i)=x(i)+h2*dx3(i);
    %
    dy3(i)=(V3(i)*sin(teta3(i)));
    y4(i)=y(i)+h2*dy3(i);
    %
    dteta3(i)=(1/V3(i))*(-g)*cos(teta3(i));
    teta4(i)=teta(i)+h2*dteta3(i);
    %
    t4(i+1)=t(i)+h2;
    %
    Cx(i)=gr1(round(V3(i)));
    ro(i)=ro0*exp(-y3(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V3(i)^2; 
    %
    dV3(i)=(R/(m0-Q*t(i))-g*sin(teta3(i))-X(i)/(m0-Q*t(i)));
    V4(i)=V(i)+h2*dV3(i);
    %%%
    %    
    Cx(i)=gr1(round(V4(i)));
    ro(i)=ro0*exp(-y(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V4(i)^2;
    %
    dV4(i)=(R/(m0-Q*t(i))-g*sin(teta4(i))-X(i)/(m0-Q*t(i)));
    %
    dteta4(i)=(1/V4(i))*(-g)*cos(teta4(i));
    %
    dx4(i)=(V4(i)*cos(teta4(i)));
    %
    dy4(i)=(V4(i)*sin(teta4(i)));
    %
    t4(i+1)=t(i)+h2;
    %
    DeltaV(i)=h2/6*(dV1(i)+2*dV2(i)+2*dV3(i)+dV4(i));
    V(i+1)=V(i)+DeltaV(i);
    Deltateta(i)=h2/6*(dteta1(i)+2*dteta2(i)+2*dteta3(i)+dteta4(i));
    teta(i+1)=teta(i)+Deltateta(i);
    Deltax(i)=h2/6*(dx1(i)+2*dx2(i)+2*dx3(i)+dx4(i));
    x(i+1)=x(i)+Deltax(i);
    Deltay(i)=h2/6*(dy1(i)+2*dy2(i)+2*dy3(i)+dy4(i));
    y(i+1)=y(i)+Deltay(i);
    t(i+1)=t(i)+h2;
end
%Расчёт траектории на пассивном участке методом Рунге-Кутта 4-го порядка
mpas=m0-mtop; 
 while y>0
    i=i+1;
    dx1(i)=(V(i)*cos(teta(i)));
    x2(i)=x(i)+h2/2*dx1(i);
    %
    dy1(i)=(V(i)*sin(teta(i)));
    y2(i)=y(i)+h2/2*dy1(i);
    %
    dteta1(i)=(1/V(i))*(-g)*cos(teta(i));
    teta2(i)=teta(i)+h2/2*dteta1(i);
    %
    t2(i+1)=t(i)+h2/2;
    %
    Cx(i)=gr1(round(V(i)));
    ro(i)=ro0*exp(-y(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V(i)^2;
    %
    dV1(i)=(-g*sin(teta(i))-X(i)/mpas);
    V2(i)=V(i)+h2/2*dV1(i);
    %%%
    dx2(i)=(V2(i)*cos(teta2(i)));
    x3(i)=x(i)+h2/2*dx2(i);
    %
    dy2(i)=(V2(i)*sin(teta2(i)));
    y3(i)=y(i)+h2/2*dy2(i);
    %
    dteta2(i)=(1/V2(i))*(-g)*cos(teta2(i));
    teta3(i)=teta(i)+h2/2*dteta2(i);
    %
    t3(i+1)=t(i)+h2/2;
    %
    Cx(i)=gr1(round(V2(i)));
    ro(i)=ro0*exp(-y2(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V2(i)^2; 
    %
    dV2(i)=(-g*sin(teta2(i))-X(i)/mpas);
    V3(i)=V(i)+h2/2*dV2(i);
    %%%
    dx3(i)=(V3(i)*cos(teta3(i)));
    x4(i)=x(i)+h2*dx3(i);
    %
    dy3(i)=(V3(i)*sin(teta3(i)));
    y4(i)=y(i)+h2*dy3(i);
    %
    dteta3(i)=(1/V3(i))*(-g)*cos(teta3(i));
    teta4(i)=teta(i)+h2*dteta3(i);
    %
    t4(i+1)=t(i)+h2;
    %
    Cx(i)=gr1(round(V3(i)));
    ro(i)=ro0*exp(-y3(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V3(i)^2; 
    %
    dV3(i)=(-g*sin(teta3(i))-X(i)/mpas);
    V4(i)=V(i)+h2*dV3(i);
    %%%
    
    Cx(i)=gr1(round(V4(i)));
    ro(i)=ro0*exp(-y(i)/9800);
    X(i)=(1/2)*Cx(i)*ro(i)*S*V4(i)^2;
    %
    dV4(i)=(-g*sin(teta4(i))-X(i)/mpas);
    %
    dteta4(i)=(1/V4(i))*(-g)*cos(teta4(i));
    %
    dx4(i)=(V4(i)*cos(teta4(i)));
    %
    dy4(i)=(V4(i)*sin(teta4(i)));
    %
    t4(i+1)=t(i)+h2;
    %%%
    DeltaV(i)=h2/6*(dV1(i)+2*dV2(i)+2*dV3(i)+dV4(i));
    V(i+1)=V(i)+DeltaV(i);
    Deltateta(i)=h2/6*(dteta1(i)+2*dteta2(i)+2*dteta3(i)+dteta4(i));
    teta(i+1)=teta(i)+Deltateta(i);
    Deltax(i)=h2/6*(dx1(i)+2*dx2(i)+2*dx3(i)+dx4(i));
    x(i+1)=x(i)+Deltax(i);
    Deltay(i)=h2/6*(dy1(i)+2*dy2(i)+2*dy3(i)+dy4(i));
    y(i+1)=y(i)+Deltay(i);
    t(i+1)=t(i)+h2;
 end  
 figure(2)
h2=0:10:x(end);
gr2=pchip(x,y,h2);
plot(x,y,h2,gr2)
xlim([0,11500]);
xticks(0:1000:11500);
grid on
xlabel('Дальность полёта,м');
ylabel('Высота полёта, м');
%
figure(3)
tetagrad=rad2deg(teta)
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
