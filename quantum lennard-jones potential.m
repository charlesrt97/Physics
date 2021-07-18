% quantum Lennard-Jones potential

function lennjonpot
E=-1.8*10^-3; %El valor de energ?a debe ser un n?mero muy cercano al esperado cuyo valor absoluto sea mayor%
nr=0 %El valor nr, checando las tablas se observa cual valor de energ?a le corresponde
l=0  %El mismo caso con el valor l, donde depende de 0,1 y 2
L=l+1/4; %N?mero cu?ntico orbital
m1=24*1836+12; %El valor de las masas
m2=m1; %Debido a que ambas masas son iguales
M=m1*m2/(m1+m2); %La obtenci?n de la masa reducida
a=4.96*10^7; %El par?metro a
b=624; %El par?metro b
p=[E*M 0 -1/2*(l^2+l+1/4) 0 0 0 M*b 0 0 0 0 0 -M*a]; %El valor de radical, para que se pueda hacer un m?todo de ra?ces
d=roots(p); %La obtenci?n de ra?ces del radical
dd=d(imag(d)==0); %Eliminar todas las ra?ces imaginarias
h=abs(dd(2))-abs(dd(1)); %Obtener el valor de h para la integral
if h>0.0001 %El if funciona para siempre obtener un valor positivo y real, donde los l?mites no sean los mismos
h=abs(dd(1))-abs(dd(2));
Ri=linspace(abs(dd(1)),abs(dd(2)),1000);
else
h=abs(dd(2))-abs(dd(3));
Ri=linspace(abs(dd(3)),abs(dd(2)),1000);
end
F=(E-a./(Ri.^12)+b./(Ri.^6)-(l^2+l+1/4)./(2*M*Ri.^2)).^(0.5); %Obtener el valor de la funci?n integral
q=trapz(F)*h/1000; %Integral por medio del trapecio
q=abs(real(q)); %Hacer que la integral s?lo tenga valores positivos
Fin= sqrt((2*M))*q/3.1415-(nr+1/2); %Sustituir la integral para obtener el valor de f(E)
while abs(Fin)>0.0001 %Ponerle un error del 10^-4, para que despu?s de ah?, el resultado de la energ?a es el correcto
E=E+0.00001*10^-3; %Irle sumando 10^-5 valores, para tener una buena exactitud
p=[E*M 0 -1/8 0 0 0 M*b 0 0 0 0 0 -M*a]; %Se repiten los pasos anteriores, obtener ra?ces, integrar, sustituir y obtener el valor f(E)
d=roots(p);
dd=d(imag(d)==0);
h=abs(dd(2))-abs(dd(1));
if h>0.0001
h=abs(dd(1))-abs(dd(2));
Ri=linspace(abs(dd(2)),abs(dd(1)),1000);
else
h=abs(dd(2))-abs(dd(3));
Ri=linspace(abs(dd(3)),abs(dd(2)),1000);
end
F=(E-a./(Ri.^12)+b./(Ri.^6)-(l^2+l+1/4)./(2*M*Ri.^2)).^(0.5);
q=trapz(F)*h/1000;
q=abs(real(q));
Fin= sqrt((2*M))*q/3.1415-(nr+1/2);
end
format long
E  %Mostrar el valor de la energ?a
h=abs(dd(1))-abs(dd(2)); %Aplicarl o mismo anteriormente, para obtener los valores donde se integr? la funci?n
if h>0.0001
r=[abs(dd(1)); abs(dd(2))]
else
r=[abs(dd(1)); abs(dd(3))]
end
%Gr?fica l=0, con los valores generar una matriz en X y Y y graficarla
x0=[6.571,6.579,6.591,6.606,6.626,6.653,6.687,6.732,6.794,6.8856,7.0525,7.7963,8.239,8.6351,9.0361,9.462,9.928,10.44,11.047,11.746,12.588,13.632];
y0=[-0.096,-0.1537,-0.2303,-.3283,-0.4500,-.5976,-0.7732,-0.9789,-1.2165,-1.488,-1.7951,-1.7951,-1.488,-1.2165,-0.9789,-.7732,-.5976,-.450,-.3283,-0.2303,-0.1537,-0.0960];
figure(1)
plot(x0,y0,'-o','markersize',4,'markerfacecolor','y')
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',16);
title ('Distancia vs Eigenenergias con numero cuantico orbital $l=0$')
xn00=x0(1):0.1:x0(22); %Darle valores a todas las eigenenerg?as, para despu?s mostrarlo en la gr?fica
xn01=x0(2):0.1:x0(21);
xn02=x0(3):0.1:x0(20);
xn03=x0(4):0.1:x0(19);
xn04=x0(5):0.1:x0(18);
xn05=x0(6):0.1:x0(17);
xn06=x0(7):0.1:x0(16);
xn07=x0(8):0.1:x0(15);
xn08=x0(9):0.1:x0(14);
xn09=x0(10):0.1:x0(13);
xn010=x0(11):0.1:x0(12);
yn00= y0(1)*xn00./xn00;
yn01= y0(2)*xn01./xn01;
yn02= y0(3)*xn02./xn02;
yn03= y0(4)*xn03./xn03;
yn04= y0(5)*xn04./xn04;
yn05= y0(6)*xn05./xn05;
yn06= y0(7)*xn06./xn06;
yn07= y0(8)*xn07./xn07;
yn08= y0(9)*xn08./xn08;
yn09= y0(10)*xn09./xn09;
yn010= y0(11)*xn010./xn010;
p=polyfit(x0,y0,10); %Aplicar una regresi?n de polinomio a la 10, para que se ajuste a los puntos y se pueda apreciar de una mejor manera la gr?fica
hold on
plot(x0,y0,'ro','markersize',4,'markerfacecolor','g') %Graficar la regresi?n al igual que los valores de las eigenerg?as
z=@(x0) polyval(p,x0);
fplot(z,[x0(1),x0(end)])
plot(xn00,yn00,'r')
plot(xn01,yn01,'r')
plot(xn02,yn02,'r')
plot(xn03,yn03,'r')
plot(xn04,yn04,'r')
plot(xn05,yn05,'r')
plot(xn06,yn06,'r')
plot(xn07,yn07,'r')
plot(xn08,yn08,'r')
plot(xn09,yn09,'r')
plot(xn010,yn010,'r')
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('Distancia de separacion')
ylabel('Energia')
grid on
set(gca,'XAxisLocation','top','YAxisLocation','left');
grid on
hold off


%Gr?fica l=1
x1=[6.5708,6.5792,6.5908,6.6061,6.6263,6.6525,6.6866,6.7318,6.7938,6.8853,7.0518,7.7976,8.2400,8.6362,9.0374,9.4633,9.9295,10.4517,11.0497,11.7505,12.5934,13.6404];
y1=[-0.0956,-0.1533,-0.2299,-.3278,-0.4495,-.5970,-0.773,-0.9782,-1.2158,-1.4872,-1.7942,-1.7942,-1.4872,-1.2158,-0.9782,-.773,-.5970,-.4495,-.3278,-0.2299,-0.1533,-0.0956];
figure(2)
plot(x1,y1,'-o','markersize',4,'markerfacecolor','r')
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
title ('Distancia vs Eigenenergias con numero cuantico orbital $l=1$')
xn10=x1(1):0.1:x1(22);
xn11=x1(2):0.1:x1(21);
xn12=x1(3):0.1:x1(20);
xn13=x1(4):0.1:x1(19);
xn14=x1(5):0.1:x1(18);
xn15=x1(6):0.1:x1(17);
xn16=x1(7):0.1:x1(16);
xn17=x1(8):0.1:x1(15);
xn18=x1(9):0.1:x1(14);
xn19=x1(10):0.1:x1(13);
xn110=x1(11):0.1:x1(12);
yn10= y1(1)*xn10./xn10;
yn11= y1(2)*xn11./xn11;
yn12= y1(3)*xn12./xn12;
yn13= y1(4)*xn13./xn13;
yn14= y1(5)*xn14./xn14;
yn15= y1(6)*xn15./xn15;
yn16= y1(7)*xn16./xn16;
yn17= y1(8)*xn17./xn17;
yn18= y1(9)*xn18./xn18;
yn19= y1(10)*xn19./xn19;
yn110= y1(11)*xn110./xn110;
p1=polyfit(x1,y1,10);
hold on
plot(x1,y1,'ro','markersize',4,'markerfacecolor','r')
z1=@(x1) polyval(p1,x1);
fplot(z1,[x1(1),x1(end)])
plot(xn10,yn10,'r')
plot(xn11,yn11,'r')
plot(xn12,yn12,'r')
plot(xn13,yn13,'r')
plot(xn14,yn14,'r')
plot(xn15,yn15,'r')
plot(xn16,yn16,'r')
plot(xn17,yn17,'r')
plot(xn18,yn18,'r')
plot(xn19,yn19,'r')
plot(xn110,yn110,'r')
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('Distancia de separacion')
ylabel('Energia')
grid on
set(gca,'XAxisLocation','top','YAxisLocation','left');
grid on
hold off

%Gr?fica l=2
x2=[6.5707,6.5791,6.5906,6.6060,6.6261,6.6523,6.6863,6.7315,6.7933,6.8847,7.0505,7.8003,8.2422,8.6386,9.0399,9.4664,9.9332,10.4563,11.056,11.7586,12.6050,13.6678];
y2=[-0.0949,-0.1524,-0.2289,-.3268,-0.4483,-.5958,-0.7712,-0.9768,-1.2143,-1.4856,-1.7926,-1.7926,-1.4856,-1.2143,-0.9768,-.7712,-.5958,-.4483,-.3268,-0.2289,-0.1524,-0.0949];
figure(3)
plot(x2,y2,'-o','markersize',4,'markerfacecolor','r')
%print -depsc myfig.eps
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
title ('Distancia vs Eigenenergias con numero cuantico orbital $l=2$')
xn20=x2(1):0.1:x2(22);
xn21=x2(2):0.1:x2(21);
xn22=x2(3):0.1:x2(20);
xn23=x2(4):0.1:x2(19);
xn24=x2(5):0.1:x2(18);
xn25=x2(6):0.1:x2(17);
xn26=x2(7):0.1:x2(16);
xn27=x2(8):0.1:x2(15);
xn28=x2(9):0.1:x2(14);
xn29=x2(10):0.1:x2(13);
xn210=x2(11):0.1:x2(12);
yn20= y2(1)*xn20./xn20;
yn21= y2(2)*xn21./xn21;
yn22= y2(3)*xn22./xn22;
yn23= y2(4)*xn23./xn23;
yn24= y2(5)*xn24./xn24;
yn25= y2(6)*xn25./xn25;
yn26= y2(7)*xn26./xn26;
yn27= y2(8)*xn27./xn27;
yn28= y2(9)*xn28./xn28;
yn29= y2(10)*xn29./xn29;
yn210= y2(11)*xn210./xn210;
p2=polyfit(x2,y2,10);
hold on
plot(x2,y2,'ro','markersize',4,'markerfacecolor','r')
z2=@(x2) polyval(p2,x2);
fplot(z2,[x2(1),x2(end)])
plot(xn20,yn20,'r')
plot(xn21,yn21,'r')
plot(xn22,yn22,'r')
plot(xn23,yn23,'r')
plot(xn24,yn24,'r')
plot(xn25,yn25,'r')
plot(xn26,yn26,'r')
plot(xn27,yn27,'r')
plot(xn28,yn28,'r')
plot(xn29,yn29,'r')
plot(xn210,yn210,'r')
set(gcf,'DefaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex','fontsize',14);
xlabel('Distancia de separacion')
ylabel('Energia')
grid on
set(gca,'XAxisLocation','top','YAxisLocation','left');
grid on
hold off
end
