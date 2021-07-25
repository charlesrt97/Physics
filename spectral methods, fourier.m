%% derivative using fourier spectral method

close all;clear all;clc
% se definen parametros L (ancho ventana) y N (numero de puntos en el dominio espacial)
L=10;
N=1000;

dx=L/N; % paso en el dominio espacial
dk=2*pi/L; % paso en el dominio de frecuencias

x=(-N/2:1:(N/2-1/2))*dx;x=x'; % vector espacial
kx=(-N/2:1:(N/2-1/2))*dk;kx=kx'; % vector frecuencias

f1=sin(cos(2*pi*x)); % se define la funcion a derivar
df1=8*pi^3*sin(2*pi*x).*(((sin(2*pi*x)).^2+1).*cos(cos(2*pi*x))-3*cos(2*pi*x).*sin(cos(2*pi*x))); % se define la tercera derivada analitica de f1

ff1=fft(f1).*fftshift((1i*kx).^3); % se realiza la transformada de fourier a f1
Desp1=real(ifft(ff1)); % se realiza la transformada inversa de fourier para cambiar al espacio orginial

figure(1)
plot(x,f1,'LineWidth',3,'color','r') % se grafica la funcion original, f1
%ylim([-1.5 2.1])
xlim([-5 5])
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
ylabel('$f_d(x)=\sin(\cos(2\pi x))$')
xlabel('$x$')

figure(2)
plot(x,df1,'LineWidth',3,'color','b') % se grafica la cuarta derivada analitica de f1
hold on
plot(x,Desp1,'k*','MarkerSize',8) % se grafica la cuarta derivada de f1, utilizando Fourier
%ylim([-1.1e6 1.1e6])
xlim([-5 5])
legend({'Analitica','Fourier'},'Location','northeast','Interpreter','latex')
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
ylabel('$f^{(3)}_d(x)$')
xlabel('$x$')

error4=norm(Desp1-df1,inf) % se calcula el error de la derivada numerica respecto a la analitica

i=find(x==0) % se encuentra la posicion donde el vector x es 0
Desp1(i) % se evalua la derivada numerica en ese punto
df1(i) % se evalua la derivada analitica en ese punto


