%% computes the derivative using Chebyshev spectral method
  
close all;clear all;clc

N=100; % se define el numero de puntos
[D,x]=cheb(N); % se utiliza la funcion cheb para definir el operador D y el vector x

f=x.^4-x.^3; % se define la funcion a derivar
df=6*x.*(2*x-1); % se define la segunda derivada analitica de la funcion f

d2=D*D; % se define el operador de segunda derivada d2
dfch=d2*f; % se opera el operador d2 a la funcion

figure(1)
plot(x,f,'LineWidth',3,'color','r') % se grafica la funcion original
ylim([-50 800])
xlim([-5 5])
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
ylabel('$f_b(x)=x^4-x^3$')
xlabel('$x$')

figure(2)
plot(x,df,'LineWidth',5,'color','b') % se grafica la segunda derivada analitica de f
hold on
plot(x,dfch,'k*','Markersize',14) % se grafica la derivada numerica de f
xlim([-5 5])
legend({'Analitica','Chebyshev'},'Location','northeast','Interpreter','latex')
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
ylabel('$f^{(2)}_b(x)$')
xlabel('$x$')

error2=norm(dfch-df,inf) % se calcula el error de la derivada numerica respecto a la analitica

i=find(abs(x)<1e-15); % se encuentra la posicion donde el vector x es 0
dfch(i) % se evalua la derivada numerica en ese punto
df(i) % % se evalua la derivada analitica en ese punto
