%% computes the derivative using Chebyshev spectral method

close all;clear all;clc

N=3000; % se define el numero de puntos
[D,x]=cheb(N); % se utiliza la funcion cheb para definir el operador D y el vector x

f=airy(x)+erf(x); % se define la funcion a derivar
df=airy(1,x)+(2*exp(-x.^2))/sqrt(pi); % se define la segunda derivada analitica de la funcion f

dfch=D*f; % se opera el operador D a la funcion

figure(1)
plot(x,f,'LineWidth',3,'color','r') % se grafica la funcion original
xlim([-10 10])
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
ylabel('$f_c(x)=Ai (x)+erf (x)$')
xlabel('$x$')

figure(2)
plot(x,df,'LineWidth',5,'color','b') % se grafica la segunda derivada analitica de f
hold on
plot(x,dfch,'k*','MarkerSize',11) % se grafica la derivada numerica de f
xlim([-10 10])
legend({'Analitica','Chebyshev'},'Location','northeast','Interpreter','latex')
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
ylabel('$f^{(1)}_c(x)$')
xlabel('$x$')

error3=norm(dfch-df,inf) % se calcula el error de la derivada numerica respecto a la analitica

i=find(abs(x)<1e-14); % se encuentra la posicion donde el vector x es 0
dfch(i) % se evalua la derivada numerica en ese punto
df(i) % % se evalua la derivada analitica en ese punto
