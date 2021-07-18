% Quantum tunneling simulation

close all;clear all;clc

n=1000; %se define el tamaño de la matriz, 1000x1000
x=linspace(-5,5,n); %se define la variable x
h=10/n; %se define el espaciamiento

%funcion de onda inicial
psi0=((10/pi)^(1/4)).*exp(-5*(x+3).^2+1i*25.*x)'; %se define la funcion de onda inicial
%Ipsi0=imag(psi0);
%Rpsi0=real(psi0);
%dens0=Ipsi0.^2+Rpsi0.^2;

% potencial
u=(5e4/7)*eye(n); %se define el potencial
u(1:494,1:494)=0; %de tal manera que sea de -0.05 a 0.05
u(506:1000,506:1000)=0;

% hamiltoniano (usando derivadas espectrales)
H=((1/2)*(-2*eye(n)+diag(ones((n-1),1),1)+diag(ones((n-1),1),-1))/(h^2)+u); %se define el hamiltoniano

%solucion y animacion
for t=0:1:22
    psit=expm(-1i*H*t*.01)*psi0; %se resuelve para obtener la funcion de onda dependiente del tiempo
    Ipsit=imag(psit);
    Rpsit=real(psit);
    dens=(abs(psit)).^2; %se define la densidad de probabiliad
    figure(1)
    plot(x,dens,'Linewidth',2); %se grafica la evolucion de la densidad de probabilidad
    hold on
    rectangle('Position',[-0.05 0 .1 12]) %se grafica el potencial
    hold off
    set(gcf,'defaulttextinterpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',12)
    xlabel('$x$'),ylabel('$|\Psi(x,t)|^2$'),axis([-4 4 0 3.5]);
    %{
    figure(2)
    plot(x,Rpsit,'Linewidth',2); %se grafica la evolucion de la parte real de la funcion de onda dep. del tiempo
    hold on
    rectangle('Position',[-0.05 -2 .1 12]) %se grafica el potencial
    hold off
    set(gcf,'defaulttextinterpreter','latex')
    set(gca,'TickLabelInterpreter','latex','fontsize',13)
    xlabel('$x$'),ylabel('$\Psi_\mathrm{real}(x,t)$'),axis([-5 5 -2 3.5]);
    %}
    format long
    pause(1e-1000000000000)
    t
end
% probabilidades
probTOTAL1=sum(dens)*h; %se asegura que la probabilidad total sea 1

%probabilidad que la onda se refleje, es una integral de 0 a n/2
R=0;
for a=1:n/2
   R=R+dens(a); 
end
probR=R*h

%probabilidad que la onda se transmita, es una integral de n/2 a n
T=0;
for a=n/2:n
   T=T+dens(a);
end
probT=T*h
probTOTAL2=probT+probR %se corrobora que la suma sea 1
