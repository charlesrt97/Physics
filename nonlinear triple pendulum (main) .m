close all; clear all; clc
%%%%%% --- PENDULO TRIPLE NO-LINEAL --- %%%%%%
%%%%%% --- Programa principal       --- %%%%%%
% Se establecen los valores de las constantes

l1=1;    % Longitud barra 1
l2=1.5;  % Longitud barra 2
l3=2;    % Longitud barra 3
m1=1;    % Valor masa 1
m2=1.5;  % Valor masa 2
m3=2;    % Valor masa 3

ts=[0 10];           % Tiempo que se dejara oscilar al pendulo triple
y0=[0 0 1 0.1 0 0.1]; % Condiciones iniciales [PosAng1 PosAng2 PosAng3 VelAng1 VelAng2 VelAng3]

[t,x]=ode45(@myeq,ts,y0); % Se utiliza ODE45 para resolver el sistema de 6 ODEs no lineales

% La matriz x guarda 6 vectores columna, los primeros 3 vectores 
% corresponden a la posicion angular de las barras 1, 2 y 3,
% respectivamente. Los ultimos 3 vectores (vector 4, 5 y 6 de la matriz x) 
% corresponden a la velocidad angular de las barras 1, 2 y 3, 
% respectivamente.
% Cada elemento del vector t corresponde al tiempo de cada
% elemento de cada vector de la matriz x.

% Se grafica la posicion angular de las tres barras (o 3 masas)
figure(1)
plot(t,x(:,1),'b','LineWidth',3) % Posicion angular de la barra 1
hold on
plot(t,x(:,2),'g','LineWidth',3) % Posicion angular de la barra 2
hold on
plot(t,x(:,3),'k','LineWidth',3) % Posicion angular de la barra 3
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
legend('$\theta_1(t)$','$\theta_2(t)$','$\theta_3(t)$','Location','northeast','Interpreter','latex')
xlabel('$t$');                   % Variable independiente, tiempo
ylabel('$\theta (t)$')           % Variable dependiente, posicion angular

% Se grafica la velocidad angular de las tres barras (o 3 masas)
figure(2)
plot(t,x(:,4),'b','LineWidth',3)  % Velocidad angular de la barra 1
hold on
plot(t,x(:,5),'g','LineWidth',3)  % Velocidad angular de la barra 2
hold on
plot(t,x(:,6),'k','LineWidth',3)  % Velocidad angular de la barra 3
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
legend('$\dot{\theta}_1(t)$','$\dot{\theta}_2(t)$','$\dot{\theta}_3(t)$','Location','northeast','Interpreter','latex')
xlabel('$t$');                    % Variable independiente, tiempo
ylabel('$\dot{\theta} (t)$')      % Variable dependiente, velocidad angular

% Se grafican los diagramas de espacio-fase, es decir,
% la velocidad angular de la barra vs. la posicion angular de la misma
figure(3)
plot(x(:,1),x(:,3),'b','LineWidth',3) % Diagrama espacio-fase barra 1
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
xlabel('$\theta_1$');                 % Posicion angular barra 1
ylabel('$\dot{\theta}_1$')            % Velocidad angular barra 2
figure(4)
plot(x(:,2),x(:,4),'g','LineWidth',3) % Diagrama espacio-fase barra 2
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
xlabel('$\theta_2$');                 % Posicion angular barra 2
ylabel('$\dot{\theta}_2$')            % Velocidad angular barra 2
xlim([-1 1])
figure(5)
plot(x(:,3),x(:,6),'k','LineWidth',3) % Diagrama espacio-fase barra 3
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
xlabel('$\theta_3$');                 % Posicion angular barra 3
ylabel('$\dot{\theta}_3$')            % Velocidad angular barra 3

% Para obtener las posiciones de las tres masas en el plano xy
x1=zeros(1,length(x)); % Se crea vector para la posicion x de la masa 1
y1=zeros(1,length(x)); % Se crea vector para la posicion y de la masa 1
x2=zeros(1,length(x)); % Se crea vector para la posicion x de la masa 2
y2=zeros(1,length(x)); % Se crea vector para la posicion y de la masa 2
x3=zeros(1,length(x)); % Se crea vector para la posicion x de la masa 3
y3=zeros(1,length(x)); % Se crea vector para la posicion y de la masa 3

% Se rellenan los vectores de posicion anteriores con la informacion de 
% los vectores de posicion angular, de la matriz x (i.e. los 3 primeros 
% vectores de la matriz x.
for t=1:length(x)
    x1(t)=l1*cos(x(t,1)); % Vector x de la masa 1
    y1(t)=l1*sin(x(t,1)); % Vector y de la masa 1
    x2(t)=l1*cos(x(t,1))+l2*cos(x(t,2)); % Vector x de la masa 2
    y2(t)=l1*sin(x(t,1))+l2*sin(x(t,2)); % Vector y de la masa 2
    x3(t)=l1*cos(x(t,1))+l2*cos(x(t,2))+l3*cos(x(t,3)); % Vector x de la masa 3
    y3(t)=l1*sin(x(t,1))+l2*sin(x(t,2))+l3*sin(x(t,3)); % Vector y de la masa 3
end

% Se grafican las posiciones de las 3 masas
figure(6)
plot(y1,-x1,'b','LineWidth',3) % Grafica posicion en plano xy, de la masa 1
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
xlabel('$x$');                 % Eje x
ylabel('$y$');                 % Eje y
figure(7)
plot(y2,-x2,'g','LineWidth',3) % Grafica posicion en plano xy, de la masa 2
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
xlabel('$x$');                 % Eje x
ylabel('$y$');                 % Eje y
figure(8)
plot(y3,-x3,'k','LineWidth',3) % Grafica posicion en plano xy, de la masa 3
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
xlabel('$x$');                 % Eje x
ylabel('$y$');                 % Eje y

% Se elabora la simulacion del pendulo triple, utilizando la informacion 
% de las posiciones de las 3 masas
s = length(x1); % Longitud de los vectores posicion
% Comandos para grabar animacion
v = VideoWriter('TriplePendulo3.avi');
v.FrameRate = 60;
open(v)
% Comienza la animacion del pendulo 
for n = 1:s
    figure(8)
    hold off
    plot([-0.5 0.5],[0 0],'LineWidth',2,'Color',[0.9 0.9 0.9]) % Simula la base del pendulo
    hold on
    plot([0 y1(n)],[0 -x1(n)],'Color',[0.9290 0.6940 0.1250],'LineWidth',3); % Simula la barra 1, linea entre punto inicial a la posicion de la masa 1
    hold on
    plot(y1(n),-x1(n),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4*m1); % Simula la masa 1, posicion de la masa 1
    hold on
    plot([y1(n) y2(n)],[-x1(n) -x2(n)],'Color',[0 0.4470 0.7410],'LineWidth',3); % Simula la barra 2, linea entre posicion masa 1 a la posicion de la masa 2
    hold on
    plot(y2(n),-x2(n),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4*m2); % Simula la masa 2, posicion de la masa 2
    hold on
    plot([y2(n) y3(n)],[-x2(n) -x3(n)],'Color',[0.4660 0.6740 0.1880],'LineWidth',3); % Simula la barra 3, linea entre posicion masa 2 a la posicion de la masa 3
    hold on
    plot(y3(n),-x3(n),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4*m3); % Simula la masa 3, posicion de la masa 3
    set(gcf,'defaulttextinterpreter','latex')
    set(gcf,'defaulttextinterpreter','latex')   % Formato latex
    set(gca,'TickLabelInterpreter','latex','fontsize',24)
    set(gca,'Color','k')                        % Fondo negro
    xlabel('$x$');  ylabel('$y$');              % Etiquetas de los ejes
    % Texto sobre datos de la animacion (masas, longitudes y condiciones iniciales)
    text(3.1,1.5,'$m_1=1$','Color','w'); % Valor masa 1
    text(3.1,1,'$m_2=1.5$','Color','w'); % Valor masa 2
    text(3.1,0.5,'$m_3=2$','Color','w'); % Valor masa 3
    text(4.5,1.5,'$l_1=1$','Color','w'); % Valor longitud barra 1
    text(4.5,1,'$l_2=1.5$','Color','w'); % Valor longitud barra 2
    text(4.5,0.5,'$l_3=2$','Color','w'); % Valor longitud barra 3
    text(2.4,-4,'Condiciones iniciales','Color','w');
    text(3.3,-4.5,'$\theta_1=0$','Color','w'); % Condicion inicial pos. ang. barra 1
    text(3.3,-5,'$\theta_2=0$','Color','w'); % Condicion inicial pos. ang. barra 2
    text(3.3,-5.5,'$\theta_3=1$','Color','w'); % Condicion inicial pos. ang. barra 3
    text(4.5,-4.5,'$\dot{\theta}_1=0.1$','Color','w'); % Condicion inicial vel. ang. barra 1
    text(4.5,-5,'$\dot{\theta}_2=0$','Color','w'); % Condicion inicial vel. ang. barra 1
    text(4.5,-5.5,'$\dot{\theta}_3=0.1$','Color','w'); % Condicion inicial vel. ang. barra 1
    axis([-6 6 -6 2]);     % Limites de los ejes
    frame = getframe(gcf); % Grabacion
    writeVideo(v,frame);   % Grabacion
end
close(v)                   % Se cierra grabacion
