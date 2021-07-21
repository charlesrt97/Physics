% Kronig - penney potential

close all;clear all;clc
% parametros iniciales
a=1;
p1=0.25;
p2=0.75;
N=12; N=N+1;
v0=-3; % se varia v0=0, -1, -3
w=200;% iteraciones

Ix=0; Iy=0; % parte en x y parte en y
l=1:N;
nx1=(1+(-1).^l.*(2*l-1))/4;

Et=zeros(3*w+1,N*N); %energia total
k=0:(1/(3*w)):1; % vector k

[nx11,ny11]=meshgrid(nx1,nx1);
n4=[zeros(N*N,1) nx11(:) ny11(:)];

%nsort=n4(:,1);
%nx=n4(:,2);
%ny=n4(:,3);

n4(:,1)=(n4(:,2)).^2+(n4(:,3)).^2; % n_sort
%nsort=n4(:,1);

[m,ind]=sort(n4(:,1));
n4=n4(ind,:);

nsort=n4(:,1);
nx=n4(:,2);
ny=n4(:,3);

h2=zeros(1,N); 
%n=N-2;

for m=1:2*N-1
    h2(m)=1i*(exp(1i*2*pi*(m-N)*p1)-exp(1i*2*pi*(m-N)*p2))./(2*pi*(m-N));
end

hv=diag(v0*(p2-p1).^2+4*nsort);

for i=1:N*N
    for j=i+1:N*N % rellena la parte de la derecha de la diagonal principal
        if nx(i)==nx(j)
            Ix=p2-p1;
        else 
            Ix=h2(nx(i)-nx(j)+N);
        end
        if ny(i)==ny(j)
            Iy=p2-p1;
        else
            Iy=h2(ny(i)-ny(j)+N);
        end
        hv(i,j)=hv(i,j)+v0*Ix*Iy; % potencial
    end
end

hv2=triu(hv,1); hv2=hv2'; % rellena solo la parte izquierda de la diagonal, incluyendo esta
hv=hv+hv2; % potencial en forma matricial

kx=0:(pi/w):(pi/a); ky=0:(pi/w):(pi/a);

% gamma a X'
m=1:N*N;
for i=1:w+1
   [VE,VA]=eig(diag((2*nx(m)+kx(1)*a/pi).^2+(2*ny(m)+ky(i)*a/pi).^2-4*((nx(m)).^2+(ny(m)).^2))+hv);
   VA=diag(VA);
   E(1,i,:)=VA;
end

% x' a M
for t=1:w+1
   [VE,VA]=eig(diag((2*nx(m)+kx(t)*a/pi).^2+(2*ny(m)+ky(w+1)*a/pi).^2-4*((nx(m)).^2+(ny(m)).^2))+hv);
   VA=diag(VA);
   E(t,w+1,:)=VA; 
end

% M a gamma
for b=w+1:-1:1
   [VE,VA]=eig(diag((2*nx(m)+kx(b)*a/pi).^2+(2*ny(m)+ky(b)*a/pi).^2-4*((nx(m)).^2+(ny(m)).^2))+hv);
   VA=diag(VA);
   E(b,b,:)=VA; 
end

% se construye el vector de energias para graficar
for g=1:w+1
    Et(g,:)=E(1,g,:);
    Et(g+w,:)=E(g,1+w,:);
    Et(g+2*w,:)=E(2+w-g,2+w-g,:);
end

k=[k(1:1+2*w) k(2+2*w:end)]; % se construye el vector k para graficar

% graficas
figure(1)
plot(k, Et(:,1),'LineWidth',4)
hold on
plot(k, Et(:,2),'LineWidth',4)
hold on
plot(k, Et(:,3),'LineWidth',4)
hold on
plot(k, Et(:,4),'LineWidth',4)
hold on
plot(k, Et(:,5),'LineWidth',4)
hold on
plot(k, Et(:,6),'LineWidth',4)
hold on
plot(k, Et(:,7),'LineWidth',4)
hold on
plot(k, Et(:,8),'LineWidth',4)
hold on
plot(k, Et(:,9),'LineWidth',4)
hold on
plot(k, Et(:,10),'LineWidth',4)
hold on
plot(k, Et(:,11),'LineWidth',4)
hold on
plot(k, Et(:,12),'LineWidth',4)
hold on
plot(k, Et(:,13),'LineWidth',4)
ylim([-2 4])
set(gca,'xtick',[])
set(gcf,'defaulttextinterpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontsize',24)
ylabel('$E/E_{ISW}$');

xticks([0 0.333 0.666 1])
xticklabels({'$\Gamma$', 'X''', 'M' , '$\Gamma$'})

%text(-0.01,       floor(Et(1,1))-0.25,'$\Gamma$','FontSize', 24);
%text(0.32,        floor(Et(1,1))-0.25,'X''',   'FontSize', 24);
%text(0.65,        floor(Et(1,1))-0.25,'M',     'FontSize', 24);
%text(k(end)-0.01, floor(Et(1,1))-0.25,'$\Gamma$','FontSize', 24);


[Kx,Ky]=meshgrid(k);

figure(2)
surf(Et(:,1),Kx,Ky)
hold on
surf(Et(:,2),Kx,Ky)
hold on
surf(Et(:,3),Kx,Ky)
hold on
surf(Et(:,4),Kx,Ky)
hold on
surf(Et(:,5),Kx,Ky)
hold on
surf(Et(:,6),Kx,Ky)
hold on
surf(Et(:,7),Kx,Ky)
hold on
surf(Et(:,8),Kx,Ky)
hold on
surf(Et(:,9),Kx,Ky)
hold on
surf(Et(:,10),Kx,Ky)
hold on
surf(Et(:,11),Kx,Ky)
hold on
surf(Et(:,12),Kx,Ky)
hold on
surf(Et(:,13),Kx,Ky)
hold on
surf(Et(:,14),Kx,Ky)
hold on
%xlim([0 6])
shading interp



%{
%}

