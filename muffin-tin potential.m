% muffin-tin potential

close all;clear all;clc
% parametros iniciales
a=1;
p1=0.25;
p2=0.75;
N=12; N=N+1;
v0=0;
w=200;% iteraciones

Ix=0; Iy=0; % parte en x y parte en y
l=1:N;
nx1=(1+(-1).^l.*(2*l-1))/4;

Et=zeros(3*w+1,N*N); %energia total
k=0:(1/(3*w)):1; % vector k

[nx11,ny11]=meshgrid(nx1,nx1);
n4=[zeros(N*N,1) nx11(:) ny11(:)];

n4(:,1)=(n4(:,2)).^2+(n4(:,3)).^2; % n_sort
%nsort=n4(:,1);

[m,ind]=sort(n4(:,1));
n4=n4(ind,:);

nsort=n4(:,1);
nx=n4(:,2);
ny=n4(:,3);

h2=zeros(1,N); 
h22=zeros(1,N);
r=0.25;

xmin=1/2-r;
xmax=1/2+r;

syms x y

for m=1:2*N-1
    %h2(m)=1i*(exp(1i*2*pi*(m-N)*p1)-exp(1i*2*pi*(m-N)*p2))./(2*pi*(m-N));
    h2(m)=exp(1i*2*pi*(m-N));
end

h3(m)=h2(m).^x;

for s=1:2*N-1
    %h2(m)=1i*(exp(1i*2*pi*(m-N)*p1)-exp(1i*2*pi*(m-N)*p2))./(2*pi*(m-N));
    h22(s)=exp(1i*2*pi*(m-N));
end

h222(s)=h22(s).^y;

ymin=1/2-sqrt(r.^2-(x-1/2).^2);
ymax=1/2+sqrt(r.^2-(x-1/2).^2);

hv=diag(pi*r.^2*v0+4*nsort);

% relllena la parte triangular de arriba
for i=1:N*N
    for j=i+1:N*N
        syms x y
        if nx(i)==nx(j)
            Ix=r;
        else 
            Ix=h3((nx(i)-nx(j)+N));
            Ix=int(Ix,[xmin xmax]);
        end
        if ny(i)==ny(j)
            Iy=r;
        else
            Iy=h222((ny(i)-ny(j)+N));
            Iy=int(Iy,[ymin ymax]);
        end
        if nx(i)==nx(j) && ny(i)==ny(j)
            hv(i,j)=hv(i,j)+pi*v0*r^2; % potencial
        else
            I=integral2(Ix*Iy,xmin,xmax,ymin,ymax);
            hv(i,j)=hv(i,j)+pi*v0*I; % potencial
            %hv(i,j)=hv(i,j);
        end
    end
end

hv2=triu(hv,1); hv2=hv2'; % rellenar solo la parte triangular de abajo de la matriz
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
ylim([0 5])
set(gca,'xtick',[])
ylabel('$E/E_{ISW}$');

xticks([0 0.333 0.666 1])
xticklabels({'$\Gamma$', 'X''', 'M' , '$\Gamma$'})
