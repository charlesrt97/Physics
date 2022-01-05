% Ito's indexing method 
% Identifies an unknown sample from its powder diffraction data

close all;clear all;clc
A=load('data2theta.txt');
lambda=1.5418;
% calculo de los valores de Q_obs
for i=1:length(A(:,3))
    Qobs(i)=4.*(sind(A(i,3)/2)).^2/(lambda).^2;
end
Q200=Qobs(3);
Q020=Qobs(9);
Q002=Qobs(12);
ar=sqrt(Q200/4); %parametro de celda reciproca, a
br=sqrt(Q020/4); %parametro de celda reciproca, b
cr=sqrt(Q002/4); %parametro de celda reciproca, c
a=1/ar; %parametro de celda a
b=1/br; %parametro de celda b
c=1/cr; %parametro de celda c

Q100=Q200/4; 
Q010=Q020/4;
Q001=Q002/4;

Q011=Q010+Q001;
Q101=Q100+Q001;
Q110=Q100+Q010;

% para alpha

w=1;
v=1;
alpha1=91;
for i=1:length(A(:,3))
    angg1(w)=Qobs(i)-Q011;
    for j=i:length(A(:,3))
       angg2(v)=Qobs(j)-Q011;
       if 0.0001 > abs(angg1(w)+angg2(v))
           ang1alpha=angg1(w);
           ang2alpha=angg2(v);
           alpha1=acos((Qobs(j)-Qobs(i))./(4*br*cr))*180/pi;
       else
           alpha2=90;
       end
       v=v+1;
    end
    w=w+1;
    if alpha1 < alpha2
        alpha=alpha1;
    else
        alpha=alpha2;
    end
end
alpha;

% para beta

k=1;
l=1;
beta1=91;
for i=1:length(A(:,3))
    angg1(k)=Qobs(i)-Q101;
    for j=i:length(A(:,3))
       angg2(l)=Qobs(j)-Q101;
       if 0.0003 > abs(angg1(k)+angg2(l))
           ang1beta=angg1(k);
           ang2beta=angg2(l);
           beta1=acos((Qobs(j)-Qobs(i))./(4*ar*cr))*180/pi;
       else
           beta2=90;
       end
       l=l+1;
    end
    k=k+1;
    if beta1 < beta2
        beta=beta1;
    else
        beta=beta2;
    end
end
%ang1beta
%ang2beta
beta;

% para gamma
t=1;
h=1;
gamma1=91;
for i=1:length(A(:,3))
    angg1(t)=Qobs(i)-Q011;
    for j=i:length(A(:,3))
       angg2(h)=Qobs(j)-Q011;
       if 0.0001 > abs(angg1(t)+angg2(h))
           ang1gamma=angg1(t);
           ang2gamma=angg2(h);
           gamma1=acos((Qobs(j)-Qobs(i))./(4*ar*br))*180/pi;
       else
           gamma2=90;
       end
       h=h+1;
    end
    t=t+1;
    if gamma1 < gamma2
        gamma=gamma1;
    else
        gamma=gamma2;
    end
end
gamma;
%calculo angulos reales
alpha2=acosd(cosd(gamma)*cosd(beta))-cosd(alpha)/(sind(gamma)*sind(beta));
beta2=acosd(cosd(gamma)*cosd(alpha))-cosd(beta)/(sind(gamma)*sind(alpha));
gamma=acosd(cosd(beta)*cosd(alpha))-cosd(gamma)/(sind(beta)*sind(alpha));

fprintf('Parametros de celda, espacio reciproco: a*=%f A, b*=%f A, c*=%f A\n',ar,br,cr)
fprintf('Parametros de celda, espacio real: a=%f A, b=%f A, c=%f A\n',a,b,c)
fprintf('Angulos reciprocos: alpha*=%f, beta*=%f, gamma*=%f\n',alpha,beta,gamma)
fprintf('Angulos reales: alpha=%f, beta=%f, gamma=%f\n',alpha2,beta2,gamma2)

%indizacion
B=load('indeex.txt');
fprintf('Linea Theta   h k l\n')
for y=1:length(B(:,3))
    Qj(y)=B(y,1).^2.*Q100+B(y,2).^2.*Q010+B(y,3).^2.*Q001;
    for f=1:length(A(:,3))
        Qobs2(f)=4.*(sind(A(f,3)/2)).^2/(lambda).^2;
        if 0.0004 > abs(Qj(y)-Qobs2(f))
            fprintf(' %d    %.3f   %d %d %d\n',A(f,2),A(f,3),B(y,1),B(y,2),B(y,3))
        %else
         %   fprintf('%d %0.3f - - -\n',A(f,2),A(f,3))
        end
    end
end
