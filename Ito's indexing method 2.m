% Ito's indexing method 
% Identifies an unknown sample from its powder diffraction data

close all;clear all;clc
% 110 100
A=load('Unknow04.dat');
lambda=1.54184;
% calculo de los valores de Q_obs
%fprintf('Intensity       Q_obs\n');
for i=1:length(A(:,2))
    Qobs(i)=4.*(sind(A(i,2)/2)).^2/(lambda).^2;
    %fprintf('%0.4f         %0.4f\n',A(i,1),Qobs(i));
end

%primero 100, segundo 110, tercero 111
Q100=Qobs(1); %condicion incial
Q110=Qobs(2); %condicion inicial
Q111=Qobs(3); %condicion inicial
Q010=Q110-Q100;

Q001=Q111-Q100-Q010;

%Q111=Q100+Q010+Q001;
Q011=Q010+Q001;
Q101=Q100+Q001;

%Q111=Q100+Q010+Q001;

ar=sqrt(Q100/1); %parametro de celda reciproca, a
br=sqrt(Q010/1); %parametro de celda reciproca, b
cr=sqrt(Q001/1); %parametro de celda reciproca, c
a=1/ar; %parametro de celda a
b=1/br; %parametro de celda b
c=1/cr; %parametro de celda c

% para alpha

w=1;
v=1;
alpha1=91;
for i=1:length(A(:,2))
    angg1(w)=Qobs(i)-Q011;
    for j=i:length(A(:,2))
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
for i=1:length(A(:,2))
    angg1(k)=Qobs(i)-Q101;
    for j=i:length(A(:,2))
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
for i=1:length(A(:,2))
    angg1(t)=Qobs(i)-Q011;
    for j=i:length(A(:,2))
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

alpha2=acosd(cosd(gamma)*cosd(beta))-cosd(alpha)/(sind(gamma)*sind(beta));
beta2=acosd(cosd(gamma)*cosd(alpha))-cosd(beta)/(sind(gamma)*sind(alpha));
gamma=acosd(cosd(beta)*cosd(alpha))-cosd(gamma)/(sind(beta)*sind(alpha));


%fprintf('Parametros de celda, espacio reciproco: a*=%f A, b*=%f A, c*=%f A\n',ar,br,cr)
%fprintf('Parametros de celda, espacio real: a=%f A, b=%f A, c=%f A\n',a,b,c)
%fprintf('Angulos reciprocos: alpha*=%f, beta*=%f, gamma*=%f\n',alpha,beta,gamma)
%fprintf('Angulos reales: alpha=%f, beta=%f, gamma=%f\n',alpha2,beta2,gamma2)

%indizacion
i=1;
in=zeros(729,3);
for h=0:8
    for k=0:8
        for l=0:8
            in(i,:)=[h k l];
            i=i+1;
            %fprintf('%d\n',l)
        end
    end
end


i=1;
%B=load('indexb.txt');
fprintf('Intensidad   2Theta    h k l      d          Q_calc     Q_obs      Diferencia\n')
for y=1:length(in(:,3))
    Qj(y)=in(y,1).^2.*Q100+in(y,2).^2.*Q010+in(y,3).^2.*Q001;
    for f=1:length(A(:,2))
        Qobs2(f)=4.*(sind(A(f,2)/2)).^2/(lambda).^2;
        if 0.00000001 > abs(Qj(y)-Qobs2(f))
            d(y)=(in(y,1)^2/a^2+in(y,2)^2/b^2+in(y,3)^2/c^2)^(-1/2);
            Qcalc(y,:)=in(y,1).^2.*ar.^2+in(y,2).^2.*br.^2+in(y,3).^2.*cr.^2+2.*in(y,3).*in(y,2).*br.*cr.*cosd(alpha)+2.*in(y,3).*in(y,1).*cr.*ar.*cosd(beta)+2.*in(y,1).*in(y,2).*ar.*br.*cosd(gamma);
            dif(y)=abs(Qobs(f)-Qcalc(y));
            mat(i,:)=[A(f,1) A(f,2) in(y,1) in(y,2) in(y,3) d(y) Qcalc(y) Qobs(f) dif(y)];
            %m=table(A(f,1),A(f,2),in(y,1),in(y,2),in(y,3),d(y),Qcalc(y),Qobs(f),dif(y))
            %fprintf('%.3f       %.3f    %d %d %d     %0.3f       %0.3f      %0.3f      %0.20f\n',A(f,1),A(f,2),in(y,1),in(y,2),in(y,3),d(y),Qcalc(y),Qobs(f),dif(y));
        %else
         %   fprintf('%0.3f - - -\n',A(f,2))
        i=i+1;
        end
    end
end
fprintf('Parametros de celda, espacio real: a=%f A, b=%f A, c=%f A, alpha= %f, beta=%f, gamma=%f \n',a,b,c,alpha,beta,gamma)
[a,b]=sort(mat(:,2));
mat2=mat(b,:);

Intensidad=mat2(:,1);
dostheta=mat2(:,2);
h=mat2(:,3);
k=mat2(:,4);
l=mat2(:,5);
d=mat2(:,6);
Q_calc=mat2(:,7);
Q_obs=mat2(:,8);
Diferencia=mat2(:,9);
tab=table(Intensidad,dostheta,h,k,l,d,Q_calc,Q_obs,Diferencia)
