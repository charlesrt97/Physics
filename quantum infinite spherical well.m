% Infinite spherical potential well

function quantumsphwell
% valores iniciales
ggrid=50;
delta = pi/ggrid;
theta = 0 : delta : pi;
phi = 0 : 2*delta : 2*pi; 
[Phi,Theta] = meshgrid(phi,theta);
g=exp(-cos(Phi+2*Theta)).*(pi-Theta).*Theta;
Nrmin=g.^2.*sin(Theta);
Nrm=1/sqrt(trapz(phi,trapz(theta,Nrmin)));
g=g*Nrm;
figure(1)
[X,Y,Z] = sph2cart(Phi,pi/2-Theta,1);
clf;
surf(X,Y,Z,g);
colorbar; 
shading interp; 
daspect([1 1 1]); 
axis tight; 
view([70 25]);
title('3D Plot of g on the sphere zone(70,25)')
hold on
figure(2)
surf(X,Y,Z,g);
colorbar; 
shading interp; 
daspect([1 1 1]); 
axis tight; 
view([40 35]);
title('3D Plot of g on the sphere zone(40,35)')

%% C?lculo angular
%Arm?nicos
l=0;
mr=0;
R=diag(10)-diag(10);
m=abs(mr);
f00=g.*conj(harmonicY(0,0,Theta,Phi)).*sin(Theta);
c00=trapz(phi,trapz(theta,f00));
R(1,1)=c00;
Prob(1,1)=norm(R(1,1)).^2;
Fun=c00*harmonicY(0,0,Theta,Phi);
j=1;
  %Rutina para hacer la expansi?n

a=(Fun-g)./g;
  a(isinf(a))=0;
  b=max((max(abs(a))));
  find(a==b);
  c=abs(a(201)); %Se encuentra el m?ximo valor absoluto y se eval?a ese valor, si ese valor es menor a la condici?n, todos los dem?s igual lo ser?n
while c>0.01
    l=l+1;
    i=1;
    for mr= -l:l
    Int=conj(harmonicY(l,mr,Theta,Phi)).*g.*sin(Theta);
    R(i,j+1)=trapz(phi,trapz(theta,Int));
    Fun=Fun+R(i,j+1).*harmonicY(l,mr,Theta,Phi);
    Prob(i,j+1)=norm(R(i,j+1)).^2;
    i=i+1;
    a=(Fun-g)./g;
    a(isinf(a))=0;
    c=abs(a(201));
    end
j=j+1;
end
R;
l;
Prob; %La probabilidad de los ?ndices donde j son los valores verticales y m son los valores horizontales
figure(3)
[X,Y,Z] = sph2cart(Phi,pi/2-Theta,1);
clf;
Fun=abs(Fun);
surf(X,Y,Z,Fun);
colorbar; 
shading interp; 
daspect([1 1 1]); 
axis tight; 
view([70 25]);
title('3D Plot of spherical harmonic on the sphere zone(70,25)')
figure(4)
surf(X,Y,Z,Fun);
colorbar; 
shading interp; 
daspect([1 1 1]); 
axis tight; 
view([40 35]);
title('3D Plot of spherical harmonic on the sphere zone(40,35)')
%% C?lculo radial
i=1;
a=1;
h=1;
M=1/2;
r_=linspace(0.0001,1,ggrid+1);
[R,Theta,Phi]=meshgrid(r_,theta,phi);
I2=(a-R.^2).^2.*g.^2.*sin(Theta).*R.^2;
B=diag(20)-diag(20);
B=besselzero(1/2,20);
for d=0:19
    B(i,:)=besselzero(d+1/2,20);
    i=i+1;
end
B;% Matriz con bessel, donde i tiene los diferentes par?metros bessel y j tiene los diferentes ceros
A=1./sqrt(trapz(r_,trapz(phi,trapz(theta,I2))));
n=1;
l=0;
mr=0;
Q=A*g.*(a-R.^2);
Pn=0;
Func=0;
Enex=0;
prob=diag(10)-diag(10);
j=1;
i=1;
for l=0:4
    for mr= -l:l
        for n=1:3
        I_alp=(a^3)*((sqrt(2./(pi*B(l+1,n))).*besselj(l+3/2,B(l+1,n))).^2)/2;
        Alph(j)=sqrt(1/I_alp);
        psi=Alph(j)*sqrt(2*a./(pi*R.*B(l+1,n))).*besselj(l+1/2,B(l+1,n).*R./a).*harmonicY(l,mr,Theta,Phi);
        Intg=conj(psi).*Q.*R.^2.*sin(Theta);
        Cn(i)=trapz(r_,trapz(phi,trapz(theta,Intg)));
        prob(i)=(abs(Cn(i))).^2;
        Pn=Pn+abs((Cn(i)))^2;
        Func=abs(Cn(i).*psi).^2;
        Func=Func+Cn(i).*psi;
        Enex=Enex+abs(Cn(i))^2.*(h^2)/(2*M*a^2)*B(l+1,n)^2;
        i=i+1;
        j=j+1;
        end
        
    end
end
format long 
Cn
prob
n=1 ;l=0;mr=0; %genera la probabilidad m?s alta de 0.64;
Enex; %valor esperado de la energ?a
prob=max(prob); %mayor probabilidad del estado
Pn;
Enpalt=(h^2)/(2*M*a^2)*B(1,1)^2; %mayor probabilidad de energ?a

%% Orbital L

   
Prob=max(max(Prob)) %La mayor probabilidad est? en el estado l=0, m=0, por lo tanto
LzY=h*m*harmonicY(0,0,Theta,Phi);

%L2

vl2=0;
for l=0:6
    for m=-l:1:l
        fj=g.*conj(harmonicY(l,m,Theta,Phi)).*sin(Theta);
        o=trapz(phi,trapz(theta,fj));
        h=(abs(o)).^2.*l*(l+1);
        vl2=vl2+h;
    end
    l=l+1/2;
end
vl2;


%% C?lculo temporal
%Datos
A=1./sqrt(trapz(r_,trapz(phi,trapz(theta,I2))));
n=1;
l=0;
mr=0;
Q=A*g.*(a-R.^2);
a=1;
[Phi,Theta] = meshgrid(phi,theta);
X=abs(X);
Y=abs(Y);
Z=abs(Z);
r_=linspace(0.0001,1,ggrid+1);
[X,Y,Z] = sph2cart(Phi,pi/2-Theta,1/2);
%Gr?fica esfera
for t=0:.1:10
   Func=0;
    i=1;
    j=1;
   for l=0:4
    for mr= -l:l
        for n=1:3
        psi=Alph(j).*sqrt(2*a./(pi*a/2.*B(l+1,n))).*besselj(l+1/2,B(l+1,n).*a/2./a).*harmonicY(l,mr,Theta,Phi);
        psi=abs(psi).^2;
       Cg(i)=abs(Cn(i)).^2;
       Func=Func+Cg(i).*psi;
        Func=Func+(psi*Cn(i)*real(exp(1i*t*B(l+1,n)^2)));
        i=i+1;
        j=j+1;
        end
     end
   end
Func2=abs(Func)^2;
normF2=Func2-min(Func2(:));
normF2=normF2./max(normF2(:));
figure(5)
surf(X,Y,Z,normF2);
colorbar; 
shading interp;
daspect([1 1 1]);
view([70 25]);
axis tight; 
pause(0.01);

figure(15)
surf(X,Y,Z,normF2);
colorbar; 
shading interp;
daspect([1 1 1]);
view(2);
axis tight; 
pause(0.01);


end

r_=linspace(0.0001,1,ggrid+1);
%Gr?fica onda
for t=0:.0001:10
    i=1;j=1;
    func1=0;
 for l=0:4
        for mr=-l:l
            for n=1:3
                m=abs(mr);
                psi2=Alph(j).*sqrt(2*a./(pi*r_.*B(l+1,n))).*besselj(l+1/2,B(l+1,n).*r_./a).*sqrt((2*l+1)/(4*pi));
                
                func=(psi2)*Cn(i)*real(exp(B(l+1,n)^2*t*(-1i)));
                func1=func1+func;
                j=j+1;
                i=1+1; 
            end
        end 
 end
func1=abs(func1).^2;
norf=func1-min(func1(:));
norf=norf./max(norf(:));
figure(6)
grid on
plot(r_,norf); xlim([0 1]); 
ylim([0 1]);
pause(0.1);  
end

end
