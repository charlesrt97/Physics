function dx=myeq(t,x)
%%%%%% --- Funcion myeq --- %%%%%%

% Se establecen los valores de las constantes
l1=1;    % Longitud barra 1
l2=1.5;  % Longitud barra 2
l3=2;    % Longitud barra 3
m1=1;    % Valor masa 1
m2=1.5;  % Valor masa 2
m3=2;    % Valor masa 3

% Estas expresiones son los elementos de la matriz M. Para mayor 
% simplicidad se hizo a=M11, d=M12, e=M13, b=M22, f=M23, c=M33.
% Para mas info ver documento PDF
a=(m1+m2+m3)*l1^2;               % M11
d=(m2+m3)*l1*l2*cos(x(1)-x(2));  % M12
e=m3*l1*l3*cos(x(1)-x(3));       % M13
b=(m2+m3)*l2^2;                  % M22
f=m3*l2*l3*cos(x(2)-x(3));       % M23
c=m3*l3^2;                       % M33
% Estas expresiones son los elementos del vector Q. Para mayor simplicidad,
% se hizo g=Q1, h=Q2, j=Q3. Para mas info ver PDF
g=-(m2+m3)*l1*l2*(x(5))^2*sin(x(1)-x(2))-m3*l1*l3*(x(6))^2*sin(x(1)-x(3))-(m1+m2+m3)*9.81*l1*sin(x(1)); % Q1
h=(m2+m3)*l1*l2*(x(4))^2*sin(x(1)-x(2))-m3*l2*l3*(x(6))^2*sin(x(2)-x(3))-(m2+m3)*9.81*l2*sin(x(2));     % Q2
j=m3*l1*l3*(x(4))^2*sin(x(1)-x(3))+m3*l2*l3*(x(5))^2*sin(x(2)-x(3))-m3*9.81*l3*sin(x(3));               % Q3
% Las tres primeras ecuaciones son las primeras derivadas temporales, y las
% ultimas tres expresiones son equivalentes a las segundas derivadas temporales.
% Para mas info ver documento PDF
dx=[x(4);x(5);x(6);(b*c*g-f^2*g-c*d*h+e*f*h-b*e*j+d*f*j)/(a*b*c-c*d^2-b*e^2+2*d*e*f-a*f^2);(e^2*h+c*(d*g-a*h)+a*f*j-e*(f*g+d*j))/(c*d^2+e*(b*e-2*d*f)+a*(-b*c+f^2));(a*f*h-d*(f*g+e*h)+d^2*j+b*(e*g-a*j))/(c*d^2+e*(b*e-2*d*f)+a*(-b*c+f^2))];
