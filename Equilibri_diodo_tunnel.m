%Analisi degli equilibri circuito con diodo tunnel 
clear all 
clc

%Definizione incognite
syms x1 x2 

%Definizione parametri fisici
C = 1; %capacità del condensatore 
L = 1; %induttanza 
R = 1; %resistenza
E = 1; %tensione di alimentazione

%Definizione del sistema 
h = @(x) 17.76*x - 103.79*x.^2 + 229.62*x.^3 - 226.31*x.^4 + 83.72*x.^5;
eq_1 = h(x1) == x2;
eq_2 = x1 + R*x2 == E;
sol = solve([eq_1, eq_2], [x1, x2]);

x1_dot = @(x1, x2) (-h(x1) + x2) / C;
x2_dot = @(x1, x2) (-x1 - R*x2 + E) / L ;
f1 = x1_dot(x1, x2);
f2 = x2_dot(x1,x2);

%Calcolo del jacobiano
J = jacobian([f1, f2], [x1, x2]);


%Ciclo FOR per il calcolo degli equilibri 
for i = 1:length(sol.x1)
    disp("-------- STEP "+ i + " ---------------")
    disp("X1 all'equilibrio: " + double(sol.x1(i)));
    disp("X2 all'equilibrio: " + double(sol.x2(i)));
    J_eq = double(subs(J,{x1, x2}, {sol.x1(i), sol.x2(i)}));
    disp("Jacobiano all'equilibrio: ")
    disp(J_eq)

    eigval = eig(J_eq);
    disp("Autovalori della matrice Jacobiana: ")
    disp(eigval)
    

    disp("--------------------------------")
end
