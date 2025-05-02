%Analisi di stabilità del sistema lineare e non lineare
%massa-molla-smorzatore
close all;
clear all;
clc;

%1°caso lineare 
%Definizione parametri 

k = 10; %coefficiente di rigidità 
m = 100; %massa
b = 50; %coefficiente di smorzamento 


%Verifica della stabilità 
A = [0 1 ; -k/m -b/m];
lambda = eig(A);
disp("Autovalori del sistema lineare: ");
disp(lambda);


%Funzione di Lyapunov 
V = @(x1,x2) 0.5*k*(x1.^2) + 0.5*m*(x2.^2);

%variabili del sistema
x1 = linspace(-10,10,100);
x2 = linspace(-10,10,100);

[X1,X2] = meshgrid(x1,x2);

Z = V(X1,X2);

%Crea il grafico 3D
figure;
surf(x1, x2, Z);  % Usa surf per un grafico 3D con superfici colorate

% Aggiungi etichette agli assi e un titolo
zlim([0,8000]);
xlabel('x1');
ylabel('x2');
zlabel('V(x1,x2)');
title('Funzione di Lyapunov per il sistema lineare');

%Definizione matrice P 
P = [0.5*k 0 ; 0 0.5*m];
%Calcolo degli autovalori di P 
lambda_min = eigs(P,1,"smallestabs");
lambda_max = eigs(P,1,"largestabs");

%Definizione di delta,eps,beta 
eps = 1;
beta = lambda_min*eps^2;
delta = sqrt(lambda_min/lambda_max);

%------------------------------------------------------
%Plot dell'evoluzione dello stato 
% Definizione del sistema
sys_lineare = @(t, x) [x(2); -(b/m)*x(2) - (k/m)*x(1)];

% Condizioni iniziali (x1, x2)
condizioni_iniziali = [
    0.25, 0;    % Condizione iniziale 1
    0.5, 0.25;    % Condizione iniziale 2
    0.5, 0;    % Condizione iniziale 3
    0, 0.5;    % Condizione iniziale 4
    -0.5, +0.5;   % Condizione iniziale 5
    0, -0.5;   % Condizione iniziale 6
    -0.25, -0.25;   % Condizione iniziale 7
    0.5, -0.25    % Condizione iniziale 8
];

% Intervallo di tempo per la simulazione
tspan = [0 15];  % Da t = 0 a t = 15

% Array per memorizzare le traiettorie
traiettorie = cell(size(condizioni_iniziali, 1), 1);

% Risoluzione del sistema per ogni condizione iniziale
for i = 1:size(condizioni_iniziali, 1)
    [t, x] = ode45(sys_lineare, tspan, condizioni_iniziali(i, :));
    traiettorie{i} = x;  % Memorizza la traiettoria
end

% Creazione del grafico delle traiettorie
figure;
hold on;
for i = 1:length(traiettorie)
    plot(traiettorie{i}(:, 1), traiettorie{i}(:, 2), 'LineWidth', 2);  % Traiettoria i-esima
end
xlabel('x1');
ylabel('x2');
xlim([-1.25,1.25]);
ylim([-1,1]);
title('Traiettorie del sistema massa-molla-smorzatore lineare');
grid on;

% Creazione delle circonferenze
theta = linspace(0, 2*pi, 100);  % Angoli per disegnare le circonferenze

%Circonferenza epsilon
x_circle_epsilon = eps * cos(theta);
y_circle_epsilon = eps * sin(theta);

%Circonferenza delta
x_circle_delta = delta * cos(theta);
y_circle_delta = delta * sin(theta);

%Elisse di semiassi delta ed epsilon 
x_ellisse = delta * cos(theta);
y_ellisse = eps * sin(theta);


% Plot delle circonferenze
plot(x_circle_epsilon, y_circle_epsilon, 'g', 'LineWidth', 1.5);  % Circonferenza epsilon
plot(x_circle_delta, y_circle_delta, 'r', 'LineWidth', 1.5);      % Circonferenza delta
plot(x_ellisse,y_ellisse,'b','LineWidth',1.5);                    % Ellisse delta-epsilon


%Plot delle condizioni iniziali 
for i = 1:size(condizioni_iniziali, 1)
    plot(condizioni_iniziali(i, 1), condizioni_iniziali(i, 2), 'ko', ...
         'MarkerSize', 5, 'LineWidth', 2, 'MarkerFaceColor', 'black');
end

hold off;

%----------------------------

%2°caso non lineare
%Definizione dei parametri 

k1 = 10; %coefficiente di rigidità lineare
k2 = 1; %coefficiente di rigidità non lineare
m = 100; %massa
b1 = 50; %coefficiente di smorzamento lineare
b2 = 5; %coefficiente di smorzamento non lineare

%Funzione di Lyapunov 
V = @(x1,x2) 0.5*m*(x2.^2) + 0.5*k1*(x1.^2) + 0.25*k2*(x1.^2);

%Calcolo intersezioni con gli assi 
A = @(C) sqrt(((-0.5*k1) + sqrt((0.25*k1^2)+k2*C))/0.5*k2);
B = @(C) sqrt(2*C/m);

%Curva di livello 
C = 5; %energia
A_calc = A(C);
B_calc = B(C);

%calcolo di eps e delta 
eps = max(A_calc,B_calc);
delta= min(A_calc,B_calc);


%Rapporto tra eps e delta per grafico
C = linspace(0,10,10);
E = A(C);
D = B(C);
eps_delta_rapp = D./E;
eps_delta_rapp(1) = eps_delta_rapp(2);

figure;

plot(C,E,'r', 'LineWidth',1);
hold on;
plot(C,D,'b', 'LineWidth',1);
hold on;
plot(C,eps_delta_rapp,'g','LineWidth',1);
grid on;
legend('eps(C)', 'delta(C)', 'delta/eps(C)'); 
title('Rapporto delta e epsilon'); 
hold off;

%------------------------------------------------------
%Plot dell'evoluzione dello stato 
% Definizione del sistema
sys_non_lineare_1 = @(t, x) [x(2); -(k1/m)*x(2)-(k2/m)*(x(1).^3)-(b1/m)*(x(2).^3)-(b2/m)*x(1)];

% Intervallo di tempo per la simulazione
tspan = [0 60];  % Da t = 0 a t = 15 

% Risoluzione del sistema per ogni condizione iniziale
for i = 1:size(condizioni_iniziali, 1)
    [t, x] = ode45(sys_non_lineare_1, tspan, condizioni_iniziali(i, :));
    traiettorie{i} = x;  % Memorizza la traiettoria
end

% Creazione del grafico delle traiettorie
figure;
hold on;
for i = 1:length(traiettorie)
    plot(traiettorie{i}(:, 1), traiettorie{i}(:, 2), 'LineWidth', 2);  % Traiettoria i-esima
end
xlabel('x1');
ylabel('x2');
xlim([-1.35,1.35]);
ylim([-1,1]);
title('Traiettorie del sistema massa-molla-smorzatore non lineare CASO 1');
grid on;

% Creazione delle circonferenze
theta = linspace(0, 2*pi, 100);  % Angoli per disegnare le circonferenze

%Circonferenza epsilon
x_circle_epsilon = eps * cos(theta);
y_circle_epsilon = eps * sin(theta);

%Circonferenza delta
x_circle_delta = delta * cos(theta);
y_circle_delta = delta * sin(theta);

%Elisse di semiassi delta ed epsilon 
x_ellisse = eps * cos(theta);
y_ellisse = delta * sin(theta);


% Plot delle circonferenze
plot(x_circle_epsilon, y_circle_epsilon, 'g', 'LineWidth', 1.5);  % Circonferenza epsilon
plot(x_circle_delta, y_circle_delta, 'r', 'LineWidth', 1.5);      % Circonferenza delta
plot(x_ellisse,y_ellisse,'b','LineWidth',1.5);                    % Ellisse delta-epsilon


%Plot delle condizioni iniziali 
for i = 1:size(condizioni_iniziali, 1)
    plot(condizioni_iniziali(i, 1), condizioni_iniziali(i, 2), 'ko', ...
         'MarkerSize', 5, 'LineWidth', 2, 'MarkerFaceColor', 'black');
end

%----------------------------

%3°caso non lineare
%Definizione dei parametri 

k1 = 10; %coefficiente di rigidità lineare
k2 = 1; %coefficiente di rigidità non lineare
m = 100; %massa
b1 = 50; %coefficiente di smorzamento lineare
b2 = 5; %coefficiente di smorzamento non lineare

%Funzione di Lyapunov 
V = @(x1,x2) 0.5*m*(x2.^2) + 0.5*k1*(x1.^2) - 0.25*k2*(x1.^4);

%-------------BEGIN PLOT FUNZIONE LYAPUNOV----------------------

%variabili del sistema
x1 = linspace(-5,5,100);
x2 = linspace(-5,5,100);

[X1,X2] = meshgrid(x1,x2);

Z = V(X1,X2);

%Crea il grafico 3D
figure;
surf(x1, x2, Z);  % Usa surf per un grafico 3D con superfici colorate

% Aggiungi etichette agli assi e un titolo
zlim([-200,1400]);
xlabel('x1');
ylabel('x2');
zlabel('V(x1,x2)');
title('Funzione di Lyapunov per il sistema non lineare');

%-----------END PLOT FUNZIONE LYAPUNOV-----------------------

%Calcolo intersezioni con gli assi 
A = sqrt(k1/k2);
B = k1/sqrt(2*k2*m);

%Curva di livello 
C = (0.25*k1^2)/k2; %energia


%calcolo di eps e delta 
eps = max(A,B);
delta= min(A,B);


%------------------------------------------------------
%Plot dell'evoluzione dello stato 
% Definizione del sistema
sys_non_lineare_2 = @(t, x) [x(2); -(k1/m)*x(1)+(k2/m)*(x(1).^3)-(b1/m)*x(2)-(b2/m)*(x(2).^3)];

% Intervallo di tempo per la simulazione
tspan = [0 60];  % Da t = 0 a t = 15 

% Condizioni iniziali (x1, x2)
condizioni_iniziali = [
    0.25, 1;    % Condizione iniziale 1
    2, 0.5;    % Condizione iniziale 2
    0.75, 0.75;    % Condizione iniziale 3
    0, 0.5;    % Condizione iniziale 4
    -2, +0.5;   % Condizione iniziale 5
    0, -0.5;   % Condizione iniziale 6
    -2, 1;   % Condizione iniziale 7
    -1, -0.25    % Condizione iniziale 8
];

% Risoluzione del sistema per ogni condizione iniziale
for i = 1:size(condizioni_iniziali, 1)
    [t, x] = ode45(sys_non_lineare_2, tspan, condizioni_iniziali(i, :));
    traiettorie{i} = x;  % Memorizza la traiettoria
end

% Creazione del grafico delle traiettorie
figure;
hold on;
for i = 1:length(traiettorie)
    plot(traiettorie{i}(:, 1), traiettorie{i}(:, 2), 'LineWidth', 2);  % Traiettoria i-esima
end
xlabel('x1');
ylabel('x2');
xlim([-5,5]);
ylim([-5,5]);
title('Traiettorie del sistema massa-molla-smorzatore non lineare CASO 2');
grid on;

% Creazione delle circonferenze
theta = linspace(0, 2*pi, 100);  % Angoli per disegnare le circonferenze

%Circonferenza epsilon
x_circle_epsilon = eps * cos(theta);
y_circle_epsilon = eps * sin(theta);

%Circonferenza delta
x_circle_delta = delta * cos(theta);
y_circle_delta = delta * sin(theta);


% Plot delle circonferenze
plot(x_circle_epsilon, y_circle_epsilon, 'g', 'LineWidth', 1.5);  % Circonferenza epsilon
plot(x_circle_delta, y_circle_delta, 'r', 'LineWidth', 1.5);      % Circonferenza delta


%Plot delle condizioni iniziali 
for i = 1:size(condizioni_iniziali, 1)
    plot(condizioni_iniziali(i, 1), condizioni_iniziali(i, 2), 'ko', ...
         'MarkerSize', 5, 'LineWidth', 2, 'MarkerFaceColor', 'black');
end

%Plot linee verticali 
xline(A,'--',"LineWidth",1.5,"Color",'r');
xline(-A,'--',"LineWidth",1.5,"Color",'r');



%Coordinate
x_eq_1 = [-5, -A, 0, A , 5];
y_eq_1 = [-1, 0, B, 0,-1];
x_eq_2 = [-5, -A, 0, A , 5];
y_eq_2 = [1, 0, -B, 0, 1];

%Creazione spline
x_interp1 = linspace(min(x_eq_1), max(x_eq_1), 100);  % 100 punti interpolati
y_interp1 = interp1(x_eq_1, y_eq_1, x_interp1, 'spline');  % Interpolazione spline
x_interp2 = linspace(min(x_eq_2), max(x_eq_2), 100);  % 100 punti interpolati
y_interp2 = interp1(x_eq_2, y_eq_2, x_interp2, 'spline');  % Interpolazione spline

%Plot curva interpolata 1
plot(x_interp1, y_interp1, 'b', 'LineWidth', 2);  % Curva interpolata


%Aggiunta del testo per evidenziare livello di energia curva 1
for i = 1:length(x_eq_1)
    text(x_eq_1(i), y_eq_1(i)+0.25, "25", ...
         'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
end

%Plot curva interpolata 2
plot(x_interp2, y_interp2, 'b', 'LineWidth', 2);  % Curva interpolata

%Aggiunta del testo per evidenziare livello di energia curva 1
for i = 1:length(x_eq_2)
    text(x_eq_2(i), y_eq_2(i)-0.25, "25", ...
         'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
end


%----------------------------

%4°caso non lineare
%Definizione dei parametri 

k1 = 10; %coefficiente di rigidità lineare
k2 = 1; %coefficiente di rigidità non lineare
m = 100; %massa
b1 = 50; %coefficiente di smorzamento lineare
b2 = 5; %coefficiente di smorzamento non lineare

%Funzione di Lyapunov 
V = @(x1,x2) 0.5*m*(x2.^2) + 0.5*k1*(x1.^2) - 0.25*k2*(x1.^4);

%-------------BEGIN PLOT FUNZIONE LYAPUNOV----------------------

%variabili del sistema
x1 = linspace(-5,5,100);
x2 = linspace(-5,5,100);

[X1,X2] = meshgrid(x1,x2);

Z = V(X1,X2);

%Crea il grafico 3D
figure;
surf(x1, x2, Z);  % Usa surf per un grafico 3D con superfici colorate

% Aggiungi etichette agli assi e un titolo
zlim([-200,1400]);
xlabel('x1');
ylabel('x2');
zlabel('V(x1,x2)');
title('Funzione di Lyapunov per il sistema non lineare');

%-----------END PLOT FUNZIONE LYAPUNOV-----------------------

%Calcolo intersezioni con gli assi 
A = sqrt(k1/k2);
B = k1/sqrt(2*k2*m);

%Curva di livello 
C = (0.25*k1^2)/k2; %energia


%calcolo di eps e delta 
eps = max(A,B);
delta= min(A,B);


%------------------------------------------------------
%Plot dell'evoluzione dello stato 
% Definizione del sistema
sys_non_lineare_3 = @(t, x) [x(2); -(k1/m)*x(1)+(k2/m)*(x(1).^3)-(b1/m)*x(2)+(b2/m)*(x(2).^3)];

% Intervallo di tempo per la simulazione
tspan = [0 60];  % Da t = 0 a t = 15 

% Condizioni iniziali (x1, x2)
condizioni_iniziali = [
    0.25, 1;    % Condizione iniziale 1
    2, 0.5;    % Condizione iniziale 2
    0.75, 0.75;    % Condizione iniziale 3
    0, 0.5;    % Condizione iniziale 4
    -2, +0.5;   % Condizione iniziale 5
    0, -0.5;   % Condizione iniziale 6
    -2, 1;   % Condizione iniziale 7
    -1, -0.25    % Condizione iniziale 8
];

% Risoluzione del sistema per ogni condizione iniziale
for i = 1:size(condizioni_iniziali, 1)
    [t, x] = ode45(sys_non_lineare_3, tspan, condizioni_iniziali(i, :));
    traiettorie{i} = x;  % Memorizza la traiettoria
end

% Creazione del grafico delle traiettorie
figure;
hold on;
for i = 1:length(traiettorie)
    plot(traiettorie{i}(:, 1), traiettorie{i}(:, 2), 'LineWidth', 2);  % Traiettoria i-esima
end
xlabel('x1');
ylabel('x2');
xlim([-5,5]);
ylim([-5,5]);
title('Traiettorie del sistema massa-molla-smorzatore non lineare CASO 3');
grid on;

% Creazione delle circonferenze
theta = linspace(0, 2*pi, 100);  % Angoli per disegnare le circonferenze

%Circonferenza epsilon
x_circle_epsilon = eps * cos(theta);
y_circle_epsilon = eps * sin(theta);

%Circonferenza delta
x_circle_delta = delta * cos(theta);
y_circle_delta = delta * sin(theta);


% Plot delle circonferenze
plot(x_circle_epsilon, y_circle_epsilon, 'g', 'LineWidth', 1.5);  % Circonferenza epsilon
plot(x_circle_delta, y_circle_delta, 'r', 'LineWidth', 1.5);      % Circonferenza delta


%Plot delle condizioni iniziali 
for i = 1:size(condizioni_iniziali, 1)
    plot(condizioni_iniziali(i, 1), condizioni_iniziali(i, 2), 'ko', ...
         'MarkerSize', 5, 'LineWidth', 2, 'MarkerFaceColor', 'black');
end


%Definizione vertici 
x_vertex = [-A, A, A, -A ,-A];
y_vertex = [A, A, -A, -A, A];

%Plot area in cui è valida la condizione di LaSalle 
plot(x_vertex, y_vertex, 'r--', 'LineWidth', 2);  % Bordi del quadrato in blu

%Coordinate
x_eq_1 = [-5, -A, 0, A , 5];
y_eq_1 = [-1, 0, B, 0,-1];
x_eq_2 = [-5, -A, 0, A , 5];
y_eq_2 = [1, 0, -B, 0, 1];

%Creazione spline
x_interp1 = linspace(min(x_eq_1), max(x_eq_1), 100);  % 100 punti interpolati
y_interp1 = interp1(x_eq_1, y_eq_1, x_interp1, 'spline');  % Interpolazione spline
x_interp2 = linspace(min(x_eq_2), max(x_eq_2), 100);  % 100 punti interpolati
y_interp2 = interp1(x_eq_2, y_eq_2, x_interp2, 'spline');  % Interpolazione spline

%Plot curva interpolata 1
plot(x_interp1, y_interp1, 'b', 'LineWidth', 2);  % Curva interpolata


%Aggiunta del testo per evidenziare livello di energia curva 1
for i = 1:length(x_eq_1)
    text(x_eq_1(i), y_eq_1(i)+0.25, "25", ...
         'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
end

%Plot curva interpolata 2
plot(x_interp2, y_interp2, 'b', 'LineWidth', 2);  % Curva interpolata

%Aggiunta del testo per evidenziare livello di energia curva 1
for i = 1:length(x_eq_2)
    text(x_eq_2(i), y_eq_2(i)-0.25, "25", ...
         'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
end


%----------------------------
