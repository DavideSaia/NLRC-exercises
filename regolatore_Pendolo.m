%Regolatore Pendolo Informazione parziale 
clear
clc 

%Parametri fisici del sistema 
m = 0.5;%kg
M = 5;%massa carrello in kg
l = 1;%lunghezza asta in m
g = 9.81; %accelerazione di gravità
j_p = m*(l^2);%momento di inerzia

%Matrici del sistema 
A = [0 1 0 0 ; ...
     0 0 -m*g/M 0; ...
     0 0 0 1 ; ... 
     0 0 (M+m)*M*g*l/M*j_p 0];
B = [ 0      0 ; ...
     1/M -m*l/M*j_p ; 
      0      0 ; 
     -m*l/M*j_p m+M/M*j_p];

C  = [ 1 0 0 0 ; 0 0 1 0 ];


S = [0 1 0 ; -1 0 0 ; 0 0 0];%matrice dell'esosistema

Q = zeros(2,1);%matrice disturbo per l'uscita  

P = zeros(4,3);%matrice disturbo per lo stato 


eig(S) %verifica H1 (S antisimmetrica)
rho_AB = rank(ctrb(A,B)); %verifica H2 (coppia (A,B) controllabile)
disp("Controllabilità del sistema per la verifica di H2, rango: " + rho_AB)


%Matrici del sistema esteso 
Q_tilde = [-Q eye(2)];
P_tilde = [P zeros(4,1)];
C_tilde = -C;

A_e = [A P ; zeros(3,4) S];
C_e = [C Q_tilde];

rho_AC = rank(obsv(A_e,C_e)); %verifica H3*
disp("Osservabilità del sistema esteso per la verifica di H3*, rango: " + rho_AC)

%PROGETTO REGOLATORE
m=2;%Numero di ingressi
n=2;%Numero di uscite
q=4;%Numero di stati
r=3;%Dim esosistema

%costruzione matrici L1,L2,N
L1 = [-A -B ; C zeros(n,m)]; %dim (n+q * n+m)
L2 = blkdiag(eye(4),zeros(n,m)); % dim(n+q * n+m)
N = [P ; -Q_tilde];% dim (n+q)*r

%IMPLEMENTAZIONE FORMULA VETTORIZZAZIONE
Z = kron(S',L2) + kron(eye(3),L1);
X = Z\N(:); 
X = reshape(X,[6 r]);

%Estrazione delle matrici GAMMA e PIGRECO
T = X(1:4,:); %pigreco
G = X(5:6,:); %gamma


%CALCOLO di K per piazzare gli autovalori
R = 1;
Q1 = eye(4);
K = -lqr(A,B,Q1,R); % autovalori

%CALCOLO DI L dell'osservatore
L = G-K*T;

%Controllore Finale
poles = -abs(randn(size(A_e,1),1)*10);
G_e = place(A_e',C_e',poles)';%piazzamento autovalori

G0 = G_e(1:4,:);
G1 = G_e(5:7,:);

F = [A-G0*C+B*K P-G0*Q_tilde+B*(G-K*T) ; -G1*C S-G1*Q_tilde];
H = [K G-K*T];



