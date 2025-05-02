%Regolatore manipolatore
clear all 
clc 


syms q1 q2 qp_1 qp_2 v1 v2

theta1 = 10.6125; 
theta2 = 0.85; 
theta3 = 2.25;
theta4 = 1.6;
theta5 = 80.9325;
theta6 = 14.7150;
F1 = 15;
F2 = 10;
lambda = 2;

M = [theta1 + 2*theta3*q2 theta2+theta3*cos(q2) ; ...
     theta2+theta3*cos(q2) theta4];
C = [-2*theta3*qp_2*sin(q2)  -theta3*qp_2*sin(q2) ; ...
      theta3*qp_1*sin(q2)              0              ];
g = [theta5*cos(q1)+theta6*cos(q1+q2) ; theta6*cos(q1+q2)];
F = [F1*qp_1 ; F2*qp_2];

x1 = [q1 ; q2];
x2 = [qp_1 ; qp_2];
v  = [v1 ; v2];

fp = [x2 ; (-inv(M) * C * x2 - inv(M)* F + inv(M) * v)];

%Calcolo jacobiano 
J = jacobian(fp,[x1 ; x2]);
J1 = jacobian(fp,v);

%qr = [pi/4 pi/4]
A_eq = double(subs(J,[x1  x2  v] ,[pi/4 0  0; pi/4  0  0 ]));
B_eq = double(subs(J1,[x1  x2  v] ,[pi/4 0  0; pi/4  0  0 ]));

%Calcolo autovalori
autovalori = eig(A_eq);
%Controllabilità
rank(ctrb(A_eq,B_eq)); %rank = 4

%LQR
A_new = A_eq + lambda/2 * eye(4);
Q = eye(4);
R = eye(2);
K = lqr(A_new,B_eq,Q,R);
kp = K(:,1:2); 
kd = K(:,3:4);

%Verifica matrice simmetriche 
eig(kp + kp');
eig(kd + kd');

%Jacobiano robot
l1 = 1.5;
l2 = 1;
f_q = [l1*cos(q1)+l2*cos(q1+q2);l1*sin(q1)+l2*sin(q1+q2)];
J = jacobian(f_q,[q1 q2]);

%Funzione di trasferimento 
tau = 1;




