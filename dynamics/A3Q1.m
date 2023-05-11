syms t h(t) theta_1(t) theta_2(t) R m g ...
    F_A_x2 F_A_y2 F_A_z2 M_A_x2 M_A_y2 M_A_z2 ...
    F_B_x2 F_B_y2 F_B_z2 M_B_x2 M_B_y2 M_B_z2 ...
    F_C_x2 F_C_y2 F_C_z2 M_C_x2 M_C_y2 M_C_z2 ...
    F_O_x2 F_O_y2 F_O_z2 M_O_x2 M_O_y2 M_O_z2
% Meaning of the variables
% F_i_x2 F_i_y2 F_i_z2 M_i_x2 M_i_y2 M_i_z2 General vector components of constraint forces and moments in frame at point "i" (i=A,B,C,O) in frame {2}
% You must use the components needed to construct the vector for each joint. 

%The following quantiites were tested in assignment 1. Here they are given for you to use
% Rotation matrices
R01 = [cos(theta_1) 0 sin(theta_1); 0 1 0; -sin(theta_1) 0 cos(theta_1)];
R12 = [cos(theta_2) sin(theta_2) 0; -sin(theta_2) cos(theta_2) 0; 0 0 1];
% Position vectors
rOC_1 = [0; -h; 0];
rCG_2 = [4*R/(3*pi); 0; 0];

% Angular velocities
w1_1 = [0; diff(theta_1, t); 0];
w21_2 = [0; 0; -diff(theta_2, t)];

%% Test 1: Enter you student ID: StudentID
StudentID = 1298022

% Test 2: "clear" or "clearvars" (Pretest)
% Do not include any form of a "clear" command in your code
R0 = [-1 0 0; 0 1 0;0 0 -1]
syms x y z
I0 = [x 0 0; 0 y 0;0 0 z]
res = R0*I0*R0

% Test 3: p_dot_2 (Pretest)
% Find the time-derivative of the linear momentum of the plate in frame {2}
R21 = R12'
rOG_2 = R21*rOC_1 + rCG_2
w2_2 = w21_2 + R21*w1_1
v_2 = diff(rOG_2,t) + cross(w2_2,rOG_2)
p_2 = m*v_2
p_dot_2 = diff(p_2,t) + cross(w2_2,p_2)

% Test 4: I_plate_2
% Find the inertia tensor of the plate about its centre of mass in frame {2}
% http://mechanicsmap.psu.edu/websites/centroidtables/centroids2D/centroids2D.html
%Iy = 0 % approximate assuming thin plate
%Iz = (pi/8 - 8/(9*pi))*R^4
%Ix = pi/8*R^4
Iy = (1/2 - 16/(9*pi^2))*m*R^2
Iz = (1/4 - 16/(9*pi^2))*m*R^2
Ix = 1/4*m*R^2
I_plate_2 = [Ix 0 0;
             0 Iy 0;
             0  0 Iz]
% Test 5: h_dot_2
% Find the time-derivative of the angular momentum of the plate about its centre of mass in frame {2}
h_2 = I_plate_2*w2_2
h_dot_2 = diff(h_2) + cross(w2_2,h_2)

% Test 6: F_A_2 (Pretest)
% Give the vector of the constraint force acting on the plate at joint A in frame {2}. Remember to use the number zero if a component will always be zero

F_A_2 = [F_A_x2;F_A_y2;F_A_z2]

% Test 7: M_A_2 (Pretest)
% Give the vector of the constraint moment acting on the plate at joint A in frame {2}. Remember to use the number zero if a component will always be zero

M_A_2 = [M_A_x2;M_A_y2;0]

% Test 8: sumF_2
% Find the sum of all forces acting on the plate in frame {2}
F_B_2 = [F_B_x2;F_B_y2;F_B_z2]
F_G_1 = [0; -m*g; 0] % gravity active in frame 1 (or 0)
sumF_2 = F_A_2 + F_B_2 + R21*F_G_1
R21

% Test 9: sumM_2
% Find the sum of all moments acting on the plate about its centre of mass in frame {2}
M_B_2 = [M_B_x2;M_B_y2;0]
rGA_2 = -rCG_2 + [0;0;R]
rGB_2 = -rCG_2 + [0;0;-R]
sumM_2 = M_B_2 + M_A_2 + cross(rGA_2,F_A_2) + cross(rGB_2,F_B_2)

