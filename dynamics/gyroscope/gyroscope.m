syms t g real
syms alpha1(t) beta1(t) gamma1(t)
syms m1 m2 m3 m4 real
syms h d1 d2 r1 r2 r3 r4 t1 t2 t3 L real
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Meaning variables% 
% m1 mas of big disc rotor (disc 1)
% m2 mass small disc on rotating axle (disc 2)
% m3 mass of small disck on rotor (disc 3)
% m4 mass of shaft (bar)
% h distance from base O to joint A
% d1 distance from joint A to start of small cylinder on rotor
% d2 distance from joint A to start of cylinder on rotating axle
% r1 radius of big disc, rotor
% r2 radius of cylinder on rotationg axle
% r3 radius of small cylinder on rotor
% r4 radius of shaft on rotating axle
% t1 thickness big disc rotor
% t2 thickness small cylinder on rotating axle
% t3 thickness small cylinder or rotor
% L length of the shaft on rotating axle
% a b c dimentions of rectangular hole
% d5 distance to centre of mass of rectangular hole along the radius
L=d1+d2

% Test 1: "clear" or "clearvars" (Pretest)
% Do not include any form of a "clear" command in your code

%% Test 1: Enter you student ID: StudentID
StudentID = 1298022


% Test 3: R10 (Pretest)
% Find the rotation matrix from the fixed frame {0} to the rotating frame {1}

R10 = [cos(alpha1), sin(alpha1), 0;
      -sin(alpha1), cos(alpha1), 0;
           0,           0,       1];

% Test 4: R21 (Pretest)
% Find the rotation matrix from the rotating frame {1} to the rotating frame {2}

R12 = [1     0          0;
       0 cos(beta1)  sin(beta1);
       0 -sin(beta1) cos(beta1)];

% Test 5: R32 (Pretest)
% Find the rotation matrix from the rotating frame {2} to the rotating frame {3}

R32 = [cos(gamma1) 0 sin(gamma1);
            0      1      0;
      -sin(gamma1) 0 cos(gamma1)];

% Test 6: w3_2 (Pretest)
% Find the absolute angular velocity of frame 3 in frame {2}
w1_1 = [0; 0; diff(alpha1)];
w21_2 = [-diff(beta1); 0; 0];
w32_2 = [0; -diff(gamma1);0];
R21 = R12.'
w3_2 = w32_2 + w21_2 + R21*w1_1;

% Test 7: rACbar_2  (Pretest)
% Find the position vector from point A to the centre of mass of the bar in frame {2}
dACbar = -d2 + L/2;
rACbar_2 = [0; -dACbar; 0];

% Test 8: rAD2_2  (Pretest)
% Find the position vector from point A to the centre of mass of disc 2 in frame {2}

rAD2_2 = [0; d2 + t2/2; 0];

% % Test 9: rAC2_2  (Test)
% Find the position vector from point A to the centre of mass (C2) of the combined structure of the rotating axle in frame {2}
rAC2_2 = (rACbar_2*m4 + rAD2_2*m2)/(m4+m2);

% Test 10: I4_2 (Pretest)
% Find the inertia tensor of the rod about its own centre of mass in frame {2}
Ix = 1/12*m4*L^2 + 1/4*m4*r4^2;
Iy = 1/2*m4*r4^2;
Iz = 1/12*m4*L^2 + 1/4*m4*r4^2;


I4_2 = [Ix 0 0;
        0 Iy 0;
        0  0 Iz];

% Test 11: I2_2 (Pretest)
% Find the inertia tensor of disc 2 about its own centre of mass in frame {2}.
Ix = 1/12*m2*t2^2 + 1/4*m2*r2^2;
Iy = 1/2*m2*r2^2;
Iz = 1/12*m2*t2^2 + 1/4*m2*r2^2;

I2_2 = [Ix 0 0;
        0 Iy 0;
        0  0 Iz];

% Test 12: I2_C2_2
% Find the inertia tensor of disc 2 about the centre of mass of the combined rotating axle (C2) in frame {2}.
r_AG2 = rAD2_2 - rAC2_2; % might be the opposite so this is exactly wrong. It's currently pointing from combined center to individual center
dx = r_AG2(1);
dy = r_AG2(2);
dz = r_AG2(3);
I_prime = m2*[dy^2+dz^2    -dx*dy    -dx*dz;
              -dx*dy    dx^2+dz^2    -dy*dz;
              -dx*dz     -dy*dz    dx^2+dy^2];

I2_C2_2 = I2_2 + I_prime;

% Test 13: I4_C2_2
% Find the inertia tensor of the rod about the centre of mass of the combined rotating axle (C2) in frame {2}.
r_AG4 = rACbar_2 - rAC2_2; % might be the opposite so this is exactly wrong.It's currently pointing from combined center to individual center
dx = r_AG4(1);
dy = r_AG4(2);
dz = r_AG4(3);
I_prime = m4*[dy^2+dz^2    -dx*dy    -dx*dz;
              -dx*dy    dx^2+dz^2    -dy*dz;
              -dx*dz     -dy*dz    dx^2+dy^2];
I4_C2_2 = I4_2 + I_prime;

% Test 14: I_RA_2 
% Find the inertia tensor of the combined rotating axle about its centre of mass {C2} in frame {2}.

I_RA_2 = I2_C2_2 + I4_C2_2;

% Test 15: rAD1_2 (Pretest)
% Find the position vector from point A to the centre of mass of disc 1 in frame {2}
dAG1 = d1 + t3+t1/2;
rAD1_2 = [0; -dAG1; 0];

% Test 16: rAD3_2 (Pretest) 
% Find the position vector from point A to the center of mass of disc 3 in frame {2}
dAG3 = d1 + t3/2;
rAD3_2 = [0; -dAG3; 0];


% Test 18: rAC1_2 (Test)
% Find the position vector from point A to the centre of mass of the combined rotor (C1) in frame {2}
R23 = R32';
rAC1_2 = (rAD1_2*m1 + rAD3_2*m3)/(m1+m3);

% Test 19: I_rotor_C1_3 (Test)
% Find the inertia tensor of the combined rotor about its centre of mass (C1) in frame {3}
% disc1
Ix1 = 1/12*m1*t1^2 + 1/4*m1*r1^2;
Iy1 = 1/2*m1*r1^2;
Iz1 = 1/12*m1*t1^2 + 1/4*m1*r1^2;
I1_3 = [Ix1 0 0;
        0 Iy1 0;
        0  0 Iz1];
r_CD1_3 =  R32*(rAD1_2 - rAC1_2)
dx = [1 0 0] * r_CD1_3;
dy = [0 1 0] * r_CD1_3;
dz = [0 0 1] * r_CD1_3;
I1_prime  = m1*[dy^2+dz^2    -dx*dy    -dx*dz;
           -dx*dy    dx^2+dz^2    -dy*dz;
           -dx*dz     -dy*dz    dx^2+dy^2];
I_rotor1 = I1_3 + I1_prime
% disc3
Ix3 = 1/12*m3*t3^2 + 1/4*m3*r3^2;
Iy3 = 1/2*m3*r3^2;
Iz3 = 1/12*m3*t3^2 + 1/4*m3*r3^2;
I3_3 = [Ix3 0 0;
        0 Iy3 0;
        0  0 Iz3];
r_CD3_3 = R32*(rAD3_2 - rAC1_2)
dx = [1 0 0] * r_CD3_3
dy = [0 1 0] * r_CD3_3
dz = [0 0 1] * r_CD3_3
I3_prime  = m3*[dy^2+dz^2    -dx*dy    -dx*dz;
           -dx*dy    dx^2+dz^2    -dy*dz;
           -dx*dz     -dy*dz    dx^2+dy^2];
I_rotor3 = I3_3 + I3_prime

%%% Here goes fucking nothing. Fingers crossed...
I_rotor_C1_3 = I_rotor1 + I_rotor3
% Test 20: p_rotor_dot_2 (Test)
% Find the time-derivative of the linear momentum of the rotor in frame {2} 
m_rot = m1+m3
rOA_1 = [0; 0; h];
rOC1_2 = R21*rOA_1 + rAC1_2
w2_2 = w21_2 + R21*w1_1
p_rotor_2 = m_rot*(diff(rOC1_2) + cross(w2_2,rOC1_2))
p_rotor_dot_2 = diff(p_rotor_2) + cross(w2_2,p_rotor_2)
% Test 21: p_RA_dot_2 (Test)
% Find the time-derivative of the linear momentum of the rotating axle in frame {2}
m_ax = m2+m4
rOC2_2 = R21*rOA_1 + rAC2_2
p_RA_2 = m_ax*(diff(rOC2_2) + cross(w2_2,rOC2_2))
p_RA_dot_2 = diff(p_RA_2) + cross(w2_2,p_RA_2)

% Test 22: h_RA_C2_dot_2 (Test)
% Find the time-derivative of the angular momentum of the rotating axle about its centre of mass in frame {2}
h_RA_C2_2 = I_RA_2*w2_2
h_RA_C2_dot_2 = diff(h_RA_C2_2) + cross(w2_2,h_RA_C2_2)

% Test 22: h_rotor_C1_dot_3 (Test)
% Find the time-derivative of the angular momentum of the combined rotor about its centre  of mass in frame {3}
w3_3 = R32*w3_2
h_rotor_C1_3 = I_rotor_C1_3*w3_3
h_rotor_C1_dot_3 = diff(h_rotor_C1_3) + cross(w3_3,h_rotor_C1_3)

h_rotor_simp = simplify(h_rotor_C1_dot_3)