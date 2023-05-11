syms m1 R vx1 vy1 vz1 wx1 wy1 wz1 % thin disk
syms m2 a b c v vx2 vy2 vz2 wx2 wy2 wz2 % right tetrahedron
syms t theta1 e F_hat
%% Meaning of the variables %%
% vx1 vy1 vz1 linear velocity components of the disc after collision in frame {1}
% wx1 wy1 wz1 angular velocity components of the disc after collision in frame {1}
% vx2 vy2 vz2 linear velocity components of the tetrahedron after collision in frame {2}
% wx2 wy2 wz2 angular velocity components of the tetrahedron after collision in frame {2}
% F_hat is the magnitude of the impulsive force
% theta (picture) was renamed to theta1 (syms)

%% Test 1: Enter you student ID: StudentID
StudentID = 1298022;

% Test 2: "clear" or "clearvars" (Pretest)
% Do not include any form of a "clear" command in your code


% Test 3: sumF1_1_hat (Pretest)
% Find the left-hand side  (LHS)  of the linear Impulse-Momentum equations for the disk in frame {1}
R21 = [cos(theta1) 0 -sin(theta1);
            0      1      0;
       sin(theta1) 0  cos(theta1)];
R12 = R21.';

sumF1_1_hat = -[0;0;F_hat]

% Test 4: deltaP_disk 
% Find the right-hand (RHS) side of the linear Impulse-Momentum equations for the disk in frame {1}
% Initially it's standing still => p_0 = 0

deltaP_disk = m1*[vx1; vy1; vz1]

% Test 5: sumMG11_1_hat (Pretest)
% Find the left-hand side of the angular Impulse-Momentum equations for the disk about its centre of mass in frame {1}
R_GC  = [R/2; R/2; 0]; % C for "contact point"
sumMG11_1_hat = cross(R_GC, sumF1_1_hat)

% Test 6: deltaH_disk
% Find the right-hand side of the angular Impulse-Momentum equations for the disk about its centre of mass in frame {1}
Ix = m1/4*R^2;
Iy = m1/4*R^2;
Iz = m1/2*R^2;
Idisk = [Ix 0 0;
         0 Iy 0;
         0  0 Iz];
w_disk = [wx1;wy1;wz1];
deltaH_disk = Idisk*w_disk

% Test 7: sumF2_2_hat_tetra (Pretest)
% Find the left-hand side of the linear Impulse-Momentum equations for the right tetrahedron in frame {2}

sumF2_2_hat_tetra = R21*(-sumF1_1_hat)

% Test 8: deltaP_tetra
% Find the right-hand side of the linear Impulse-Momentum equations for the tetrahedron in frame {2}
% in frame 2, initial momentum of tetrahedron is 0 (when all of a sudden this rude disk bumps into it)
% so deltaP is final - initial = final - 0 = final
deltaP_tetra = m2*([vx2;vy2;vz2] - [v;0;0])

% Test 9: sumMG22_2_hat
% Find the left-hand side of the angular Impulse-Momentum equations for the tetrahedron about its centre of mass in frame {2}
r_Gpoint = [3*a/4;-b/4;-c/4] % center of mass to point,(bring y and z to x axis, then continue a - a/4 in x direction)
%sumMG22_2_hat = R21*(-sumMG11_1_hat)  % didn't think about this carefully, just did what I did for force
sumMG22_2_hat = cross(r_Gpoint, sumF2_2_hat_tetra)
% Test 10: deltaH_tetra
% Find the right-hand side of the angular Impulse-Momentum equations for the tetrahedron about its centre of mass in frame {2}.
Ix = 3/80*m2*(b^2+c^2);
Iy = 3/80*m2*(a^2+c^2);
Iz = 3/80*m2*(a^2+b^2);
Ixy = 1/80*m2*a*b;  % these are positive (from ed discussion and trying the operation on a known result)
Ixz = 1/80*m2*a*c;
Iyz = 1/80*m2*b*c;
I_tetra = [Ix Ixy Ixz;
           Ixy Iy Iyz;
           Ixz Iyz Iz]

w_tetra = [wx2; wy2; wz2]  % initial w is zero so final - initial = delta

deltaH_tetra = I_tetra*w_tetra

%%%%%%%%%%%%%%%%%% IMPOTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For technical reasons, you must not redefine the variables given at the beginning. You must create new variables for the next questions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hint: define equations using declarations like 'eqn = LHS == RHS' and then use solve() to assign values to your new variables.

%Test 11: rOG1_dot_plus_1
% Determine the absolute linear velocity of the centre of mass of the disc after the collision in frame {1}.
% coeff restitution equation
z_tilde = [0 0 1]*R12  % R21 is R12 transpose 0 0 1 gets the 3rd column of that which is what we need I think???
v_before = z_tilde*[v;0;0] %denominator: approaching velocity rotated to frame 1
v_after = z_tilde*[vx2;vy2;vz2] - vz1 % numerator: subtract 2 from 1 to get to with respect to 1 whih is what we need I think????

rest_eqn = e == v_after/v_before

% force-momentum disk (frame 1)
deltaP_disk_new = m1*[vx1; vy1; vz1]
sumF1_1_hat_new = -[0;0;F_hat]
Fd_eqn = sumF1_1_hat_new  == deltaP_disk_new



% force-momentum tetra (frame 2) rotate to frame 1  (by not rotating to 2 as above)
%deltaP_tetra = m2*R12*[vx2;vy2;vz2]   % % subtract v????
deltaP_tetra_new = m2*R12*([vx2;vy2;vz2] - [v;0;0])   % % subtract v????
sumF2_2_hat_tetra_new = -sumF1_1_hat_new
Ft_eqn = sumF2_2_hat_tetra_new == deltaP_tetra_new

% moment-h disk  (frame 1)
sumMG11_1_hat_new = cross(R_GC, sumF1_1_hat)
deltaH_disk_new = Idisk*w_disk
Md_eqn = sumMG11_1_hat_new == deltaH_disk_new

% moment-h tetra  (frame 2) rotate to frame 1 (by not rotating to 2 as above)
sumMG22_2_hat_new = R12*cross(r_Gpoint, sumF2_2_hat_tetra)
deltaH_tetra_new = R12*I_tetra*w_tetra
Mt_eqn = sumMG22_2_hat_new == deltaH_tetra_new

sol = solve([Fd_eqn;Ft_eqn;Md_eqn;Mt_eqn;rest_eqn],[vx1,vy1,vz1,vx2,vy2,vz2,wx1,wy1,wz1,wx2,wy2,wz2,F_hat])


rOG1_dot_plus_1 = [sol.vx1;sol.vy1;sol.vz1]
%Test 12: w1_plus_1
% Determine the absolute angular velocity of the disc after the collision in frame {1}.

w1_plus_1 = [sol.wx1;sol.wy1;sol.wz1]

%Test 13: rOG2_dot_plus_2
% Determine the absolute linear velocity of the centre of mass of the tetrahedron after the collision in frame {2}.

rOG2_dot_plus_2 = [sol.vx2;sol.vy2;sol.vz2]

%Test 14: w2_plus_2
% Determine the absolute angular velocity the tetrahedron after the collision in frame {2}.

w2_plus_2 = [sol.wx2;sol.wy2;sol.wz2]

%Test 15: ImpulsiveForce
% Determine the magnitude of the impulsive force.

ImpulsiveForce= sol.F_hat




