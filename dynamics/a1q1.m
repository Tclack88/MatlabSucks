syms t theta_1(t) theta_2(t) theta_dot_3 L_1 L_2 L_3 rho
%% Test 1: Enter you student ID: StudentID
StudentID = 1298022

% Test 2: "clear" or "clearvars" (Pretest)
% Do not include any form of a "clear" command in your code


% Test 3: R10 (Pretest)
% Find the rotation matrix that brings a vector from frame {0} into frame {1}
R10 = [cos(theta_1) sin(theta_1) 0;
       -sin(theta_1) cos(theta_1)  0;
           0                0     1]

% Test 4: R21 (Pretest)
% Find the rotation matrix that brings a vector from frame {1} into frame {2}
R21 = [cos(theta_2) 0 sin(theta_2);
           0        1        0;
       -sin(theta_2) 0 cos(theta_2)]

% Test 5: R32
% Find the rotation matrix that brings a vector from frame {2} into frame {3}
theta_3 = theta_dot_3*t
R32 = [1        0         0;
       0 cos(theta_3) sin(theta_3);
       0 -sin(theta_3) cos(theta_3)]

% Test 6: R30
% Find the rotation matrix that brings a vector from frame {0} into frame {3}
R30 = R32*R21*R10

% Test 7: rOA_2 (Pretest)
% Find the position vector from point O to point A in frame {2}
rOA_2 = R21*[0;L_1;0]

% Test 8: rAB_2 (Pretest)
% Find the position vector from point A to point B in frame {2}
rAB_2 = [L_2; 0; 0]

% Test 9: rBC_3
% Find the position vector from point B to point C in frame {3}
rBC_3 = [0; 0; -L_3] % or is it [0;L_3;0] ?

% Test 10: rOC_2
% Find the position vector from point O to point C in frame {2}
rOC_2 = rOA_2 + rAB_2 + rBC_3

% Test 11: w1_1 (Pretest)
% Find the absolute angular velocity of the link that frame {1} is attached to, and express that vector in frame {1}
w1_1 = [0;0;diff(theta_1)]

% Test 12: w21_2 (Pretest)
% Find the angular velocity of the link that frame {2} is attached to, relative to the link that frame {1} is attached to, and express that vector in frame {2}
w21_2 = [0;-diff(theta_2);0]

% Test 13: w32_3
% Find the angular velocity of the link that frame {3} is attached to, relative to the link that frame {2} is attached to, and express that vector in frame {3}
w32_3 = [theta_dot_3;0;0]

% Test 14: w3_3
% Find the absolute angular velocity of the link that frame {3} is attached to, and express that vector in frame {3}
w3_3 = R32*w21_2 + R32*R21*w1_1

% Test 15: rOC_dot_2 (Pretest)
% Find the absolute linear velocity of point C in frame {2}
w2_2 = w21_2 + R21*w1_1
rOC_dot_2 = diff(rOC_2) + cross(w2_2,rOC_2)
% Test 16: rOC_ddot_2
% Find the absolute linear acceleration of point C in frame {2}


% Test 17: w3_dot_2
% Find the absolute angular acceleration of the link that frame {3} is attached to, and express that vector in frame {2}


