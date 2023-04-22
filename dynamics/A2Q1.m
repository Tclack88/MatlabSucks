syms theta(t) R a b c g rho 
%% Test 1: Enter you student ID: StudentID
StudentID = 1298022

% Test 2: "clear" or "clearvars" (Pretest)
% Do not include any form of a "clear" command in your code


% Test 3: R01 (Pretest)
% Find the rotation matrix from the rotating frame {1} to the fixed frame {0} 

R01 = [1       0         0;
       0 cos(theta) -sin(theta);
       0 sin(theta) cos(theta)]

% Test 4: ang_vel (Pretest)
% Find the absolute angular velocity of the eccentric rotating mass (ERM) in frame {1}

ang_vel = [diff(theta); 0; 0];

% Test 5: mass_semicyl (Pretest)
% Find the mass of the semicylinder part of the ERM

mass_semicyl = rho*(a*pi*R^2)/2

% Test 6: mass_cuboid
% Find the mass of the cuboid part of the ERM

mass_cuboid = rho*a*b*c

% Test 7: CoM_semicyl (Pretest)
% Find the position vector from the fixed point O to the centre of mass (CoM) of the semicylinder part in frame {1}

CoM_semicyl = [a/2;0;-4*R/(3*pi)]

% Test 8: CoM_cuboid
% Find the position vector from O to the CoM of the cuboid part in frame {1}

CoM_cuboid = [a/2;0;c/2]

% Test 9: CoM_total
% Find the position vector from O to the CoM of the entire ERM (i.e., both the semicylinder part and the cuboid part together) in frame {1}
% probably wrong. I need to do all parts separately I think. x should definitely be at a/2
% currently it's giving (a*pi*R^2)/4 + (a*b*c)/2)/((pi*R^2)/2 + b*c). I'm unsure why it's not reducing, 
% setting isequal(left,right) is giving False which seems like a problem
A_cube = b*c*a*rho
A_semicyl = .5*pi*R^2*a*rho
CoM_total = (CoM_semicyl*A_semicyl + CoM_cuboid*A_cube)/(A_semicyl + A_cube)    %%%%%% I have A, but is is masss?
% isequal(((a*pi*R^2)/4 + (a*b*c)/2)/((pi*R^2)/2 + b*c),a/2)
% Test 10: lin_mom_dot
% Find the time-derivative of the linear momentum of the ERM in frame {1}
ERM_mass = mass_semicyl + mass_cuboid
%lin_mom = ERM_mass*diff(CoM_total)   % this was my first thought
lin_mom = ERM_mass*cross(ang_vel,CoM_total) % but this makes more sense (only y component since frame 1 is moving withi t wxr should only give a y component)
lin_mom_dot = diff(lin_mom) + cross(ang_vel, lin_mom)

% Test 11: IG_semicyl (Pretest)
% Find the inertia tensor of the semicylinder part about its CoM in frame {1}. 
% from table: https://www.efunda.com/math/solids/solids_display.cfm?SolidName=HalfCircularCylinder
IxG_cyl = (1/2 - 16/(9*pi^2))*mass_semicyl*R^2;
IyG_cyl = mass_semicyl*((1/4 - 16/(9*pi^2))*R^2 + a^2/12);
IzG_cyl = mass_semicyl*(R^2/4 + a^2/12)
IG_semicyl = [IxG_cyl  0      0;
               0   IyG_cyl    0;
               0       0   IzG_cyl]


% Test 12: IG_cuboid
% Find the inertia tensor of the cuboid part about its CoM in frame {1}
IxG_cub = mass_cuboid*(b^2 + c^2)/12;
IyG_cub = mass_cuboid*(c^2 + a^2)/12;
IzG_cub = mass_cuboid*(b^2 + a^2)/12;

IG_cuboid =  [IxG_cub  0      0;
               0   IyG_cub    0;
               0       0   IzG_cub]

% Test 13: IO_semicyl (Pretest)
% Find the inertia tensor of the semicylinder part about the fixed point O in frame {1}
% from table: https://www.efunda.com/math/solids/solids_display.cfm?SolidName=HalfCircularCylinder
%IxO_cyl = (mass_semicyl*R^2)/2;
%IyO_cyl = mass_semicyl*(R^2/4 + a^2/3);
%IzO_cyl = mass_semicyl*(R^2/4 + a^2/3);
%IO_semicyl= [IxO_cyl  0      0;
%              0    IyO_cyl    0;
%              0       0   IzO_cyl]
dx = CoM_semicyl(1)
dy = CoM_semicyl(2)
dz = CoM_semicyl(3)

I_prime_cyl = mass_semicyl*[dy^2+dz^2 -dx*dy  -dx*dz;
                            -dx*dy  dx^2+dz^2 -dy*dz;
                            -dx*dz    -dy*dz  dx^2+dy^2];

IO_semicyl= IG_semicyl + I_prime_cyl

% Test 14: IO_cuboid
% Find the inertia tensor of the cuboid part about the fixed point O in frame {1}
dx = CoM_cuboid(1)
dy = CoM_cuboid(2)
dz = CoM_cuboid(3)

I_prime_cuboid = mass_cuboid*[dy^2+dz^2 -dx*dy  -dx*dz;
                         -dx*dy  dx^2+dz^2 -dy*dz;
                         -dx*dz    -dy*dz  dx^2+dy^2];
IO_cuboid = IG_cuboid + I_prime_cuboid
% Test 15: IO_total
% Find the inertia tensor of the entire ERM about O in frame {1}

IO_total = IO_semicyl + IO_cuboid

% Test 16: ang_mom_dot
% Find the time-derivative of the angular momentum of the ERM about O in frame {1}
ang_mom = IO_total*ang_vel
ang_mom_dot = diff(ang_mom) + cross(ang_vel,ang_mom)

% Test 17: vibration_force_0
% Find the force that the mass exerts on the body of the motor at point O IN THE FIXED FRAME {0}

%r_OG = (CoM_semicyl*mass_semicyl + CoM_cuboid*mass_cuboid)/(mass_semicyl + mass_cuboid)  % same as CoM_total (except Matlab fucking sucks and won't cancel out common products in the numerator and denominator)

rOG_0 = R01*CoM_total
pOG_0 = ERM_mass*(diff(rOG_0) + cross(ang_vel,rOG_0))
vibration_force_0 = diff(pOG_0) + cross(ang_vel,pOG_0)

% Test 18: vibration_moment_0
% Find the moment that the mass exerts on the body of the motor about point O in frame {0}
%r_OG = (CoM_semicyl + CoM_cuboid)/(mass_semicyl + mass_cuboid)
vibration_moment_0 = cross(rOG_0,vibration_force_0)

