global h1 h2 h3 alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS
syms t bigr smallr bigt smallt holeBigr holeSmallr pullyr pullyt mPullyAndBig rho axlel dNearSmall dNearPully dSmallEnd dBigEnd;
syms g Fr3x Fr3y Fr3z Mr3x Mr3z Fr2x Fr2y Fr2z Mr2y
syms alpha(t) beta(t) gamma(t);
syms alpha_ddotS alpha_dotS alphaS beta_ddotS beta_dotS betaS gamma_ddotS gamma_dotS gammaS
R10 = [cos(alpha(t)),sin(alpha(t)),0;-sin(alpha(t)),cos(alpha(t)),0;0,0,1];
R21 = [1,0,0;0,cos(beta(t)),-sin(beta(t));0,sin(beta(t)),cos(beta(t))];
R32 = [cos(gamma(t)),0,sin(gamma(t));0,1,0;-sin(gamma(t)),0,cos(gamma(t))];
bigrA = 0.127;%m (A stands for actual value)
smallrA = 0.035;%m
holeBigrA = 0.0093/2;%m
bigtA = 0.0221;%m
smalltA = 0.0313;%m
holeSmallrA = 0.0121/2;%m
pullyrA = 0.0582/2;%m
pullytA = 0.022;%m
mPullyAndBigA = 1.747;%kg
axlelA = 0.48;%m
rhoA = 2700;%kg/m3
gA = 9.8;%m/s2

%d = displacement along y1 from axle's connection to base (point O)
%B = scaled value in image of equipment.
dSmallEndB = -51;
dBigEndB = 35;
dNearPullyB = 20;
dNearSmallB = -41;
dAccelB = 10;
dScale = axlelA / (-dSmallEndB + dBigEndB);
dBigEndA = dScale * dBigEndB;
dSmallEndA = dScale * dSmallEndB;
dNearPullyA = dScale * dNearPullyB;
dNearSmallA = dScale * dNearSmallB;
dAccelA = dScale * dAccelB;

%counterrA = 0.07/2;%m
%holeCounterrA = 0.0121/2;%m
%countertA = 0.0313;%m
%(axlelA * pi * (holeCounterrA^2) + (countertA * pi * ((counterrA^2)-(holeCounterrA^2))))*2700
%(axlelA * pi * (holeBigrA^2) + (countertA * pi * ((counterrA^2)-(holeCounterrA^2))))*2700

pullyV = (pullyr ^ 2 - holeBigr ^ 2) * pi * pullyt;
bigV = (bigr ^ 2 - holeBigr ^ 2) * pi * bigt;
mPullyNoHole = pi * (pullyr ^ 2) * pullyt * mPullyAndBig / (bigV + pullyV);
mBigNoHole = pi * (bigr ^ 2) * bigt * mPullyAndBig / (bigV + pullyV);
mBigHole = pi * (holeBigr ^ 2) * bigt * mPullyAndBig / (bigV + pullyV);
mPullyHole = pi * (holeBigr ^ 2) * pullyt * mPullyAndBig / (bigV + pullyV);
mPully = pullyV * mPullyAndBig / (bigV + pullyV);
mBig = bigV * mPullyAndBig / (bigV + pullyV);
IPullyNoHoleCom_2 = mPullyNoHole * [(3*pullyr^2+pullyt^2)/12,0,0;0,pullyr^2/2,0;0,0,(3*pullyr^2+pullyt^2)/12];
IPullyHoleCom_2 = mPullyHole * [(3*holeBigr^2+pullyt^2)/12,0,0;0,holeBigr^2/2,0;0,0,(3*holeBigr^2+pullyt^2)/12];
IBigNoHoleCom_2 = mBigNoHole * [(3*bigr^2+bigt^2)/12,0,0;0,bigr^2/2,0;0,0,(3*bigr^2+bigt^2)/12];
IBigHoleCom_2 = mBigHole * [(3*holeBigr^2+bigt^2)/12,0,0;0,holeBigr^2/2,0;0,0,(3*holeBigr^2+bigt^2)/12];
IPullyCom_2 = IPullyNoHoleCom_2 - IPullyHoleCom_2;
IBigCom_2 = IBigNoHoleCom_2 - IBigHoleCom_2;
deltaYPully = (mBig/mPullyAndBig) * ((bigt+pullyt)/2);
deltaYBig = (mPully/mPullyAndBig) * ((bigt+pullyt)/2);
%Note Co1 means center of object 1 (pully and large disk)
IPully_Co1_2 = IPullyCom_2 + mPully * [deltaYPully^2,0,0;0,0,0;0,0,deltaYPully^2];
IBig_Co1_2 = IBigCom_2 + mBig * [deltaYBig^2,0,0;0,0,0;0,0,deltaYBig^2];
I1_Co1_2 = IBig_Co1_2 + IPully_Co1_2;

smallV = (smallr ^ 2 - holeSmallr ^ 2) * pi * smallt;
axleV = (holeSmallr ^ 2) * pi * axlel;
mSmall = rho * smallV;
mAxle = rho * axleV;
mSmallNoHole = rho * (smallr ^ 2) * pi * smallt;
mSmallHole = rho * (holeSmallr ^ 2) * pi * smallt;
ISmallNoHoleCom_2 = mSmallNoHole * [(3*smallr^2+smallt^2)/12,0,0;0,smallr^2/2,0;0,0,(3*smallr^2+smallt^2)/12];
ISmallHoleCom_2 = mSmallHole * [(3*holeSmallr^2+smallt^2)/12,0,0;0,holeSmallr^2/2,0;0,0,(3*holeSmallr^2+smallt^2)/12];
IAxleCom_2 = mAxle * [(3*holeSmallr^2+axlel^2)/12,0,0;0,holeSmallr^2/2,0;0,0,(3*holeSmallr^2+axlel^2)/12];
ISmallCom_2 = ISmallNoHoleCom_2 - ISmallHoleCom_2;
dAxleCom = (dBigEnd + dSmallEnd)/2;
dSmallCom = dNearSmall - smallt/2;
%Note Co2 means center of object 2 (axle and small disk)
dCo2 = (dSmallCom * mSmall + dAxleCom * mAxle)/(mAxle + mSmall);
deltaYAxle = abs(dAxleCom - dCo2);
deltaYSmall = abs(dSmallCom - dCo2);
ISmall_Co2_2 = ISmallCom_2 + mSmall * [deltaYSmall^2,0,0;0,0,0;0,0,deltaYSmall^2];
IAxle_Co2_2 = IAxleCom_2 + mAxle * [deltaYAxle^2,0,0;0,0,0;0,0,deltaYAxle^2];
I2_Co2_2 = IAxle_Co2_2 + ISmall_Co2_2;

rOCo1_2 = [0;dNearPully + pullyt/2 + deltaYPully;0];
rOCo2_2 = [0;dCo2;0];
w1_1 = [0;0;diff(alpha(t))];
w21_2 = [diff(beta(t));0;0];
w32_2 = [0;diff(gamma(t));0];
w1_2 = R21 * w1_1;
w2_2 = w21_2 + w1_2;
w3_2 = w32_2 + w2_2;
w3_3 = R32 * w3_2;
w2_dot_2 = diff(w2_2,t) + cross(w2_2,w2_2);
w2_ddot_2 = diff(w2_dot_2,t) + cross(w2_2,w2_dot_2);
w3_dot_2 = diff(w3_2,t) + cross(w2_2,w3_2);
w3_ddot_2 = diff(w3_dot_2,t) + cross(w2_2,w3_dot_2);
rOCo1_dot_2 = diff(rOCo1_2,t) + cross(w2_2,rOCo1_2);
rOCo1_ddot_2 = diff(rOCo1_dot_2,t) + cross(w2_2,rOCo1_dot_2);
rOCo2_dot_2 = diff(rOCo2_2,t) + cross(w2_2,rOCo2_2);
rOCo2_ddot_2 = diff(rOCo2_dot_2,t) + cross(w2_2,rOCo2_dot_2);
hCo1_2 = I1_Co1_2 * w3_2;
hCo2_2 = I2_Co2_2 * w2_2;
hCo1_dot_2 = diff(hCo1_2,t) + cross(w2_2,hCo1_2);
hCo2_dot_2 = diff(hCo2_2,t) + cross(w2_2,hCo2_2);
pCo2_2 = (mAxle+mSmall)*rOCo2_dot_2;
pCo1_2 = mPullyAndBig*rOCo1_dot_2;
pCo1_dot_2 = diff(pCo1_2,t) + cross(w2_2,pCo1_2);
pCo2_dot_2 = diff(pCo2_2,t) + cross(w2_2,pCo2_2);

Fg2_1 = [0;0;-(mAxle+mSmall)*g];
Fr2_2 = [Fr2x;Fr2y;Fr2z];
Fr3_2 = [Fr3x;Fr3y;Fr3z];
Mr2_1 = [0;Mr2y;0];
Mr3_2 = [Mr3x;0;Mr3z];
Fg3_1 = [0;0;-mPullyAndBig*g];
F2net_2 = Fr2_2 + R21 * Fg2_1 - Fr3_2;
F3net_2 = Fr3_2 + R21 * Fg3_1;
M3net_2 = Mr3_2;
M2net_2 = R21 * Mr2_1 + cross(rOCo2_2,R21*Fg2_1) - cross(rOCo1_2,Fr3_2) - Mr3_2;
left = [F2net_2;F3net_2;M3net_2;M2net_2];
right = [pCo2_dot_2;pCo1_dot_2;hCo1_dot_2;hCo2_dot_2];
eqns = left==right;
eqns = subs(eqns,pullyt,pullytA);
eqns = subs(eqns,pullyr,pullyrA);
eqns = subs(eqns,bigt,bigtA);
eqns = subs(eqns,smallt,smalltA);
eqns = subs(eqns,bigr,bigrA);
eqns = subs(eqns,smallr,smallrA);
eqns = subs(eqns,holeBigr,holeBigrA);
eqns = subs(eqns,holeSmallr,holeSmallrA);
eqns = subs(eqns,rho,rhoA);
eqns = subs(eqns,axlel,axlelA);
eqns = subs(eqns,dBigEnd,dBigEndA);
eqns = subs(eqns,dSmallEnd,dSmallEndA);
eqns = subs(eqns,dNearPully,dNearPullyA);
eqns = subs(eqns,dNearSmall,dNearSmallA);
eqns = subs(eqns,mPullyAndBig,mPullyAndBigA);
eqns = subs(eqns,g,gA);
eqns = subs(eqns,diff(diff(alpha,t),t),alpha_ddotS);
eqns = subs(eqns,diff(alpha,t),alpha_dotS);
eqns = subs(eqns,alpha,alphaS);
eqns = subs(eqns,diff(diff(beta,t),t),beta_ddotS);
eqns = subs(eqns,diff(beta,t),beta_dotS);
eqns = subs(eqns,beta,betaS);
eqns = subs(eqns,diff(diff(gamma,t),t),gamma_ddotS);
eqns = subs(eqns,diff(gamma,t),gamma_dotS);
eqns = subs(eqns,gamma,gammaS);
% eqns = subs(eqns,gammaS,1);
% eqns = subs(eqns,alphaS,1);
% eqns = subs(eqns,betaS,1);
%eqns = subs(eqns,gamma_dotS,1);
%eqns = subs(eqns,alpha_dotS,1);
%eqns = subs(eqns,beta_dotS,1);
%eqns = subs(subs(subs(eqns,alpha(t),alpha2),beta(t),beta2),gamma(t),gamma2)
%simplify(eqns)
%eqns = vpa(eqns);
eqns;
S = solve(eqns,[Fr3x,Fr3y,Fr3z,Fr2x,Fr2y,Fr2z,Mr3x,Mr3z,Mr2y,alpha_ddotS,beta_ddotS,gamma_ddotS]);
h1=S.alpha_ddotS
h2=S.beta_ddotS
h3=S.gamma_ddotS
%tx = linspace(0,30,1000);
%k = [0.0 0.0 0.0 0.0 0.0 4.0];
%[tx,values] = ode45(@next,[0 30],k);
%plot(values(:,0),tx);

% function out = next(t1,in1,in2,in3,in4,in5,in6)
%    global h1 h2 h3 alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS
%    out = [in4;in5;in6;
%        subs(h1,{alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS}, {in1 in2 in3 in4 in5 in6});
%        subs(h2,{alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS}, {in1 in2 in3 in4 in5 in6});
%        subs(h3,{alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS}, {in1 in2 in3 in4 in5 in6})];
% end


function out = next(unused, in)
   global h1 h2 h3 alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS
   out = vpa([in(4);in(5);in(6);
       subs(h1,{alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS}, {in(1) in(2) in(3) in(4) in(5) in(6)});
       subs(h2,{alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS}, {in(1) in(2) in(3) in(4) in(5) in(6)});
       subs(h3,{alphaS betaS gammaS alpha_dotS beta_dotS gamma_dotS}, {in(1) in(2) in(3) in(4) in(5) in(6)})])
   class(out)
end

