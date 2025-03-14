function system_dynamic_construction()
% Brief: System dynamics for 1 DOF Torsion System
% Details: See Section IV.B of Fusion Kalman/UFIR Filter for State Estimation
% With Uncertain Parameters and Noise Statistics
%    None
% 
% Syntax:  
%     system_dynamic_construction()
% 
% Inputs:
%    None
% 
% Outputs:
%    None
% 
% Example: 
%    None
% 
% See also: None

% Author:                          ShileiLi
% Email:                           slidk@connect.ust.hk
% Created:                         13-Mar-2025 12:09:43
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%




bm = 0.015;
bt = 0.004;
Jm = 5.4538e-4;
Jt = 2.1841e-4;
Ks = 1.0;

A = [...
  0         0         1         0; ...
  0         0         0         1; ...
 -Ks/Jm   +Ks/Jm   -bm/Jm       0; ...
 +Ks/Jt   -Ks/Jt       0     -bt/Jt ...
];

B1 = [0; 0; 1/Jm; 0];

B2 = [0; 0; 0; 1/Jt];


% (If you need an output matrix, for example measuring both angles:)
C = [1 0 0 0;
     0 1 0 0];

D = [0; 0];

%

dt = 0.01;
sysc = ss(A,B1,C,D);
sysd1 = c2d(sysc, dt, 'zoh');  % zero-order hold

sysc = ss(A,B2,C,D);
sysd2 = c2d(sysc, dt, 'zoh');  % zero-order hold

Ad = sysd1.A;
Bd1 = sysd1.B;
Cd = sysd1.C;
Dd = sysd1.D;

Bd2 = sysd2.B;


sys.Ad=Ad;
sys.Bd1=Bd1;
sys.Bd2=Bd2;
sys.Cd=Cd;
sys.Dd=Dd;

save sys.mat sys



end
