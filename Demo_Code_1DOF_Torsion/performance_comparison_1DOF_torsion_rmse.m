function performance_comparison_1DOF_torsion_rmse()

% Brief: 1DOF Torsion Simulation
% Details:
%    None
% 
% Syntax:  
%     performance_comparison_1DOF_torsion_rmse()
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
% Created:                         13-Mar-2025 14:52:52
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%
clear all

addpath(genpath(pwd));  % 将当前目录加入路径

load('sys.mat')

%
dt=0.01;
F=sys.Ad;
G1=sys.Bd1;
G2=sys.Bd2;
H=sys.Cd;

% kd-dob 
Fa=[1 zeros(1,4);       % state transfer matrix
    G2,F];
Ga=[0;G1];              % input matrix
Ha=[zeros(2,1),H];      % observation matrix 

Vwd=0.01;
Vwn=0.01;
Vvn=0.5;


Q=eye(4);
Qa=diag([Vwd;Vwn*diag(Q)]);

Ra=Vvn*eye(2);              % measurement covariance
x0a=[0;0;0;0;0];        % initial guess
P0a=eye(5);             % initial covariance

tsim = 10;
len=tsim/dt;            % step number
t    = 0:dt:tsim-dt;          % time

kfa.dt=dt;
kfa.F=Fa; 
kfa.H=Ha;
kfa.Ga=Ga; 
kfa.Q=Qa;
kfa.R=Ra;
kfa.x0=x0a;
kfa.P0=P0a;
kfa.len=len;
kfa.n=5;
kfa.m=2;

%%
% ground truth and measurement generation
uk=0*ones(1,len);
[state,z,d]=measurement_generation(kfa,uk);


% 
%% 
% kf - dob - imm 
T=[0.98,0.02;0.5,0.5];
kfa.eta=0;              
Qa=diag([Vwd*exp(kfa.eta);Vwd*diag(Q)]);
kfa.Q=Qa;
kfa.eta=4;
Qsise=diag([Vwd*exp(kfa.eta);Vwd*diag(Q)]);
kfa.Qsise=Qsise;
kffimm=kf_dob_forward_imm(kfa,uk,z,T);
% mkmckf-dob
kfa.sigma_p=[2.5 10^8 10^8 10^8 10^8]';
kfa.sigma_r=[10^8 10^8]';
kfa.eta=0;              
Qa=diag([Vwd*exp(kfa.eta);Vwd*diag(Q)]);
kfa.Q=Qa;
mkckff=mkc_dob_forward(kfa,uk,z);


%% kf-dob
j=0;
for num=-3:0.5:20
j=j+1;
kfa.eta=num;              % assign different disturbance noise coefficient
Qa=diag([Vwd*exp(kfa.eta);Vwd*diag(Q)]);
kfa.Q=Qa;
kffa{j}=kf_dob_forward(kfa,uk,z);
end





%% calculate the RMSE of the state and the disturbance 
xtrue=[d;state];
mkcdob=mkckff.statef;
immdob=kffimm.statef;
mkcdobe=(mkcdob-xtrue)';
immdobe=(immdob-xtrue)';
mkcdobe_rms=rms(mkcdobe);
immdobe_rms=rms(immdobe);

for cnt=1:j
    kfdob{cnt}=kffa{cnt}.statef;
    kfdobe{cnt}=(kfdob{cnt}-xtrue)';
    kfdobe_rms(cnt,:)=rms(kfdobe{cnt});
end

% extract the disturbance error

logeta=((1:j)-4)/2; % note that we start from eta=exp(-3)
figure
box on
hold on
plot(logeta,kfdobe_rms(:,1),'color','blue','LineWidth',2)
plot(logeta,mkcdobe_rms(1)*ones(j,1),'color','black','LineWidth',1)
plot(logeta,immdobe_rms(1)*ones(j,1),'color','red','LineWidth',1)
xlabel('$\log(\eta)$','interpreter','latex')
ylabel('RMSE of $d$','interpreter','latex')
xlim([0,10])
lg=legend('KF-DOB ($\eta$)','MKCKF-DOB','IMMKF-DOB','interpreter','latex','NumColumns',3);
set(lg,'fontsize',14)
set(gca,'fontsize',16)
set(gcf,'Position',[100 100 700 600]);


end