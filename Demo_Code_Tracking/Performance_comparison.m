function Performance_comparison()

% Brief: Demo code for Bias-Variance Trade-off in Kalman Filter-Based Disturbance Observers
% Details: The disturbance estimate of KF-DOB, MKCKF-DOB, and IMMKF-DOB.
%    None
% 
% Syntax:  
%     Performance_comparison()
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
% Created:                         12-Mar-2025 11:54:48
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%


% This function is to demonstrate the bias-variance effects of 
% different algorithms
clear all
addpath(genpath(pwd));  % 将当前目录加入路径

% linear tracking example with unknown input
dt=0.1;
R=0.01*[10 0;0,2];
G=[0.5 * dt ^2;
     dt];
F=[1 dt;
     0 1];
H=[1 0; 0 1];
tsim = 300;
len=tsim/dt;
t    = 0:dt:tsim-dt;          % time

% generate measurement and true state
vk=sqrt(R)*randn(2,len);
uk=sqrt(0.5)*randn(1,len);
[state,z,u]=measurement_generation(len,vk,uk);
statetrue.true=[u;state];

%% sise
D_nom=0.5;
Q=D_nom*G*G';
%Q=diag(diag(Q));
x0=[0;0];
P0 = diag([1;1]);
kf.Q=Q;
kf.R=R;
kf.F=F;
kf.H=H;
kf.G=G;
kf.x0=x0;
kf.P0=P0;
kf.len=len;
kf.n=2;
kf.m=2;
u_n=0*randn(1, len); 
ksise=kf_sise_forward(kf,u_n,z);


%% kf-dob
x0a=[0;0;0]; % d + x 
P0a = diag([1;1;1]);
Fa=[1 0 0;
    G,F];
Ha=[[0;0],H];
kfa.D=1;
Qa=[D_nom*exp(kfa.D-1) 0 0;
    [0;0], Q]; % nominal covairance is 0.01
Ga=[1;G];    % for noise 
kfa.Q=Qa;
kfa.R=R;
kfa.F=Fa;
kfa.H=Ha;
kfa.G=G; % 
kfa.Ga=Ga; 
kfa.x0=x0a;
kfa.P0=P0a;
kfa.len=len;
kfa.n=3;
kfa.m=2;

% kf-dob
kffa=kf_dob_forward(kfa,u_n,z);
% kf - dob - imm 
T=[0.98,0.02;0.5,0.5];
kfa.DD=5;
Qsise=[D_nom*exp(kfa.DD) 0 0;
    [0;0], Q];
kfa.Qsise=Qsise;
kffimm=kf_dob_forward_imm(kfa,u_n,z,T);
% mkmckf-dob
kfa.sigma_p=[3 10^8 10^8]';
kfa.sigma_r=[10^8 10^8]';
Qa=[D_nom 0 0;  % 
    [0;0], Q];
kfa.Q=Qa;
mkckff=mkc_dob_forward(kfa,u_n,z);

% figure
% subplot(2,1,1)
% hold on
% box on
% plot(t,u,'color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'Marker','diamond','MarkerIndices',1:20:length(t))
% %plot(t,dkf,'color',[0.8500 0.3250 0.0980],'LineWidth',0.8,'LineStyle','-')
% hold on
% plot(t,ksise.dist,'color','blue','LineWidth',1,'LineStyle','-')
% plot(t,kffa.statef(1,:),'color','black','LineWidth',1,'LineStyle','-')
% plot(t,mkckff.statef(1,:),'color','red','LineWidth',1.5,'LineStyle','-','Marker','pentagram','MarkerIndices',1:15:length(t))
% legend('True','SISE','KF-DOB','MKCKF-DOB','interpreter','latex','Orientation','horizontal')
% set(gca,'fontsize',16)
% xlim([110,150])
% set(gcf,'Position',[100 100 700 600]);
% xticks([])
% ylabel('dist','interpreter','latex')
% subplot(2,1,2)
% hold on
% box on
% %plot(t,dkf,'color',[0.8500 0.3250 0.0980],'LineWidth',0.8,'LineStyle','-')
% hold on
% plot(t,ksise.dist-u,'color','blue','LineWidth',1,'LineStyle','-')
% plot(t,kffa.statef(1,:)-u,'color','black','LineWidth',1,'LineStyle','-')
% plot(t,mkckff.statef(1,:)-u,'color','red','LineWidth',1.5,'LineStyle','-','Marker','pentagram','MarkerIndices',1:15:length(t))
% %legend('sise','akf-dob','mkckf-dob')
% set(gca,'fontsize',16)
% xlim([117,142])
% set(gcf,'Position',[100 100 700 600]);
% xlabel('time (s)','interpreter','latex')
% ylabel('dist error','interpreter','latex')
% 
% 
% 
% 
% figure
% subplot(2,1,1)
% hold on
% box on
% plot(t,u,'color',[0.4660 0.6740 0.1880],'LineWidth',1.5,'Marker','diamond','MarkerIndices',1:20:length(t))
% plot(t,ksise.dist,'color','blue','LineWidth',1,'LineStyle','-')
% plot(t,kffa.statef(1,:),'color','black','LineWidth',1,'LineStyle','-')
% hold on
% plot(t,kffimm.statef(1,:),'color','r','LineWidth',1.5,'LineStyle','-')
% legend('True','SISE','KF-DOB','IMMKF-DOB','interpreter','latex','Orientation','horizontal')
% set(gca,'fontsize',16)
% xticks([])
% ylabel('dist','interpreter','latex')
% xlim([117,143]);
% subplot(2,1,2)
% hold on
% box on
% plot(t,ksise.dist-u,'color','blue','LineWidth',1,'LineStyle','-')
% plot(t,kffa.statef(1,:)-u,'color','black','LineWidth',1,'LineStyle','-')
% hold on
% plot(t,kffimm.statef(1,:)-u,'color','r','LineWidth',1.5,'LineStyle','-')
% %legend('SISE','AKF','IMMKF','interpreter','latex')
% set(gca,'fontsize',16)
% xlabel('time (s)','interpreter','latex')
% ylabel('dist error','interpreter','latex')
% set(gcf,'Position',[100 100 700 600]);
% xlim([117,143])


%
figure
hold on
box on
plot(t,ksise.dist-u,'color','blue','LineWidth',1.5,'LineStyle','-')
plot(t,kffa.statef(1,:)-u,'color','black','LineWidth',1.5,'LineStyle','-')
plot(t,mkckff.statef(1,:)-u,'color','g','LineWidth',1.5,'LineStyle','-','Marker','pentagram','MarkerIndices',1:15:length(t))
hold on
plot(t,kffimm.statef(1,:)-u,'color','r','LineWidth',1.0,'LineStyle','-','Marker','square','MarkerIndices',1:15:length(t))
legend('SISE','KF-DOB','MKCKF-DOB','IMMKF-DOB','interpreter','latex','NumColumns',2)
set(gca,'fontsize',16)
%xlabel('time (s)','interpreter','latex')
%ylabel('disturbance error','interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
xlim([117,143])



% figure
% hold on
% box on
% plot(t,kffimm.MU(1,:),'LineWidth',2,'color','black')
% plot(t,kffimm.MU(2,:),'LineWidth',2,'color','red')
% xlim([110,150])
% set(gca,'fontsize',16)
% xlabel('time (s)','interpreter','latex')
% ylabel('$\mu$','interpreter','latex')
% set(gcf,'Position',[100 100 700 600]);
% legend('$\mu_1\left(\eta=1\right)$','$\mu_2 \left(\eta=\exp(10)\right)$','interpreter','latex')
% ylim([-0.1,1.1])

end
