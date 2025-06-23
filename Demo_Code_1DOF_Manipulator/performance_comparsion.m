function performance_comparsion()
% Brief: Performance Comparison of Different Estimators
% Details:
%    None
% 
% Syntax:  
%     performance_comparsion()
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
% Created:                         14-Mar-2025 11:25:56
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%



clear all
KF=load('KF_error.mat');
MK=load('MKC_error.mat');
IMM=load('IMM_error.mat');

t=0:0.01:0.01*(1000-1);
theta=15*sin(2*pi*0.2*t);
theta=theta';

figure
box on
x1=subplot(3,1,1);
plot(t,KF.error.d(:,1),'red','LineWidth',1)
hold on
plot(t,MK.error.d(:,1),'black','LineWidth',1.2)
plot(t,IMM.error.d(:,1),'m','LineWidth',1.2,'Linestyle','--')
ylabel('$d$ error','interpreter','latex','FontSize',20)
xlabel([])
xticks([])
set(gca,'FontSize',18)
xlim([3,7])

x2=subplot(3,1,2);
box on
plot(t,KF.error.v(:,1),'red','LineWidth',1)
hold on
plot(t,MK.error.v(:,1),'black','LineWidth',1.2)
plot(t,IMM.error.v(:,1),'m','LineWidth',1.2,'Linestyle','--')
ylabel('$\dot{\theta}$ error','interpreter','latex','FontSize',20)
set(gca,'FontSize',18)
xlabel([])
xlim([3,7])

x3=subplot(3,1,3);
box on
plot(t,KF.error.p(:,1),'red','LineWidth',1)
hold on
plot(t,MK.error.p(:,1),'black','LineWidth',1.2)
plot(t,IMM.error.p(:,1),'m','LineWidth',1.2,'Linestyle','--')
ylabel('${\theta}$ error','interpreter','latex','FontSize',20)
set(gca,'FontSize',16)
lg=legend('KF-DOB','MKCKF-DOB','IMMKF-DOB','interpreter','latex','FontSize',16,'orientation','horizontal');
set(lg,'fontsize',12)
xlabel('time/s', 'interpreter','latex')
linkaxes([x1,x2],'x')
xlim([3,7])

set(gcf,'Position',[100 100 700 600]);

%
figure
hold on
plot(t,theta-KF.error.p_a(:,1),'color','red','LineWidth',1.2)
plot(t,theta-MK.error.p_a(:,1),'color','black','LineWidth',1.2)
plot(t,theta-IMM.error.p_a(:,1),'color','m','LineWidth',1)
ylabel('tracking error ($\deg$)','interpreter','latex','FontSize',20)
set(gca,'FontSize',18)
lg=legend('KF-DOB','MKCKF-DOB','IMMKF-DOB','interpreter','latex','FontSize',16,'orientation','horizontal');
set(lg,'fontsize',12)
xlabel('time/s', 'interpreter','latex')
box off
xlim([3,7])
set(gcf,'position',[100 100 800 600])
box on


end