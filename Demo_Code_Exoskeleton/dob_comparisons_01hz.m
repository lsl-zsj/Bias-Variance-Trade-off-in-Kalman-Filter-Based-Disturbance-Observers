function dob_comparisons_01hz()

% Brief: One line description of what the function or class performs
% Details:
%    None
% 
% Syntax:  
%     dob_comparisons_01hz()
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
% Created:                         14-Mar-2025 17:11:26
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%

clear all
clc
addpath(genpath(pwd));

% f=0.1 cc=8 vardtimes=1
dob{1}=load('10hznodob');     % 1 nodob
dob{2}=load('10hzmkekf_1');   % 2 mkcekf
dob{3}=load('10hzimm_1');     % 3 immkf
dob{4}=load('10hzekf_1');     % 4 ekf

%% smoother bandwidth
smc = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',0.2,'DesignMethod','butter');

%% performance comparison 3002:end
index=3002:33002;
dt=0.001;
t=0:dt:dt*(length(index)-1);
dobd=[];
dobm=[];
for i=1:4
    dobd{i}.q=dob{i}.rlt.q(index,:); % q
    dobd{i}.dq=dob{i}.rlt.dq(index,:); % dq
    dobd{i}.ddq=dob{i}.rlt.ddq(index,:); % ddq
    dobd{i}.qe=dob{i}.rlt.qe(index,:); % qe
    dobd{i}.dqe=dob{i}.rlt.dqe(index,:); % dqe
    dobd{i}.ndob=dob{i}.rlt.dist(index,3:4); % ndoe
    dobd{i}.mkcekf=dob{i}.rlt.mkekfxk(index,:); % mkcekf
    dobd{i}.immekf=dob{i}.rlt.immxk(index,:); % immkf
    dobd{i}.cekf=dob{i}.rlt.cxk(index,:); % cekf
    dobd{i}.ekf=dob{i}.rlt.ekfxk(index,:); % ekf
    dobd{i}.mu=dob{i}.rlt.mu(index,:); % mu
    dobd{i}.tor=dob{i}.rlt.tor(index,:); % tor
    %% 
    torhips = filtfilt(smc,dobd{i}.tor(:,1));
    torknees = filtfilt(smc,dobd{i}.tor(:,2));
    dobd{i}.tors=[torhips,torknees];
    %
    dobd{i}.cmd=dob{i}.rlt.cmd(index,:); % tor
    torhips = filtfilt(smc,dobd{i}.cmd(:,1));
    torknees = filtfilt(smc,dobd{i}.cmd(:,2));
    dobd{i}.cmds=[torhips,torknees];
end



%% hip and knee native disturbance
dq=dobd{3}.dq;
ddq(:,1)=gradient(dq(:,1))*1000;
ddq(:,2)=gradient(dq(:,2))*1000;
naivedob=native_disturbance_observer(dobd{3}.q,dobd{3}.dq,ddq,dobd{3}.tor);
%% angle
% figure
% x1=subplot(2,1,1);
% box on
% % plot(t,10*sign(dobd{3}.dq(:,1)));
% plot(t,dobd{3}.q(:,1),'LineWidth',1,'Color','red');
% legend('hip','intepreter','latex','Orientation','horizontal')
% set(gca,'fontsize',16)
% ylabel('rad','Interpreter','latex')
% xticks([])
% x2=subplot(2,1,2);
% box on
% hold on
% % plot(t,10*sign(dobd{3}.dq(:,2)));
% plot(t,dobd{3}.q(:,2),'LineWidth',1,'Color','red');
% legend('knee','intepreter','latex','Orientation','horizontal')
% set(gca,'fontsize',16)
% linkaxes([x1,x2],'x')
% xlabel('time (s)','Interpreter','latex')
% ylabel('rad','Interpreter','latex')
% set(gcf,'Position',[100 100 700 600]);
% xlim([5,25])
% 
% figure
% x1=subplot(2,1,1);
% box on
% hold on
% %plot(t,dobd{1}.qe(:,1),'LineWidth',1,'Color','black');
% plot(t,dobd{2}.qe(:,1),'LineWidth',1,'Color','black');
% plot(t,dobd{3}.qe(:,1),'LineWidth',1,'Color','red');
% plot(t,dobd{4}.qe(:,1),'LineWidth',1,'Color','blue');
% set(gca,'fontsize',16)
% legend('MKCEKF','IMMEKF','EKF','intepreter','latex','Orientation','horizontal')
% xticks([])
% x2=subplot(2,1,2);
% box on
% hold on
% %plot(t,dobd{1}.qe(:,2),'LineWidth',1,'Color','black');
% plot(t,dobd{2}.qe(:,2),'LineWidth',1.5,'Color','black');
% plot(t,dobd{3}.qe(:,2),'LineWidth',1.5,'Color','red');
% plot(t,dobd{4}.qe(:,2),'LineWidth',1,'Color','blue');
% set(gca,'fontsize',16)
% linkaxes([x1,x2],'x')
% xlabel('time (s)','Interpreter','latex')
% set(gcf,'Position',[100 100 700 600]);
% xlim([5,25])


% zoom-in 
figure
x1=subplot(2,1,1);
box on
hold on
% plot(t,10*sign(dobd{3}.dq(:,1)));
plot(t(1:2:end),naivedob(1:2:end,1),'LineWidth',0.4,'Color',[0 0.4470 0.7410],'LineStyle','--');
plot(t,dobd{3}.mkcekf(:,1),'LineWidth',1.5,'Color','black');
plot(t,dobd{3}.immekf(:,1),'LineWidth',1.5,'Color','red');
plot(t,dobd{3}.ekf(:,1),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
set(gca,'fontsize',16)
lg=legend('RDOB','MKCEKF','IMMEKF','EKF','Interpreter','latex','NumColumns',4);
ylabel('Nm','Interpreter','latex')
xticks([])
x2=subplot(2,1,2);
box on
hold on
% plot(t,10*sign(dobd{3}.dq(:,2)));
plot(t(1:2:end),naivedob(1:2:end,2),'LineWidth',0.4,'Color',[0 0.4470 0.7410],'LineStyle','--');
plot(t,dobd{3}.mkcekf(:,2),'LineWidth',1.5,'Color','black');
plot(t,dobd{3}.immekf(:,2),'LineWidth',1.5,'Color','red');
plot(t,dobd{3}.ekf(:,2),'LineWidth',1,'Color',[0.4940 0.1840 0.5560]);
set(gca,'fontsize',16)
linkaxes([x1,x2],'x')
xlabel('time (s)','Interpreter','latex')
ylabel('Nm','Interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
xlim([7.9,11])

%% position error plot
% 1 nodob
% 2 mkcekf
% 3 immkf
% 4 cekf
% 5 ndob
figure
x1=subplot(2,1,1);
box on
hold on
%plot(t,dobd{1}.qe(:,1),'LineWidth',1,'Color','black');
plot(t,dobd{2}.qe(:,1),'LineWidth',1.5,'Color','black');
plot(t,dobd{3}.qe(:,1),'LineWidth',1.5,'Color','red');
plot(t,dobd{4}.qe(:,1),'LineWidth',1,'Color','blue');
set(gca,'fontsize',16)
legend('MKCEKF','IMMEKF','EKF','Interpreter','latex','NumColumns',5)
xticks([])
ylabel('rad','Interpreter','latex')
x2=subplot(2,1,2);
box on
hold on
%plot(t,dobd{1}.qe(:,2),'LineWidth',1,'Color','black');
plot(t,dobd{2}.qe(:,2),'LineWidth',1.5,'Color','black');
plot(t,dobd{3}.qe(:,2),'LineWidth',1.5,'Color','red');
plot(t,dobd{4}.qe(:,2),'LineWidth',1,'Color','blue');
set(gca,'fontsize',16)
linkaxes([x1,x2],'x')
xlabel('time (s)','Interpreter','latex')
ylabel('rad','Interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
xlim([7.9,11])





end