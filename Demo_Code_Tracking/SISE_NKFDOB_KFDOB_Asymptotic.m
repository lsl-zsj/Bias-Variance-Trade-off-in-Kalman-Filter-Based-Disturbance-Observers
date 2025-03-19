function SISE_NKFDOB_KFDOB_Asymptotic()
% Brief: Demo code for Bias-Variance Trade-off in Kalman Filter-Based Disturbance Observers
% Details: The performances of different estimators.
%    None
% 
% Syntax:  
%     bias_variance_effect_yes1()
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
%% linear tracking example with unknown input
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
%
D_nom=0.5;
% generate measurement and true state
vk=sqrt(R)*randn(2,len);     % measurement
uk=sqrt(D_nom)*randn(1,len); % process noise 
[state,z,u]=measurement_generation(len,vk,uk);
%
%% kf instance
Q=D_nom*G*G';
%Q=diag(diag(Q));
x0=[0;0];
P0 = 10*diag([1;1]);
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
kf.D=0;
%% sise
u_n=zeros(1,len);
ksise=kf_sise_forward(kf,u_n,z);
samp=21;
%% nkfdob
for i=1:1:samp
kf.D=D_nom*exp((i-1));
kff=nkf_forward(kf,u_n,z);
nkfdob{i}=kff;
end
%% akfdob
x0a=[0;0;0]; 
P0a = 10^(10)*diag([1;1;1]);
Fa=[1 0 0;
    G,F];
Ha=[[0;0],H];
Ga=[1;G];    % for noise 
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
kfa.D=1;
for i=1:1:samp
kfa.D=D_nom*exp((i-1));
kfa.Q=[kfa.D 0 0;
    [0;0], Q];
% kf-dob
kffa=kf_dob_forward(kfa,u_n,z);
akfdob{i}=kffa;
end

%% identity property
figure
hold on
box on
plot(t,u,'color','black','LineWidth',1.5)
plot(t,ksise.dist,'blue','LineWidth',0.8,'Marker','square','MarkerIndices',1:10:length(t),'MarkerSize',10)
%plot(t,pkfdob{end}.dist,'blue','LineWidth',0.6)
plot(t,nkfdob{end}.dist,'color',[0.8500 0.3250 0.0980],'LineWidth',0.6,'Marker','+','MarkerIndices',1:10:length(t),'MarkerSize',5)
plot(t,akfdob{end}.statef(1,:),'Color','red','LineWidth',0.6,'Marker','o','MarkerIndices',1:10:length(t),'MarkerSize',5)
legend('True','SISE','PKF-DOB','KF-DOB','interpreter','latex')
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('disturbance','interpreter','latex')
xlim([117,143]);
set(gcf,'Position',[100 100 700 600]);


figure
hold on
box on
plot(t,u,'color',[0.8500 0.3250 0.0980],'LineWidth',1)
plot(t,ksise.dist,'m','LineWidth',0.6)
%plot(t,pkfdob{end}.dist,'blue','LineWidth',0.6)
plot(t,nkfdob{end}.dist,'red','LineWidth',0.6)
plot(t,akfdob{end}.statef(1,:),'black','LineWidth',0.6)
legend('True','SISE','PKF-DOB','KF-DOB','interpreter','latex')
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('disturbance','interpreter','latex')
xlim([115,125]);
ylim([-60,60]);
set(gcf,'Position',[100 100 700 600]);

%% 
figure
x1=subplot(2,1,1);
box on
hold on
%plot(t,z(1,:),'red','LineWidth',1.0,'LineStyle',':')
plot(t,ksise.statef(1,:),'m','LineWidth',1.0,'LineStyle','-','Marker','o','MarkerIndices',1:30:length(t),'MarkerSize',10)
%plot(t,pkfdob{end}.statef(1,:),'blue','LineWidth',1.0,'LineStyle','-','Marker','pentagram','MarkerIndices',1:30:length(t),'MarkerSize',10)
plot(t,nkfdob{end}.statef(1,:),'red','LineWidth',1.0,'LineStyle','-','Marker','square','MarkerIndices',1:30:length(t),'MarkerSize',10)
plot(t,akfdob{end}.statef(2,:),'black','LineWidth',1.0,'LineStyle','-','Marker','+','MarkerIndices',1:30:length(t),'MarkerSize',10)
%legend('sise','pkf','nkf','akf','interpreter','latex')
set(gca,'fontsize',16)
xticks([])
ylabel('pos (m)','interpreter','latex')
x2=subplot(2,1,2);
box on
hold on
%plot(t,z(2,:),'red','LineWidth',1.0,'LineStyle',':')
plot(t,ksise.statef(2,:),'m','LineWidth',1.0,'LineStyle','-','Marker','o','MarkerIndices',1:30:length(t),'MarkerSize',10)
%plot(t,pkfdob{end}.statef(2,:),'blue','LineWidth',1.0,'LineStyle','-','Marker','pentagram','MarkerIndices',1:30:length(t),'MarkerSize',10)
plot(t,nkfdob{end}.statef(2,:),'red','LineWidth',1.0,'LineStyle','-','Marker','square','MarkerIndices',1:30:length(t),'MarkerSize',10)
plot(t,akfdob{end}.statef(3,:),'black','LineWidth',1.0,'LineStyle','-','Marker','+','MarkerIndices',1:30:length(t),'MarkerSize',10)
legend('SISE','NKF-DOB','KF-DOB','interpreter','latex','orientation','horizontal')
set(gca,'fontsize',16)
linkaxes([x1,x2],'x')
xlabel('time (s)','interpreter','latex')
ylabel('vel (m/s)','interpreter','latex')
xlim([110,150]);
set(gcf,'Position',[100 100 700 600]);

%% 
figure
x1=subplot(2,1,1);
box on
hold on
plot(t,reshape(ksise.covf(1,1,:),[len,1]),'m','LineWidth',1,'Marker','o','MarkerIndices',1:200:length(t),'MarkerSize',10)
%plot(t,reshape(pkfdob{end}.covf(1,1,:),[len,1]),'blue','LineWidth',1,'Marker','pentagram','MarkerIndices',1:200:length(t),'MarkerSize',10)
plot(t,reshape(nkfdob{end}.covf(1,1,:),[len,1]),'red','LineWidth',1,'Marker','square','MarkerIndices',1:200:length(t),'MarkerSize',10)
plot(t,reshape(akfdob{end}.covf(2,2,:),[len,1]),'black','LineWidth',1,'Marker','+','MarkerIndices',1:200:length(t),'MarkerSize',10)
xticks([])
ylabel('','interpreter','latex')
set(gca,'fontsize',16)
legend('SISE $P_{k|k}(x_1)$','NKF-DOB $P_{k|k}(x_1)$','KF-DOB $P_{k|k}(x_1)$','interpreter','latex')
x2=subplot(2,1,2);
box on
hold on
plot(t,reshape(ksise.covf(2,2,:),[len,1]),'color','m','LineWidth',1,'Marker','o','MarkerIndices',1:200:length(t),'MarkerSize',10,'LineStyle','--')
%plot(t,reshape(pkfdob{end}.covf(2,2,:),[len,1]),'color','blue','LineWidth',1,'Marker','pentagram','MarkerIndices',1:200:length(t),'MarkerSize',10,'LineStyle','--')
plot(t,reshape(nkfdob{end}.covf(2,2,:),[len,1]),'color','red','LineWidth',1,'Marker','square','MarkerIndices',1:200:length(t),'MarkerSize',10,'LineStyle','--')
plot(t,reshape(akfdob{end}.covf(3,3,:),[len,1]),'color','black','LineWidth',1,'Marker','+','MarkerIndices',1:200:length(t),'MarkerSize',10,'LineStyle','--')
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('','interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
legend('SISE $P_{k|k}(x_2)$','NKF-DOB $P_{k|k}(x_2)$','KF-DOB $P_{k|k}(x_2)$','interpreter','latex')
yrange=[0.01996,0.02000];
ylim(yrange)
linkaxes([x1,x2],'x')
set(gcf,'Position',[100 100 700 600]);
xlim([0,150])
%% trade-off property
%% pkfdob
% 
% %% almost linear decay with the growth of D
% for i=1:1:samp
% seg=1322:1359;
% pkfdob{i}.dbias=pkfdob{i}.dist(seg)-40*ones(1,length(seg));
% pkfdob{i}.dbiasmean=mean(pkfdob{i}.dbias);
% biasm(i)=abs(pkfdob{i}.dbiasmean);
% end
% figure
% plot(log(exp((1:1:samp)-1)),log(biasm),'LineWidth',2)

% combination of pkf and nkf: suprising result




%% nkfdob
figure
hold on
box on
plot(t,u,'color',[0.8500 0.3250 0.0980],'LineWidth',2)
%plot(t,ksise.dist,'m','LineWidth',0.6)
plot(t,nkfdob{1}.dist,'blue','LineWidth',0.6)
plot(t,nkfdob{2}.dist,'red','LineWidth',0.6)
plot(t,nkfdob{3}.dist,'m','LineWidth',0.6)
plot(t,nkfdob{4}.dist,'black','LineWidth',0.6)
plot(t,nkfdob{5}.dist,'color',[0 0.4470 0.7410],'LineWidth',0.6)
plot(t,nkfdob{end}.dist,'color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
legend('True','$\eta=0$','$\eta=1$','$\eta=2$','$\eta=3$','$\eta=4$','$\eta=20$','interpreter','latex','fontsize',13)
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('disturbance','interpreter','latex')
xlim([117,143]);
set(gcf,'Position',[100 100 700 600]);
title('NKF-DOB','interpreter','latex')

figure
x1=subplot(2,1,1);
box on
hold on
plot(t,nkfdob{1}.statef(1,:)-state(1,:),'blue','LineWidth',0.6)
plot(t,nkfdob{2}.statef(1,:)-state(1,:),'red','LineWidth',0.6)
plot(t,nkfdob{3}.statef(1,:)-state(1,:),'m','LineWidth',0.6)
plot(t,nkfdob{4}.statef(1,:)-state(1,:),'black','LineWidth',0.6)
plot(t,nkfdob{5}.statef(1,:)-state(1,:),'color',[0 0.4470 0.7410],'LineWidth',0.6)
plot(t,nkfdob{end}.statef(1,:)-state(1,:),'color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
xticks([])
ylabel('pos error','interpreter','latex')
set(gca,'fontsize',16)
title('NKF-DOB','interpreter','latex')
%legend('pkf $\eta=1$','pkf $\eta=2$','pkf $\eta=3$','pkf $\eta=4$','pkf $\eta=5$','pkf $\eta=21$','interpreter','latex')
x2=subplot(2,1,2);
box on
hold on
plot(t,nkfdob{1}.statef(2,:)-state(2,:),'blue','LineWidth',0.6)
plot(t,nkfdob{2}.statef(2,:)-state(2,:),'red','LineWidth',0.6)
plot(t,nkfdob{3}.statef(2,:)-state(2,:),'m','LineWidth',0.6)
plot(t,nkfdob{4}.statef(2,:)-state(2,:),'black','LineWidth',0.6)
plot(t,nkfdob{5}.statef(2,:)-state(2,:),'color',[0 0.4470 0.7410],'LineWidth',0.6)
plot(t,nkfdob{end}.statef(2,:)-state(2,:),'color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('vel error','interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
legend('$\eta=0$','$\eta=1$','$\eta=2$','$\eta=3$','$\eta=4$','$\eta=20$','interpreter','latex','NumColumns',3)
linkaxes([x1,x2],'x')
set(gcf,'Position',[100 100 700 600]);
xlim([117,143]);

%% akfdob
figure
hold on
box on
plot(t,u,'color',[0.8500 0.3250 0.0980],'LineWidth',1.5,'color','black')
%plot(t,ksise.dist,'m','LineWidth',0.6)
plot(t,akfdob{1}.statef(1,:),'blue','LineWidth',0.5)
plot(t,akfdob{2}.statef(1,:),'red','LineWidth',0.5)
plot(t,akfdob{3}.statef(1,:),'m','LineWidth',0.5)
plot(t,akfdob{4}.statef(1,:),'cyan','LineWidth',0.5)
plot(t,akfdob{5}.statef(1,:),'color',[0 0.4470 0.7410],'LineWidth',0.5)
plot(t,akfdob{end}.statef(1,:),'color',[0.8500 0.3250 0.0980],'LineWidth',0.5)
legend('True',' $\eta=0$',' $\eta=1$',' $\eta=2$',' $\eta=3$',' $\eta=4$',' $\eta=20$','interpreter','latex','fontsize',13)
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('disturbance','interpreter','latex')
xlim([117,143]);
set(gcf,'Position',[100 100 700 600]);
title('KF-DOB','interpreter','latex')

figure
x1=subplot(2,1,1);
box on
hold on
plot(t,akfdob{1}.statef(2,:)-state(1,:),'blue','LineWidth',0.6)
plot(t,akfdob{2}.statef(2,:)-state(1,:),'red','LineWidth',0.6)
plot(t,akfdob{3}.statef(2,:)-state(1,:),'m','LineWidth',0.6)
plot(t,akfdob{4}.statef(2,:)-state(1,:),'black','LineWidth',0.6)
plot(t,akfdob{5}.statef(2,:)-state(1,:),'color',[0 0.4470 0.7410],'LineWidth',0.6)
plot(t,akfdob{end}.statef(2,:)-state(1,:),'color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
xticks([])
ylabel('pos error','interpreter','latex')
set(gca,'fontsize',16)
title('KF-DOB','interpreter','latex')
%legend('pkf $\eta=1$','pkf $\eta=2$','pkf $\eta=3$','pkf $\eta=4$','pkf $\eta=5$','pkf $\eta=21$','interpreter','latex')
x2=subplot(2,1,2);
box on
hold on
plot(t,akfdob{1}.statef(3,:)-state(2,:),'blue','LineWidth',0.6)
plot(t,akfdob{2}.statef(3,:)-state(2,:),'red','LineWidth',0.6)
plot(t,akfdob{3}.statef(3,:)-state(2,:),'m','LineWidth',0.6)
plot(t,akfdob{4}.statef(3,:)-state(2,:),'black','LineWidth',0.6)
plot(t,akfdob{5}.statef(3,:)-state(2,:),'color',[0 0.4470 0.7410],'LineWidth',0.6)
plot(t,akfdob{end}.statef(3,:)-state(2,:),'color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('vel error','interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
legend(' $\eta=0$',' $\eta=1$',' $\eta=2$',' $\eta=3$',' $\eta=4$',' $\eta=20$','interpreter','latex','NumColumns',3)
linkaxes([x1,x2],'x')
set(gcf,'Position',[100 100 700 600]);
xlim([117,143]);


end