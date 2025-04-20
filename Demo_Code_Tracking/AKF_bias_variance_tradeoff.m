function AKF_bias_variance_tradeoff()
% Brief: Demo code for Bias-Variance Trade-off in Kalman Filter-Based Disturbance Observers
% Details: Bias-variance tradeoff effect of the Augmented Kalman filter
%    None
% 
% Syntax:  
%     AKF_bias_variance_tradeoff()
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
% Created:                         12-Mar-2025 11:50:31
% Version history revision notes:
%                                  None
% Implementation In Matlab R2023a
% Copyright © 2025 ShileiLi.All Rights Reserved.
%

%% bias-variance tradeoff visualization
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
D_nom=0.5;

%% kfdob
Q=G*G';
%Q=diag(diag(Q));
x0a=[0;0;0]; 
P0a = 10*diag([1;1;1]);
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

%% KF-DOB
mcsc=100;
for i=1:1:mcsc
% ground truth generation
vk=sqrt(R)*randn(2,len);     % measurement
uk=sqrt(D_nom)*randn(1,len); % process noise 
[state,z,u]=measurement_generation(len,vk,uk);
statetrue{i}.true=[u;state];
%
u_n=zeros(1,3000);
% kf-dob exp(0)
kfa.D=exp(0)*D_nom;
kfa.Q=[kfa.D 0 0;
    [0;0], Q];
kffa=kf_dob_forward(kfa,u_n,z);
akfdob1{i}=kffa;
% kf-dob exp(20)
kfa.D=exp(20)*D_nom;
kfa.Q=[kfa.D 0 0;
    [0;0], Q];
kffa=kf_dob_forward(kfa,u_n,z);
akfdob2{i}=kffa;

% kf-dob exp(2)
kfa.D=exp(2)*D_nom;
kfa.Q=[kfa.D 0 0;
    [0;0], Q];
kffa=kf_dob_forward(kfa,u_n,z);
akfdob3{i}=kffa;
end

%
for i=1:3000
    %% exp(0)
    statef=zeros(3,mcsc);
    for j=1:mcsc
        statef(:,j)=akfdob1{j}.statef(:,i);
        trues(:,j)=statetrue{j}.true(:,i);
    end
    bias=mean(statef-trues,2);
    statefmean1(i,:)=bias;
    statefstd1(i,:)=std(statef-trues,0,2);
    %% exp(20)
    statef=zeros(3,mcsc);
    for j=1:mcsc
        statef(:,j)=akfdob2{j}.statef(:,i);
    end
    bias=mean(statef-trues,2);
    statefmean2(i,:)=bias;
    statefstd2(i,:)=std(statef-trues,0,2);
    %% exp(2)
    statef=zeros(3,mcsc);
    for j=1:mcsc
        statef(:,j)=akfdob3{j}.statef(:,i);
    end
    bias=mean(statef-trues,2);
    statefmean3(i,:)=bias;
    statefstd3(i,:)=std(statef-trues,0,2);
end

statefl1=statefmean1(:,1)-3*statefstd1(:,1);
statefh1=statefmean1(:,1)+3*statefstd1(:,1);

statefl2=statefmean2(:,1)-3*statefstd2(:,1);
statefh2=statefmean2(:,1)+3*statefstd2(:,1);

statefl3=statefmean3(:,1)-3*statefstd3(:,1);
statefh3=statefmean3(:,1)+3*statefstd3(:,1);

t1=[t,fliplr(t)];
inbew_s1=[statefl1',fliplr(statefh1')];
inbew_s2=[statefl2',fliplr(statefh2')];
inbew_s3=[statefl3',fliplr(statefh3')];



%%%%%%%%%

figure
box on
hold on
patch(t1,inbew_s1,'m','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,statefmean1(:,1),'LineWidth',0.6,'color','black');
%patch(t1,inbew_s2,'g','FaceAlpha',0.4,'Edgecolor','none')
%legend('off')
%plot(t,statefmean2(:,1),'LineWidth',0.6,'color','red');
patch(t1,inbew_s3,'g','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,statefmean3(:,1),'LineWidth',0.6,'color','red');
%xlim([116,125.9]);
xlim([127,131]);
ylim([-10,25]);
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('disturbance error','interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
lg=legend('$3\sigma$ $\eta=0$','${b}_{d,k}$ $\eta=0$','$3\sigma$ $\eta=2$','${b}_{d,k}$ $\eta=2$','interpreter','latex'...
    ,'fontsize',13);
%set(lg,'fontsize',16)
title('KF-DOB','interpreter','latex')
%

figure
box on
hold on
patch(t1,inbew_s1,'m','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,statefmean1(:,1),'LineWidth',0.6,'color','black');
%patch(t1,inbew_s2,'g','FaceAlpha',0.4,'Edgecolor','none')
%legend('off')
%plot(t,statefmean2(:,1),'LineWidth',0.6,'color','red');
patch(t1,inbew_s2,'c','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,statefmean2(:,1),'LineWidth',0.6,'color','blue');
%xlim([116,125.9]);
xlim([127,131]);
ylim([-10,25]);
set(gca,'fontsize',16)
xlabel('time (s)','interpreter','latex')
ylabel('disturbance error','interpreter','latex')
set(gcf,'Position',[100 100 700 600]);
legend('$3\sigma$ $\eta=0$','${b}_{d,k}$ $\eta=0$','$3\sigma$ $\eta=20$','${b}_{d,k}$ $\eta=20$','interpreter','latex'...
    ,'fontsize',13)
title('KF-DOB','interpreter','latex')



figure
box on
hold on
patch(t1,inbew_s1,'m','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,statefmean1(:,1),'LineWidth',0.6,'color','black');
%
patch(t1,inbew_s3,'cyan','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,statefmean3(:,1),'LineWidth',0.6,'color','red');
%
patch(t1,inbew_s2,'g','FaceAlpha',0.2,'Edgecolor','none')
legend('off')
plot(t,statefmean2(:,1),'LineWidth',0.6,'color','blue');

xlim([127,131.9]);
ylim([-10,30]);
lg=legend('$3\sigma_{d,k}$ $\eta=0$','${b}_{d,k}$ $\eta=0$','$3\sigma_{d,k}$ $\eta=2$','${b}_{d,k}$ $\eta=2$','$3\sigma_{d,k}$ $\eta=20$','${b}_{d,k}$ $\eta=20$','interpreter','latex'...
    ,'fontsize',13);
set(lg,'fontsize',16)
title('KF-DOB','interpreter','latex')
set(gca,'fontsize',16)
set(gcf,'Position',[100 100 700 600]);
xlabel('time (s)','interpreter','latex')
ylabel('value','interpreter','latex')


end