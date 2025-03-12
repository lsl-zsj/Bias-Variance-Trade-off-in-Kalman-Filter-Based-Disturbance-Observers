function Bias_variance_effect()

% Brief: Demo code for Bias-Variance Trade-off in Kalman Filter-Based Disturbance Observers
% Details: bias-variance visualization of  the KF-DOB, MKCKF-DOB, and IMMKF-DOB
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
%% kfdob
D_nom=0.5;
Q=D_nom*G*G';
%Q=diag(diag(Q));
x0a=[0;0;0]; 
P0a = 10*diag([1;1;1]);
Ga=[1;G];    % for noise 
kfa.Ga=Ga;
Fa=[1 0 0;
    G,F];
Ha=[[0;0],H];
kfa.R=R;
kfa.F=Fa;
kfa.H=Ha;
kfa.G=G; % 
kfa.x0=x0a;
kfa.P0=P0a;
kfa.len=len;
kfa.n=3;
kfa.m=2;
kfa.D=1;

%% KF-DOB
mcsc=20;
for i=1:1:mcsc
vk=sqrt(R)*randn(2,len);
uk=sqrt(0.5)*randn(1,len);
[state,z,u]=measurement_generation(len,vk,uk);
statetrue{i}.true=[u;state];
% kf-dob exp(0) item=1:20
itemlen=23;
for item=1:itemlen
    kfa.Q=[D_nom*exp((item-1)/2) 0 0;
    [0;0], Q]; % nominal covairance is 0.01
    u_n=zeros(1,len);
    kffa=kf_dob_forward(kfa,u_n,z);
    akfdob{i,item}=kffa;
end

%% mkckf 
for item=1:itemlen
    kfa.Q=[D_nom*exp((item-1)/2) 0 0;
    [0;0], Q];
    u_n=zeros(1,len);
    sigmad=3+item; %%%%%%%%%%%%%%%%%%%note
    kfa.sigma_p=[sigmad 10^8 10^8]';
    kfa.sigma_r=[10^8 10^8]';

    mkckf=mkc_dob_forward(kfa,u_n,z);
    mkckfdob{i,item}=mkckf;
end
%%
% immkf
T=[0.98,0.02;0.5,0.5];
for item=1:itemlen
    kfa.Q=[D_nom 0 0;
    [0;0], Q];
    kfa.Qsise=[D_nom*exp((item-1)/2) 0 0;
    [0;0], Q];
    kffimm=kf_dob_forward_imm(kfa,u_n,z,T);
    immkfdob{i,item}=kffimm;
end

end

% statistics
for i=1:3000
    %% exp(0)
    statef=zeros(3,mcsc);
    for item=1:itemlen
        for j=1:mcsc
            statef(:,j)=akfdob{j,item}.statef(:,i);
            trues(:,j)=statetrue{j}.true(:,i);
            statemkc(:,j)=mkckfdob{j,item}.statef(:,i);
            stateimm(:,j)=immkfdob{j,item}.statef(:,i);
        end
    % 
    biasvar{item}.bias(i,:)=mean(statef-trues,2);
    biasvar{item}.std(i,:)=std(statef-trues,0,2);
    biasvarmkc{item}.bias(i,:)=mean(statemkc-trues,2);
    biasvarmkc{item}.std(i,:)=std(statemkc-trues,0,2);
    biasvarimm{item}.bias(i,:)=mean(stateimm-trues,2);
    biasvarimm{item}.std(i,:)=std(stateimm-trues,0,2);
    end
end



t1=[t,fliplr(t)];
statefl1=biasvar{1}.bias(:,1)-3*biasvar{1}.std(:,1);
statefh1=biasvar{1}.bias(:,1)+3*biasvar{1}.std(:,1);
inbew_s1=[statefl1',fliplr(statefh1')]; 

mkcstatefl1=biasvarmkc{1}.bias(:,1)-3*biasvarmkc{1}.std(:,1);
mkcstatefh1=biasvarmkc{1}.bias(:,1)+3*biasvarmkc{1}.std(:,1);
mkcinbew_s1=[mkcstatefl1',fliplr(mkcstatefh1')]; 

immstatefl1=biasvarimm{9}.bias(:,1)-3*biasvarimm{9}.std(:,1);
immstatefh1=biasvarimm{9}.bias(:,1)+3*biasvarimm{9}.std(:,1);
imminbew_s1=[immstatefl1',fliplr(immstatefh1')]; 


figure
box on
hold on
patch(t1,inbew_s1,'m','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,biasvar{1}.bias(:,1),'LineWidth',0.6,'color','black');

patch(t1,mkcinbew_s1,'g','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,biasvarmkc{1}.bias(:,1),'LineWidth',0.6,'color','red');

patch(t1,imminbew_s1,'cyan','FaceAlpha',0.4,'Edgecolor','none')
legend('off')
plot(t,biasvarimm{9}.bias(:,1),'LineWidth',0.6,'color','blue');
lg=legend('KF-DOB $3\sigma_{d,k}$','KF-DOB $b_{d,k}$',...
    'MKCKF-DOB $3\sigma_{d,k}$','MKCKF-DOB $b_{d,k}$',...
    'IMMKF-DOB $3\sigma_{d,k}$','IMMKF-DOB $b_{d,k}$','interpreter','latex','fontsize',13);
set(gca,'fontsize',16)
xlim([127,131.9]);
set(gcf,'Position',[100 100 700 600]);



% bias-variance tradeoff plot
seg=1265:1315;
seg=1270:1320;
lenseg=length(seg);
for item=1:itemlen
% kf-dob
biasvar{item}.bias_seg=mean(biasvar{item}.bias(seg,:));
biasvar{item}.std_seg=mean(biasvar{item}.std(seg,:));
bias_seg_square(item,:)=biasvar{item}.bias_seg.^2;
var_seg(item,:)=biasvar{item}.std_seg.^2;
% mkckf-dob
biasvarmkc{item}.bias_seg=mean(biasvarmkc{item}.bias(seg,:));
biasvarmkc{item}.std_seg=mean(biasvarmkc{item}.std(seg,:));
bias_seg_squaremkc(item,:)=biasvarmkc{item}.bias_seg.^2;
var_segmkc(item,:)=biasvarmkc{item}.std_seg.^2;
% immkf-dob
biasvarimm{item}.bias_seg=mean(biasvarimm{item}.bias(seg,:));
biasvarimm{item}.std_seg=mean(biasvarimm{item}.std(seg,:));
bias_seg_squareimm(item,:)=biasvarimm{item}.bias_seg.^2;
var_segimm(item,:)=biasvarimm{item}.std_seg.^2;
end

ETA=((1:itemlen)-1)/2;
figure
hold on
box on
plot(ETA,bias_seg_square(:,1),'Color','red','LineWidth',1.0)
plot(ETA,var_seg(:,1),'Color','black','LineWidth',1.0)
legend('bias$^{2}$','variance','interpreter','latex')
xlabel('$\log(\eta)$','interpreter','latex')
ylabel('value','interpreter','latex')
title('Bais-variance of $\hat{d}_k$ in KF-DOB','interpreter','latex')
set(gca,'fontsize',16)
set(gcf,'Position',[100 100 700 600]);
xlim([0,(10)])
% 
% % mkc
% ETA=((1:itemlen)-1)/2;
% figure
% hold on
% box on
% plot(ETA,bias_seg_square(:,1),'Color','red','LineWidth',1.0)
% plot(ETA,var_seg(:,1),'Color','black','LineWidth',1.0)
% plot(ETA(1:end),bias_seg_squaremkc(1:end,1),'Color','cyan','LineWidth',1.0)
% plot(ETA(1:end),var_segmkc(1:end,1),'Color','m','LineWidth',1.0)
% plot(ETA,bias_seg_squareimm(:,1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.0)
% plot(ETA,var_segimm(:,1),'Color',[0.4660 0.6740 0.1880],'LineWidth',1.0)
% legend('KF bias$^{2}$','KF variance','MKC bias$^{2}$','MKC variance','IMM bias$^{2}$','IMM variance','interpreter','latex')
% xlabel('$\log(\eta)$','interpreter','latex')
% ylabel('value','interpreter','latex')
% title('Bais-variance of $\hat{d}_k$ in KF-DOB','interpreter','latex')
% set(gca,'fontsize',16)
% set(gcf,'Position',[100 100 700 600]);


% 

% visualization
% figure 
% hold on
% box on
% %plot(ETA,bias_seg_square(:,1),'Color','red','LineWidth',1.0,'LineStyle','--')
% %plot(ETA,var_seg(:,1),'Color','black','LineWidth',1.0,'LineStyle','--')
% plot(ETA,bias_seg_square(:,1)+var_seg(:,1),'Color','red','LineWidth',2.0,'LineStyle','--')
% plot(ETA,bias_seg_squaremkc(:,1)+var_segmkc(:,1),'Color','red','LineWidth',1.0)
% plot(ETA,bias_seg_squareimm(:,1)+var_segimm(:,1),'Color',[0.4660 0.6740 0.1880],'LineWidth',1.0)
% 
% %



% %
% figure 
% hold on
% box on
% plot(ETA,bias_seg_square(:,1),'Color','red','LineWidth',1.0,'LineStyle','--')
% plot(ETA,var_seg(:,1),'Color','black','LineWidth',1.0,'LineStyle','--')
% plot(ETA,bias_seg_square(:,1)+var_seg(:,1),'Color','black','LineWidth',2.0,'LineStyle','-')
% scatter(ETA(1:end),bias_seg_squaremkc(1:end,1)+var_segmkc(1:end,1),'SizeData', 40,'Color','red','LineWidth',2.0,'Marker','+')
% scatter(ETA(1:end),bias_seg_squareimm(1:end,1)+var_segimm(1:end,1),'SizeData', 40,'Color',[0.4660 0.6740 0.1880],'LineWidth',2.0,'Marker','square')
% lg=legend('KF-DOB bias$^{2}$','KF-DOB variance','KF-DOB bias$^{2}$+variance',...
%     'MKCKF-DOB bias$^{2}$+variance','IMMKF-DOB bias$^{2}$+variance','interpreter','latex');
% set(lg,'fontsize',12)
% set(gca,'fontsize',16)
% xlabel('$\log(\eta)$','interpreter','latex')
% ylabel('value','interpreter','latex')
% ylim([0,5.5])
% set(gcf,'Position',[100 100 700 600]);


end