function pos_error_plot()

clear all
clc
addpath(genpath(pwd));

for i=1:6
    for j=1:3
        name1=num2str(i*10);
        if(j==1)
        name2=num2str(j);
        else
        name2=num2str((j-1)*5);
        end
        name=[name1,'_',name2,'.mat'];
        rms{i,j}=load(name);
    end
end
%rms.ij=['hz;','var'];
%rms.channel=['nodob;','mkcekf;','immkf;','cekf;','cekf;','ndob'];

%% rms plot
for i=1:6
    % hip
    yh_f(:,i)=[rms{i,1}.dobm{1}.qe_rms(1);rms{i,1}.dobm{2}.qe_rms(1);rms{i,1}.dobm{3}.qe_rms(1);rms{i,1}.dobm{4}.qe_rms(1);rms{i,1}.dobm{5}.qe_rms(1);rms{i,1}.dobm{6}.qe_rms(1)];
    
    yh_f5(:,i)=[rms{i,2}.dobm{1}.qe_rms(1);rms{i,2}.dobm{2}.qe_rms(1);rms{i,2}.dobm{3}.qe_rms(1);rms{i,2}.dobm{4}.qe_rms(1);rms{i,2}.dobm{5}.qe_rms(1);rms{i,2}.dobm{6}.qe_rms(1)];

    yh_f10(:,i)=[rms{i,3}.dobm{1}.qe_rms(1);rms{i,3}.dobm{2}.qe_rms(1);rms{i,3}.dobm{3}.qe_rms(1);rms{i,3}.dobm{4}.qe_rms(1);rms{i,3}.dobm{5}.qe_rms(1);rms{i,3}.dobm{6}.qe_rms(1)];
    
    % knee
    yk_f(:,i)=[rms{i,1}.dobm{1}.qe_rms(2);rms{i,1}.dobm{2}.qe_rms(2);rms{i,1}.dobm{3}.qe_rms(2);rms{i,1}.dobm{4}.qe_rms(2);rms{i,1}.dobm{5}.qe_rms(2);rms{i,1}.dobm{6}.qe_rms(2)];
    
    yk_f5(:,i)=[rms{i,2}.dobm{1}.qe_rms(2);rms{i,2}.dobm{2}.qe_rms(2);rms{i,2}.dobm{3}.qe_rms(2);rms{i,2}.dobm{4}.qe_rms(2);rms{i,2}.dobm{5}.qe_rms(2);rms{i,2}.dobm{6}.qe_rms(2)];

    yk_f10(:,i)=[rms{i,3}.dobm{1}.qe_rms(2);rms{i,3}.dobm{2}.qe_rms(2);rms{i,3}.dobm{3}.qe_rms(2);rms{i,3}.dobm{4}.qe_rms(2);rms{i,3}.dobm{5}.qe_rms(2);rms{i,3}.dobm{6}.qe_rms(2)];

    % dq
    yh_df(:,i)=[rms{i,1}.dobm{1}.dqe_rms(1);rms{i,1}.dobm{2}.dqe_rms(1);rms{i,1}.dobm{3}.dqe_rms(1);rms{i,1}.dobm{4}.dqe_rms(1);rms{i,1}.dobm{5}.dqe_rms(1);rms{i,1}.dobm{6}.dqe_rms(1)];
    yk_df(:,i)=[rms{i,1}.dobm{1}.dqe_rms(2);rms{i,1}.dobm{2}.dqe_rms(2);rms{i,1}.dobm{3}.dqe_rms(2);rms{i,1}.dobm{4}.dqe_rms(2);rms{i,1}.dobm{5}.dqe_rms(2);rms{i,1}.dobm{6}.dqe_rms(2)];

    % tors
    yh_tors(:,i)=[rms{i,1}.dobm{1}.tors_rms(1);rms{i,1}.dobm{2}.tors_rms(1);rms{i,1}.dobm{3}.tors_rms(1);rms{i,1}.dobm{4}.tors_rms(1);rms{i,1}.dobm{5}.tors_rms(1);rms{i,1}.dobm{6}.tors_rms(1)];
    yk_tors(:,i)=[rms{i,1}.dobm{1}.tors_rms(2);rms{i,1}.dobm{2}.tors_rms(2);rms{i,1}.dobm{3}.tors_rms(2);rms{i,1}.dobm{4}.tors_rms(2);rms{i,1}.dobm{5}.tors_rms(2);rms{i,1}.dobm{6}.tors_rms(2)];


    % snr
    yh_snr(:,i)=[rms{i,1}.dobm{1}.snr(1);rms{i,1}.dobm{2}.snr(1);rms{i,1}.dobm{3}.snr(1);rms{i,1}.dobm{4}.snr(1);rms{i,1}.dobm{5}.snr(1);rms{i,1}.dobm{6}.snr(1)];
    yk_snr(:,i)=[rms{i,1}.dobm{1}.snr(2);rms{i,1}.dobm{2}.snr(2);rms{i,1}.dobm{3}.snr(2);rms{i,1}.dobm{4}.snr(2);rms{i,1}.dobm{5}.snr(2);rms{i,1}.dobm{6}.snr(2)];
    
end


%% rmse display
% Format the matrix elements
accuracy=5;
formattedMatrix = sprintf(['%0.' num2str(accuracy) 'f\t'], yh_f);

% Display the formatted matrix
disp(formattedMatrix);


yh_f'*1000;

yk_f'*1000;

%%  plot

fgait=0.1:0.1:0.6;
% q error
figure
x1=subplot(2,1,1);
hold on
box on
%plot(yh_f(1,:))
plot(fgait,yh_f(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yh_f(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yh_f(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
xticks([])
ylabel('hip error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
legend('MKCEKF','IMMEKF','EKF','Orientation','horizontal')
x2=subplot(2,1,2);
hold on
box on
%plot(yk_f(1,:))
plot(fgait,yk_f(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yk_f(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yk_f(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
linkaxes([x1,x2],'x')
xticks([0.1 0.2 0.3 0.4 0.5 0.6])
xlabel('frequency (hz)','Interpreter','latex')
ylabel('knee error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
set(gcf,'position',[100,100,800,700])


% snr
figure
x1=subplot(2,1,1);
hold on
box on
%plot(yh_f(1,:))
plot(fgait,yh_snr(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_snr(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yh_snr(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_snr(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yh_snr(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
xticks([])
ylabel('SNR (db)','Interpreter','latex')
set(gca,'Fontsize',16)
legend('MKCEKF','IMMEKF','EKF','Orientation','horizontal')
x2=subplot(2,1,2);
hold on
box on
%plot(yk_f(1,:))
plot(fgait,yk_snr(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_snr(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yk_snr(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_snr(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yk_snr(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
linkaxes([x1,x2],'x')
xticks([0.1 0.2 0.3 0.4 0.5 0.6])
xlabel('frequency (hz)','Interpreter','latex')
ylabel('SNR (db)','Interpreter','latex')
set(gca,'Fontsize',16)
set(gcf,'position',[100,100,800,700])

%

figure
hold on
box on
%plot(yh_f(1,:))
plot(fgait,yh_f(2,:)+yk_f(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f(3,:)+yk_f(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yh_f(4,:)+yk_f(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f(5,:)+yk_f(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yh_f(6,:)+yk_f(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
xticks([0.1 0.2 0.3 0.4 0.5 0.6])
ylabel('hip error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
legend('MKCEKF','IMMEKF','EKF','Orientation','horizontal')
xlabel('frequency (hz)','Interpreter','latex')
ylabel('total angle error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
set(gcf,'position',[100,100,800,700])

% 

figure
hold on
box on
%plot(yh_f(1,:))
plot(fgait,yh_tors(2,:)+yk_tors(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_tors(3,:)+yk_tors(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_tors(4,:)+yk_tors(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_tors(5,:)+yk_tors(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_tors(6,:)+yk_tors(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
xticks([0.1 0.2 0.3 0.4 0.5 0.6])
ylabel('hip error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
legend('MKCEKF','IMMEKF','CEKF','EKF','NDOB','Orientation','horizontal')
xlabel('frequency (hz)','Interpreter','latex')
ylabel('total angle error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
set(gcf,'position',[100,100,800,700])


%% var=5
figure
x1=subplot(2,1,1);
hold on
box on
%plot(yh_f(1,:))
plot(fgait,yh_f5(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f5(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f5(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f5(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f5(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
xticks([])
ylabel('hip error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
legend('MKCEKF','IMMEKF','CEKF','EKF','NDOB','Orientation','horizontal')
x2=subplot(2,1,2);
hold on
box on
%plot(yk_f(1,:))
plot(fgait,yk_f5(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f5(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f5(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f5(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f5(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
linkaxes([x1,x2],'x')
xticks([0.1 0.2 0.3 0.4 0.5 0.6])
xlabel('frequency (hz)','Interpreter','latex')
ylabel('knee error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
set(gcf,'position',[100,100,800,700])


%% var=10
figure
x1=subplot(2,1,1);
hold on
box on
%plot(yh_f(1,:))
plot(fgait,yh_f10(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f10(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f10(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yh_f10(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yh_f10(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
xticks([])
ylabel('hip error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
legend('MKCEKF','IMMEKF','CEKF','EKF','NDOB','Orientation','horizontal')
x2=subplot(2,1,2);
hold on
box on
%plot(yk_f(1,:))
plot(fgait,yk_f10(2,:),'Marker','+','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f10(3,:),'Marker','square','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f10(4,:),'Marker','diamond','MarkerSize',10,'LineWidth',1.5)
plot(fgait,yk_f10(5,:),'Marker','pentagram','MarkerSize',10,'LineWidth',1.5)
%plot(fgait,yk_f10(6,:),'Marker','o','MarkerSize',10,'LineWidth',1.5)
linkaxes([x1,x2],'x')
xticks([0.1 0.2 0.3 0.4 0.5 0.6])
xlabel('frequency (hz)','Interpreter','latex')
ylabel('knee error (rad)','Interpreter','latex')
set(gca,'Fontsize',16)
set(gcf,'position',[100,100,800,700])

%



end
