function [state,zz,uu]=measurement_generation(len,vk,uk)

amp=30;
width=40;
start1=1200;
start2=1300;
dt=0.1;
R=0.01*[10 0;0,2];
G=[0.5 * dt ^2;
     dt];
F=[1 dt;
     0 1];
H=[1 0; 0 1];
x    = zeros(2, len);
u   = zeros(1, len);          % control input
z    = zeros(2, len);         % measurement at hz
x0 = [0;0];
u0=0;
for i = 1:len
    % process disturbance generation
    if(i>=start1&&i<=(start1+2*width))
        if(mod(floor((i-start1)/(width)),2)==0)
            sign=1;
        else
            sign=-1;
        end
        if(mod(i-start1,width)<=width/2)
           u(:,i) = sign*amp;   
        else
           u(:,i) = -sign*amp; 
        end
    elseif(i>=start2&&i<=(start2+2*width))
        if(mod(floor((i-start2)/(width)),2)==0)
            sign=1;
        else
            sign=-1;
        end
        if(mod(i-start2,width)<=width/2)
           u(:,i) = -sign*amp;   
        else
           u(:,i) = sign*amp; 
        end
    else
    u(:,i) = 0;    
    end
    %u_n(:,i) = randn(1); % this actually is the sequence of 0 to N-1
    % real state generation
    if(i==1)
    x(:,i) =F*x0+G*u0;
    else
    %x(:,i) = F * x(:,i-1) + G*u(:,i)+  G* u_n(:,i);   % Generate truth  u(:,i) actually is u(:,i-1)
    x(:,i) = F * x(:,i-1) + G*(u(:,i)+uk(:,i));   % Generate truth  u(:,i) actually is u(:,i-1)
    end
    % measurement
    z(:,i) = H*x(:,i) +   vk(:,i);      
end
state=x;
zz=z;
uu=u+uk;
end