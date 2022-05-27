%同时计算太阳能集热器最佳倾斜角以及方位角
%由于计算方法,采用系数转化过程存在一定误差。
%% 加载数据
clc
clear
load data1
%%  处理透明度经验公式法
%由于所选测量点海拔47m，而且处于中纬度地区，已经选好了相关修正系数
a0=0.4237-0.00821*36;
a1=0.5055+0.00595*6.6*6.5;
k=0.2711+0.01858*2.5*2.5;
alpha=alpha*pi/180;
taub=a0+a1*exp(-k./cos(pi/2-alpha));
taud=0.271-0.293*taub;

n=size(tau,1);
for i=1:n
    if tau(i)==0
       tau(i)=0;
       taub(i)=0;
       taud(i)=0;
    elseif alpha(i)<0
            tau(i)=0;
            taub(i)=0;
            taud(i)=0;
        else
        tau(i)=1;
    end
end
%% 其他处理
phi=37.4;
phi=repmat(phi,8760,1);
%角度转化
phi=phi*pi/180;
gammas=gammas*pi/180;
omega=omega*pi/180;
%计算赤纬角
delta=chiweijiao(gammas,alpha,omega);
%跟踪以及记录每次计算的结果
i=1;
trace=zeros(91,2);
%% 计算辐射
for gamma=-180:1:180
    for beta=0:1:90
        gamma1=repmat(gamma,8760,1);
        beta1=repmat(beta,8760,1);
        gamma1=gamma1*pi/180;
        beta1=beta1*pi./180;
        theta=rushejiao(delta,phi,beta1,gamma1,omega);
        G=1353*(1+0.033.*cos(360.*n/365)).*cos(theta).*(taub+taud.*(1+cos(beta))./2+(taub+taud).*0.2.*(1-cos(beta))./2);
        %无阳光时强制归0
        for j=1:8760
            if tau(j)==0
                G(j)=0;
            end
        end
        G=sum(G);
        trace(i,1)=gamma;
        trace(i,2)=beta;
        trace(i,3)=G;
        i=i+1;
    end
end
%% 选取计算结果的最大值，将方位角和其倾斜角展示
 [G,i]=max(trace(:,3));
 %显示方位角和倾斜角的最佳组合
 trace(i,1)
 trace(i,2)
 %展示方位角为15,30,45,60时，随着倾斜角的变化，吸收太阳辐射量
 k=0:90;
 k=k';
 n1=(15+180)*91+1;
 n2=(15+180)*91+91;
 m=trace(n1:n2,3);
 plot(k,m,'r')
 hold on
 n1=(30+180)*91+1;
 n2=(30+180)*91+91;
 m1=trace(n1:n2,3);
 plot(k,m1,'b')
 n1=(45+180)*91+1;
 n2=(45+180)*91+91;
 m2=trace(n1:n2,3);
 plot(k,m2,'g')
 n1=(60+180)*91+1;
 n2=(60+180)*91+91;
 m3=trace(n1:n2,3);
 plot(k,m3,'y')
 hold off
 lengend('15','30','45','60');


 
 