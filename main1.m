%计算太阳能集热器最佳倾斜角（方位角为0）
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
%% 太阳能集热器吸收的辐射能计算
    for beta=0:1:90
        beta1=repmat(beta,8760,1);
        beta1=beta1*pi./180;
        thetaz=pi/2-alpha;
        Rb=(cos(phi-beta1).*cos(delta).*cos(omega)+sin(phi-beta1).*sin(delta))./(cos(phi).*cos(delta).*cos(omega)+sin(phi).*sin(delta));
        G=1353*(1+0.033.*cos(360.*n/365)).*cos(thetaz).*(taub.*Rb+taud.*(1+cos(beta1))./2+(taub+taud).*0.2.*(1-cos(beta1))./2);
        %无阳光时系数强制归0
        for j=1:8760
            if tau(j)==0
                G(j)=0;
            end
        end
        G=sum(G);
        trace(i,1)=beta;
        trace(i,2)=G;
        i=i+1;
    end

%% 选取计算结果的最大值，将方位角和其倾斜角展示
[G,i]=max(trace(:,2));
plot(trace,'b-')
xlabel('集热器的倾斜角')
ylabel('单位面积每小时全年吸收辐射能的近似值');