% ԭʼ��
%����̫���ܼ����������б���Լ���λ��
%% ��������
clc
clear
load data1
%% ����͸��������
% n=size(tau,1);
% taub=tau;
% taud=tau;
% for i=1:n
%     if taub(i)==0
%        taub(i)=0;
%     else
%         taub(i)=(tau(i)-0.271)/0.707;
%     end
%     if taud(i)==0
%         taud(i)=0;
%     else
%         taud(i)=0.271-0.293*taub(i);
%     end
% end
%%  ����͸���Ⱦ��鹫ʽ��
%������ѡ�����㺣��47m�����Ҵ�����γ�ȵ������Ѿ�ѡ�����������ϵ��
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
%% ��������
phi=37.4;
phi=repmat(phi,8760,1);
%�Ƕ�ת��
phi=phi*pi/180;
gammas=gammas*pi/180;
omega=omega*pi/180;
%�����γ��
delta=chiweijiao(gammas,alpha,omega);
%�����Լ���¼ÿ�μ���Ľ��
i=1;
% trace=zeros(32851,3);
trace=zeros(91,2);
%% δ��ɣ�������ֱ�����
% for gamma=-180:1:180
%     for beta=0:1:90
%         gamma1=repmat(gamma,8760,1);
%         beta1=repmat(beta,8760,1);
%         gamma1=gamma1*pi/180;
%         beta1=beta1*pi./180;
%         theta=rushejiao(delta,phi,beta1,gamma1,omega);
%         G=1353*(1+0.033.*cos(360.*n/365)).*cos(theta).*(taub+taud.*(1+cos(beta))./2+tau.*0.2.*(1-cos(beta))./2);
%         G=sum(G);
%         trace(i,1)=gamma;
%         trace(i,2)=beta;
%         trace(i,3)=G;
%         i=i+1;
%     end
% end
%% ��gammaΪ0 �ļ���
    for beta=0:1:90
        beta1=repmat(beta,8760,1);
        beta1=beta1*pi./180;
        thetaz=pi/2-alpha;
        Rb=(cos(phi-beta1).*cos(delta).*cos(omega)+sin(phi-beta1).*sin(delta))./(cos(phi).*cos(delta).*cos(omega)+sin(phi).*sin(delta));
        ceshi=(cos(phi-beta1).*cos(delta).*cos(omega)+sin(phi-beta1).*sin(delta));
%         for j=1:8760
%             if Rb(j)<0
%                 Rb(j)=0;
%             end
%         end 
        G=1353*(1+0.033.*cos(360.*n/365)).*cos(thetaz).*(taub.*Rb+taud.*(1+cos(beta1))./2+(taub+taud).*0.2.*(1-cos(beta1))./2);
        %������ʱǿ�ƹ�0
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

%% ѡȡ�����������ֵ������λ�Ǻ�����б��չʾ
[G,i]=max(trace(:,2));
plot(trace)