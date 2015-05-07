clear all;
close all;
%% IMMת�Ʒ�����һ���ģ������ν�������һ���ǿ���ģ��ת�Ƹ��ʽ����˲�Ԥ�����ڶ������������
%--------------------------Ԥ�����-------------------------
T=1;para.T=T;
X0=[0 10 0 0 0 0]';para.X0=X0;
e=[5 1 0.5 5 1 0.5]';para.e=e;
P0=diag([5 1 0.5 5 1 0.5],0);para.P0=P0;
H=[1 0 0 0 0 0; 0 0 0 1 0 0];para.H=H;
R=diag([1 1],0);para.R=R;


%--------------------------ģ����Ϣ-------------------------
%CV
F_cv=@(X)([1 T 0 0 0 0; 
      0 1 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 1 T 0;
      0 0 0 0 1 0;
      0 0 0 0 0 0]*X);
G_cv=[0.5*T^2 T 1 0 0 0;
      0 0 0 0.5*T^2 T 1]';
Q_cv=diag([0.005 0.005],0);%Q��ϵͳ����

%CA
F_ca=@(X)([1 T 0.5*T^2 0 0 0;
      0 1 T 0 0 0;
      0 0 1 0 0 0;
      0 0 0 1 T 0.5*T^2;
      0 0 0 0 1 T;
      0 0 0 0 0 1]*X);
G_ca=[0.5*T^2 T 1 0 0 0;
      0 0 0 0.5*T^2 T 1]';
Q_ca=diag([0.005 0.005],0);

%CT

F_ct=@(w)(@(X)([X(1)+(sin(w*T)/w)*X(2)+((cos(w*T)-1)/w)*X(5);
      cos(w*T)*X(2)-sin(w*T)*X(5);
      -w*X(5);
      (-(cos(w*T)-1)/w)*X(2)+X(4)+(sin(w*T)/w)*X(5);
      sin(w*T)*X(2)+cos(w*T)*X(5);
      w*X(2)]));
 G_ct=[0.5*T^2 T 1 0 0 0;
      0 0 0 0.5*T^2 T 1]';
 Q_ct=diag([0.005 0.005],0);
  
%---------------------------�켣��Ϣ------------------------
i=1;
para.period{i}.F=F_cv;
para.period{i}.F0=[0 10 5 0 10 5]';
para.period{i}.time=80;
para.period{i}.G=G_cv;
para.period{i}.Q=Q_cv;
i=i+1;
para.period{i}.F=F_ct(0.1);
para.period{i}.F0=[0 0 0 0 0 0]';
para.period{i}.time=80;
para.period{i}.G=G_ct;
para.period{i}.Q=Q_ct;
% i=i+1;
% para.period{i}.F=F_ct(-0.1);
% para.period{i}.F0=[0 0 0 0 0 0]';
% para.period{i}.time=30;
% para.period{i}.G=G_ct;
% para.period{i}.Q=Q_ct;
% i=i+1;
% para.period{i}.F=F_ca;
% para.period{i}.F0=[0 0 2 0 0 2]';
% para.period{i}.time=30;
% para.period{i}.G=G_ca;
% para.period{i}.Q=Q_ca;



%% ---------------------����X-----------------------
[X]=create_x(para);
para.X=X;


%% ---------------------����Z-----------------------
[Z]=create_z(para);
para.Z=Z;


%% ---------------------�˲�------------------------
[X_e]=filter_ca(para);
para.X_e=X_e;
tic;
[X_e_imm,u_ca]=filter_imm(para);
para.X_e_imm=X_e_imm;
toc;
%% ---------------------��ͼ------------------------
figure(1);
subplot(3,2,1);
title('x����λ��');
plot(X(1,:));
hold on;plot(Z(1,:),'g');
hold on;plot(X_e(1,:),'r');
hold on;plot(X_e_imm(1,:),'c');

subplot(3,2,2);
plot(X(4,:));
title('y����λ��');
hold on;plot(Z(2,:),'g');
hold on;plot(X_e(4,:),'r');
hold on;plot(X_e_imm(4,:),'c');

subplot(3,2,3);
plot(X(2,:));
title('x�����ٶ�');
hold on;plot(X_e(2,:),'r');
hold on;plot(X_e_imm(2,:),'c');

subplot(3,2,4);
title('y�����ٶ�');
plot(X(5,:));
hold on;plot(X_e(5,:),'r');
hold on;plot(X_e_imm(5,:),'c');

subplot(3,2,5);
plot(X(3,:));
hold on;plot(X_e(3,:),'r');
hold on;plot(X_e_imm(3,:),'c');

subplot(3,2,6);
plot(X(6,:));
hold on;plot(X_e(6,:),'r');
hold on;plot(X_e_imm(6,:),'c');

figure(2);
title('�켣')
plot(X(1,:),X(4,:));
hold on;plot(Z(1,:),Z(2,:),'g');
hold on;plot(X_e(1,:),X_e(4,:),'r');
hold on;plot(X_e_imm(1,:),X_e_imm(4,:),'c');
axis equal;

figure(3);
plot(u_ca);


