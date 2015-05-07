function [X_e,u_ca]=filter_imm(para)

X0=para.X0;%��ʼ��ֵ
P0=para.P0;%��ʼЭ����
e=para.e;%��ʼ�۲����
Z=para.Z;%�۲�
H=para.H;%�۲����
T=para.T;%��������
R=para.R;%�۲�����
X_e(:,1)=X0+e;X_e_ca(:,1)=X_e(:,1);X_e_ct(:,1)=X_e(:,1);
P(:,:,1)=P0;P_ca(:,:,1)=P(:,:,1);P_ct(:,:,1)=P(:,:,1);
Q=diag([1 1],0);
G=[0.5*T^2 T 1 0 0 0;
   0 0 0 0.5*T^2 T 1]';
F_ca=[1 T 0.5*T^2 0 0 0;
      0 1 T 0 0 0;
      0 0 1 0 0 0;
      0 0 0 1 T 0.5*T^2;
      0 0 0 0 1 T;
      0 0 0 0 0 1];
w=0.1;%�ɱ����
F_ct=[1 sin(w*T)/w 0 0 (cos(w*T)-1)/w 0;
      0 cos(w*T) 0 0 -sin(w*T) 0;
      0 0 0 0 -w 0;
      0 (-(cos(w*T)-1)/w) 0 1 (sin(w*T)/w) 0;
      0 sin(w*T) 0 0 cos(w*T) 0;
      0 w 0 0 0 0];
d(:,1)=H*e;d_ca(:,1)=d(:,1);d_ct(:,1)=d(:,1); %Ԥ�����
S(:,:,1)=d(:,1)*d(:,1)';S_ca(:,:,1)=S(:,:,1);S_ct(:,:,1)=S(:,:,1);
K_ca(:,:,1)=zeros(6,2);K_ct(:,:,1)=zeros(6,2);
trans=[0.9 0.1;0.1 0.9];%1Ϊca��2Ϊct ״̬ת�ƾ���
u_ca(1)=0.5;u_ct(1)=0.5;%����
for i=2:size(Z,2)
    %% ���뽻��
    c_ca=trans(1,1)*u_ca(i-1)+trans(2,1)*u_ct(i-1);%����ĸ��ʣ���ʽ4-6�ķ�ĸ����
    X_e_ca_pre(:,i)=(trans(1,1)*u_ca(i-1)*X_e_ca(:,i-1)+trans(2,1)*u_ct(i-1)*X_e_ct(:,i-1))/c_ca;
    c_ct=trans(1,2)*u_ca(i-1)+trans(2,2)*u_ct(i-1);
    X_e_ct_pre(:,i)=(trans(1,2)*u_ca(i-1)*X_e_ca(:,i-1)+trans(2,2)*u_ct(i-1)*X_e_ct(:,i-1))/c_ct;    
    
    %��ʽ4-7
    P_ca_pre(:,:,i)=(trans(1,1)*u_ca(i-1)*(P_ca(:,:,i-1)+(X_e_ca(:,i-1)-X_e_ca_pre(:,i))*(X_e_ca(:,i-1)-X_e_ca_pre(:,i))')+...
    trans(2,1)*u_ct(i-1)*(P_ct(:,:,i-1)+(X_e_ct(:,i-1)-X_e_ca_pre(:,i))*(X_e_ct(:,i-1)-X_e_ca_pre(:,i))'))/c_ca;
    P_ct_pre(:,:,i)=(trans(1,2)*u_ca(i-1)*(P_ca(:,:,i-1)+(X_e_ca(:,i-1)-X_e_ct_pre(:,i))*(X_e_ca(:,i-1)-X_e_ct_pre(:,i))')+...
    trans(2,1)*u_ct(i-1)*(P_ct(:,:,i-1)+(X_e_ct(:,i-1)-X_e_ct_pre(:,i))*(X_e_ct(:,i-1)-X_e_ct_pre(:,i))'))/c_ct;    
    
    %% �˲�����
    X_e_ca(:,i)=F_ca*X_e_ca_pre(:,i);
    P_ca(:,:,i)=F_ca*P_ca_pre(:,:,i)*F_ca'+G*Q*G';
    d_ca(:,i)=Z(:,i)-H*X_e_ca(:,i);
    S_ca(:,:,i)=H*P_ca(:,:,i)*H'+R;
    K_ca(:,:,i)=P_ca(:,:,i)*H'*inv(S_ca(:,:,i));
    X_e_ca(:,i)= X_e_ca(:,i)+K_ca(:,:,i)*d_ca(:,i);
    P_ca(:,:,i)=P_ca(:,:,i)-K_ca(:,:,i)*H*P_ca(:,:,i);
    
    X_e_ct(:,i)=F_ct*X_e_ct_pre(:,i);
    P_ct(:,:,i)=F_ct*P_ct_pre(:,:,i)*F_ct'+G*Q*G';
    d_ct(:,i)=Z(:,i)-H*X_e_ct(:,i);
    S_ct(:,:,i)=H*P_ct(:,:,i)*H'+R;
    K_ct(:,:,i)=P_ct(:,:,i)*H'*inv(S_ct(:,:,i));
    X_e_ct(:,i)= X_e_ct(:,i)+K_ct(:,:,i)*d_ct(:,i);
    P_ct(:,:,i)=P_ct(:,:,i)-K_ct(:,:,i)*H*P_ct(:,:,i);
    
    %% ģ�͸��ʸ���
    % ��ʽ4-17
    p_ca=(1/((2*pi)*sqrt(det(S_ca(:,:,i)))))*exp(-d_ca(:,i)'*inv(S_ca(:,:,i))*d_ca(:,i)/2);
    p_ct=(1/((2*pi)*sqrt(det(S_ct(:,:,i)))))*exp(-d_ct(:,i)'*inv(S_ct(:,:,i))*d_ct(:,i)/2);
    c=p_ca*c_ca+p_ct*c_ct; %��ʽ4-18
    u_ca(i)=p_ca*c_ca/c;
    u_ct(i)=p_ct*c_ct/c;
    
    X_e(:,i)=u_ca(i)*X_e_ca(:,i)+u_ct(i)*X_e_ct(:,i);
    P(:,:,i)=u_ca(i)*(P_ca(:,:,i)+(X_e_ca(:,i)-X_e(:,i))*(X_e_ca(:,i)-X_e(:,i))')+...
        u_ct(i)*(P_ct(:,:,i)+(X_e_ct(:,i)-X_e(:,i))*(X_e_ct(:,i)-X_e(:,i))');
    
end

end