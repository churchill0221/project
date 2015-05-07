function [X_e]=filter_ca(para)
X0=para.X0;
P0=para.P0;
e=para.e;
Z=para.Z;
H=para.H;
T=para.T;
R=para.R;
X_e(:,1)=X0+e;
P(:,:,1)=P0;
Q=diag([1 1],0);
G=[0.5*T^2 T 1 0 0 0;
      0 0 0 0.5*T^2 T 1]';
F=[1 T 0.5*T^2 0 0 0;
      0 1 T 0 0 0;
      0 0 1 0 0 0;
      0 0 0 1 T 0.5*T^2;
      0 0 0 0 1 T;
      0 0 0 0 0 1];
d(:,1)=H*e;
S(:,:,1)=d(:,1)*d(:,1)';
K(:,:,1)=zeros(6,2);
for i=2:size(Z,2)
    X_e(:,i)=F*X_e(:,i-1);
    P(:,:,i)=F*P(:,:,i-1)*F'+G*Q*G';
    d(:,i)=Z(:,i)-H*X_e(:,i);
    S(:,:,i)=H*P(:,:,i)*H'+R;
    K(:,:,i)=P(:,:,i)*H'*inv(S(:,:,i));
    X_e(:,i)= X_e(:,i)+K(:,:,i)*d(:,i);
    P(:,:,i)=P(:,:,i)-K(:,:,i)*H*P(:,:,i);
end

end
