function [Z]=create_z(para)
X=para.X;
H=para.H;
R=para.R;
Z=H*X+mvnrnd([0 0],R,size(X,2))';
end
