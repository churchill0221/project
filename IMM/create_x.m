function [X]=create_x(para)
X0=para.X0;
period_num=size(para.period,2);
for i=1:period_num
   F0{i}=para.period{i}.F0;
   F{i}=para.period{i}.F;
   time{i}=para.period{i}.time;
   Q{i}=para.period{i}.Q;
   G{i}=para.period{i}.G;
end

count=0;
X(:,1)=X0;
j=0;
for i=1:period_num
    X(:,j+1)=X(:,j+1)+F0{i};
    for j=(count+1):(count+time{i})
        X(:,j+1)=F{i}(X(:,j))+G{i}*mvnrnd(zeros(size(Q{i},1),1),Q{i},1)';
    end
    count=j;
end

end