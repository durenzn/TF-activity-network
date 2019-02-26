function [TFtN,NW]=MutualNet1(LA,TFid,arfa)
TF=LA(TFid,:);
G=MI(LA',TF');
for i=1:size(G,2)
    G(TFid(i),i)=mean(G(:,i));
end
beta=norminv(1-arfa/2);
NW=zeros(size(TF,1),size(LA,1));
for j=1:size(TF,1)
    Z(:,j)=zscore(G(:,j));
    NW(j, find(abs(Z(:,j))>beta))=1;
    % NW(j,TFid(j))=1;
end
a=[];
for i=1:size(TF,1)
    q=find(NW(i,:)~=0);
    q=q';
    p=[i*ones(size(q,1),1),q];
    a=[a;p];
end
r=(a(:,1)-1)*size(LA,1)+a(:,2);
a(:,3)=abs(GG(r));
Q=1-normcdf(abs(Z));
a(:,4)=Q(r);
a(:,1)=TFid(a(:,1));
W=a(:,1)==a(:,2);
a(W,4)=eps;
a(find(a(:,4)==0),4)=eps;
TFtN=a;
