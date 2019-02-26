function [nor_a, nor_p]=normalization(a,p)
[N,L]=size(a);
B=abs(a);
B1=B~=0;
N1=sum(B1);
NF=sum(B)./N1;
NF(find(NF==0))=1;
NF1=NF.^-1;
for k=1:L
   nor_a(:,k)=a(:,k)*NF1(1,k);
   nor_p(k,:)=p(k,:)*NF(1,k);
end

