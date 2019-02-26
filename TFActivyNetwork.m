function [NW ccmi Modulator AA TFA BB TFtN MTFtN MTFNet TFtNet]=TFActivyNetwork(LA,TFid,ppi,arfa,lambda)
tic
[TFtN,NW]=MutualNet(LA,TFid,arfa);
toc
disp('Transcriptional Regulatory Network was reconstructed....')
A=LA';
ccmi=Modulator1(A,NW,TFid,ppi);
Modulator=zeros(size(ccmi,2),size(LA,1));
for i=1:size(ccmi,2)
    if L0(ccmi{1,i})~=0
        S=ccmi{1,i};
        a=[];
        ab=find(ppi(:,1)==i);
        for j=1:size(S,2)
            if L0(S(:,j))>=size(S,1)/10
                a=[a j];
            end
        end
        if size(a,2)~=0
             Modulator(i,ppi(ab(a),2))=1;
        end
    end
end
sizeA2=size(LA,2);
mm=floor(sizeA2*0.35);
R=[];
for i=1:size(TFid,1)
a=find(NW(i,:)~=0);
b=find(Modulator(i,:)~=0);
c=ppi(find(ppi(:,1)==i),2);
if min(size(a,2),size(b,2))>0
TF=LA(TFid(i),:);
t=LA(a,:);
for j=1:size(b,2)
[dd ff]=sort(LA(b(1,j),:));
Z11=zscore(TF(1,ff(1:mm)));
Z12=zscore(t(:,ff(1:mm))');
Z21=zscore(TF(1,ff(sizeA2-mm+1:sizeA2)));
Z22=zscore(t(:,ff(sizeA2-mm+1:sizeA2))');
mu1=Z11*Z12/(mm-1);
mu2=Z21*Z22/(mm-1);
d=find(c==b(1,j));
Q=ccmi{1,i}(:,d(1))>0;
Q1=[b(1,j)*ones(size(a,2),1),TFid(i)*ones(size(a,2),1),a',(mu1-mu2)',Q];
R=[R;Q1];
end
end
end
MTFtN=R;
toc
disp('Modulational Regulatory Network was reconstructed....')
[AA TFA BB] = ITFA(LA,TFid,Modulator,NW,lambda);
M=Modulator;
kid=[];
for j=1:size(M,2)
	if L0(M(:,j))>0
	kid=[kid j];
	end
end

kkid=zeros(1,size(M,2));
for i=1:size(kid,2)
kkid(1,kid(i))=i;
end

Q=find(NW'~=0);
[Q1,Q2]=find(NW~=0);
TFtNet=[TFid(Q1),Q2,AA(Q)];

[Q1,Q2]=find(M~=0);
BQ2=kkid(Q2)';
Q=(BQ2-1)*size(M,1)+Q1;
MTFNet=[Q2,TFid(Q1),BB(Q)];
toc
disp('TF activity Network was reconstructed')