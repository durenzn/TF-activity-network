function ccmi=Modulator1(A,NW,id,ppi)
TF=A(:,id)';
sizeA2=size(A,2);
sizeA1=size(A,1);
mm=floor(size(A,1)*0.35);
for i=1:size(NW,1)
    cmi=[];
    aa=find(NW(i,:)~=0);
    a=A(:,aa);
    if size(a,2)>1
        ab=find(ppi(:,1)==i);
        if size(ab,1)~=0 
        for j=1:size(aa,2)
            if abs(pearsoncoeff(TF(i,:)',a(:,j)))<1-10^(-4)
            for l=1:100
                         index=randperm(sizeA1);
                         [dd ff]=sort(index);
                         mu1=pearsoncoeff(TF(i,ff(1:mm))',a(ff(1:mm),j));
                         mu2=pearsoncoeff(TF(i,ff(sizeA1-mm+1:sizeA1))',a(ff(sizeA1-mm+1:sizeA1),j));
                         cc(l)=mu1-mu2;
            end
            maxcc=max(cc);
            mincc=min(cc);
                for k=1:size(ab,1)
                    [dd ff]=sort(A(:,ppi(ab(k,1),2)));
                    mu1=pearsoncoeff(TF(i,ff(1:mm))',a(ff(1:mm),j));
                    mu2=pearsoncoeff(TF(i,ff(sizeA1-mm+1:sizeA1))',a(ff(sizeA1-mm+1:sizeA1),j));
                    cmi(j,k)=mu1-mu2;

                    if cmi(j,k)>=max(cc) || cmi(j,k)<=min(cc)
                        cmi(j,k)=i;
                        else cmi(j,k)=0;
                    end
                end
            else cmi(j,:)=zeros(1,size(ab,1));
            end
        end
        else cmi=0;
        end
    else cmi=0;
    end
    ccmi{i}=cmi;
end
