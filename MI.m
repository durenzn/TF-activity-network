function MI  = MI(X,Y)
L = size(X,1);
MI = zeros(size(X,2),size(Y,2));
sigmax=std(X);
sigmay=std(Y);
h1=1.06*sigmax*L^(-0.2);
h2=1.06*sigmay*L^(-0.2);
h1square = h1.^2;
h2square = h2.^2;
X=X';
Y=Y';
H1=diag(1./h1square);
H2=diag(1./h2square);
h = waitbar(0,'Computing Mutual information...');
for i=1:L
    tmpx = X - repmat(X(:,i),1,L);
    tmpy = Y - repmat(Y(:,i),1,L);
    tmpx = exp(-H1*(tmpx.^2)/2);
    tmpy = exp(-H2*(tmpy.^2)/2);
    tmp1x = sum(tmpx,2);
    tmp1y = sum(tmpy,2);
    
    
    tmp2 = tmpx*tmpy';
    tmp3=tmp1x*tmp1y';
    tmp2=tmp2./tmp3;
    MI = MI + log(tmp2);
    clear tmp2

% %   The following commented line does the same job as lines 16~22
%     MIs = MIs + log((tmp*tmp')./(tmp1*tmp1'));
    per = i /L;
    waitbar(per, h ,sprintf('%2.0f%%',per*100))
end
MI = MI/L + log(L);
return
close(h)