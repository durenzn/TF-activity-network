function k=L0(X)      
[m,n]=size(X);
[p q]=size(find(X==0));
k=m*n-p*q;