function r=pearsoncoeff(X,Y)
if length(X) ~= length(Y)  
    error('������ֵ���е�ά�������');  
    return;  
end  
  
fenzi =  length(X)*sum(X .* Y) - (sum(X) * sum(Y)) ;  
fenmu = sqrt((length(X)*sum(X .^2) - sum(X)^2 ) * (length(X)*sum(Y .^2) - sum(Y)^2 ));
r = fenzi / fenmu;    
end %����myPearson����  
