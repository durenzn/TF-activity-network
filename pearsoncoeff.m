function r=pearsoncoeff(X,Y)
if length(X) ~= length(Y)  
    error('两个数值数列的维数不相等');  
    return;  
end  
  
fenzi =  length(X)*sum(X .* Y) - (sum(X) * sum(Y)) ;  
fenmu = sqrt((length(X)*sum(X .^2) - sum(X)^2 ) * (length(X)*sum(Y .^2) - sum(Y)^2 ));
r = fenzi / fenmu;    
end %函数myPearson结束  
