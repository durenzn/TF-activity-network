function [ A , TFA , B ] = GGM( T , A0 , TF , B0 , K ,lambda )
% Guassian Graphical Model to inferring TF activity additional with TF-gene
% regulatory and kinase-TF activation.
%s
% -----------------
%   Parameters
% -----------------
%     T----[n,m]  expression level of n genes in m time points
%     A0--[n,t]    connectivity matrix TF-gene, t TFs
%     TF--[t,m]   basic level of TF
%     B0--[t,k]    connectivity matrix kinase-TF, k kinases
%     K---[k,m]   concentration of kinases

% options=optimset( 'quadprog' );
% options=optimset( options , 'Display' , 'off' , 'LargeScale' , 'off' );
[ n , m ] = size( T );
[ t , k ] = size( B0 );

N = 100;
i = 1; e = 1;
rmsd = zeros( N , 1 );
A = A0;
B=B0;
A(find(A==1))=randn(size(find(A==1)));
B(find(B==1))=randn(size(find(B==1)));

TFA = zeros( t , m );

while( i < N  )
    % iteration until e not change much or steps more than given N
    i = i + 1;
    %     disp(e)
    %     TFA = [ eye( t ) ; A ] \ [ B * K + TF; T ];
    TFA1 = TF + B * K;
    TFA2 = pinv(A)*T;
    % for k=1:60
    %    TFA2(:,k) = lsqr(A,T(:,k),[],5);
    % end
    TFA_old = TFA;
    TFA = lambda*TFA1 + (1-lambda)*TFA2;
    % E step
    
    %    B = ( TFA - TF ) / K;
    
    for j = 1 : t
        % calculate activating strengh TF by TF
        index= find( B( j , : ) ~= 0 );
        K_nnz = K( index , : );
        if size( K_nnz , 1 ) == 0
            continue
        end
        TFAAA=TFA( j , : ) - TF( j , : ) ;
        B( j , index ) = TFAAA*pinv(K_nnz);
    end
    % M step: 1. B
    for j = 1 : n
        % calculate regulatory strengh gene by gene
        index= find( A( j , : ) ~= 0 );
        TFA_nnz = TFA( index , : );
        A( j , index ) = T( j , : )*pinv(TFA_nnz);
    end
    
    [A, TFA]=Normalization(A,TFA);
    
    for j = 1 : t
        % normalize B
        index= find( B( j , : ) ~= 0 );
        K_nnz = K( index , : );
        if size( K_nnz , 1 ) == 0
            continue
        end
        B( j , index ) = ( TFA( j , : ) - TF( j , : ) )*pinv(K_nnz);
    end
    
     err1 = A * TFA - T;
     err2 = TF + B * K - TFA;
     %%disp(sum(err1(:).^2)/m/n);
     %%disp(sum(err2(:).^2)/t/m);
%     rmsd_error( i ) = ((1-lambda)*sum(err1(:).^2)/m/n+lambda*sum(err2(:).^2)/t/m)^0.5;
err = TFA - TFA_old;
%     rmsd( i ) = ( sum( sum( abs(  ).^2 ) ) / t / m ) ^ 0.5;
rmsd( i ) = sum(err(:).^2)/m/t;
end