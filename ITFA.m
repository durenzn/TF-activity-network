function [AA TFA BB] = ITFA(A,sss,Modulator,NW,lambda)
kid=[];
for j=1:size(Modulator,2)
	if L0(Modulator(:,j))>0
	kid=[kid j];
	end
end
T=A;
TF=A(sss,:);
K=A(kid,:);
A0=NW';
B0=Modulator(:,kid);
[ AA , TFA , BB ] = GGM( T , A0 , TF , B0 , K ,lambda );

