load A.mat;
load ppi.mat;
load TFid.mat;
[NW ccmi Modulator AA TFA BB TFtN MTFtN MTFNet TFtNet]=TFActivyNetwork(A,TFid,ppi,10^(-6),0.5);
disp('Exporting TF-target network....')
load name.mat;
a=[name(TFtN(:,1)),name(TFtN(:,2))];
 for i=1:size(TFtN,1)
     a{i,3}=TFtN(i,3);
     a{i,4}=TFtN(i,4);
 end
title{1,1}='TF';
title{1,2}='Target';
title{1,3}='Mutual information';
title{1,4}='P-value';
TFtN=[title;a];
filename='TF-target network.txt';
 fid=fopen(filename,'wt');
A=TFtN;
 [m,n]=size(A);
 for j=1:4
     if j==4 fprintf(fid,'%s\n',A{1,j});
     else fprintf(fid,'%s\t',A{1,j});
     end
 end
 for i=2:m
     for j=1:n
         if j==4
             fprintf(fid,'%g\n',A{i,j});
         else if j==3 
             fprintf(fid,'%g\t',A{i,j});
         else  fprintf(fid,'%s\t',A{i,j});
         end
end
     end
 end
 fclose(fid);
 disp('Exporting Modulator-TF-target network....')
 a=[name(MTFtN(:,1)),name(MTFtN(:,2)),name(MTFtN(:,3))];
 for i=1:size(MTFtN,1)
     a{i,4}=MTFtN(i,4);
     a{i,5}=MTFtN(i,5);
 end
title{1,1}='Modulator';
title{1,2}='TF';
title{1,3}='Target';
title{1,4}='Conditional mutual information';
title{1,5}='significance or not(arfa=0.01)';
MTFtN=[title;a];
A=MTFtN;
filename='Modulator-TF-target network.txt';
 fid=fopen(filename,'wt');
 [m,n]=size(A);
 for j=1:5
     if j==5 fprintf(fid,'%s\n',A{1,j});
     else fprintf(fid,'%s\t',A{1,j});
     end
 end
 for i=2:m
     for j=1:n
         if j==5
             fprintf(fid,'%g\n',A{i,j});
         else if j==4
             fprintf(fid,'%g\t',A{i,j});
         else  fprintf(fid,'%s\t',A{i,j});
         end
end
     end
 end
 fclose(fid);
 disp('Exporting TF activity network....')
 TFActivity=[name(MTFNet(:,1)),name(MTFNet(:,2));name(TFtNet(:,1)),name(TFtNet(:,2))];
 for i=1:size(TFActivity,1)
     if i<size(MTFNet,1)+1
         TFActivity{i,3}=sign(MTFNet(i,3));
     else TFActivity{i,3}=sign(TFtNet(i-size(MTFNet,1),3));
     end
 end
 A=TFActivity;
filename='TF activity network.txt';
 fid=fopen(filename,'wt');
 [m,n]=size(A);
  for i=1:m
     for j=1:n
         if j==3
             fprintf(fid,'%g\n',A{i,j});
         else  fprintf(fid,'%s\t',A{i,j});

end
     end
  end
  fclose(fid);