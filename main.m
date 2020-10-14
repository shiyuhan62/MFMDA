
          Vp=find(B1()==1);
        Vn=find(B1()==0);
         MatPredict=zeros(383,495);
       Ip=crossvalind('Kfold',numel(Vp),5);
        In=crossvalind('Kfold',numel(Vn),5);
         for I=1:5
           vp=Ip==I;
           vn=In==I;
              matDT=B1;
             matDT(Vp(vp))=0; 
      
     recMatrix=MFMDA(matDT,kd,km,100,2000,0.1,1,1);
            V=[Vn(vn);Vp(vp)];
             MatPredict(V)=recMatrix(V);
           end
     [AUC,AUPR,Acc,Sen,Spe,Pre]=ROCcompute(MatPredict(),B1(),1);  

        
       
       
  