function score=MFMDA(Y,m_mat_new,d_mat_new,k,iterate,lamda, lamda_m,lamda_d)


[rows,cols] = size(Y);
rand('state',1);
A=abs(rand(rows,k));  
rand('state',1);
B=abs(rand(cols,k));
I=Y;
diag_m = diag(sum(m_mat_new,2));
diag_d = diag(sum(d_mat_new,2));
L_m =diag_m -m_mat_new; 
L_d =diag_d -d_mat_new;


fid = fopen( 'RunResult.txt','wt+');
for step=1:iterate
        YB = (I.*Y)*B;
        
        ABB = I.*(A*B')*B; 
        
        if lamda > 0 && lamda_m >0
            SA = m_mat_new*A;
            DA = diag_m*A;            
            YB = YB + lamda_m*SA;
            ABB = ABB + lamda*A + lamda_m*DA;
        end
        A = A.*(YB./ABB);
        YA = (I'.*Y')*A;
        BAA = ((B*A').*I')*A;
        if lamda > 0 && lamda_d >0
            SB = d_mat_new*B;
            DB = diag_d*B;
            YA = YA + lamda_d*SB;
            BAA = BAA + lamda*B + lamda_d*DB;
        end
        B = B.*(YA./BAA);
        
        dY = Y-A*B';
        obj_NMF = sum(sum(dY.^2));
        ABF = sum(sum(A.^2))+sum(sum(B.^2));
        ALA = sum(sum((A'*L_m).*A'));
        BLB = sum(sum((B'*L_d).*B'));
        obj = obj_NMF+lamda*ABF+lamda_m*ALA+lamda_d*BLB;
        error = mean(mean(abs(Y-A*B')))/mean(mean(Y));      
       
        fprintf('step=%d  obj=%d  error=%d\n',step, obj, error);   
        if error< 10^(-5)
            fprintf('step=%d\n',step);
            break;
        end            
       
end
score=A*B';
end
