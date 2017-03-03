function hr_im  = Superresolution( lr_im, par )

hr_im     =   imresize( lr_im, par.scale, 'bicubic' );
[h  w]   =   size(hr_im);

y         =   lr_im;
lamada    =   par.lamada;
BTY       =   par.B'*y(:);
BTB       =   par.B'*par.B;
Tau       =   zeros(0);
flag      =   0;
for k   =  1:1
    
    if par.method==2
        A            =   Compute_AR_Matrix( hr_im, par );
        N            =   Compute_NLM_Matrix( hr_im, 5, par );
        ATA          =   A'*A*par.lam_AR;
        NTN          =   N'*N*par.lam_NL;
    elseif par.method==1
        A            =   Compute_AR_Matrix( hr_im, par );
        ATA          =   A'*A*par.lam_AR;
    end
    
    [PCA_D, D0, cls_idx, s_idx, seg]    =   Set_PCA_idx( hr_im, par, par.Codeword, flag );
    [Arr  Wei]    =    find_blks( hr_im, par );
    S             =    @(x) LPCA_IS_fast( x, cls_idx, PCA_D, par, s_idx, seg, D0, flag, Tau, Arr, Wei );
    f             =    hr_im;

    for  iter = 1 : par.nIter
              
        f_pre    =  f;
          
        if (mod(iter, 150) == 0) 
            if  (iter>=700)   flag = 1;  end
            [PCA_D, D0, cls_idx, s_idx, seg]   =  Set_PCA_idx( f, par, par.Codeword, flag );
            
            if ( iter>= 700 )
                Tau     =   Cal_Sparsity_Parameters( f, cls_idx, PCA_D, par, s_idx, seg, D0, Arr, Wei );
            end
            S           =   @(x) LPCA_IS_fast( x, cls_idx, PCA_D, par, s_idx, seg, D0, flag, Tau, Arr, Wei );
                        
            
            if (iter== 300 || iter==600)
                if par.method==2
                    A            =   Compute_AR_Matrix( f, par );
                    ATA          =   A'*A*par.lam_AR;            
                    N            =   Compute_NLM_Matrix( f, 5, par );                
                    NTN          =   N'*N*par.lam_NL;
                    [Arr  Wei]   =    find_blks( f, par );
                    S            =    @(x) LPCA_IS_fast( x, cls_idx, PCA_D, par, s_idx, seg, D0, flag, Tau, Arr, Wei );
                elseif par.method==1
                    A            =   Compute_AR_Matrix( f, par );
                    ATA          =   A'*A*par.lam_AR;            
                end                
            end            
        end        
        
        f        =   f_pre(:);
        for i = 1:5
            f        =   f + lamada.*(BTY - BTB*f);
        end
        if par.method==2
            f         =  f  - ATA*f_pre(:) - NTN*f_pre(:);      
        elseif par.method==1
            f         =  f  - ATA*f_pre(:);
        end
        f    =  S( reshape(f, h,w) );
        
        
        if (mod(iter, 40) == 0)
            
            dif       =  mean((f(:)-f_pre(:)).^2);
            if (dif<par.eps && iter>=800) 
                break; 
            end            
        end

    end
    hr_im   =  f;
    fprintf('\n');    
end
function  Tau1    =   Cal_Sparsity_Parameters( im, PCA_idx, PCA_D, par, s_idx, seg, A, Arr, Wei )
b        =  par.win;
s        =  par.step;
b2       =  b*b;
[h  w]   =  size(im);

N       =  h-b+1;
M       =  w-b+1;
r       =  [1:s:N];
r       =  [r r(end)+1:N];
c       =  [1:s:M];
c       =  [c c(end)+1:M];
X0      =  zeros(b*b, N*M);
X_m     =  zeros(b*b,length(r)*length(c),'single');

N       =  length(r);
M       =  length(c);
L       =  N*M;

% For the Y component
k    =  0;
for i  = 1:b
    for j  = 1:b
        k        =  k+1;        
        blk      =  im(i:end-b+i,j:end-b+j);
        X0(k,:)  =  blk(:)';
    end
end
% Compute the mean blks
idx      =   s_idx(seg(1)+1:seg(2));
set      =   1:size(X_m,2);
set(idx) =   [];
for i = 1:par.nblk
   v            =  Wei(i,set);
   X_m(:,set)   =  X_m(:,set) + X0(:, Arr(i, set)) .*v(ones(b2,1), :);
end


Cu0         =   zeros(b2, L, 'single' );
coe         =   A*X0(:, idx);
Cu0(:,idx)  =   abs(coe);

% for  k  =  1 : length(idx)
%     i          =   idx(k);
% %     coe        =   A*(X0(:, Arr(:, i)) - repmat(X_m(:, i), 1, par.nblk) );
%     coe        =   A*X0(:, Arr(:, i));
%     Cu0(:,i)   =   sqrt( mean(coe.^2, 2) );    
% end

set        =   1:L;
set(idx)   =   [];
L          =   length(set);

for  k  =  1 : L
    i         =   set(k);
    cls       =   PCA_idx(i);
    P         =   reshape(PCA_D(:, cls), b2, b2);
    
    coe       =   P*(X0(:, Arr(:, i)) - repmat(X_m(:, i), 1, par.nblk) );
    Cu0(:,i)  =   sqrt( mean(coe.^2, 2) );    
end
e        =   0.35;
Tau1     =   par.c1./(abs(Cu0) + e);    
