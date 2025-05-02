function Pdot_vec = DRE_vec(P,A,B,Q,R)
    
    n = size(A,1);
    P = reshape(P,n,n);
    Pdot = -A'*P-P*A-Q+P*B*inv(R)*B'*P;
    Pdot_vec = Pdot(:);

end