function KE=AXTH4_N(k0,coordinate)
quadorder = 2;
[W,Q] = quadrature(quadorder,'GAUSS', 2 );
KE=zeros(4,4);
for i=1:size(Q,1)
    [~,dNdxi] = lagrange_basis('Q4',Q(i,:));  % element shape functions
    J0 = coordinate'*dNdxi;                   % element Jacobian matrix
    invJ0 = inv(J0);
    dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
    B = zeros(2,4);
    B(1,1:4)  = dNdx(:,1)';
    B(2,1:4)  = dNdx(:,2)';
    KE = KE + B'*(eye(2)*k0)*B*W(i)*det(J0);
end