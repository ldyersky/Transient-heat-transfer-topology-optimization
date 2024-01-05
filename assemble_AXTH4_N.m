function [iK,jK,sK]=assemble_AXTH4_N(element,node,k)
numdofs=1;
numke=4*4;
iK=zeros(numke*size(element,1),1);
jK=zeros(numke*size(element,1),1);
sK=zeros(numke*size(element,1),1);

for i=1:size(element,1)
    dofs=ones(numdofs,1)*element(i,:)*numdofs-(numdofs-1:-1:0)'*ones(1,size(element,2));
    dofs=dofs(:);
    Ke=AXTH4_N(k(i),node(element(i,:),:));
    [ii,jj]=meshgrid(dofs,dofs);
    iK((i-1)*numke+1:i*numke)=jj(:);
    jK((i-1)*numke+1:i*numke)=ii(:);
    sK((i-1)*numke+1:i*numke)=Ke(:);
end