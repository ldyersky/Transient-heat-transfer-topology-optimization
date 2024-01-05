
% clc;
clear;

%% Input

node_info = csvread('node6400.csv');
node = node_info(:,1:2)*0.001;

element = csvread('element6400.csv');

%%

NumOfNode = length(node(:,1));

designdomain = 1:length(element(:,1));

Fixnodes = 6561;

FixDofs = Fixnodes;

AllDofs = 1:NumOfNode;

freedofs = setdiff(AllDofs,FixDofs);

%freedofs
FixTemp = 100;
%%
% node_h = [Fixnodes1 Fixnodes];


loop = 0; 
change = 1;
k = 40;
rho = 7850;
cp = 60;

penalK = 2;
penalC = 1.1;
volfrac = 0.25;
xPhys = volfrac*ones(length(element(:,1)),1);
xPhys(1) = 1;
xPhys(6400) = 1;
xold1 = zeros(6400,1);
xold2 = xold1;
upp = ones(6400,1);
low = zeros(6400,1);

Q = zeros(length(AllDofs)); 
Q(1,1) = 1;
Q = Q(freedofs,freedofs);
kxx = k*ones(length(element(:,1)),1);
[iH,jH,sH0]=assemble_AXTH4_N(element,node,kxx);
Se = (node(element(1,1),1) - node(element(1,2),1))^2 + .......
    (node(element(1,1),2) - node(element(1,2),2))^2;
[iC,jC,sC0]=assemble_AXTH4_C(element,rho,cp,Se);

tic
for jj = 1:1
    loop = loop+1;
    sH=sH0.*reshape(repmat(xPhys.^penalK,1,16)',[],1);
    sC=sC0.*reshape(repmat(xPhys.^penalC,1,16)',[],1);

    K=sparse(iH,jH,sH);
    C=sparse(iC,jC,sC);

    K = (K+K')/2;
    C = (C+C')/2;

    Kc = K(freedofs,freedofs);
    Cc = C(freedofs,freedofs);
    K12 = K(freedofs,1);
    A = -Cc\Kc;
    Q = zeros(length(AllDofs)); 
    Q(1,1) = 1;
    Q = Q(freedofs,freedofs);
    %%%%%%%%%%%%%%%%%%%%%%%
    tic
    x0 = zeros(1,6560);
    func = @(t,x) A * x - Cc\(100*K12) ;
    [t,x] = ode15s(func,0:0.001:15,x0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Kc = full(Kc);
%     Cc = full(Cc);
%     [VX,~] = eig(Kc);
    VX = pca(x);
    toc
    VX = VX(:,1:10);
    VX = VX';
    VKc = VX*Kc*VX';
    VKc = (VKc+VKc')/2;
    VCc = VX*Cc*VX';
    VCc = (VCc+VCc')/2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = -100*ones(size(Kc,1),1);
    VA = -VCc\VKc;
    VQ = VX*Q*VX';
    VKc_inv = inv(VKc);
    VKc_inv = (VKc_inv+VKc_inv')/2;   
    VCc_inv = inv(VCc);
    VCc_inv = (VCc_inv+VCc_inv')/2;    
    VK12 = VX*K12;
    VP = lyap(VA',VQ);
    T0 = VX*T;
    J = T0'*VP*T0
    numofelem = length(element(:,1));

    for i = 1:1
        dxPhys = zeros(numofelem,1);   
        dxPhys(50) = xPhys(50);
        sH =sH0.*reshape(repmat(penalK*dxPhys.^(penalK-1),1,16)',[],1);
        sC =sC0.*reshape(repmat(penalC*dxPhys.^(penalC-1),1,16)',[],1);
   
        dK=sparse(iH,jH,sH);
        dC=sparse(iC,jC,sC);

        dK = (dK+dK')/2;
        dC = (dC+dC')/2;

        dKc = dK(freedofs,freedofs);
        dCc = dC(freedofs,freedofs);

        dVKc = VX*dKc*VX';
        dVCc = VX*dCc*VX';
        dVKc = (dVKc+dVKc')/2;
        dVCc = (dVCc+dVCc')/2;
        dCc_inv = -VCc\dVCc/VCc; 
        dVA = - dCc_inv * VKc - VCc\dVKc;
%         dVA = - dVCc\VKc - VCc\dVKc;
        dVP = lyap(VA',dVA'*VP + VP*dVA);
        dJp(i) = T0'*dVP*T0;
    end
end
toc
tic
for jj = 1:1
    sH=sH0.*reshape(repmat(xPhys.^penalK,1,16)',[],1);
    sC=sC0.*reshape(repmat(xPhys.^penalC,1,16)',[],1);

    K=sparse(iH,jH,sH);
    C=sparse(iC,jC,sC);

    K = (K+K')/2;
    C = (C+C')/2;

    Kc = K(freedofs,freedofs);
    Cc = C(freedofs,freedofs);
    T = -100*ones(size(Kc,1),1);
    A = -Cc\Kc;
    P = lyap(A',Q);
    J2 = T'*P*T
    numofelem = length(element(:,1));

    for i = 1:1
        dxPhys = zeros(numofelem,1);   
        dxPhys(50) = xPhys(50);
        sH =sH0.*reshape(repmat(penalK*dxPhys.^(penalK-1),1,16)',[],1);
        sC =sC0.*reshape(repmat(penalC*dxPhys.^(penalC-1),1,16)',[],1);
   
        dK=sparse(iH,jH,sH);
        dC=sparse(iC,jC,sC);

        dK = (dK+dK')/2;
        dC = (dC+dC')/2;

        dKc = dK(freedofs,freedofs);
        dCc = dC(freedofs,freedofs);

        dCc_inv = -Cc\dCc/Cc;

        dA = - dCc_inv * Kc - Cc\dKc;
        dQ2 = dA'*P + P*dA ;

        dP2 = lyap(A',dQ2);

        dJ2(i) = T'*dP2*T;

    end
end
toc

dlta = 0.000001;
xPhys(50) = xPhys(50) + dlta;


sH2=sH0.*reshape(repmat(xPhys.^penalK,1,16)',[],1);
sC2=sC0.*reshape(repmat(xPhys.^penalC,1,16)',[],1);

K2=sparse(iH,jH,sH2);
C2=sparse(iC,jC,sC2);

K2 = (K2+K2')/2;
C2 = (C2+C2')/2;

Kc2 = K2(freedofs,freedofs);
Cc2 = C2(freedofs,freedofs);
T = -100*ones(size(Kc,1),1);
A2 = -Cc2\Kc2;
P1 = lyap(A2',Q);
J3 = T'*P1*T
(J3-J2)/dlta