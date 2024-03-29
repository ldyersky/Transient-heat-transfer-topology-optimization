
% clc;
clear;

%% Input

node_info = csvread('node6400.csv');
node = node_info(:,1:2)*0.001;

element = csvread('element6400.csv');

%%

NumOfNode = length(node(:,1));

designdomain = 1:length(element(:,1));

Fixnodes = 1;

FixDofs = Fixnodes;

AllDofs = 1:NumOfNode;

freedofs = setdiff(AllDofs,FixDofs);

FixTemp = 100;

%%
[HH,Hs] = matrix2DH(80,80,1.5);
%%
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
kk = xPhys'*xPhys;

kxx = k*ones(length(element(:,1)),1);
[iH,jH,sH0]=assemble_AXTH4_N(element,node,kxx);
Se = (node(element(1,1),1) - node(element(1,2),1))^2 + .......
    (node(element(1,1),2) - node(element(1,2),2))^2;
[iC,jC,sC0]=assemble_AXTH4_C(element,rho,cp,Se);

tic
for jj = 1:100
    loop = loop+1;
    sH=sH0.*reshape(repmat(xPhys.^penalK,1,16)',[],1);
    sC=sC0.*reshape(repmat(xPhys.^penalC,1,16)',[],1);
    K=sparse(iH,jH,sH);
    C=sparse(iC,jC,sC);

    K = (K+K')/2;
    C = (C+C')/2;
    Kc = K(freedofs,freedofs);
    Cc = C(freedofs,freedofs);
    K12 = K(freedofs,Fixnodes);

    Kc_inv = inv(Kc);
    Kc_inv = (Kc_inv+Kc_inv')/2;
    Cc_inv = inv(Cc);
    Cc_inv = (Cc_inv+Cc_inv')/2;

    T = 100;
    A = -Cc\Kc;
    x0 = zeros(1,6560);
    func = @(t,x) A * x - Cc\(K12*T);
    [t,x] = ode15s(func,0:0.001:20,x0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VX = pca(x);
    VX = VX(:,1:50);
    VX = VX';
    VKc = VX*Kc*VX';
    VKc = (VKc+VKc')/2;
    VCc = VX*Cc*VX';
    VCc = (VCc+VCc')/2;
    VVX = pinv(VX);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = zeros(length(AllDofs));
    Q(6561,6561) = 1;
    Q = Q(freedofs,freedofs);
    VQ = VX*Q*VX';
    VKc_inv = inv(VKc);
    VKc_inv = (VKc_inv+VKc_inv')/2;   
    VCc_inv = inv(VCc);
    VCc_inv = (VCc_inv+VCc_inv')/2;    

    VA = -VCc\VKc;
    VP = lyap(VA',VQ);

    T0 = -VKc\VX*K12*T;
    J = T0'*VP*T0
    JJ(jj) = J;

    numofelem = length(element(:,1));
    dT = VX*K12*T;

    for i = 1:6400
        dxPhys = zeros(numofelem,1);   
        dxPhys(i) = xPhys(i);
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
        dVP = lyap(VA',dVA'*VP + VP*dVA);
        dKc_inv = -VKc\dVKc/VKc;
        dKc_inv = (dKc_inv+dKc_inv')/2;
        dK12 = dK(freedofs,Fixnodes);
        dT0 = -dKc_inv*dT-VKc\VX*dK12*T;
        dJ(i) = T0'*dVP*T0 + dT0'*VP*T0 +T0'*VP*dT0;
    end

    xval  = xPhys;
    xmin = 0.01*ones(6400,1);
    xmax = ones(6400,1);
    f0val = J;
    df0dx = HH*(dJ'./Hs);
    fval = zeros(2,1);
    fval(1,:) = sum(xval)-volfrac*6400;
    fval(2,:) = -xval'*xval+kk;
    dfdx = zeros(2,6400);
    dfdx(1,:) = ones(1,6400);
    dfdx(2,:) = -2*xval';
    a = zeros(2,1);
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low1,upp1] = ...
    mmasub(2, 6400, loop, xval, xmin, xmax, xold1, xold2, ...
    f0val,df0dx,fval,dfdx,low,upp,1,a,1000,0);
    xold2 = xold1;
    xold1 = xPhys;
    low = low1;
    upp = upp1;
    change = max(abs(xmma(:)-xPhys(:)));
    xsum = sum(xmma(:));
    xPhys = (HH*xmma)./Hs;
    xPhys(1) = 1;
    xPhys(6400) = 1;
    kk = xmma'*xmma;
    fprintf(' It.:%5i   Obj.:%11.4f ch.:%7.3f\n',loop,mean(xmma(:)),change);
    colormap(gray); imagesc(-reshape(xmma,80,80)); axis equal; axis tight; axis off;pause(1e-6);
end
toc
% csvwrite('11.csv',xPhys);
% csvwrite('11JJ0.25.csv',JJ);