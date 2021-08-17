clc;
clear;
addpath ..;
% addpath ..\mma;
% addpath ..\ipopt;
load FEMmodel.mat;
global inputstrs
%%
%Select the filter type, ft=2 Helmholtz,ft=3 Helmholtz + Heaviside
ft=3;
%%
%Compute the elastic modulus and the yielding stress of every element
%Case 1: uniform elastic modulus and yielding stress
% E=E*ones(size(element,1),1);stressp=818*ones(size(element,1),1);
%Case 2: temperature dependent elastic modulus and yielding stress
% load nonuniformtempnofix.mat
% load nonuniformtempnofix_zenggao.mat
aaa=element';Te=mean(T(aaa))';
[E,stressp]=Inconel718(Te);
%%
% Scale factor for objective and constraint
strainenergyscale=1e8;volscale=1e2;strsscale=1e3;
%%
%Small vlue used to set the minimum allowed value the physical properties
smallfac=1e-6;
E0=E*smallfac;%Elastic moudulus
kheat=k;kheat0=kheat*smallfac;%Heat conduction factor
rho0=rho*smallfac;%Density
%%
%Penalty factor for the interpolation of various physical properties
%Elastic modulus
penal=3;%SIMP
Re=16;%RAMP
%density
penalrho=1;
%heat conduction
Rc=8;
%%
%parameters for stress constraint
q=0.5;
pnorm=8;%the exponent for P-norm stress
% stressp=818;%maximum allowable stress
%%
%Initialize the design variables and save it to *.mat file
sss = 0;
deltax=0.5e-6;
variablenum=length(designdomain);
if sss == 1
    testi = randi(variablenum);
    x=0.2+(1-0.2)*rand(variablenum,1);
    xint0=x;
    save mindu.mat x testi;
elseif sss==2
    load mindu.mat
    x(testi)=x(testi)+deltax;
else
    x=ones(variablenum,1)*volfrac;
    x=ones(variablenum,1)*0.5;
    x=0.2+(1-0.2)*rand(variablenum,1);
%     load StrainEnergyStress07.mat xint0
%     x=xint0;
    xint0=x;
end
%%
%Assemble the matrix related to Helmholtz filter
filterRadius=6; % Filter radius must be larger than 2.5.
filterType=1;
if filterType==1
    %1:Isotropic;
    k0=filterRadius^2/12*eye(2);
else
    %2:anisotropic.
    k0=diag(filterRadius).^2/12;
end
[KF,TF,VE]=assemble_filter_Q4(element(designdomain,:),node,k0);
KF=KF(freedofs_F,freedofs_F);
LF=ichol(KF);
TF=TF(freedofs_F,:);
%%
%Obtain the physical density via Helmholtz filter and/or the Heaviside
%projection
beta=2;
xPhys=ones(size(element,1),1);
if ft==1
    xPhys(designdomain)=x;
elseif ft==2
    xTilde=(TF'*(KF\(TF*x(:))))./VE(:);
    xPhys(designdomain) =xTilde;
elseif ft == 3
    xTilde = (TF'*(KF\(TF*x(:))))./VE(:);
    xTilde(xTilde<0)=1e-3;
    eta=0.5;
    xPhys(designdomain)=(tanh(beta.*eta)+tanh(beta.*(xTilde-eta)))./(tanh(beta.*eta)+tanh(beta.*(1-eta)));
end
%%
%Display the material distribution
input.node=node;
input.element=element;
input.Output_DensityField=1;
input.isovalue=0;
input.Output_DeformationField=0;
input.x=xPhys;
inputstrs=input;inputtemp=input;
figure(1)
postprocessing(input,designdomainvertex);
%%
%Obtain the dof for every element
edofMat=repmat(element,[1,2])*2-repmat([ones(1,4),zeros(1,4)],[size(element,1),1]);
edofMat=edofMat(:,[1,5,2,6,3,7,4,8]);
%%
%Other parameter setting
low=0;upp=1;
xmin=1e-3;xmax=1;
xold1=x;xold2=xold1;
change=1;
loop=0;
loopbeta=0;
objlist=[];
rholist=[];
kkkscale=1;
%%
% xPhys=ones(size(element,1),1);
% load MassStress03.mat rholist
% xPhys=rholist(:,end-1);
while change > 0.001 && loop<500
    loop = loop+1;
    loopbeta=loopbeta+1;
    %%
    %The interpolation of elastic modulus used in the stiffness matrix
    %Case 1: SIMP, parameter penal
%         Ek=E0+(E-E0)*xPhys.^penal;
%         dEk=penal*(E-E0)*xPhys.^(penal-1);
%         dkr=dEk./Ek;
    %Case 2: RAMP, parameter Re
    Ek= E0+(E-E0).*xPhys./(1+Re*(1-xPhys));
    dEk=(E-E0).*(1+Re)./(1+Re*(1-xPhys)).^2;
    dkr=dEk./Ek;
    %%
    %The interpolation of equivalent heat load
    %Case 1: The heat load interpolated using the elastic modulus of
    %stiffness
    thdependbool=0;
    Eo = Ek;
    dEo = dEk;
    Betao = Eo*alpha;
    dBetao = dEo*alpha;
    diffBeta2toEi=(2*Betao.*dBetao.*Eo-dEo.*Betao.^2)./(Eo).^2;
    %Case 2: The heat load interpolated using independent elastic modulus
    %     thdependbool=1;
    %     Rb=8;
    %     Eo = Ek;
    %     dEo = dEk;
    %     Betao = xPhys./(1+Rb*(1-xPhys))*E*alpha;
    %     alphaInterpolated=Betao./Eo;
    %     dBetao = (1+Rb)./(1+Rb*(1-xPhys)).^2*E*alpha;
    %     diffBeta2toEi=(2*Betao.*dBetao.*Eo-dEo.*Betao.^2)./(Eo).^2;
    %%
    %The interpolation of density
    %Case 1: SIMP, parameter penal
    RhoF=rho0+(rho-rho0).*xPhys.^penalrho;
    dfr=penalrho*(1-smallfac)*xPhys.^(penalrho-1)./(smallfac+(1-smallfac)*xPhys.^penalrho);
    %Case 2: RAMP, parameter Re
    %%
    %The interpolation of the heat conduction factor k
    %SIMP
    %     Kheat=kheat0+(kheat-kheat0)*xPhys.^penal;
    %     dKheat=penal*(1-smallfac)*xPhys.^(penal-1)./(smallfac+(1-smallfac)*xPhys.^penal);
    %RAMP
    Kheat=kheat0+(kheat-kheat0)*xPhys./(1+Rc*(1-xPhys));
    dKheat=(kheat-kheat0)*(1+Rc)./(1+Rc*(1-xPhys)).^2;
    dKheat=dKheat./Kheat;
    %%
    %The interpolation of elastic modulus used in the calculatio of
    %stress,q-p method
    Estrs=E.*xPhys.^q;
    dEstrs=q*E.*xPhys.^(q-1);
    %     Estrs=E*xPhys.^penal./(xPhys.^q);
    %     dEstrs=E*(penal-q)*xPhys.^(penal-q-1);
    %%
    %Temperature field 1: Compute the nodal temperature by solving heat conduction equation
    [iH,jH,sH]=assemble_AXTH4_mex(element,node,Kheat);
    H=sparse(iH,jH,sH);
%     Fth=H(freedofs_Th,Ls)*ones(length(Ls),1)*Lth+H(freedofs_Th,Rs)*ones(length(Rs),1)*Rth;
%     T=zeros(size(H,1),1);
%     T(Ls,:)=Lth;T(Rs,:)=Rth;
%     T(freedofs_Th,1)=H(freedofs_Th,freedofs_Th)\-Fth;
    %Temperature field 2: Uniform temperature
    %     Rth=500;
    %         Rth=T_ref;
    %         T=ones(size(node,1),1)*Rth;
    %load nonuniformtempnofix.mat
    %load nonuniformtempnofix.mat
    inputtemp.x=T;
    figure(2);
    postprocessing_stress(inputtemp)
    %%
    %Assemble global stiffness matrix
    [iK,jK,sK]=assemble_AXQ4_mex(element,node,Ek,nu*ones(size(xPhys)));
    K=sparse(iK,jK,sK);K=(K+K')/2;
    %%
    %Assemble inertial load
    [iF,jF,sF,V]=assemble_AXQ4_inertia_mex(element,node,RhoF,w);
    force_inertia=sparse(iF,jF,sF,size(node,1)*2,1);
    %%
    %Compute the equivalent thermal load and initial thermal strain energy
    if thdependbool==0
        %[force_t,initialThermalStrainEnergy,BDC_t,dFp2_t]=assemble_AXQ4_ThForceandThStrainEnergy_mex(element,node,Ek,dEk,nu*ones(size(xPhys)),alpha*ones(size(xPhys)),T-T_ref);
        [force_t,initialThermalStrainEnergy,BDC_t,dFp2_t,sE]=assemble_AXQ4_ThForceandThStrainEnergyConduction_mex(element,node,Ek,dEk,nu*ones(size(xPhys)),alpha*ones(size(xPhys)),T-T_ref,Betao.^2./Eo);
    else
        [force_t,initialThermalStrainEnergy,BDC_t,dFp2_t]=assemble_AXQ4_ThForceandThStrainEnergy_mex(element,node,Ek,dBetao./alphaInterpolated,nu*ones(size(xPhys)),alphaInterpolated,T-T_ref);
    end
    BDC_t=sparse(BDC_t.i,BDC_t.j,BDC_t.s);
    %Computation method 1:
    force_t=sparse(force_t.i,force_t.j,force_t.s);
    %Computation method 2:
    %force_t=BDC_t*(T-T_ref);
    %%
    %Compute total load
    Fm=force1+force_inertia;Fth=force_t;
    F=Fm+Fth;
    %%
    %Solve displacement
    U=zeros(size(F));Um=zeros(size(F));Uth=zeros(size(F));
    Utemp=K(freedofs,freedofs)\[Fm(freedofs,1) Fth(freedofs,1)];
    Um(freedofs)=Utemp(:,1);Uth(freedofs)=Utemp(:,2);
    U(freedofs)=Um(freedofs)+Uth(freedofs);
    %%
    %Compute the strain energy and its sensitivity w.r.t physical desntity
    objlist(loop)=(0.5*(K*U-2*Fth)'*U+sum(Betao.^2./Eo.*initialThermalStrainEnergy))/strainenergyscale*2;
    
    %dse=sensitivity_AXQ4_inertia_ThermalStrainEnergy(Um,edofMat,sF,dfr)+sensitivity_AXQ4_ThermalStrainEnergy(Um,Uth,edofMat,sK,dkr)+...
    %    sensitivity_AXQ4_ThForce_ThermalStrainEnergy(element,Uth,T,BDC_t,dFp2_t,sH,H,freedofs_Th,dKheat)+initialThermalStrainEnergy.*diffBeta2toEi;
    
    A=-Uth'*BDC_t;lambda=zeros(size(BDC_t,2),1);lambda(freedofs_Th)=H(freedofs_Th,freedofs_Th)\A(freedofs_Th)';
    lambda1=zeros(size(BDC_t,2),1);
    lambda1(freedofs_Th)=-H(freedofs_Th,freedofs_Th)\sE(freedofs_Th);
    dse=sensitivity_AXQ4_inertia_ThermalStrainEnergy_mex(Um,edofMat,sF,dfr)+sensitivity_AXQ4_ThermalStrainEnergy_mex(Um,Uth,edofMat,sK,dkr)+...
        sensitivity_AXQ4_ThForceMex_ThermalStrainEnergy_mex(element,Uth,T,lambda,dFp2_t,sH-sH,dKheat)+...
        sensitivity_AXQ4_InitStrainE_ThermalStrainEnergy_mex(initialThermalStrainEnergy,diffBeta2toEi,element,T,sH-sH,lambda1,dKheat);
    %     dse=initialThermalStrainEnergy.*diffBeta2toEi;sensitivity_AXQ4_ThForceMex_ThermalStrainEnergy_mex(element,Uth,T,lambda,dFp2_t,sH,dKheat)+...
    %     dse=sensitivity_AXQ4_inertia_ThermalStrainEnergy(Um,edofMat,sF,dfr)+sensitivity_AXQ4_ThermalStrainEnergy(Um,Uth,edofMat,sK,dkr)+...
    %         sensitivity_AXQ4_ThForce_ThermalStrainEnergy(element,Uth,T,BDC_t,dFp2_t,sH-sH,H,freedofs_Th,dKheat);
    dse=dse(designdomain)/strainenergyscale*2;
    %%
    %Compute volume constraint value and its sensitivity w.r.t physical density
    cv(loop)=(xPhys(designdomain)'*V(designdomain)/sum(V(designdomain))-volfrac)*volscale;
    dcvdx=V(designdomain)/sum(V(designdomain))*volscale;
    %%
    %     [maxstress,dmaxstress,maxstress_exact]=maxstress_calc_AXQ4TH(Estrs,dEstrs,nu*ones(size(xPhys)),alpha*ones(size(xPhys)),node,element,U,...
    %         T,T_ref,K,sK,sF,edofMat,freedofs,dkr,dfr,pnorm,stressp,1:size(element,1),BDC_t,dFp2_t,sH,H,freedofs_Th,dKheat);
    [maxstress,dmaxstress,maxstress_exact]=maxstress_calc_AXQ4TH_MatDependent(Estrs,dEstrs,nu*ones(size(xPhys)),alpha*ones(size(xPhys)),node,element,U,...
        T,T_ref,K,sK,sF,edofMat,freedofs,dkr,dfr,pnorm,stressp,1:size(element,1),BDC_t,dFp2_t,sH-sH,H,freedofs_Th,dKheat);
    %Stablizing corection paramter
    %if loop==1
    %    cp(loop)=maxstress_exact/maxstress;
    %else
    %    cp(loop)=cp(loop-1)+0.1*(maxstress_exact/maxstress-cp(loop-1));
    %end
    k=maxstress_exact/maxstress;
    %Case 2: Used only in the verification of the sensitivity of stress
    %k=0.8;
    %Compute the stress constraint value and its sensitivity w.r.t physical
    %density
    cs(loop)=(maxstress*k-1)*strsscale;
    dmaxstress=dmaxstress(designdomain)*k*strsscale;
    %%
    %Compute the sensitivity of objective and constrint w.r.t design
    %variables by the chain rule
    if ft == 1
        dcvdx(:)=dcvdx(:);
        dmaxstress(:) = dmaxstress(:) ;
    elseif ft == 2
%         dc(:) = TF'*(KF\(TF*(dc(:)./VE(:))));
        dcvdx(:) = TF'*(KF\(TF*(dcvdx(:)./VE(:))));
        dmaxstress(:) = TF'*(KF\(TF*(dmaxstress(:)./VE(:))));
        dse(:) = TF'*(KF\(TF*(dse(:)./VE(:))));
    elseif ft == 3
        eta=0.5;
        dx = -(beta.*(tanh(beta.*(eta - xTilde)).^2 - 1))./(tanh(beta.*eta) - tanh(beta.*(eta - 1)));
%         dc = TF'*(KF\(TF*(dx(:).*dc(:)./VE(:))));
        dcvdx = TF'*(KF\(TF*(dx(:).*dcvdx(:)./VE(:))));
        dmaxstress = TF'*(KF\(TF*(dx(:).*dmaxstress(:)./VE(:))));
        dse(:) = TF'*(KF\(TF*(dx(:).*dse(:)./VE(:))));
    end
    %%
    %Display the validation result of the sensitivities of strain energy and the stress
    if sss == 1
        value1o=cv(loop);value2o=cs(loop);value3o=objlist(loop);ana1=dcvdx(testi);ana2=dmaxstress(testi);ana3=dse(testi);
        save("sss1.mat","value1o","value2o","value3o","ana1","ana2","ana3");
    elseif sss==2
        load sss1.mat;
        value1e=cv(loop);value2e=cs(loop);value3e=objlist(loop);
        fdm1=(value1e-value1o)/deltax;fdm2=(value2e-value2o)/deltax;fdm3=(value3e-value3o)/deltax;
        disp([' 初始： ' sprintf('%22.16f ',value1o)  sprintf('%22.16f ',value2o) sprintf('%22.16f ',value3o)]);
        disp([' 摄动： ' sprintf('%22.16f ',value1e)  sprintf('%22.16f ',value2e) sprintf('%22.16f ',value3e)]);
        disp([' 差分： ' sprintf('%22.16f ',fdm1)  sprintf('%22.16f ',fdm2) sprintf('%22.16f ',fdm3)]);
        disp([' 解析： ' sprintf('%22.16f ',ana1)  sprintf('%22.16f ',ana2) sprintf('%22.16f ',ana3)]);
        disp([' 误差： ' sprintf('%22.16f ',(ana1-fdm1)/fdm1*100)  sprintf('%22.16f ',(ana2-fdm2)/fdm2*100) sprintf('%22.16f ',(ana3-fdm3)/fdm3*100)]);
    else
%         break;
    end
    %%
    %Construct the vector of constraints and its sensitivity
    v=[cv(loop);cs(loop)];
    dv=[dcvdx(:)';dmaxstress(:)'];
%     v=[cs(loop)];
%     dv=[dmaxstress(:)'];
    %%
    %Prepare the parameters for MMA solver
    m=length(v);n=variablenum;
    ccc=1000*ones(m,1);a0 = 1;a = zeros(m,1);d = ones(m,1);
    %%
    %Obtain new design variables via MMA solver
    [xnew,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,x,xmin,xmax,xold1,xold2, ...
        objlist(loop),dse,v,dv,low,upp,a0,a,ccc,d);
    xold2=xold1;xold1=x;
    if loop==1
        x=xnew;
    else
        x=xold1+0.3*(xnew-xold1);
    end
    %%
    %Calculate the change of design varibles
    change=max(abs(x-xold1));
    chageintercal=10;
    if loop>chageintercal
        change=max( abs( objlist(end-chageintercal:end)-mean(objlist(end-chageintercal:end)) ) )/mean(objlist(end-chageintercal:end));
    end
    %%
    %Obatin the physical densities via Helmholtz filtering and/or Heaviside
    %projection
    if ft == 1
        xPhys(designdomain) = x;
    elseif ft == 2
        xPhys(designdomain) = TF'*(KF\(TF*x(:)))./VE(:);
    elseif ft == 3
        xTilde=(TF'*(KF\(TF*x(:))))./VE(:);
        xTilde(xTilde<0)=1e-3;
        eta=0.5;
        xPhys(designdomain)=(tanh(beta.*eta)+tanh(beta.*(xTilde-eta)))./(tanh(beta.*eta)+tanh(beta.*(1-eta)));
    else
        xPhys(designdomain) = x;
    end
    %%
    %Save the density of every iteration
    rholist=[rholist,xPhys];
    %%
    %Diplay the objective and constraint values in the commmand line
    disp([' It.: ' sprintf('%4i',loop) ' StrainE.: ' sprintf('%10.4f',objlist(loop)) ' VolFrac.: ' sprintf('%10.4f',cv(loop)/volscale) ' StrsCons.: ' sprintf('%10.4f',cs(loop)/strsscale)...
        ' ch.: ' sprintf('%6.3f',change )]);
    %%
    %Display the density field
    input.node=node;
    input.element=element;
    input.Output_DensityField=1;
    input.isovalue=0;
    input.Output_DeformationField=0;
    input.x=xPhys;
    figure(1)
    postprocessing(input,designdomainvertex);
    %%
    %Adjust the beta parameter for the Heaviside projection
    if ft == 3 && beta < 32 && (loopbeta >= 50 || change <= 0.001)
        beta = 1.5*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
end
save result.mat;