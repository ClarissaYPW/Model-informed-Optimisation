function [Data] = kineticgeneratorflexy(Trange,ini, mon, cta,pol, tmax, resiconc,CReact_SM, In0,t2)

ctaconc=resiconc;
T=Trange;

%% Section 1. Rate constants for kinetics     
%% initiator information
%find the decomposition rate constants for each condition systems 
halflife = ini.hl; % half life time in seconds = 10 h at 44 degs 
R= 8.314; % Gas constant (J/mol/K)
kT1 = log(2)/halflife;  % use general half life equation 
% finding the preexponential factor, A for the above data using
% Arrhenius.

T1=ini.T1+273.15;
Xp1 = exp(-ini.Ea/(R*T1)); %find exponential segment of arrhenius (1) equation 
preXpA = kT1/Xp1; %find preexponential factor, A.
T2a = T; 
T2=T2a+273.15;
Xp2 = exp(-ini.Ea./(R.*T2)); %find exponential segment of arrhenius (2) equation.
kT2 = (preXpA).*Xp2; %solve for kT2.
 
%Termination rate constant assuming full control
%% monomer information       

kt=mon.At.*exp(-(mon.Eat./(R*T2))); %initiator termination rate constant            
kd=kT2; 
kp3 = mon.Ap.*exp(-(mon.Eap./(R.*T2)));
kct=0.25.*kt;
%%cta information

Aadd =cta.Aadd;
 
kadd = cta.Aadd.*exp(-(cta.Eaadd./(R.*T2)));
kdprime= 2.*ini.f.*kd;
kfrag = kadd/cta.K
% symbolic solvation 
syms I(t) k I0
eqnx= -diff(I,t)== k.*I;
condx = I(0)==I0;
I(t)=dsolve(eqnx, condx);
initime= matlabFunction(I);
Ini = initime(t2,In0, kdprime);

I=Ini;
dt = numel(t2);

%% Section 5. Finding the changes in starting materials involved in the RAFT Equilibria%% 
% Macromolecules 2023, 56, 4, 1581–1591

rI= (2.*kd.*ini.f.*In0).*(exp(-(kd.*t2)));
num = numel(In0);
%solve ODE symbolically to retain the time frame used
syms T(t) kx ky kz ko ri T0
eqn= diff(T,t)== -((ri)./((kx.*ky)./((kz.*ko).*T)+2)); 
cond = T(0)==T0;
T(t)=dsolve(eqn, cond);
ctagent= matlabFunction(T);
      
for ii=1:dt
CTA (1,ii)= ctagent(t2(1,ii),resiconc,kadd,kt,kfrag,kct,rI(1,ii));

end  
for ii=1:dt
R(1,ii) = sqrt(((kd*kfrag)/(kct*kadd)).*(1+((kt*kfrag)/(2*resiconc*kct*kadd))))./(sqrt(((1+((((1+(((kt*kfrag)/(2*resiconc*kct*kadd)))).^(2)))*((resiconc/(ini.f*In0))-1))*exp(-kd.*t2(1,ii))))));
end
    

R2 = R(:,2)';

   

%% Section 6. Finding monomer conversion via pseudo first order kinetic equatons 
kp2 = kp3.*R2;  
% finds the kprime constant =kp2=k[R]
%first order ode solved symbolically due to speed. 
%first iteration of monomer conversion
    syms M(t) kr M0 
    eqn9 = -diff(M,t) == kr.*M;
    cond9 = M(0)==M0;
    M(t) = dsolve(eqn9,cond9);
    monomer = matlabFunction(M);
     for ii=1:dt 
    Monomer(1,ii) = monomer(t2(1,ii),CReact_SM,kp2); % substitute initial valu
    M1(1,ii) = Monomer(1,ii)./CReact_SM;
    
        end 
     conver = 1- M1;
     %plot Monomer Conversion 

%% Section 7. Chainlength dependent termination and its effect on conversion and MWD 
% finding chain length dependent termination effects on MWD 
%adjust kt to include a dependence
L = ((CReact_SM-Monomer)./resiconc);
alphas= mon.alphas;
Lc = mon.Lc;
alphaL = mon.alphaL;
ktL1=zeros(1,dt);
ktL2=zeros(1,dt);
% two schemes observed for chainlength dependence. 
% L is the DP, kt is the termination of initiator rdicals, alpha s is the
% powerlaw for low molecular weight chains 
% Lc is the cross over DP between small and long chains 

for ii=1:dt
if L(1,ii)<Lc 
ktL1(1,ii) = kt.*L(1,ii).^(-alphas);
elseif L(1,ii)>Lc 
ktL2(1,ii)= kt.*(L(1,ii).^(-alphas)).*(Lc.^(-alphas+alphaL));
end 
ktL(1,ii)=ktL1(1,ii)+ktL2(1,ii); 
end 
kct=0.25*ktL;
%livingness and its relationship to initiator concentration.
fc=0.8;

      for ii= 1:dt
%                  R(1,ii) = sqrt(rI(1,ii)./(kt(1,ii).*(1+2.*CTA(1,ii).*((kct(1,ii).*kadd)./(ktL(1,ii).*kfrag)))));         
R(1,ii) = (sqrt(((kd*kfrag)/(kct(1,ii)*kadd))*(1+((ktL(1,ii)*kfrag)/(2*resiconc*kct(1,ii)*kadd)))))/(sqrt(((1+((((1+(((ktL(1,ii)*kfrag)/(2*resiconc*kct(1,ii)*kadd))))^2))*((resiconc/(ini.f*In0))-1))*exp(-kd*t2(1,ii))))));
      end     

    R2 = R(:,2);
    kp2 = kp3.*R2;
    Monomer = monomer(t2,CReact_SM,kp2); 
    M1 = Monomer./CReact_SM;
          
    
     conver = 1- M1;

%Finding average molecular weight using equation seen in RAFT
%Review(Macromolecules 2017, 50, 19, 7433–7447)
    fc=0.8; %?? is the coupling factor, defined as the ratio of polymers that terminated via coupling. ?_?=1, 100 % termination by coupling and ?_?= 0 is 100%coupling via disproportionation, this is dependent on monomer, acrylamides tend to terminate by coupling so a value of 0.8 is assumed 
  
        for ii=1:dt
    Mn(1,ii) = ((CReact_SM.*(conver(1,ii)).*(mon.Mr))./(resiconc+((2*ini.f).*(In0)).*(1-exp(-kd.*t2(1,ii))).*(1-(fc./2))))+cta.Mr;
        end 

    
    Mna= Mn(1,dt);% selects the molecular weight at a certain time 
    
c=conver;
XDP = 0:1:pol.DP*4;
aveDP = ((Mna-cta.Mr)./mon.Mr);
kmadd =0;
ini.kd = kd;
p=1;
c = conver(1,dt);
[Disp3]=Dispersity_Calc_ap(aveDP,p,conver,kadd,kfrag,kp3,resiconc,t2, Monomer,kt, cta.phi,ini);
Mp=aveDP*mon.Mr;
Disp=Disp3';
[Data] = [c  Disp];
