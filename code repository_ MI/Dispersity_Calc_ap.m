function [Disp3]=Dispersity_Calc_ap(aveDP,p,conver,kadd,kfrag,kp3,resiconc,t2, Monomer,kt, phi, ini)
%% this script calculates the molar mass dispersity from the simulated conversion 
c=conver(:,end); 
tmax = t2(:,end)

xp = kt; %get values of ktL at iteration p as defined in main code
cx = c;
kctr =(kadd*phi);
kctri = (kfrag/phi)-kfrag;
Ctr=kctr/(kp3);
kctri = (kfrag/phi)-kfrag;
f= ini.f
kd = ini.kd
monco= Monomer(1,1);
for i=1:numel(aveDP) 
   Disp(i)=1+(resiconc./(monco*c(i)))+((2./c(i))-1).*(1./Ctr);%(2-c(i)/c(i))))(2-c(i))/c(i)
    
end 
Init = 2*f*kd*exp(-(kd.*xp)/(kp3));
syms R ktr ktri CTAr ri kter kp M
eqn1 = 1*ri +ktri*CTAr*R - ktr*CTAr*R - 2*kter*R^2 == 0;
Rad = solve(eqn1,R); %solve quadratic equation above 
Rad = Rad(2); %positive solution only 
ra = matlabFunction(Rad);

 Kraft= kctr./(kctri);
  for i=1:numel(aveDP)
   rad(i) = ((ra(resiconc, xp(i),kctr, kctri,Init(i))));
   a1(i) = (xp(p)*monco)./(4*kp3*resiconc);
  a(i) =a1(i).*rad(i).*cx(i);
  end 
  for i=1:numel(aveDP)

  
Disp3(i) = Disp(i)+a(i)+0.103;%% add 0.10 to account for online broadening 

  end 

for i=1:numel(aveDP)
    
if Disp3(i)>=2
    Disp3(i)=2;
end 


end 

for i=1:numel(aveDP) 
  
if c(i)> 0.50
    Disp3(i) = Disp3(i)+0.07;%% accounts for broadening due to RTD

end 
end 

