function [Final] =  Initialisationsteadstate(Mon,cta, ini, sol, pol, nc)
%% code to calculate changes in relative concentrations of reagents and radicals%%
%% Section 1. tackles finding the initial concentrations of reagents and reservoir concentration%%
%Reactor volume and residence time are changable and  allow the
%concentrations to be calculated by telling us a volume of monomer to start
%the calculation
ReactV =pol.ReactV; %take input frsom reactor volume as ReactV(mL)
%%collect data from app relevant fields 

m_Mr =Mon.Mr;
m_dens=Mon.dense;

CTA_Mr =cta.Mr ;

 Ini_Mr=ini.Mr;
 
  
 dense_sol=sol.dense;
  m_V = 4*ReactV;%monomer volume  is 4 x reactor volume(arbitrary value
  m_mass = m_V; %find mass of monomer (g)m_dens/
  m_moles = (m_mass)/(m_Mr); %find moles of monomer (mol)
 
DP = nc(:,4);
defCTAIni = nc(:,1);
wtpercent=nc(:,5);
RTim = nc(:,3)
 Temp = nc(:,2)
 
%calculate concentration of CTA
  CTA_moles= (m_moles./DP); % find CTA moles (mol)
  CTA_mass=CTA_moles.*CTA_Mr;
  ratio_CTA=1;
 molefactor=ratio_CTA./defCTAIni;
Iniesmol = CTA_moles./molefactor;
Iniemass= Iniesmol.*Ini_Mr;
wtpercent = 30; %put wt percentage here 
m_solute = m_mass+CTA_mass+Iniemass; % sum of all solute masses (g)
wt=wtpercent/100;
m_solvent=((100/wtpercent)*m_solute)-m_solute;
%dense_sol = 0.999; %methanol(g/mL);
V_solvent = (m_solvent/dense_sol)/1000 %(mol/dm3)
Ini0=Iniesmol./V_solvent; %concentration of initiator in reactor 
CTAconc= CTA_moles./V_solvent;
Monconc = m_moles./V_solvent;

%% set time matrix

st=zeros(numel(RTim),1);
       for i=1:numel(RTim)
            dt = 199;      % 30 s intervals
            dt2(i) =((RTim(i))/dt);  %vector containing 267 elements
       end 
      % dt2=dt2';
%       
       for i= 1:numel(RTim)
           
time(i,:)= st(i):dt2(i):RTim(i);
       end 

       ti=time.*60;

% set the conditions

        
            
            
        Dat=zeros(numel(RTim),2);
      for j=1:numel(RTim)
        
         resiconc=CTAconc(j);
         In0=Ini0(j)
         CReact_SM=Monconc(j);
         Trange=Temp(j);
          t2=ti(j,:)
          tmax=RTim(j);
          
  [Data]= kineticgeneratorflexy(Trange,ini, Mon, cta, pol, tmax, resiconc,CReact_SM, In0,t2);  
     Dat(j,:)= Data;
  
      end 
     
Final = [nc Dat];    
     end 
 