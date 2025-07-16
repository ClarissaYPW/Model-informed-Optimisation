function varargout = DigitalaugParetomultix(varargin)
%% Author Clarissa Wilding, University of Leeds (C.Y.Wilding@leeds.ac.uk)
%% Serial commands are sent via the PolymerSynthesis_multi.m to CO2 GUI 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DigitalaugParetomultix_OpeningFcn, ...
                   'gui_OutputFcn',  @DigitalaugParetomultix_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Paretomultix is made visible.
    function DigitalaugParetomultix_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Paretomultix (see VARARGIN)

% Choose default command line output for Paretomultix
handles.output = hObject;

if ~(nargin == 4)
    error('Incorrect number of input arguments')
else
    % stashes away the data
    PolymerGUI = varargin{1};
end
%Save PolymerGUI handles to callback to later

setappdata(hObject,'PolymerGUI',PolymerGUI)
ReactVol=get(PolymerGUI.reactorvolume,'String');
SteadyState=get(PolymerGUI.steadystate,'String');
deswt=get(PolymerGUI.deswt,'String'); %get desired wt %
fs1 =get(PolymerGUI.fs1,'String');%get feedstockconc wt %
fs2 =get(PolymerGUI.fs2,'String'); 
fs3 =get(PolymerGUI.fs3,'String');
DilF =get(PolymerGUI.DilF,'String');
reswt =get(PolymerGUI.reswt,'String');
tDP=get(PolymerGUI.tDP,'String');
monmass =get(PolymerGUI.monmass,'String');
CTAmass =get(PolymerGUI.CTAmass,'String');
initmass =get(PolymerGUI.initmass,'String');


setappdata(hObject,'fs1',fs1)
setappdata(hObject,'fs2',fs2)
setappdata(hObject,'fs3',fs3)
setappdata(hObject,'tDP',tDP)
setappdata(hObject,'reswt',reswt)
setappdata(hObject,'monmass',monmass)
setappdata(hObject,'CTAmass',CTAmass)
setappdata(hObject,'initmass', initmass)
setappdata(hObject,'DilF',DilF)
% Choose default command line output for KineticsGUI
handles.output = hObject;
set(handles.ReactVol, 'String',ReactVol)
set(handles.SteadyVol, 'String',SteadyState)

%% Get All connected serial objects and commands from CO2GUI
%get serial object info from main CO2 GUI
%PolymerGUI string contains the Main CO2 GUI data 
objectConfig=getappdata(PolymerGUI.MainGUI,'objectConfig');
objectTypes=getappdata(PolymerGUI.MainGUI,'objectTypes');

%get the serial object names from this
names = {objectConfig.name};
%set empty value for linked objects
setappdata(hObject,'linked',zeros(8,8));
% gets names of pumps, min and max values, dx values from Pareto GUI for
% all 8 attachments 
for m = 1:8
    %get the type of object from the config setting
    type = objectConfig(m).type;
    class = objectTypes(type).class;
    nametag = sprintf('nametag%d', m); %define handle for name field to change
    checkbox = sprintf('checkbox%d', m); %define handle for checkbox to change
    min = sprintf('min%d', m); %define handle for min to change
    max = sprintf('max%d', m); %define handle for max to change
    dx = sprintf('dx%d', m); %define handle for dx to change
    serial = sprintf('serialObject%d', m); %define handle for serialobject to change
    set(handles.(nametag),'String',names{m}); % Changes Object Names
    % if it's not connected, disable the box
    if PolymerGUI.connectedObjects(m)
        set(handles.(nametag), 'Enable', 'on')
        set(handles.(checkbox), 'Enable', 'on')
        set(handles.(min), 'Enable', 'on')
        set(handles.(max), 'Enable', 'on')
        set(handles.(dx), 'Enable', 'on')        
    end
    if PolymerGUI.connectedObjects(m)
        if strcmp(class, 'loop')
        set(handles.(nametag), 'Enable', 'on')
        set(handles.(checkbox), 'Enable', 'Off')
        set(handles.(min), 'Enable', 'Off')
        set(handles.(max), 'Enable', 'Off')
        set(handles.(dx), 'Enable', 'Off')
        end
     %gets connected pumps and names 
       connected(m).name = names{m};
       connected(m).type = objectConfig(m).type;
       connected(m).sobjects = getappdata(PolymerGUI.MainGUI,serial);
       connected(m).commands = objectTypes(type);
       if strcmp(names{m},'GPC')
       connected(m).pumptype = 'GPC';
       else
       connected(m).pumptype = 'pump';    
       end
    end
end

axes(handles.axes1)
matlabImage = imread('C:\Users\Spinsolve\Documents\MATLAB\Badger\PolymerSynthesis\Rig.jpg');
image(matlabImage)
axis off
axis image

 setappdata(hObject,'connected', connected);

 handles.output = hObject;
%  setappdata(handles.output,'PolymerGUI', PolymerGUI);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Paretomultix wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DigitalaugParetomultix_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Start.
function Start_Callback(hObject, ~, handles)
% hObject    handle to Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Start
set(handles.Start,'Enable','off');
lb=[]; %empty numeric array for lb
ub=[]; %empty numeric array for ub
dx=[]; %empty numeric array for dx
connected = getappdata(handles.output, 'connected');
NumberObjects = 0;
PumpColumn = []  %empty numeric array for Pump Column
PumpObjects =[]  %empty numeric array for PumpObjects
TempObjects = []; %empty numeric array for TempObjects
TempColumn = []; %empty numeric array for TempColumn
RTColumn = []; %empty numeric array for Residence Time 
RTObjects=[]; %empty numeric array for Residence Time 
%% collects monomer type from UI
montype = get(handles.montype,'Value'); 
% kinetic rate parameters 
switch montype
    case 2
        %DMAm
       Mon.Eap = 14.1e3; %(J mol)
       Mon.Ap= 11e6;
       Mon.Mr= 99.13;  %(g/mol)
       Mon.dense= 0.962;%(g/mL)
    case 3
       %tBuAm
        Mon.Eap = 1.955e3; %(J/mol)
       Mon.Ap= 8.68e6;
       Mon.Mr= 127.13;  %(g/mol)
       Mon.dense= 1.2;%(g/mL)
    case 4
        %nBuA
       Mon.Eap = 17.9e3; %(J mol)
       Mon.Ap= 2.24e7;
       Mon.Mr= 128.17;  %(g/mol)
       Mon.dense= 0.894; %(g/mL)
         case 5
        %HEAm
       Mon.Eap = 14.1e3; %(J mol)
       Mon.Ap= 11e6;
       Mon.Mr= 115.13;  %(g/mol)
       Mon.dense= 1.111; %(g/mL)
end 
setappdata(handles.output,'Mon',Mon)
%% collects CTA type from UI inputs
ctatype = get(handles.ctatype,'Value') 
switch ctatype
    case 2
        %BM1640
       cta.Mr= 238.4; %(g/mol)
       cta.phi=0.67;
       cta.Ctr = 100.9;
       cta.K=300;
    case 3
        %BM1429
        cta.Mr= 254.4; %(g/mol)
        cta.phi=0.67;
       
       cta.K=1000;
end 
setappdata(handles.output,'cta',cta);
%% collect initiator type from UI
initype = get(handles.initype,'Value');
switch initype
    case 2
        %AIBN
     ini.Mr = 164.1;%(g/mol)
     ini.f =0.42;
     ini.T1=65; %degsC
     ini.Ea=132400;%J/mol
    case 3 
        %VA044
      ini.Mr = 323.3;%(g/mol)
     ini.f =0.6;
     ini.T1=44; %(degsC)
     ini.Ea=108200;%(J/mol)
     case 4
        %ACVA
      ini.Mr = 280.8;%(g/mol)
     ini.f =0.3;
     ini.T1=69; %(degsC)
     ini.Ea=132800;%(J/mol)
end 
    etappdata(handles.output,'ini',ini);
 soltype = get(handles.soltype,'Value') ;
 
 switch soltype
     %water
     case 2
         sol.dense=0.999;
     %methanol
     case 3
         sol.dense=0.792; %(g/mL)
         %dioxane  
     case 4
         sol.dense =1.03; %g/mL       
 end 
    setappdata(handles.output,'sol',sol);    
%% check which variables are selected for optimisation     
for i=1:8
    check = sprintf('checkbox%i',i); % collects all ticked checkboxes
    optimise = get(handles.(check),'Value'); %get check values 
    if optimise %if optimise? is checked
        NumberObjects = NumberObjects+1  
        if connected(i).type == 4
            TempObjects = horzcat(TempObjects,connected(i));
            TempColumn = horzcat(TempColumn,NumberObjects);
        elseif connected(i).type == 7
            PumpObjects = horzcat(PumpObjects,connected(i));
            PumpColumn = horzcat(PumpColumn,NumberObjects);
        end    
    end 
end 

pcol=numel(PumpColumn);
tcol=numel(TempColumn);
Colnum=pcol+tcol;
%% get optimisation bounds (defines reaction space limits)
    for i=1:9
        check = sprintf('checkbox%i',i); % collects all ticked checkboxes
         optimise = get(handles.(check),'Value'); %get check values 
      
        if optimise %if optimise? is checked
         if get(handles.checkbox9,'Value')==1 % if RT checkbox 
             RTColumn = horzcat(Colnum+1);
             NumberObjects=NumberObjects+1;
         end   
        %gets handles of min(i) where i is 1:8 and turns it into a double
        min= str2double(get(handles.(sprintf('min%i',i)),'String'));
        %gets handles of max(i) where i is 1:8 and turns it into a double
        max= str2double(get(handles.(sprintf('max%i',i)),'String'));
        %gets handles of dx(i) where i is 1:8 and turns it into a double
        step = str2double(get(handles.(sprintf('dx%i',i)),'String'));
        %concanate lb and min to form a lower bount
        lb=horzcat(lb,min);
        ub=horzcat(ub,max);
        dx=horzcat(dx,step); 
        end
    end

%in the string 'Limits' define lb, ub and dx 
Limits.lb = lb;
Limits.ub = ub;
Limits.dx = dx;
Columns.TempColumn = TempColumn; %set Columns.TempColumn string
Columns.PumpColumn = PumpColumn; %set Columns.PumpColumn string
Columns.RTColumn = RTColumn;%set Columns.RTColumn string 
isrtpresent = isfield(Columns,'RTColumn');
setappdata(handles.output,'Columns',Columns)%attach Columns string to app data

setappdata(handles.output,'Limits',Limits)%attach Limits string to app data
setappdata(handles.output,'PumpObjects',PumpObjects)%attach PumpObjects string to app data
setappdata(handles.output,'TempObjects',TempObjects)%attach TempObjects string to app data
setappdata(handles.output,'RTObjects',RTObjects)%attach TempObjects string to app data

%% conduct LHS in-silico
NumberExperiments = str2double(get(handles.training, 'string')); %convert string from training data points to Number of Experiments 
LHS = lhsdesign(NumberExperiments,NumberObjects); %perform an LHS experiment where NumberExperiments and NumberObjects are defined from the GUI
for i=1:numel(lb)
Conditions(:,i)  = LHS(:,i).*(ub(:,i)-lb(:,i))+lb(:,i); %find the conditions from multiplying all rows from colums 1:number of variables
Conditions(:,i) = round(Conditions(:,i)./dx(:,i))*dx(:,i); %rounds the conditions to the nearest integer value
end
Conditions = sortrows(Conditions,TempColumn);
LHCconditions = Conditions
setappdata(handles.output,'Conditions',Conditions)
%ask if constant values are logical, if nx==1 then forms a matrix of constants  
nc=zeros(NumberExperiments,5);
askforvariables(handles,NumberExperiments)
nc= getappdata(handles.output, 'nc')
%% This is where the model is called 
[Final] =  Initialisationsteadstate(handles)
% input objectives simulated by model into TSEMO
setappdata(handles.output,'Final',Final)
SaveLHC(handles)
%% format digital LHC into objective functions 
     ObjectiveFunctions = [-log(Final(:,4)) log(Final(:,7))];
     opt = TSEMO_options;                      
     opt.NoOfBachSequential = 5;
    numexp = opt.NoOfBachSequential;
    trainLims = Limits.dx.*3
  %increase the dx limits by 3 * increase spread   
    [NewConditions] = TSEMO_V2_SOFR(Conditions,ObjectiveFunctions,Limits.lb,Limits.ub,opt);
   for i= 1:numel(Limits.lb)
         NewConditions(:,i) = round(NewConditions(:,i)./Limits.dx(:,i))*Limits.dx(:,i);
   end
   comprows = unique(NewConditions, 'rows');
   OldConditions = comprows;
   
   % ask if rows are unique, if not then it will reiterate 
   if numel(comprows(:,1))< numexp
    opt.NoOfBachSequential = numexp - numel(comprows(:,1));
    [NewConditions] = TSEMO_V2_SOFR(Conditions,ObjectiveFunctions,Limits.lb,Limits.ub,opt);
   for i= 1:numel(Limits.lb)
   NewConditions(:,i) = round(NewConditions(:,i)./Limits.dx(:,i))*Limits.dx(:,i);
   end
   NewConditions = vertcat(OldConditions,NewConditions);
   comprows = unique(NewConditions, 'rows');
   else
   end    
Conditions = NewConditions;
setappdata(handles.output,'Conditions',Conditions)
NumberExperiments = numexp;
%assign experimental conditions to String 'Experiments' 
for i=1:NumberExperiments % for 1 to the number of trainging experiments defined from the GUI
   Experiments(i).Temperature = Conditions(i,TempColumn); %Experiments.Temperature is the 1st column in Conditions 
   Experiments(i).RT = Conditions(i,RTColumn);%Experiments.RT is the 2nd column in Conditions 
   Experiments(i).Ratios = Conditions(i,PumpColumn);
   setappdata(handles.output,'Experiments',Experiments)  
   ReactVol = str2double(get(handles.ReactVol,'String'));
 % Asks if Residence time exists in the Column Structure
if isrtpresent==1 %if rt exists gets the LHS generated RTs and finds the total flowrat
    Experiments(i).Flowrate = ReactVol./Experiments(i).RT;
   
elseif isrtpresent==0 %if rt doesnt exist gets the static RT from GUI and finds the total flowrate
  RTStatic = str2double(get(handles.RTStatic,'String')) ;
   Experiments(i).Flowrate = ReactVol./RTStatic;
  %setappdata(handles.output,'Experiments',Experiments)   
end
%% Calculate Flow rates
relativeFlowrate(handles,i)
Experiments = getappdata(handles.output,'Experiments');
flowratecalc(handles,i)
 Experiments = getappdata(handles.output,'Experiments');
end

SteadyVol = str2double(get(handles.SteadyVol, 'String')); %number of reactor volumes to reach steady state
Iteration = 0;
% Sets flow rate 
for i = 1:NumberExperiments
    StartDelay = Experiments(i).RT*60*SteadyVol;
    ReactTimer(i) = timer('StartDelay', StartDelay ,...
        'TimerFcn', {@gpcrun,handles,i,Experiments(i).RT},...
        'ExecutionMode', 'singleShot',...
        'ObjectVisibility', 'off',...
        'BusyMode', 'drop','Name',sprintf('Reaction%c', i),'UserData',Conditions);
end
setappdata(handles.output,'ReactTimer',ReactTimer)
setappdata(handles.output,'Experiments',Experiments)
Iteration = 0;
ChangeExperiment(handles,Iteration)

    function askforvariables(handles, NumberExperiments)
%asks the GUI what variables are selected and generates a matrix of
%conditions n x 5 Column 1 is DP; Column 2 is Initiator e;, Column 3 is
%wt%solvent dilution; Column 4 is Temperature; column 5 is RT
PumpObjects = getappdata(handles.output,'PumpObjects')
Conditions= getappdata(handles.output,'Conditions')
Columns = getappdata(handles.output,'Columns')

connected = getappdata(handles.output,'connected')
PumpO=[];
for j = 1:numel(connected)
if connected(j).type == 7
            PumpO = horzcat(PumpO,connected(j));
end
end 
   for j= 1:numel(PumpO)
      tarDP =get(handles.tarDP,'String')
 ndp=isempty(tarDP)

if ndp == 0
    tarDP=str2double(tarDP);
    newConditions.targetDP= repmat(tarDP,NumberExperiments,1);
    nc(:,1)=newConditions.targetDP;
else 
    if strcmp(PumpO(j).name,'Monomer') == 1
    nc(:,1) = Conditions(:,j);
    else 
    end 
end 
   end 
   for j= 1:numel(PumpO)
initeq =get(handles.initequ,'String');
ni=isempty(initeq);
if ni == 0
    initeq=str2double(initeq);
    newConditions.initeq= repmat(initeq,NumberExperiments,1);
    nc(:,2) = newConditions.initeq;
else 
 if strcmp(PumpO(j).name,'Initiator') == 1
    nc(:,2)=Conditions(:,1);
 else 
end 
end
   end
    for j= 1:numel(PumpO)
wtpa = get(handles.wtpa,'String')
nwt= isempty(wtpa)
if nwt==0
   wtpa=str2double(wtpa);
    newConditions.wtpa= repmat(wtpa,NumberExperiments,1);  
    nc(:,3) = newConditions.wtpa;
else
    if strcmp(PumpO(j).name,'Solvent')==1
    nc(:,3) = Conditions(:,j);
    else
    end       
end 
    end
  
     for j= 1:numel(PumpO)
tem=get(handles.SetT,'String');
nT = isempty(tem);
if nT==0
   temp=str2double(tem);
    newConditions.temp= repmat(temp,NumberExperiments,1); 
     nc(:,4) = newConditions.temp;
else
   
    nc(:,4) = Conditions(:,Columns.TempColumn);
end 
     end 
      for j= 1:numel(PumpO)
ResT = get(handles.RTStatic,'String');
nrt= isempty(ResT);
if nrt==0
   ResT=str2double(ResT);
    newConditions.ResT= repmat(ResT,NumberExperiments,1);  
    nc(:,5) = newConditions.ResT;
else
    
     nc(:,5) = Conditions(:,Columns.RTColumn);
end 
      end 
      
  
setappdata(handles.output,'nc',nc)

function relativeFlowrate(handles,i,NumberExperiments)

PumpObjects=getappdata(handles.output,'PumpObjects');% get pump objects
connected = getappdata(handles.output,'connected') %ask if connected
PumpO=[];
% if object type = 7 then this is a JASCO HPLC pump, this will find the
% number of pumps active which will enable assignment of a flowrate and
% ensure the correct flowrate calculation is conducted
for j = 1:numel(connected)
if connected(j).type == 7
            PumpO = horzcat(PumpO,connected(j));
end
end 
Conditions=getappdata(handles.output,'Conditions');
fs1 = getappdata(handles.output,'fs1'); % feed stock conc of monomer(M) from Polymer synthesis gui
fs2 = getappdata(handles.output,'fs2'); % feed stock conc of CTA (M) from Polymer synthesis gui
fs3 = getappdata(handles.output,'fs3'); % feed stock conc of Initiator(M) from Polymer synthesis gui
reswt = getappdata(handles.output,'reswt'); % wt % of reservoir from Polymer synthesis gui
tDP = getappdata(handles.output,'tDP'); % target degree of polymerisation from Polymer synthesis gui
monmass=getappdata(handles.output,'monmass'); % mass of monomer  from Polymer synthesis gui
CTAmass=getappdata(handles.output,'CTAmass'); % mass of CTA  from Polymer synthesis gui
initmass=getappdata(handles.output,'initmass');% mass of initiator from Polymer synthesis gui
fs1=str2double(fs1); % convert feedstock conc of monomer to double 
fs2=str2double(fs2); % convert feedstock conc of CTA to double 
fs3=str2double(fs3);% convert feedstock conc of initiator to double 
reswt=str2double(reswt);% convert wt %  to double 
tDP=str2double(tDP); % % convert target degree of polymerisation to double 
monmass=str2double(monmass); % convert mass of monomer to double 
CTAmass=str2double(CTAmass);% convert mass of CTA to double 
initmass=str2double(initmass)% convert mass of initiator to double 
Experiments = getappdata(handles.output,'Experiments') % get suggested experiments 
Conditions = getappdata(handles.output,'Conditions') % get the suggested conditions
initequ=get(handles.initequ,'String');
 initequ= str2double(initequ);
 %work out relative ratios
for k=1:numel(PumpObjects)
   P_name = strcat('Pumpname_',PumpObjects(k).name);
    eval(sprintf('%s = %g', P_name, k))
end 
Experiments = getappdata(handles.output,'Experiments');
% if a constant initiator equivalence does not exist 
if  isnan(initequ)==1
 %so the number associated with the pumpname_ should correspond to the column number in Experiments.Ratios
 inieq= Conditions(:,Pumpname_Initiator); 
 CTAini(i) = inieq(i)./(fs3./fs2); % ratio of LHS Initeq/ feedstock conc
else
    initequ=get(handles.initequ,'String');
     initequ= str2double(initequ); 
    initequ = repmat(initequ,numel(Experiments),1);
 CTAini(i) = initequ(i)./(fs3./fs2); %if iniitiator isnt ticked takes the static value inputed on the GUI
end 

Experiments(i).CTAIni = CTAini(i);
setappdata(handles.output,'Experiments',Experiments) 
%if Solvent pump ticked, takes the range calculated by the LHS/ suggested
%bt the algorithm 
wtpa= get(handles.wtpa,'String');
 wtpa=str2double(wtpa);
if isnan(wtpa)==1
    deswt = Conditions(:,Pumpname_Solvent);
    tmasssolute = monmass+CTAmass+initmass;
    reswt = reswt/100;
    ressol = ((100/reswt).*tmasssolute)- tmasssolute;
    deswt(i) = deswt(i)/100;
    dessol(i) = ((100/deswt(i)).*tmasssolute)- tmasssolute;
    Ratiowt(i) = ressol./dessol(i);
else 
    %if solvent pump isnt ticked it takes the static value inputted into
    %the GUI
     tmasssolute = monmass+CTAmass+initmass;
 wtpa= get(handles.wtpa,'String');
 wtpa=str2double(wtpa);
 wtpa=repmat(wtpa,i,1);
 deswt= wtpa;
 ressol = ((100/reswt).*tmasssolute)- tmasssolute;
     dessol(i) = ((100/deswt(i)).*tmasssolute)- tmasssolute;
    Ratiowt(i) = ressol./dessol(i);
end 
Experiments(i).Ratiowt=Ratiowt(i);
setappdata(handles.output,'Experiments',Experiments)
 
 tarDP= get(handles.tarDP,'String');
    tarDP =str2double(tarDP);
if isnan( tarDP)==1
    Degpol = Conditions(:,Pumpname_Monomer);
    DPrat(i) = (fs1/fs2)./Degpol(i);
else 
   
    tarDP= get(handles.tarDP,'String');
    tarDP =str2double(tarDP);
    tarDP = repmat(tarDP,numel(Experiments),1);
   
    DPrat(i) = tarDP(i)/tDP;
end 

Experiments(i).DPrat=DPrat(i);
setappdata(handles.output,'Experiments',Experiments)

function flowratecalc(handles,i)
% this function calculates flowrates depending on the inlet concentrations
% wt % etc
Experiments = getappdata(handles.output,'Experiments');
connected = getappdata(handles.output,'connected'); % get connected objects
pnum = numel(connected);
% if object type = 7 then this is a JASCO HPLC pump, this will find the
% number of pumps active which will enable assignment of a flowrate and
% ensure the correct flowrate calculation is conducted
for j=1:pnum
   if connected(j).type == 7 
       pumpnump(1,j) = horzcat(connected(j));
   end 
end 
 activepumps = numel(pumpnump)-1;% if GPC is online -1 as this pump has a constant flowrate
 ReactV = str2double(get(handles.ReactVol,'String')); % get reactor volume 

 for ii=1:numel(activepumps)
 if  activepumps == 3 && strcmp(pumpnump(ii).name,'M/CTA')==1  
     totalFRate(i)=ReactV./ Experiments(i).RT;
        Experiments(i).flowrate1= (totalFRate(i).*Experiments(i).Ratiowt);
       Experiments(i).flowrate2= Experiments(i).CTAIni.*(totalFRate(i)-Experiments(i).flowrate1);
        Experiments(i).flowrate3 = totalFRate(i)-(Experiments(i).flowrate1+Experiments(i).flowrate2);
  if Experiments(i).flowrate3 <= 0 
           Experiments(i).flowrate3 =0;
         end 
 end 
 end   

setappdata(handles.output,'pumpnump',pumpnump)
setappdata(handles.output,'Experiments',Experiments)




function ChangeExperiment(handles,Iteration)
    Experiments = getappdata(handles.output,'Experiments');
Iteration = Iteration+1;
MnData =  [0 0 0 0];
NMRData = [0 0];
try
MnData =  getappdata(handles.output,'MnData');
NMRData = getappdata(handles.output,'NMRData');
MatchData = MnData(:,1)== NMRData(:,1);
% This is where you define your objective functions to be passed to TS-EMO
% code later on
ObjectiveFunctions = [-log(NMRData(MatchData,2)) log(MnData(MatchData,4))];
GraphData = [NMRData(MatchData,2) MnData(MatchData,4)];
GraphConditions = getappdata(handles.output, 'Conditions');
UpdateGraph2(handles,GraphConditions)
UpdateGraph(handles,GraphData)
catch Exception
end
%% asks if the current iteration is greater than the number of experiments 
if Iteration>numel(Experiments)
    Columns =getappdata(handles.output,'Columns');
    Conditions = getappdata(handles.output,'Conditions');
    Limits = getappdata(handles.output,'Limits');
    opt = TSEMO_options;                      
    opt.NoOfBachSequential = str2double(get(handles.predictnumber,'String'));      
  
   % % while loop so there is definitely a new experiment (i.e. if unique
    % line cancels it out
    NewExp = 0;
    while NewExp==0
    [NewConditions] = TSEMO_V2_SOFR(Conditions,ObjectiveFunctions,Limits.lb,Limits.ub,opt);
    for i= 1:numel(Limits.lb)
    NewConditions(:,i) = round(NewConditions(:,i)./Limits.dx(:,i))*Limits.dx(:,i);
    end
    % % The below commented out for batches of 1
    %NewConditions = unique(NewConditions,'rows');
    %NewConditions = sortrows(NewConditions,Columns.TempColumn);
    
    OldConditions = Conditions;
    Conditions = vertcat(Conditions,NewConditions);
    Conditions = unique(Conditions,'stable','rows');
    NewExp = size(Conditions)-size(OldConditions);
    end
    
    % % Need this line if batch size > 1, to account for potential
    % % eliminations at second unique line, leaving NewExp>1, but fewer exps
    %NewConditions = Conditions(size(OldConditions,1)+1:size(OldConditions,1)+NewExp,:);
     setappdata(handles.output,'Conditions',Conditions)
    for i=1:NewExp(1)
       
        Experiments(end+1).RT = NewConditions(i,Columns.RTColumn);
        Experiments(end).Temperature = NewConditions(i,Columns.TempColumn)
        ReactVol = str2double(get(handles.ReactVol,'String'));
        Experiments(end).Flowrate = ReactVol./Experiments(end).RT;
       setappdata(handles.output,'Experiments',Experiments)
        setappdata(handles.output,'NewConditions',NewConditions)
         relativeFlowrate(handles,i)
        flowratecalc(handles,i)
    end 
%     set the multipump FRs
       for  j=Iteration:numel(Experiments)
       relativeFlowrate(handles,j)
       Experiments=getappdata(handles.output,'Experiments');
        flowratecalc(handles,j)
        Experiments=getappdata(handles.output,'Experiments');
       end 
     setappdata(handles.output,'Experiments',Experiments)
       
 
  SteadyVol = str2double(get(handles.SteadyVol, 'String'));
%   Experiments=getappdata(handles.output,'Experiments');
    Experiments=getappdata(handles.output,'Experiments');
 
  ReactTimer = getappdata(handles.output,'ReactTimer');
  for i = Iteration:numel(Experiments)

    StartDelay = Experiments(i).RT*60*SteadyVol;
    ReactTimer(i) = timer('StartDelay', StartDelay ,...
        'TimerFcn', {@gpcrun,handles,i,Experiments(i).RT},...
        'ExecutionMode', 'singleShot',...
        'ObjectVisibility', 'off',...
        'BusyMode', 'drop','Name',sprintf('Reaction%c', i),'UserData',Conditions);
  end
  setappdata(handles.output,'ReactTimer',ReactTimer)
  SetTemp(handles, Experiments(Iteration),Iteration)
else
 SetTemp(handles, Experiments(Iteration),Iteration)
end

set(handles.ExpNo, 'String', Iteration)
set(handles.CurrRT, 'String', Experiments(Iteration).RT)
set(handles.CurrTemp, 'String', Experiments(Iteration).Temperature)
set(handles.Flowr1, 'String',Experiments(Iteration).flowrate1)
set(handles.Flowr2, 'String',Experiments(Iteration).flowrate2)
set(handles.Flowr3, 'String',Experiments(Iteration).flowrate3)

function SetTemp(handles,Experiment,Iteration)
PolymerGUI= getappdata(handles.output,'PolymerGUI');
connected = getappdata(PolymerGUI.PolymerGUI,'connected');
TempObjects = connected(5);
SetT = str2double(get(handles.SetT, 'String'));
SetT(isnan(SetT(:,1)),:) = [];

nT = isempty(SetT);
if nT ==0
    SetT = repmat(SetT, numel(Experiment),1)
   Experiment.Temperature=SetT;
else 
end 
    
   for i =1:numel(TempObjects)
   tempobjwritesettemp(TempObjects(i).sobjects,Experiment.Temperature)
   end

CheckTemp = timer('Period', 10 ,...
        'TimerFcn', {@StableTemp,handles,TempObjects,Experiment,Iteration},...
        'ExecutionMode', 'fixedSpacing',...
        'ObjectVisibility', 'off',...
        'BusyMode', 'drop','Name','Check Temperature');
setappdata(handles.output,'CheckTemp',CheckTemp)    
start(CheckTemp)

function StableTemp(timerobj,~,handles,TempObjects,Experiment,Iteration)
    for i =1:numel(TempObjects)
        currentTemp(i) = tempobjcurrenttemp(TempObjects(i).sobjects);
    end
    if (Experiment.Temperature-0.5<currentTemp) && (currentTemp<Experiment.Temperature+0.5)
        stop(timerobj)
        delete(timerobj)
        
        WriteFlow(handles,Experiment,Iteration)
    end
    
function WriteFlow(handles,Experiment,Iteration)
    PolymerGUI= getappdata(handles.output,'PolymerGUI');
    connected = getappdata(PolymerGUI.PolymerGUI,'connected');
%     for m=1:8
%          serial = sprintf('serialObject%d', m);
% connected(m).sobjects= getappdata(PolymerGUI.MainGUI,serial);
%     end 
%   
for i=1:numel(connected)
    if connected(i).type == 7
        switch connected(i).pumptype{1} 
            case 'CTA'
                feval(str2func(connected(i).commands.(['writeFlowCommand'])), connected(i).sobjects, Experiment(1).flowrate1);
                feval(str2func(connected(i).commands.('startCommand')), connected(i).sobjects);
            case 'Initiator'
                feval(str2func(connected(i).commands.(['writeFlowCommand'])), connected(i).sobjects, Experiment(1).flowrate2 );
                feval(str2func(connected(i).commands.('startCommand')), connected(i).sobjects);
            case 'Solvent'
                feval(str2func(connected(i).commands.(['writeFlowCommand'])), connected(i).sobjects, Experiment(1).flowrate3);
                feval(str2func(connected(i).commands.('startCommand')), connected(i).sobjects);
        end
    end
end
       ReactTimer = getappdata(handles.output,'ReactTimer');
       start(ReactTimer(Iteration))
      

function SlowFlow(handles)
   PolymerGUI= getappdata(handles.output,'PolymerGUI');
    connected = getappdata(PolymerGUI.PolymerGUI,'connected');
%     for m=1:8
%          serial = sprintf('serialObject%d', m);
% connected(m).sobjects= getappdata(PolymerGUI.PolymerGUI,serial);
%     end    
     for i=1:numel(connected)
        if connected(i).type == 7
        switch connected(i).pumptype{1}
            case 'CTA'
                feval(str2func(connected(i).commands.(['writeFlowCommand'])), connected(i).sobjects, 0.05);
                feval(str2func(connected(i).commands.('startCommand')), connected(i).sobjects);
            case 'Initiator' 
            feval(str2func(connected(i).commands.(['writeFlowCommand'])), connected(i).sobjects, 0.05);
                feval(str2func(connected(i).commands.('startCommand')), connected(i).sobjects);
            case 'Solvent' 
            feval(str2func(connected(i).commands.(['writeFlowCommand'])), connected(i).sobjects, 0.11);
                feval(str2func(connected(i).commands.('startCommand')), connected(i).sobjects);
        end
        end 
    end 
  %guidata(handles.output,handles)       
  function Stopflow(handles)
PolymerGUI= getappdata(handles.output,'PolymerGUI');
connected = getappdata(PolymerGUI.PolymerGUI,'connected');
for i=1:numel(connected)
    if connected(i).type == 7
        switch connected(i).pumptype{1}
                    
                    case 'CTA'
                        feval(str2func(connected(i).commands.('stopCommand')), connected(i).sobjects);
                    case 'Initiator'
                        feval(str2func(connected(i).commands.('stopCommand')), connected(i).sobjects);
                  
                    case 'Solvent'
                        feval(str2func(connected(i).commands.('stopCommand')), connected(i).sobjects);
        end
    end
end

function gpcrun(timerobj,~,handles,iteration,RetentionTime)   
PolymerGUI = getappdata(handles.output,'PolymerGUI');
connected = getappdata(PolymerGUI.PolymerGUI,'connected');
stop(timerobj)
delete(timerobj);
%search connected objects for switching valve 
for i=1:numel(connected)
    if connected(i).type == 14
        feval(str2func('sampleloopobjsample'), connected(i).sobjects);
    end 
end
% set up picotimer 
picodata = getappdata(PolymerGUI.MainGUI, 'picoData'); %get already collected picodata
picoDataTicker = getappdata(PolymerGUI.MainGUI, 'picoDataTicker');
picoDataTicker = picoDataTicker-1;
picodata = picodata(1:picoDataTicker, :);
StartTime = size(picodata);
SlowFlow(handles)

%Start NMR Scan
NMRCalibrate(handles,iteration)
% Collect GPC data
Collect = timer('StartDelay', 350,...
    'TimerFcn', {@gpccollect, handles, RetentionTime},...
    'ExecutionMode', 'singleshot',...
    'ObjectVisibility', 'off',...
    'BusyMode', 'queue', 'Name', 'GPC Collect','UserData',StartTime);
setappdata(handles.output,'collect',Collect)
start(Collect)

function gpccollect(timerobj,~,handles,RetentionTime)
% handles    structure with handles and user data (see GUIDATA)
%% Stop Timers and Slow GPC pump flow rate 
PolymerGUI = getappdata(handles.output,'PolymerGUI');
% connected = getappdata(handles.output,'connected');
StartTime = get(timerobj,'UserData');
stop(timerobj)
delete(timerobj);
%% Collect RI response from PicoLogger
picodata = getappdata(PolymerGUI.MainGUI, 'picoData'); %get already collected picodata
picoDataTicker = getappdata(PolymerGUI.MainGUI, 'picoDataTicker');
picoDataTicker = picoDataTicker-1;
picodata = picodata(1:picoDataTicker, :);
EndTime = picoDataTicker;
gpctime = picodata(StartTime:EndTime,1);
gpctime = gpctime-gpctime(1,1);
gpctime = gpctime.*(24*3600);
gpcdata = picodata(StartTime:EndTime,2);
gpcdata = horzcat(gpctime,gpcdata);
%Store GPC Data acquired
chromatograms = getappdata(handles.output,'Chromatograms');
sx = size(chromatograms);
sy = size(gpcdata);
a = max(sx(1),sy(1));
chromatograms = [[chromatograms;zeros(abs([a 0]-sx))],[gpcdata;zeros(abs([a,0]-sy))]];
setappdata(handles.output,'Chromatograms',chromatograms)
    %processes MWdata
    [Mn ,Mp, PD] = GPCAnalysis(gpcdata);
    MnData = getappdata(handles.output,'MnData');
    if isempty(MnData)
        MnData = [RetentionTime Mn Mp PD];
    else
        NewData = [RetentionTime Mn Mp PD];
        MnData = vertcat(MnData,NewData);
    end
setappdata(handles.output,'MnData',MnData)
SaveData(handles)



function WriteNMR(NMRPort,Command)
fclose(NMRPort);
fopen(NMRPort);
fwrite(NMRPort,Command)

function NMRCalibrate(handles,iteration)
PolymerGUI = getappdata(handles.output,'PolymerGUI');
NMRPort = getappdata(PolymerGUI.PolymerGUI,'NMRPort');
NMRShimCommand = getappdata(PolymerGUI.PolymerGUI,'NMRShimCommand');
WriteNMR(NMRPort,NMRShimCommand)

StartScan = timer('StartDelay', 180,...
    'TimerFcn', {@NMRScan, handles,iteration},...
    'ExecutionMode', 'singleShot',...
    'ObjectVisibility', 'off',...
    'Name','Start NMR Scan',...
    'BusyMode', 'queue');
start(StartScan)

function NMRScan(timerobj,~,handles,iteration)
% handles    structure with handles and user data (see GUIDATA)
%Run Scan on NMR
delete(timerobj)
%Run Scan on NMR
PolymerGUI = getappdata(handles.output,'PolymerGUI');
NMRPort = getappdata(PolymerGUI.PolymerGUI,'NMRPort');
FileLocation = getappdata(PolymerGUI.PolymerGUI,'NMRSaveDirectory');
CurrentDirectory = dir(FileLocation);
NMRScanCommand = getappdata(PolymerGUI.PolymerGUI,'NMRScanCommand');
WriteNMR(NMRPort,NMRScanCommand{1})

MestReProcessing = timer('StartDelay', 70,...
    'TimerFcn', {@MestReProcess, handles, FileLocation, CurrentDirectory, iteration},...
    'ExecutionMode', 'singleShot',...
    'ObjectVisibility', 'off',...
    'Name','RunMestReProcessing',...
    'BusyMode', 'drop');
setappdata(handles.output,'MestReProcessing',MestReProcessing)
start(MestReProcessing)

function MestReProcess(~,~,handles,FileLocation,OldDirectory,iteration)
% handles    structure with handles and user data (see GUIDATA)

%% Process NMR 
CurrentDirectory = dir(FileLocation);
ExportedSpectraPath = {};
i=0;
oldfiles = numel(OldDirectory);
newfiles = numel(CurrentDirectory);
if newfiles>oldfiles
    for newFile = oldfiles+1:newfiles
         SpectraLocation = strcat(FileLocation,'/',CurrentDirectory(newFile).name); %concanate current directory strings as name 
        mkdir(SpectraLocation,'Enhanced')% make a folder in the string SpectraLocation and call it Enhanced
        NoEnhanced=strcat(SpectraLocation,'/data.1d');%concantates the spectralocation to the location of the data.1d folder
        NewEnhanced=strcat(SpectraLocation,'/Enhanced'); %concantates the spectralocation to the location of the enhanced data.1d folder
        NoEnhanced2=strcat(SpectraLocation,'/acqu.par'); %to be able to access data.1d files the acqu.par is required toenable Mnova to open the spectra
       NoEnhanced3=strcat(SpectraLocation,'/processing.script');
        copyfile(NoEnhanced,NewEnhanced);% make a copy of data.1d and put it in Enhance
        copyfile(NoEnhanced2,NewEnhanced)
        copyfile(NoEnhanced3,NewEnhanced)
        TemplateID = fopen('AutoProcess.txt');
        MestReScriptTemplate = fscanf(TemplateID,'%c');
        fclose(TemplateID);
        MestReScriptTemplate = strrep(MestReScriptTemplate,FileLocation,SpectraLocation);
        ScriptID = fopen('AutoProcess.qs','w');
        fprintf(ScriptID,'%c',MestReScriptTemplate);
        fclose(ScriptID);
        system('"C:\Program Files\Mestrelab Research S.L\MestReNova\MestReNova.exe" "C:\Users\Spinsolve\Documents\MATLAB\Badger\PolymerSynthesis\AutoProcess.qs" -sf "AutoProcess"')
        i=i+1;
        ExportedSpectraPath{1} = sprintf('%c',SpectraLocation,'/exported.txt');
    end
    
SpectraProcessTimer = timer('StartDelay', 30,...
    'TimerFcn', {@ProcessNMR, handles, ExportedSpectraPath, iteration},...
    'ExecutionMode', 'singleShot',...
    'ObjectVisibility', 'off',...
    'Name','ProcessExportedMestReSpectra',...
    'BusyMode', 'drop');
setappdata(handles.output,'SpectraProcessTimer',SpectraProcessTimer)
start(SpectraProcessTimer)  
MestReProcessing = getappdata(handles.output,'MestReProcessing');
stop(MestReProcessing)
delete(MestReProcessing)
end

function ProcessNMR(timeobj,~,handles,ExportedSpectraPath,iteration)
% handles    structure with handles and user data (see GUIDATA)
stop(timeobj)
delete(timeobj)
for i=1:numel(ExportedSpectraPath)
    Spectra = dlmread(ExportedSpectraPath{i});
    Conversion = getNMRConversion(Spectra);
    Conversion = Conversion;
end
Experiments = getappdata(handles.output,'Experiments');
SampleTime = Experiments(iteration).RT;
NMRData = [SampleTime Conversion];
OldNMRData = getappdata(handles.output,'NMRData');
NMRData = vertcat(OldNMRData,NMRData);
setappdata(handles.output,'NMRData',NMRData)
%% Next Experiment
Experiments = getappdata(handles.output,'Experiments');
if iteration == numel(Experiments)
WaitTimer = timer('StartDelay', 300,...
    'TimerFcn', {@GPCWait, handles, iteration},...
    'ExecutionMode', 'singleShot',...
    'ObjectVisibility', 'off',...
    'Name','GPC Wait',...
    'BusyMode', 'drop');
start(WaitTimer)
else
ChangeExperiment(handles,iteration)
end

function GPCWait(timerobj,~,handles,iteration)
delete(timerobj)
ChangeExperiment(handles,iteration)

function SaveData(handles)
PolymerGUI = getappdata(handles.output,'PolymerGUI');
FileLocation =getappdata(PolymerGUI.PolymerGUI,'NMRSaveDirectory');
idcs   = strfind(FileLocation,'/');
newdir = FileLocation(1:idcs(end)-1);
SaveLocation = strcat(newdir,'/ExportedParetoData.xlsx');
chromatograms = getappdata(handles.output,'Chromatograms');
MnData =  getappdata(handles.output,'MnData');
NMRData = getappdata(handles.output,'NMRData');
Conditions = getappdata(handles.output,'Conditions');
xlswrite(SaveLocation,Conditions,'Experiments');
if ~isempty(NMRData)
xlswrite(SaveLocation,NMRData,'NMRData');
end
if ~isempty(chromatograms)
xlswrite(SaveLocation,chromatograms,'GPC chromatograms');
end
if ~isempty(MnData)
    xlswrite(SaveLocation,MnData,'MnData');
end


function UpdateGraph(handles, Data)
NumberofInputs = size(Data)
set(handles.objectax,'Xlim',[0 100]);
set(handles.objectax,'Ylim',[1 2]);
displabel = horzcat('Dispersity, ', char(208));
convlabel = 'Conversion (%)';
xlabel(handles.objectax, 'Conversion');
ylabel(handles.objectax, displabel);

if NumberofInputs(2) == 2
scatter(handles.objectax, Data(:,1)*100,Data(:,2),50,'filled');
else
end
function UpdateGraph2(handles, Data)
Numberofconds = size(Data)
Limits = getappdata(handles.output, 'Limits')
condname ={}
%%get axis labels 
for i=1:9
        check = sprintf('checkbox%i',i); % collects all ticked checkboxes
         optimise = get(handles.(check),'Value'); %get check values 
      
 handles.nametag9 = handles.ResT;

      if optimise %if optimise? is checked
         
 %gets handles of min(i) where i is 1:9 and turns it into a double
      condname{i} = get(handles.(sprintf('nametag%i',i)),'String');
       
    
      %  end
      end
end
for i=1:9
    check = sprintf('checkbox%i',i); % collects all ticked checkboxes
        
    truevals(i) = get(handles.(check),'Value') 
end 
truev = find(truevals);
for i =1:numel(truev)
     names{i}= horzcat(condname{truev(i)});
end 

if Numberofconds(2) == 3
  scatter3(handles.condax, Data(:,2),Data(:,1),Data(:,3),'filled')
  zlabel(handles.condax,names{3});
  ylabel(handles.condax,names{1});
  xlabel(handles.condax,names{2});
set(handles.condax,'Xlim',[Limits.lb(2), Limits.ub(2)]);
set(handles.condax,'Ylim',[Limits.lb(1), Limits.ub(1)]);
set(handles.condax,'Zlim',[Limits.lb(3), Limits.ub(3)]);

elseif Numberofconds(2) == 2
     scatter(handles.condax, Data(:,1),Data(:,2),'filled')
  ylabel(handles.condax,names{2});
  xlabel(handles.condax,names{1});
  set(handles.condax,'Xlim',[Limits.lb(1), Limits.ub(1)]);
set(handles.condax,'Ylim',[Limits.lb(2), Limits.ub(2)]);
else
end
function SaveLHC(handles)
    %this preps the data from the model to use the experimental Save Data
    %function as LHC is run in bulk 
       Final=getappdata(handles.output,'Final'); 
       
PolymerGUI = getappdata(handles.output,'PolymerGUI');
FileLocation =getappdata(PolymerGUI.PolymerGUI,'NMRSaveDirectory');
idcs   = strfind(FileLocation,'/');
newdir = FileLocation(1:idcs(end)-1);
SaveLocation = strcat(newdir,'/ExportedParetoData.xlsx');
xlswrite(SaveLocation,Final,'SimulatedLHC')

% --- Executes on button press in SaveData.
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SaveData(handles)


function SteadyVol_Callback(hObject, eventdata, handles)
% hObject    handle to SteadyVol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SteadyVol as text
%        str2double(get(hObject,'String')) returns contents of SteadyVol as a double


% --- Executes during object creation, after setting all properties.
function SteadyVol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SteadyVol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ReactVol_Callback(hObject, eventdata, handles)
% hObject    handle to ReactVol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ReactVol as text
%        str2double(get(hObject,'String')) returns contents of ReactVol as a double


% --- Executes during object creation, after setting all properties.
function ReactVol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ReactVol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
checkbox1Status= get(hObject,'Value');
if checkbox1Status == 1
    set(handles.tarDP, 'Enable', 'off')
elseif checkbox1Status == 0
    set(handles.tarDP, 'Enable', 'on')
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
checkbox3Status= get(hObject,'Value');
if checkbox3Status == 1
    set(handles.initequ, 'Enable', 'off')
elseif checkbox3Status == 0
    set(handles.initequ, 'Enable', 'on')
end

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
checkbox4Status = get(hObject,'Value');
if checkbox4Status == 1
    set(handles.wtpa, 'Enable', 'off')
elseif checkbox4Status == 0
    set(handles.wtpa, 'Enable', 'on')
end

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8



function min1_Callback(hObject, eventdata, handles)
% hObject    handle to min1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min1 as text
%        str2double(get(hObject,'String')) returns contents of min1 as a double


% --- Executes during object creation, after setting all properties.
function min1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min2_Callback(hObject, eventdata, handles)
% hObject    handle to min2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min2 as text
%        str2double(get(hObject,'String')) returns contents of min2 as a double


% --- Executes during object creation, after setting all properties.
function min2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min3_Callback(hObject, eventdata, handles)
% hObject    handle to min3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min3 as text
%        str2double(get(hObject,'String')) returns contents of min3 as a double


% --- Executes during object creation, after setting all properties.
function min3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min4_Callback(hObject, eventdata, handles)
% hObject    handle to min4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min4 as text
%        str2double(get(hObject,'String')) returns contents of min4 as a double


% --- Executes during object creation, after setting all properties.
function min4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min5_Callback(hObject, eventdata, handles)
% hObject    handle to min5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min5 as text
%        str2double(get(hObject,'String')) returns contents of min5 as a double


% --- Executes during object creation, after setting all properties.
function min5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min6_Callback(hObject, eventdata, handles)
% hObject    handle to min6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min6 as text
%        str2double(get(hObject,'String')) returns contents of min6 as a double


% --- Executes during object creation, after setting all properties.
function min6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min7_Callback(hObject, eventdata, handles)
% hObject    handle to min7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min7 as text
%        str2double(get(hObject,'String')) returns contents of min7 as a double


% --- Executes during object creation, after setting all properties.
function min7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function min8_Callback(hObject, eventdata, handles)
% hObject    handle to min8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min8 as text
%        str2double(get(hObject,'String')) returns contents of min8 as a double


% --- Executes during object creation, after setting all properties.
function min8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max1_Callback(hObject, eventdata, handles)
% hObject    handle to max1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max1 as text
%        str2double(get(hObject,'String')) returns contents of max1 as a double


% --- Executes during object creation, after setting all properties.
function max1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max2_Callback(hObject, eventdata, handles)
% hObject    handle to max2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max2 as text
%        str2double(get(hObject,'String')) returns contents of max2 as a double


% --- Executes during object creation, after setting all properties.
function max2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max3_Callback(hObject, eventdata, handles)
% hObject    handle to max3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max3 as text
%        str2double(get(hObject,'String')) returns contents of max3 as a double


% --- Executes during object creation, after setting all properties.
function max3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max4_Callback(hObject, eventdata, handles)
% hObject    handle to max4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max4 as text
%        str2double(get(hObject,'String')) returns contents of max4 as a double


% --- Executes during object creation, after setting all properties.
function max4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max5_Callback(hObject, eventdata, handles)
% hObject    handle to max5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max5 as text
%        str2double(get(hObject,'String')) returns contents of max5 as a double


% --- Executes during object creation, after setting all properties.
function max5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max6_Callback(hObject, eventdata, handles)
% hObject    handle to max6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max6 as text
%        str2double(get(hObject,'String')) returns contents of max6 as a double


% --- Executes during object creation, after setting all properties.
function max6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max7_Callback(hObject, eventdata, handles)
% hObject    handle to max7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max7 as text
%        str2double(get(hObject,'String')) returns contents of max7 as a double


% --- Executes during object creation, after setting all properties.
function max7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max8_Callback(hObject, eventdata, handles)
% hObject    handle to max8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max8 as text
%        str2double(get(hObject,'String')) returns contents of max8 as a double


% --- Executes during object creation, after setting all properties.
function max8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx1_Callback(hObject, eventdata, handles)
% hObject    handle to dx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx1 as text
%        str2double(get(hObject,'String')) returns contents of dx1 as a double


% --- Executes during object creation, after setting all properties.
function dx1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx2_Callback(hObject, eventdata, handles)
% hObject    handle to dx2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx2 as text
%        str2double(get(hObject,'String')) returns contents of dx2 as a double


% --- Executes during object creation, after setting all properties.
function dx2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx3_Callback(hObject, eventdata, handles)
% hObject    handle to dx3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx3 as text
%        str2double(get(hObject,'String')) returns contents of dx3 as a double


% --- Executes during object creation, after setting all properties.
function dx3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx4_Callback(hObject, eventdata, handles)
% hObject    handle to dx4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx4 as text
%        str2double(get(hObject,'String')) returns contents of dx4 as a double


% --- Executes during object creation, after setting all properties.
function dx4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx5_Callback(hObject, eventdata, handles)
% hObject    handle to dx5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx5 as text
%        str2double(get(hObject,'String')) returns contents of dx5 as a double


% --- Executes during object creation, after setting all properties.
function dx5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx6_Callback(hObject, eventdata, handles)
% hObject    handle to dx6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx6 as text
%        str2double(get(hObject,'String')) returns contents of dx6 as a double


% --- Executes during object creation, after setting all properties.
function dx6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx7_Callback(hObject, eventdata, handles)
% hObject    handle to dx7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx7 as text
%        str2double(get(hObject,'String')) returns contents of dx7 as a double


% --- Executes during object creation, after setting all properties.
function dx7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx8_Callback(hObject, eventdata, handles)
% hObject    handle to dx8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx8 as text
%        str2double(get(hObject,'String')) returns contents of dx8 as a double


% --- Executes during object creation, after setting all properties.
function dx8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function predictnumber_Callback(hObject, eventdata, handles)
% hObject    handle to predictnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of predictnumber as text
%        str2double(get(hObject,'String')) returns contents of predictnumber as a double


% --- Executes during object creation, after setting all properties.
function predictnumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to predictnumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function training_Callback(hObject, eventdata, handles)
% hObject    handle to training (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of training as text
%        str2double(get(hObject,'String')) returns contents of training as a double


% --- Executes during object creation, after setting all properties.
function training_CreateFcn(hObject, eventdata, handles)
% hObject    handle to training (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ExpNo_Callback(hObject, eventdata, handles)
% hObject    handle to ExpNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExpNo as text
%        str2double(get(hObject,'String')) returns contents of ExpNo as a double


% --- Executes during object creation, after setting all properties.
function ExpNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExpNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrRT_Callback(hObject, eventdata, handles)
% hObject    handle to CurrRT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrRT as text
%        str2double(get(hObject,'String')) returns contents of CurrRT as a double


% --- Executes during object creation, after setting all properties.
function CurrRT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrRT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function CurrTemp_Callback(hObject, eventdata, handles)
% hObject    handle to CurrTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CurrTemp as text
%        str2double(get(hObject,'String')) returns contents of CurrTemp as a double


% --- Executes during object creation, after setting all properties.
function CurrTemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
% hObject    handle to Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Start,'Enable','on');
ReactTimer = getappdata(handles.output,'ReactTimer');
delete(ReactTimer)
SaveData(handles)
MnData = [];
Chromatograms = [];
NMRData = [];
setappdata(handles.output,'NMRData',NMRData)
setappdata(handles.output,'MnData',MnData)
setappdata(handles.output,'Chromatograms',Chromatograms)
%Stopflow(handles)
delete(timerfindall)


function Flowr1_Callback(hObject, eventdata, handles)
% hObject    handle to Flowr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Flowr1 as text
%        str2double(get(hObject,'String')) returns contents of Flowr1 as a double


% --- Executes during object creation, after setting all properties.
function Flowr1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Flowr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Flowr2_Callback(hObject, eventdata, handles)
% hObject    handle to Flowr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Flowr2 as text
%        str2double(get(hObject,'String')) returns contents of Flowr2 as a double


% --- Executes during object creation, after setting all properties.
function Flowr2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Flowr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Flowr3_Callback(hObject, eventdata, handles)
% hObject    handle to Flowr3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Flowr3 as text
%        str2double(get(hObject,'String')) returns contents of Flowr3 as a double


% --- Executes during object creation, after setting all properties.
function Flowr3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Flowr3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in monconc.
function monconc_Callback(hObject, eventdata, handles)
% hObject    handle to monconc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of monconc
newValue = get(hObject,'Value');
if newValue
    % defines the fields
    enableProperty = 'on';
    status = 'Ready to connect';
else
    enableProperty = 'off';
    status = 'Disabled';
end

set([   handles.nametag1;...
        handles.min1;...
        handles.max1;... 
        handles.dx1], 'Enable', enableProperty)


% --- Executes on button press in ctaconc.
function ctaconc_Callback(hObject, eventdata, handles)
% hObject    handle to ctaconc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ctaconc


% --- Executes on button press in iniconc.
function iniconc_Callback(hObject, eventdata, handles)
% hObject    handle to iniconc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iniconc


% --- Executes on button press in RTcb.
function RTcb_Callback(hObject, eventdata, handles)
% hObject    handle to RTcb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RTcb


% --- Executes on button press in tempcb.
function tempcb_Callback(hObject, eventdata, handles)
% hObject    handle to tempcb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tempcb


% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox22



function min9_Callback(hObject, eventdata, handles)
% hObject    handle to min9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min9 as text
%        str2double(get(hObject,'String')) returns contents of min9 as a double


% --- Executes during object creation, after setting all properties.
function min9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max9_Callback(hObject, eventdata, handles)
% hObject    handle to max9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max9 as text
%        str2double(get(hObject,'String')) returns contents of max9 as a double


% --- Executes during object creation, after setting all properties.
function max9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dx9_Callback(hObject, eventdata, handles)
% hObject    handle to dx9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dx9 as text
%        str2double(get(hObject,'String')) returns contents of dx9 as a double


% --- Executes during object creation, after setting all properties.
function dx9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9
checkboxStatus = get(hObject,'Value');
if checkboxStatus == 1
    set(handles.RTStatic, 'Enable', 'off')
elseif checkboxStatus == 0
    set(handles.RTStatic, 'Enable', 'on')
end


function RTStatic_Callback(hObject, eventdata, handles)
% hObject    handle to RTStatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RTStatic as text
%        str2double(get(hObject,'String')) returns contents of RTStatic as a double


% --- Executes during object creation, after setting all properties.
function RTStatic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RTStatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wtpa_Callback(hObject, eventdata, handles)
% hObject    handle to wtpa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wtpa as text
%        str2double(get(hObject,'String')) returns contents of wtpa as a double

% --- Executes during object creation, after setting all properties.
function wtpa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wtpa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tarDP_Callback(hObject, eventdata, handles)
% hObject    handle to tarDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tarDP as text
%        str2double(get(hObject,'String')) returns contents of tarDP as a double

% --- Executes during object creation, after setting all properties.
function tarDP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tarDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function initequ_Callback(hObject, eventdata, handles)
% hObject    handle to initequ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initequ as text
%        str2double(get(hObject,'String')) returns contents of initequ as a double


% --- Executes during object creation, after setting all properties.
function initequ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initequ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%function to get data from polymer synthesis GUI



function SetT_Callback(hObject, eventdata, handles)
% hObject    handle to SetT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SetT as text
%        str2double(get(hObject,'String')) returns contents of SetT as a double


% --- Executes during object creation, after setting all properties.
function SetT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SetT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in montype.
function montype_Callback(hObject, eventdata, handles)
% hObject    handle to montype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns montype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from montype
 montype=get(hObject,'Value')

% --- Executes during object creation, after setting all properties.
function montype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to montype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ctatype.
function ctatype_Callback(hObject, eventdata, handles)
% hObject    handle to ctatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ctatype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ctatype

 ctatype=get(hObject,'Value')

% --- Executes during object creation, after setting all properties.
function ctatype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ctatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in initype.
function initype_Callback(hObject, eventdata, handles)
% hObject    handle to initype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns initype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from initype

 initype=get(hObject,'Value')

% --- Executes during object creation, after setting all properties.
function initype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in soltype.
function soltype_Callback(hObject, eventdata, handles)
% hObject    handle to soltype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns soltype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from soltype

soltype=get(hObject,'Value')

% --- Executes during object creation, after setting all properties.
function soltype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to soltype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function objectax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to objectax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate objectax


% --- Executes during object creation, after setting all properties.
function condax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate condax
