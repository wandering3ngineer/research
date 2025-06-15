clear;
clc;
clf;
cla;

%% ------------------------------------------------------------------------
% Get the input data
%--------------------------------------------------------------------------
%samplefilename='data/sample-minus1-hyst';
%samplefilename='data/sample-0-hyst';
%samplefilename='data/sample-1-hyst';
%samplefilename='data/sample-2-hyst';
%samplefilename='data/sample-4-hyst';
%samplefilename='data/sample-7-hyst';
%samplefilename='data/sample-7.1-hyst';
%samplefilename='data/sample-7.2-hyst';
%samplefilename='data/sample-7.3-hyst';
%samplefilename='data/sample-7.4-hyst';
%samplefilename='data/sample-9-hyst';
%samplefilename='data/sample-9.5-hyst';
%samplefilename='data/sample-10-hyst';
%samplefilename='data/sample-10.1-hyst';
%samplefilename='data/sample-10.2-hyst';
%samplefilename='data/sample-10.3-hyst';
%samplefilename='data/sample-11.1-hyst';
%samplefilename='data/sample-11.2-hyst';
%samplefilename='data/sample-11.3-hyst';
%samplefilename='data/sample-11.4-hyst';
%samplefilename='data/sample-12.1-hyst';
%samplefilename='data/sample-12.2-hyst';
%samplefilename='data/sample-12.3-hyst';
%samplefilename='data/sample-13.1-hyst';
%samplefilename='data/sample-13.2-hyst';
%samplefilename='data/sample-13.3-hyst';
%samplefilename='data/sample-13.4-hyst';
samplefilename='data/sample-17.2-hyst';
%samplefilename='data/sample-18.3-hyst';
%samplefilename='data/sample-18.4-hyst';
videofilename=strcat(samplefilename,'-anim.avi'); 
fighystfilename=strcat(samplefilename,'-hyst.fig');
figoptfilename=strcat(samplefilename,'-opt.fig');
outputfilename=strcat(samplefilename,'-out.csv');
samplefilename=strcat(samplefilename,'-in.csv');
RAWDATA = csvread(samplefilename,1,0);

% Constant or measured data inputs
T_=RAWDATA(:,1);    % Temperature [K]
H_=RAWDATA(:,2);    % Applied field [Oe]
M_=RAWDATA(:,3);    % Total Magnetization [emu/g]
ma_=RAWDATA(:,4);   % Mass [mg]
pd_=RAWDATA(:,5);   % Material Density  [g/cm3]
A_=RAWDATA(:,6);    % Lattice parameter [nm]
Ms_=RAWDATA(:,7);   % Saturation Magnetization [emu/g]
Mr_=RAWDATA(:,8);   % Remanent Magnetization [emu/g]
Hc_=RAWDATA(:,9);   % Coercivity [Oe]
Hs_=RAWDATA(:,10);   % Applied field at saturation [Oe]
Hex_=RAWDATA(:,11); % Exchange bias field [Oe]
u0_=RAWDATA(:,12);  % Permeability of free space []
kb_=RAWDATA(:,13);  % The boltzmann constant [ergs/K]
Na_=RAWDATA(:,14);  % Avagadro's Const. [/mol]
Mw_=RAWDATA(:,15);  % Molecular weight [(g/mol)/f.u.]
mb_=RAWDATA(:,16);  % Magnetic moment of Bohr magnetron [erg/G]
eT_=RAWDATA(:,17);  % Temperature Error [K]
eH_=RAWDATA(:,18);  % Applied field Error [Oe]
eM_=RAWDATA(:,19);  % Total Magnetization Error [emu/g]
ema_=RAWDATA(:,20); % Mass Error [mg]
epd_=RAWDATA(:,21); % Material Density Error [g/cm3]
eA_=RAWDATA(:,22);  % Lattice parameter Error [g/cm3]
eMs_=RAWDATA(:,23); % Saturation Magnetization Error [emu/g]
eMr_=RAWDATA(:,24); % Remanent Magnetization Error [emu/g]
eHc_=RAWDATA(:,25); % Coercivity Error [Oe]
eHs_=RAWDATA(:,26); % Applied field at Saturation Error [Oe]
eHex_=RAWDATA(:,27);% Exchange bias field Error [Oe]
eu0_=RAWDATA(:,28); % The error on the permeability of free space []
ekb_=RAWDATA(:,29); % The boltzmann constant Error [ergs/k]
eNa_=RAWDATA(:,30); % Avagadro's Const. error [/mol]
eMw_=RAWDATA(:,31); % Molecular weight error [g/mol]
emb_=RAWDATA(:,32); % Magnetic moment of Bohr magnetron error [erg/G]

%% ------------------------------------------------------------------------
% Convert the data to SI units
%--------------------------------------------------------------------------
T=T_;                         % Temperature [K]
H=(10^3)/(4*pi()).*H_;        % Applied Field [A/m]
M=(10^3).*pd_.*M_;            % Total Magnetization [A/m]
ma=(1e-6).*ma_;               % Mass [Kg]
pd=(10^3).*pd_;               % Material Density  [Kg/m3]
A=(10^-9).*A_;                % Lattice Parameter [m]
Ms=(10^3).*pd_.*Ms_;          % Saturation Magnetization [A/m]
Mr=(10^3).*pd_.*Mr_;          % Remanent Magnetization [A/m]
Hc=((10^3)/(4*pi())).*Hc_;    % Coercivitiy [A/m]
Hs=((10^3)/(4*pi())).*Hs_;    % Applied field at saturation [A/m]
Hex=((10^3)/(4*pi())).*Hex_;  % Exchange bias field [A/m]
u0=(4*pi()*(10^-7)).*u0_;     % Permeability of free space [T.m/A]
kb=(10^-7).*kb_;              % The boltzmann constant [J/K]
Na=Na_;                       % Avagadro's const [/mol]
Mw=(10^-3).*Mw_;              % Molecular weight [(Kg/mol)/f.u.]
mb=(10^-3).*mb_;              % Magnetic Moment of bohr magnetron [J/T]
eT=eT_;                       % Temperature Error [K]
eH=((10^3)/(4*pi())).*eH_;    % Applied field Error [A/m]
eM=(10^3).*abs(pd_.*M_).*hypot((epd_./pd_),(eM_./M_));       % Magnetization Field Error [A/m]
ema=(1e-6).*ema_;             % Mass Error [Kg]
epd=(10^3).*epd_;             % Material Density Error [Kg/m3]
eA=(10^-9).*eA_;              % Lattice Parameter error [m]
eMs=(10^3).*abs(pd_.*Ms_).*hypot((epd_./pd_),(eMs_./Ms_));   % Saturation Magnetization Error [A/m]
eMr=(10^3).*abs(pd_.*Mr_).*hypot((epd_./pd_),(eMr_./Mr_));   % Remanent Magnetization Error [A/m]
eHc=((10^3)/(4*pi())).*eHc_;  % Coercivity Error [A/m]
eHs=((10^3)/(4*pi())).*eHs_;  % Applied field at Saturation Error [A/m]
eHex=((10^3)/(4*pi())).*eHex_;% Exchange bias field error [A/m]
eu0=(4*pi()*(10^(-7))).*eu0_; % The error on the permeability of free space [T.m/A]
ekb=(10^(-7)).*ekb_;          % The boltzmann constant Error [J/K]
eNa=eNa_;                     % Avagadro's const [/mol]
eMw=(10^-3).*eMw_;            % Molecular weight error [(Kg/mol)/f.u.]
emb=(10^-3).*emb_;            % Magnetic moment of bohr magnetron error [J/T]

%% ------------------------------------------------------------------------
% Produce an initial guess for our parameters {n,m,al,Xaf} and set an
% initial value for the learning rate (or tuning parameter / step size)
%--------------------------------------------------------------------------
n0=(Na.*pd)./Mw;              % Num. Dipole moments in a formula unit [f.u. (dipoles./m3)]
m0=Ms./n0;                    % Magnetic moment [(J/T) /f.u.]
al0=100;                      % Weiss parameter []
c0=0.5;                       % Coefficient of proportionality []
k0=Hc;                        % Coercivity coefficienty [A/m]
Xaf0=pd*(1.44e-6);            % Antiferromagnetic Susceptibility []

% Construct parameter guess. 
% Note: n  => decreasing value lowers the position of loop tips vertically
%       m  => decreasing value makes loop less sharp and more curved &
%             decreasing value lowers the position of loop tips vertically
%       a  => decreasing value makes loop less sharp and more curved
%       Ms => decreasing value lowers the position of loop tips vertically
%       al => decreasing value narrows hysteresis loop and increases
%             softness of curvature of the loop knee
%       Xaf=> increasing value increases the slope angle at saturation 
%       c  => increasing value sharpens the bend (knee) of the curve as it 
%             slopes down from saturation
%       k  => increasing value widens the hysteresis curve

% Some adjustments need to be made to the initial guess so that Man falls
% in the right ball park 
debug=0;
%n0=n0/1400; m0=m0*1600; al0=al0/1400; Xaf0=Xaf0; c0=c0; k0=k0;        % for sample -1
%n0=n0/1400; m0=m0*1600; al0=al0/1400; Xaf0=Xaf0; c0=c0; k0=k0;        % for sample 0
%n0=n0*0.00013; m0=m0*8000; al0=al0*0.00045; Xaf0=Xaf0; c0=c0; k0=k0;  % for sample 1
n0=n0*0.00013; m0=m0*8000; al0=al0*0.0001; Xaf0=Xaf0; c0=c0; k0=k0;    % for sample 2,4,7,7.1,7.2,7.3,7.4,9,9.5,10,10.1,10.2,10.3,11.1,11.2,11.3,11.4
P0=[n0(1),m0(1),al0(1),Xaf0(1),c0(1),k0(1)];

% Construct some meta data for our parameters
PDESC={'Num. dipole moments per unit Vol.', 'n'  ,'[1/m^3]' ; ...
       'Magnetic moment'                  , 'm'  ,'[J/T]'   ; ...
       'Weiss mean field coefficient'     , 'al' ,'[]'      ; ...
       'Antiferromagnetic susceptibility' , 'Xaf','[]'      ; ...
       'Coefficient of proportionality'   , 'c'  ,'[]'      ; ...
       'Coefficient of coercivity'        , 'k'  ,'[A/m]'};

% Define the boundaries within which the above design variables
% (parameters) must absolutely fall. Use an empty LB or UB array if 
% no boundaries
nlb=1e22;  nub=1e29;
mlb=1e-29; mub=1e-19;
allb=0;    alub=1;
Xaflb=0;   Xafub=1e-1;
clb=0;     cub=1;
klb=0;     kub=12000;
PLB=[nlb; mlb; allb; Xaflb; clb; klb];
PUB=[nub; mub; alub; Xafub; cub; kub];    
    
% Construct the data variable and associated error variable
D=[H, M, T, u0, kb];
eD=[eH, eM, eT, eu0, ekb];

%% ------------------------------------------------------------------------
% Generate an artificial hysteresis loop 
% The hysteresis curve can be divided into 5 sections, 
%  Initial, Dec +M, Dec -M, Inc -M, Inc +M
%--------------------------------------------------------------------------
Nd=250;  secs=5;
Hts=max(H);
Htinit=linspace(0,Hts,round(Nd/secs))'; % Must always begin at zero
Htdec=linspace(Hts,-Hts,round(2*Nd/secs))'; Htdec=Htdec(2:end);
Htinc=linspace(-Hts,Hts,round(2*Nd/secs))'; Htinc=Htinc(2:end);
Ht=[Htinit; Htdec; Htinc];
% The number of elements will change to deal with rounding issues and to 
% remove duplicate points when transitioning from initial to decreasing to
% increasing sections of the hysteresis curve.
Nd=numel(Ht);                           

%% ------------------------------------------------------------------------
% Setup video writer
%--------------------------------------------------------------------------
video=videofilename;
videofileObject = VideoWriter(video);
videofileObject.FrameRate=10;
videofileObject.Quality=50;
VID=videofileObject;

%% ------------------------------------------------------------------------
% Create the constraints and objective function
%--------------------------------------------------------------------------
objfun=@(D_,P_) hysteretic(Ht,D_,P_,PLB,PUB,VID);
confun=@(P_) hystereticConstraints(P_);

% Number of iterations for optimization, number of expected outputs and the
% optimization algorithm to use.
ITER=40; ALG='levenberg-marquardt';

Nm=15;  % Number of montecarlo datasets to generate for estimating error
Nmd=40; % Number of data points in each dataset.

%% ------------------------------------------------------------------------
% Open the video file for animation.
%--------------------------------------------------------------------------
open(videofileObject);

%% ------------------------------------------------------------------------
% Use this debug mode to manually adjust the parameters to the right range
% by repeatedly runnin this program
%--------------------------------------------------------------------------
if (debug==1)
    % Compute the magnetization for the full hysteresis curve considering
    % so that parameters can be manually adjusted and their effect seen
    [Mt,~,O]=objfun(D,P0);
    O(1)
    figure(3); hold on;
    plot(H,M,'r*');
    plot(Ht,Mt,'b-');
    
       
    % Do the same as above but this time for a minor hysteresis loop.
    % To do this, we need to redefine our Ht as going up to a smaller value
    % and then recompute our Mt value.
    
%     % Generate a set of H values for the minor hystersis loop
%     Nd1=250;  secs1=5;
%     Hminor=100;                   % max applied field in orsteds [Oe]
%     Hts1=(10^3)/(4*pi()).*(100);  % Max value of H-field 100 [Oe]
%     Htinit1=linspace(0,Hts1,round(Nd1/secs1))'; % Must always begin at zero
%     Htdec1=linspace(Hts1,-Hts1,round(2*Nd1/secs1))'; Htdec1=Htdec1(2:end);
%     Htinc1=linspace(-Hts1,Hts1,round(2*Nd1/secs1))'; Htinc1=Htinc1(2:end);
%     Ht1=[Htinit1; Htdec1; Htinc1];
%     % The number of elements will change to deal with rounding issues and to 
%     % remove duplicate points when transitioning from initial to decreasing to
%     % increasing sections of the hysteresis curve.
%     Nd1=numel(Ht1); 
%     
%     % Recreate the objective function with this new set of H values
%     objfun1=@(D_,P_) hysteretic(Ht1,D_,P_,PLB,PUB,VID);
%     
%     % Compute this magnetization for the minor hysteresis loop
%     [Mt1,~,O1]=objfun1(D,P0);
%     O1(1)
%     figure(4); hold on;
%     plot(H,M,'r*');
%     plot(Ht1,Mt1,'b-');
%     hold off;
    break;
end

for j=1:Nm     
    %% --------------------------------------------------------------------
    % Generate a monte carlo dataset over which to run the algorithm. 
    % Compute the parameters n,m,al,Xaf,c,k using this new data set. 
    % Repeating the optimization computation of our parameters repeatedly 
    % for randomly regenerated data sets will yield a distribution with a 
    % mean and standard deviation - the standard deviation will become 
    % the error. 
    %----------------------------------------------------------------------
    [Dm_] = montecarloError(D, eD, Nmd);
    Dm(j,:,:)=Dm_;
    
    % ---------------------------------------------------------------------
    % Do parameter estimation
    %----------------------------------------------------------------------
    [P, S, O]=paramEstimate(objfun, confun, Dm_, P0, ITER, ALG);

    % Grab the number of iterations actually completed and the dense
    % magnetization output Mt, for dense input Ht - so we can graph our
    % hysteresis curves and grab our optimal parameter values
    I=numel(S); [Mt,~,~]=objfun(Dm_,P(I,:));

    %% --------------------------------------------------------------------
    % Collect the optimal values for each montecarlo data set
    %----------------------------------------------------------------------
    Pm(j,:,:)=P;
    Sm(j,:)=S;
    Om(j,:,:)=O;
    Mtm(j,:)=Mt;       
end

%% ------------------------------------------------------------------------
% Close the video file and save it.
%--------------------------------------------------------------------------
close(videofileObject);

%% ------------------------------------------------------------------------
% Compute mean output values
%--------------------------------------------------------------------------
Pavg=squeeze(mean(Pm,1));
Savg=squeeze(mean(Sm,1));
Oavg=squeeze(mean(Om,1));
Mtavg=squeeze(mean(Mtm,1)');

% Assign our average as the actual output
P=Pavg;
S=Savg;
O=Oavg;
Mt=Mtavg;

% Extract the individual variables from the optimized parameter array.
n=P(I,1); m=P(I,2); al=P(I,3); Xaf=P(I,4); c=P(I,5); k=P(I,6);
W=O(I,1); Mdec_r=O(I,2);  Minc_r=O(I,3);  Mexr=O(I,4);  Mr_sym=O(I,5);  Msp=O(I,6);  Msn=O(I,7);  Ms_avg=O(I,8); Ms_nm=O(I,9); 
          Hdec_c=O(I,10); Hinc_c=O(I,11); Hexc=O(I,12); Hc_sym=O(I,13); Hsp=O(I,14); Hsn=O(I,15); Hs_avg=O(I,16);

% Compute the error on the parameters {n,m,al,Xaf,c,k}
eP=squeeze(std(Pm,1,1));
en=eP(I,1); em=eP(I,2); eal=eP(I,3); eXaf=eP(I,4); ec=eP(I,5); ek=eP(I,6);

% Compute the error on the output and SSD
eS=squeeze(std(Sm,1,1)); 
eO=squeeze(std(Om,1,1)); 
eS=eS(I);
eW=eO(I,1); eMdec_r=eO(I,2);  eMinc_r=eO(I,3);  eMexr=eO(I,4);  eMr_sym=eO(I,5);  eMsp=eO(I,6);  eMsn=eO(I,7);  eMs_avg=eO(I,8); eMs_nm=eO(I,9);
            eHdec_c=eO(I,10); eHinc_c=eO(I,11); eHexc=eO(I,12); eHc_sym=eO(I,13); eHsp=eO(I,14); eHsn=eO(I,15); eHs_avg=eO(I,16);

% Compute the error on the output magnetization
eMt=squeeze(std(Mtm,1,1))';

%% ------------------------------------------------------------------------
% Compute output ERRORS (standard deviation) and some other useful values
%--------------------------------------------------------------------------
nfu=(((A.^3).*pd)./Mw).*Na;         % Num Formula units in unit cell [f.u./unitcell]
g=m./mb;                            % G-factor bohr magnetrons/f.u. [1/f.u.]
v=ma./pd;                           % Volume of material used in measurement [m^3]
Wm=W/pd(1);                         % Power loss per unit mass [W/(Kg.Hz)]

% Error on Num. Formula units [f.u.]
enuf1=(3.*abs(A.^3).*(eA./abs(A))).^2;
enuf2=(epd./abs(pd)).^2+(eNa./abs(Na)).^2;
enuf3=(eMw./abs(Mw)).^2;
enuf=abs(nfu).*sqrt(enuf1+enuf2+enuf3);

% Error on G-Factor - bohr magnetrons/f.u. [1/f.u.]
eg=abs(g).*sqrt((em./abs(m)).^2+(emb./abs(mb)).^2);                           

% Error on volume estimate [m^3]
ev=abs(v).*sqrt((ema./abs(ma)).^2 + (epd./abs(pd)).^2);                           

% Error on the Power loss per unit mass 
eWm = abs(Wm)*sqrt((epd(1)/abs(pd(1)))^2 + (eW/abs(W))^2);

% Susceptibilities
X=M./H;
Xt=Mt./Ht;
Xmax=max(abs(X));
Xtmax=max(abs(Xt));
[val ind]=min(abs(H));
Xinit=X(ind);
Xtinit=Xt(2);

%% ------------------------------------------------------------------------
% Compute Minor loop Power Loss
%--------------------------------------------------------------------------
% Generate minor hysteresis loop
% To do this, we need to redefine our Ht as going up to a smaller value
% and then recompute our Mt value.

% The max applied field from the remanent field produced via 
% the core material is given as Hr=Mr/Xi where Xi is the initial 
% susceptibility and Mr is the remeanent magnetization of the core material

% We can thus compute the initial permeability from 
% our parameterized magnetization curve
HtsM=Mr_sym/Xinit;

% Generate a set of H values for the minor hystersis loop
NdM=250;  secsM=5;
HtinitM=linspace(0,HtsM,round(NdM/secsM))'; % Must always begin at zero
HtdecM=linspace(HtsM,-HtsM,round(2*NdM/secsM))'; HtdecM=HtdecM(2:end);
HtincM=linspace(-HtsM,HtsM,round(2*NdM/secsM))'; HtincM=HtincM(2:end);
HtM=[HtinitM; HtdecM; HtincM];

% The number of elements will change to deal with rounding issues and to 
% remove duplicate points when transitioning from initial to decreasing to
% increasing sections of the hysteresis curve.
NdM=numel(HtM);

% Recreate the objective function with this new set of H values
objfun1=@(D_,P_) hysteretic(HtM,D_,P_,PLB,PUB,VID);

% Compute this magnetization for the minor hysteresis loop
[MtM,~,OM]=objfun1(D,P0);

WM=OM(1);
eWM=eW;
WMm=OM(1)/pd(1);
eWMm=eWm;

%% ------------------------------------------------------------------------
% Generate csv file containing output results
%--------------------------------------------------------------------------
fid = fopen(outputfilename,'w');

% Optimization Algorithm details
fprintf(fid,'OPTIMIZATION VARIABLE OUTPUTS\n');
fprintf(fid,'%s Optimization,%s\n',ALG,'Monte-Carlo Error Estimation');
fprintf(fid,'Description, Variable, Value (Mean), Unit\n');
fprintf(fid,'%s,%s,%0.5g,%s\n','Iterations','I',I,'[]');
fprintf(fid,'%s,%s,%0.5g,%s\n','Monte Carlo num. datasets'       ,'Nm' ,Nm, '[]');
fprintf(fid,'%s,%s,%0.5g,%s\n','Monte Carlo dataset size'        ,'Nmd',Nmd,'[]');
fprintf(fid,'%s,%s,%0.5g,%s\n','Densified (Ht Mt) dataset size'  ,'Nd' ,Nd, '[]');

% Write optimized Parameter values
fprintf(fid,'\n');
fprintf(fid,'OPTIMIZED PARAMETERS\n');
fprintf(fid,'Parameter, Variable, Value (Mean), Error (Std), Unit, Init Guess, Opt Lower Bound, Opt Upper Bound\n');
for i=1:6
    fprintf(fid,'%s,%s,%0.5g,%0.5g,%s,%0.5g,%0.5g,%0.5g\n',PDESC{i,1},PDESC{i,2},P(I,i),eP(I,i),PDESC{i,3},P0(i),PLB(i),PUB(i));
end
 
% Write down other optimal output values
fprintf(fid,'\n');
fprintf(fid,'OTHER OUTPUTS (AFTER OPTIMIZATION)\n');
fprintf(fid,'Description, Variable, Value (Mean), Error (Std), Value (Graphical), Error (Graphical), Unit\n');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Sum of squared difference'              ,'S'      ,S(I)   ,eS                       ,'[(A/m)^2]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Sample Mass'                            ,'ma'     ,ma(1)  ,ema(1)                   ,'[Kg]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Vol. of sample used in SQUID'           ,'v'      ,v(1)   ,ev(1)                    ,'[m^3]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Sample Density'                         ,'pd'     ,pd(1)  ,epd(1)                   ,'[Kg/m^3]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Lattice parameter (per unit cell)'      ,'A'      ,A(1)   ,eA(1)                    ,'[m/unit cell]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Num. of formula units in unit cell'     ,'nfu'    ,nfu(1) ,enuf(1)                  ,'[f.u./unit cell]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Num. bohr magnetrons per f.u.'          ,'g'      ,g(1)   ,eg(1)                    ,'[bohr.mag/f.u.]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Hysteresis Power Loss (DC)'             ,'W'      ,W      ,eW                       ,'[W/(m^3.Hz)]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Hysteresis Power Loss (DC)'             ,'Wm'     ,Wm     ,eWm                      ,'[W/(Kg.Hz)]');
fprintf(fid,'%s,%s,%0.5g,     ,     ,     ,%s\n','Hysteresis Power Loss Minor Loop (DC)'  ,'WM'     ,WM                               ,'[W/(m^3.Hz)]');
fprintf(fid,'%s,%s,%0.5g,     ,     ,     ,%s\n','Hysteresis Power Loss Minor Loop (DC)'  ,'WMm'    ,WMm                              ,'[W/(Kg.Hz)]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Decreasing curve M-intercept'           ,'Mdec_r' ,Mdec_r ,eMdec_r                  ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Increasing curve M-intercept'           ,'Minc_r' ,Minc_r ,eMinc_r                  ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','M exchange bias'                        ,'Mexr'   ,Mexr   ,eMexr                    ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,%0.5g,%0.5g,%s\n','Remanence (exchange compensated)'       ,'Mr_sym' ,Mr_sym ,eMr_sym , Mr(1) , eMr(1) ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Saturation Magnetization + (tip avg)'   ,'Msp'    ,Msp    ,eMsp                     ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Saturation Magnetization - (tip avg)'   ,'Msn'    ,Msn    ,eMsn                     ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,%0.5g,%0.5g,%s\n','Avg Saturation Magnetization (tip avg)' ,'Ms_avg' ,Ms_avg ,eMs_avg , Ms(1) , eMs(1) ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,%0.5g,%0.5g,%s\n','Avg Saturation Magnetization (from n*m)','Ms_nm'  ,Ms_nm  ,eMs_nm  , Ms(1) , eMs(1) ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Decreasing curve H-intercept'           ,'Hdec_c' ,Hdec_c ,eHdec_c                  ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Increasing curve H-intercept'           ,'Hinc_c' ,Hinc_c ,eHinc_c                  ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,%0.5g,%0.5g,%s\n','H exchange bias'                        ,'Hexc'   ,Hexc   ,eHexc   , Hex(1), eHex(1),'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,%0.5g,%0.5g,%s\n','Coercivity (exchange compensated)'      ,'Hc_sym' ,Hc_sym ,eHc_sym , Hc(1) , eHc(1) ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Applied Field at Saturation + (tip avg)','Hsp'    ,Hsp    ,eHsp                     ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,     ,     ,%s\n','Applied Field at Saturation - (tip avg)','Hsn'    ,Hsn    ,eHsn                     ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,%0.5g,%0.5g,%0.5g,%s\n','Avg App. field at Saturation (tip avg)' ,'Hs_avg' ,Hs_avg ,eHs_avg , Hs(1) , eHs(1) ,'[A/m]');
fprintf(fid,'%s,%s,%0.5g,     ,%0.5g,     ,%s\n','Max. Susceptibility'                    ,'Xtmax'  ,Xtmax  ,          Xmax           ,'[]');
fprintf(fid,'%s,%s,%0.5g,     ,%0.5g,     ,%s\n','Init. Susceptibility'                   ,'Xtinit' ,Xtinit ,          Xinit          ,'[]');

% Add RAW INPUT data set here
fprintf(fid,'\n');
fprintf(fid,'RAW INPUT DATA\n');
fprintf(fid,' ,H [A/m], M [A/m], T [K], u0 [T.m/A], kb [J/K], eH [A/m], eM [A/m], eT [K], eu0 [T.m/A], ekb [J/K]\n');
% Print ith element of all data sets
Dgrp=[D, eD];
for i=1:size(Dgrp,1)
    fprintf(fid,' ,');
    fprintf(fid,'%0.5g,',Dgrp(i,:));
    fprintf(fid,'\n');
end

% Add montecarlo INPUT data sets here
fprintf(fid,'\n');
fprintf(fid,'MONTECARLO INPUT DATA SETS\n');
% Construct the labels for all the data sets
str1='';
for i=1:Nm
    str1=sprintf('%sHm %0.5g [A/m], Mm %0.5g [A/m], Tm %0.5g [K], u0m %0.5g [T.m/A], kbm %0.5g [J/K], ,',str1,i,i,i,i,i);
end
fprintf(fid,' ,%s\n',str1);
% Print ith element of all data sets
for j=1:size(Dm,2)
    fprintf(fid,' ,');
    for i=1:size(Dm,1)      
        for k=1:size(Dm,3)
            fprintf(fid,'%0.5g,',Dm(i,j,k));
        end
        fprintf(fid,' ,');
    end
    fprintf(fid,'\n');
end

% Add montecarlo OUTPUT data sets here
fprintf(fid,'\n');
fprintf(fid,'MONTECARLO OUTPUT DATA SETS\n');
% Grab the average output, standard deviation and all the individual 
% datasets and assemble them to display the stats for the montecarlo
% datasets
Mtgrp=[Ht,Mtavg,eMt,Mtm'];
str1=sprintf('Dataset Mt %0.5g,',(4:size(Mtgrp,2))-3);
str2=sprintf('Ht [A/m], Mt Avg [A/m], Mt Std.Err [A/m], %s',str1);
fprintf(fid,' ,%s',str2);
fprintf(fid,'\n');
% Print ith element of all data sets
for i=1:size(Mtgrp,1)  
    fprintf(fid,' ,');
    fprintf(fid,'%0.5g,',Mtgrp(i,:));
    fprintf(fid,'\n');
end

% Add the raw data for the optimization process here
fprintf(fid,'\n');
fprintf(fid,'OPTIMIZATION RAW DATA\n');
% Combine all values into a single matrix for display
for i=1:Nm
    % Print the header 
    fprintf(fid,' , Nm [], I [], n [1/m^3], m [J/T], al [], Xaf [], c [], k [A/m], S [(A/m)^2], W [W/(m^3.Hz)], Mdec_r [A/m], Minc_r [A/m], Mexr [A/m], Mr_sym [A/m], Msp [A/m], Msn [A/m], Ms_avg [A/m], Ms_nm [A/m], Hdec_c [A/m], Hinc_c [A/m], Hexc [A/m], Hc_sym [A/m], Hsp [A/m], Hsn [A/m], Hs_avg [A/m]\n');
    
    % Group the display values
    OPTgrpi=[(1:I)',squeeze(Pm(i,:,:)),Sm(i,:)',squeeze(Om(i,:,:))];
    for j=1:size(OPTgrpi,1)
        fprintf(fid,' , %0.5g,',i);
        for k=1:size(OPTgrpi,2)
            fprintf(fid,'%0.5g,',OPTgrpi(j,k));
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

% Close the file
fclose(fid);

%% ------------------------------------------------------------------------
% PLOT THINGS!
%--------------------------------------------------------------------------
% Save the current hysteresis figure after opening up toolbar and 
% menubar again ... the last run of the montecarlo data set
set(gcf,'MenuBar','figure');
set(gcf,'ToolBar','figure');
saveas(gcf,fighystfilename,'fig');

figure(2); clf(2);
% Draw parameter value n
subplot(3,3,1);
plot(1:I,P(1:I,1),'b-');
lbl=sprintf('Min parameter value: n=%0.5g',P(I,1));
title(lbl);
xlabel('iteration'); ylabel('n');

% Draw parameter value m
subplot(3,3,2);
plot(1:I,P(1:I,2),'b-');
lbl=sprintf('Min parameter value: m=%0.5g',P(I,2));
title(lbl);
xlabel('iteration'); ylabel('m');

% Draw parameter value al
subplot(3,3,3);
plot(1:I,P(1:I,3),'b-');
lbl=sprintf('Min parameter value: al=%0.5g',P(I,3));
title(lbl);
xlabel('iteration'); ylabel('al');

% Draw parameter value Xaf
subplot(3,3,4);
lbl=sprintf('Min parameter value: Xaf=%0.5g',P(I,4));
plot(1:I,P(1:I,4),'b-');
title(lbl);
xlabel('iteration'); ylabel('Xaf');

% Draw parameter value c
subplot(3,3,5);
lbl=sprintf('Min parameter value: c=%0.5g',P(I,5));
plot(1:I,P(1:I,5),'b-');
title(lbl);
xlabel('iteration'); ylabel('c');

% Draw parameter value k
subplot(3,3,6);
lbl=sprintf('Min parameter value: k=%0.5g',P(I,6));
plot(1:I,P(1:I,6),'b-');
title(lbl);
xlabel('iteration'); ylabel('k');

% Draw the change in the hysteresis power loss value
subplot(3,3,7);
lbl=sprintf('Optimized Power Loss Value: W=%0.5g',O(I,1));
plot(1:I,O(1:I,1),'b-');
title(lbl);
xlabel('iteration'); ylabel('log10(powerloss)');

% Draw the change in the SSD value
subplot(3,3,8);
lbl=sprintf('Min SSD value: SSD=%0.5g',S(I));
plot(1:I,log10(S(1:I)),'b-');
title(lbl);
xlabel('iteration'); ylabel('log10(error)');

% Draw the actual data points along with the line of best fit
subplot(3,3,9);
hold on;
plot(Ht,Mt,'b-');
plot(H,M,'r*');
hold off;
title('MH curve at last iteration');
xlabel('H'); ylabel('M');

% Save the optimization figure of the last run of the montecarlo data set
saveas(gcf,figoptfilename,'fig');
