function [M, S, O]=hysteretic(H,D,P,PLB,PUB,VID)
    % Extract the measured data variables
    Hm=D(:,1); Mm=D(:,2);
        
    % Extract the constants
    u0=D(:,4); kb=D(:,5); T=D(:,3);
    u0=u0(1); kb=kb(1); T=T(1);
    
    % Extract the parameters
    n=P(1); m=P(2); al=P(3); Xaf=P(4); c=P(5); k=P(6);
        
    % Compute a and Ms values
    a=(m.*u0)./(kb.*T);
    Ms=n.*m;
        
    %% COMPUTING THE ANHYSTERETIC MAGNETIZATION FOR OUR SIMULATED LOCATIONS
    % Initialize my magnetization output array
    Man=zeros(size(H));

    % Iterate through step by step and build the hysteresis curve output.
    for i=1:numel(H)
        % Estimate starting magnetization starting from the initial
        % magnetization curve at the origin (value 0). To estimate Heff, we
        % use Mest as an estimate of the mean field contribution of the
        % magnetization by taking the previous calculated value as a guess.
        if(i>1)
            Mest=Man(i-1);
        else
            Mest=0;        
        end
        H(i)=H(i);

        % Compute the effective applied field Heff at increment i
        Heff=H(i)+al.*Mest;
        z=a*Heff;

        % Compute the value of Man. In the event that z=0, then we have a
        % limit condition such that lim(z=0) of (coth(z)-1/z) = 0 - this
        % limit itself was determined numerically by approaching from z=0+
        % and z=0-.
        if (z==0)
            Man(i)=Xaf.*Heff;
        else
            Man(i)=Ms.*(coth(z)-1./z) + Xaf.*Heff;
        end   
    end
    
    %% DEFINING THE DIFFERENTIAL EQUATION dMIRR/dH
    
    % Compute the value of Mirr by solving the differential equation
    % numerically. To do this lets create a nested function that describes
    % our differential equation sub problem.
    function [dMirrDH]=irreversible(h,H,Man,Mirr,del,delbar,k,al)
        interpolation='linear';
        
        % Do log log quadratic interpolation for more accuracy
        if(strcmp(interpolation,'loglogquadratic'))
            % Grab the 3 nearest neighbours (KNN stores the index location of
            % these nearest neighbours).
            KNN=knnsearch(H,h,'K',3);

            % Temporarily transform our data into log-log form such that 
            % Hl=log10(H-min(H)+s) and Manl=log10(Man-min(Man)+s). For 
            % log computations to be valid everywhere we need to shift our
            % data into the domain of natural numbers (all positive). This
            % shift should be reversed to extract our original values back.
            % Note that s=1 is a small value added before the log to ensure
            % that the value being logged is above zero since log10(0)=NaN.
            s=1;
            Hmin=min(H);
            Hl=log10(H-Hmin+s);
            Manmin=min(Man);
            Manl=log10(Man-Manmin+s);
            hl=log10(h-Hmin+s);

            % Use the (Hl,Manl) arrays to interpolate value (hl,Manl_)
            pl=polyfit(Hl(KNN),Manl(KNN),2);
            Hl_=hl;
            H_=10.^(Hl_) + Hmin - s;
            Manl_=polyval(pl,hl);
            Man_=10.^(Manl_) + Manmin - s;
            
        % Do linear interpolation for more speed
        elseif(strcmp(interpolation,'linear'))
            % Use the (H,Man) arrays to interpolate values H_=h and Man_=Man(h)
            % This is done by computing a linear interpolation between 
            % those values in H that are close to h. This is basically nearest
            % neighbours interpolation. 
            KNN=knnsearch(H,h,'K',2);
            p=polyfit(H(KNN),Man(KNN),1);  % Do linear fitting around 3 nearest neighbours
            H_=h;                          % Set our current analysis H point
            Man_=polyval(p,h);             % Set our current analysis M point
        end

        % DEBUG THE INTERPOLATION PROCESS. This will plot the interpolation
        % steps and ensure that the 
        debug=0;
        if(debug)
            if(strcmp(interpolation,'loglogquadratic'))
                fitHdensel=linspace(min(Hl(KNN)),max(Hl(KNN)),5);
                fitMandensel=polyval(pl,fitHdensel); 
                fitHdense=10.^(fitHdensel) + Hmin - s;
                fitMandense=10.^(fitMandensel) + Manmin - s;
            elseif(strcmp(interpolation,'linear'))
                fitHdense=linspace(min(H(KNN)),max(H(KNN)),5);
                fitMandense=polyval(p,fitHdense); 
            end
            figure(1); hold on; plot(H,Man,'r*'); xlabel('H'); ylabel('Man');
            plot(H_,Man_,'bo');        
            plot(fitHdense,fitMandense,'g-');
            hold off;
        end

        % The differential equation
        dMirrDH=(delbar.*(Man_-Mirr))./(k.*del-al.*(Man_-Mirr));
    end

    %% COMPUTING THE CURVE DIRECTIONALITY AND INDICES OF EACH DIRECTION
    
    % Compute the derivative dH/dt, dHm/dt for every data point H and Hm
    % Compute del to ensure that we know whether the H, Hm field is being
    % decreased or increased over time. Also find the index at which the 
    % change in direction happens. Note that the first saturation field 
    % value on either end of the hysteresis curve we will assume is either 
    % decreasing at the top half of the graph or increasing at the bottom 
    % half of the graph
    dHdt=diff(H);                     dHmdt=diff(Hm);
    delH=sign(dHdt);                  delHm=sign(dHmdt);
    delH(end+1)=delH(end);            delHm(end+1)=delHm(end);
    indchangeH=find(diff(delH)~=0);   indchangeHm=find(diff(delHm)~=0);
    indchangeH=[0;indchangeH];        indchangeHm=[0;indchangeHm];
    indchangeH=[indchangeH;numel(H)]; indchangeHm=[indchangeHm;numel(Hm)];
    
    % Combine the two data sets for processing below
    indchange={indchangeH,indchangeHm};
    del={delH,delHm};
    
    % We need to solve the differential equation for decreasing H, dH/dt<0
    % and for increasing H, dH/dt>=0. These must be separately solved to
    % fit with existing solvers - i.e. solvers only accept monotonically
    % increasing or decreasing spans. We should divide the indexing by
    % whether we are on the initial curve, the decreasing curve or the
    % increasing curve. At the end of this for loop we will have 
    %     ind{j}(i,1:2)=[ind{1} ind{2}]
    %         with ind(1)=[indices{...}  curveDir]
    %                     [indices{....} curveDir]
    %                     [indices{..}   curveDir]
    %              ind(2)=[indices{...}  curveDir]
    %                     [indices{....} curveDir]
    for j=1:numel(indchange)
        iend=numel(indchange{j})-1;
        for i=1:iend
            rng=((indchange{j}(i)+1):indchange{j}(i+1))';
            ind{j}(i,1)={rng};

            % label whether initial, decreasing or increasing curves
            indDec=find(del{j}==-1);
            indInc=find(del{j}==1);

            % Are the indices for the increasing, decreasing or initial 
            % part of the curve? Note that if we begin with a increasing
            % curve it will be assumed to be the initial curve. If we begin
            % with the decreasing curve no such assumption will be made -
            % not all cases are handled the order is assumed to be:
            % {initial, decreasing, increasing} or {decreasing, increasing}
            condDec=strfind(indDec',rng');
            condInc=strfind(indInc',rng');
            
            if (isempty(condDec) && ~isempty(condInc) && i==1)
                ind{j}(i,2)={'initial'};
            elseif (isempty(condDec) && ~isempty(condInc))
                ind{j}(i,2)={'increasing'};
            elseif (~isempty(condDec) && isempty(condInc))
                ind{j}(i,2)={'decreasing'};
            else
                ind{j}(i,2)={'unknown'};
            end
        end   
    end       

    %% GRABING THE DIFFERENT CURVE SECTIONS
    % Making the indices simpler to understand and use
    indInit=ind{1}{1,1};
    indDec=ind{1}{2,1};
    indInc=ind{1}{3,1};
    indmInit=[];
    indmDec=ind{2}{1,1};
    indmInc=ind{2}{2,1};  
    
    %% COMPUTING MIRR
    % The function associated with the dMirr/dH differential equation
    odefunInit=@(h,Mirr) irreversible(h,H(indInit),Man(indInit),Mirr,+1,1,k,al);  
    odefunDec=@(h,Mirr) irreversible(h,H(indDec),Man(indDec),Mirr,-1,1,k,al);  
    odefunInc=@(h,Mirr) irreversible(h,H(indInc),Man(indInc),Mirr,+1,1,k,al);  
    
    % Calculate the integral for every point in H. Note that we use an 
    % initial guess from either the zero condition from which our initial 
    % curve starts or our estimate of magnetic saturation as calculated 
    % from our parameters, Ms=n*m. This guess would be +Ms for the
    % decreasing part of the hysteresis curve and -Ms for increasing part.
    MirrInit0=Man(indInit); MirrInit0=MirrInit0(1);
    MirrDec0=Man(indDec);   MirrDec0=MirrDec0(1);
    MirrInc0=Man(indInc);   MirrInc0=MirrInc0(1);
    [~,MirrInit]=ode15s(odefunInit, H(indInit), MirrInit0);
    [~,MirrDec]=ode15s(odefunDec, H(indDec), MirrDec0);
    [~,MirrInc]=ode15s(odefunInc, H(indInc), MirrInc0);
    Mirr=[MirrInit;MirrDec;MirrInc];

    %% FIND MREV AND THE FINAL MANGETIZATION
    % Output the anhysteretic curve as the estimate of the actual
    % magnetization curve.
    Mrev=c.*(Man-Mirr);
    M=Mrev+Mirr;
    
    %% GRABING THE DIFFERENT CURVE SECTIONS
    % Grab the initial (I), decreasing (D), increasing (I), unknown (U)
    % curves for nearest neighbour interpolation        
    Hinit=H(indInit);  Hminit=[];           Haninit=H(indInit);  
    Hdec=H(indDec);    Hmdec=Hm(indmDec);   Handec=H(indDec);
    Hinc=H(indInc);    Hminc=Hm(indmInc);   Haninc=H(indInc);
    Minit=M(indInit);  Mminit=[];           Maninit=Man(indInit);
    Mdec=M(indDec);    Mmdec=Mm(indmDec);   Mandec=Man(indDec);
    Minc=M(indInc);    Mminc=Mm(indmInc);   Maninc=Man(indInc); 

    %% COMPUTE THE DC POWER LOSS FROM HYSTERESIS
    % Take the area under the decreasing curve and subtract from area under 
    % the increasing curve to get the power loss in watts [W] for a single 
    % hysteresis cycle only 
    Adec=trapz(Hdec,u0.*(Hdec+Mdec));
    Ainc=trapz(Hinc,u0.*(Hinc+Minc));
    W=abs((Adec-Ainc));
    
    %% PUT MEASURED DATA ONTOP OF SIMULATED DATA, EXTRACT COMMON POINTS
    % For the every field location in Hm (measured) we need to find the 
    % interpolated value of M (simulated) - this value is called Mt. 
    % This way we have two data sets of the same size: (Hm,Mm) and (Hm,Mt) 
    % with which we can compute the residual error for SSD. Two types of 
    % interpolation can be used:
    % (1) quadratic (after converting graph to
    %     log-log domain - this is because we know that the hysteresis
    %     roughly approximates logarithmic behaviour)
    % (2) linear (this is done in the normal domain near the origin as
    %     the origin is likely to be linear. Given the negative value
    %     for the field as a result of coercivity, the log computations
    %     in this region are also undefined
    % *For curves in the first quadrant log-log calculations aren't a
    %  problem. However, for other regions the data must be transformed
    %  to complete the interpolation. For this reason, at the moment only
    %  linear interpolation is used - this is valid as long as our 
    %  simulated H and M data set is sufficiently dense to justify this
    %  linear approximation. The log-log quadratic approach is valid for
    %  sparse simulation data and should be implemented in these cases. 
    interpolation='linear';
    
    % For higher accuracy do log-log quadratic interpolation
    % NOT WORKING YET!!!!! Use Linear
    if (strcmp(interpolation,'loglogquadratic'))
        %KNNinit=knnsearch(Hinit,Hminit,'K',3); % Throws err b/c Hminit empty. 
        KNNdec=knnsearch(Hdec,Hmdec,'K',3);
        KNNinc=knnsearch(Hinc,Hminc,'K',3);

        % Temporarily transform our data into log-log form such that 
        % Hdecl=log10(Hdec-min(Hdec)+sdec)
        % Hincl=log10(Hinc-min(Hinc)+sinc)
        % Mdecl=log10(Mdec-min(Mdec)+sdec)
        % Mincl=log10(Minc-min(Minc)+sinc)
        % For log computations to be valid everywhere we need to shift our
        % data into the domain of natural numbers (all positive). This
        % shift should be reversed to extract our original values back.
        % Note that sdec=1/sinc=1 is a small value added before the log to 
        % ensure that the value being logged is above zero since 
        % log10(0)=NaN.  
        
        % Do the transformation to log domain
        sdec=1;
        Hdecmin=min(Hdec);
        Hdecl=log10(H-Hdecmin+sdec);
        Mdecmin=min(Mdec);
        Mdecl=log10(Mdec-Mdecmin+sdec);
        Hmdecl=log10(Hmdec-Hdecmin+sdec);
        
        % Interpolate
        Mtdecl=zeros(size(Mmdec));
        for i=1:size(KNNdec,1)
            pdecl=polyfit(Hdecl(KNNdec(i,:)),Mdecl(KNNdec(i,:)),2);
            Mtdecl(i)=polyval(pdecl,Hmdecl(i));
        end
        % Do the inverse transformation to the normal domain.
        Mtdec=10.^(Mtdecl) + Mdecmin - sdec;
        
        % Do the transoformation to log domain
        sinc=1;
        Hincmin=min(Hinc);
        Hincl=log10(H-Hincmin+sinc);
        Mincmin=min(Minc);
        Mincl=log10(Minc-Mincmin+sinc);
        Hmincl=log10(Hminc-Hincmin+sinc);
        
        % Interpolate 
        Mtincl=zeros(size(Mminc));
        for i=1:size(KNNinc,1)
            pincl=polyfit(Hincl(KNNinc(i,:)),Mincl(KNNinc(i,:)),2);
            Mtincl(i)=polyval(pincl,Hmincl(i));
        end
        Mtinc=10.^(Mtincl) + Mincmin - sinc;
        
        % Combine the increasing and decreasing curves.
        Mt=[Mtdec;Mtinc];
        
    % For higher speed do linear interpolation
    elseif (strcmp(interpolation,'linear'))
        %KNNinit=knnsearch(Hinit,Hminit,'K',2); % Throws err b/c Hminit empty. 
        KNNdec=knnsearch(Hdec,Hmdec,'K',2);
        KNNinc=knnsearch(Hinc,Hminc,'K',2);

        % Interpolate linear
        Mtdec=zeros(size(Mmdec));
        for i=1:size(KNNdec,1)
            pdec=polyfit(Hdec(KNNdec(i,:)),Mdec(KNNdec(i,:)),1);
            Mtdec(i)=polyval(pdec,Hmdec(i));
        end
        Mtinc=zeros(size(Mminc));
        for i=1:size(KNNinc,1)
            pinc=polyfit(Hinc(KNNinc(i,:)),Minc(KNNinc(i,:)),1);
            Mtinc(i)=polyval(pinc,Hminc(i));
        end
        Mt=[Mtdec;Mtinc];
    end
    
    %% COMPUTE THE SSD FUNCTION VALUE BEFORE APPLYING PENALTIES
    % Compute the sum of square difference from the residuals. Note that
    % since the measured data (Mm) is sparse in comparison to the data 
    % generated as output (M) we need to create a new data set sampled at
    % measurement points (Hm).    
    r = (Mm - Mt); 
    S=sum(r.^2);
    
    %% COMPUTE THE MAGNETIC SATURATION VALUES Msp, Msn, Hsp, Hsn
    % Use the decreasing curve as reference to find the magnetization,
    % Minci values in the increasing curve at the particular Hdeci values 
    % of the decreasing curve. Use interpolation to find the Magnetization 
    % Minci from the two nearest neighbour Hinc1i Hinc2i at Hdeci.
    KNNs=knnsearch(Hinc,Hdec,'K',2);
    psfun=@(i) polyfit(Hinc(KNNs(i,:)),Minc(KNNs(i,:)),1);
    ps=arrayfun(psfun,(1:size(KNNs,1))','UniformOutput',0);
    Mincfun=@(i) polyval(ps{i},Hdec(i));
    Minc_new=arrayfun(Mincfun,(1:numel(Hdec))','UniformOutput', 0);
    Hinc_new=Hdec;
    Minc_new=cell2mat(Minc_new);
    
    % Grab the loop tip of the decreasing curve and increasing curve in
    % first quadrant and third quadrant. Find the magnetization value and 
    % H value at which the values are within <=0.1% difference of each 
    % other.    
    sthreshold=0.001;
    percentdiff=(Mdec-Minc_new)./((Mdec+Minc_new)./2);
    diffindp=find(percentdiff<=sthreshold & percentdiff>=0);
    diffindn=find(percentdiff>=-sthreshold & percentdiff<=0);
    Msp=Mdec(diffindp); Msn=Mdec(diffindn);
    Hsp=Hdec(diffindp); Hsn=Hdec(diffindn);
    Msp=mean(Msp);      Msn=mean(Msn);
    Hsp=mean(Hsp);      Hsn=mean(Hsn);
    
    %% PRODUCE A PENALTY FUNCTION FOR INEQUALITY CONSTRAINTS
    
    % We have the following basic boundary constraints in our problem:
    %  nlb   <=  n  <= nub   ===>>  nlb-n<=0,     n-nub<=0
    %  mlb   <=  m  <= mub   ===>>  mlb-m<=0,     m-mub<=0
    %  allb  <= al  <= alub  ===>>  allb-al<=0,   al-alub<=0
    %  Xaflb <= Xaf <= Xafub ===>>  Xaflb-Xaf<=0, Xaf-Xafub<=0
    %  clb   <=  c  <= cub   ===>>  clb-c<=0,     c-cub<=0
    %  klb   <=  k  <= kub   ===>>  klb-k<=0,     k-kub<=0
    nlb=PLB(1);     nub=PUB(1);
    mlb=PLB(2);     mub=PUB(2);
    allb=PLB(3);    alub=PUB(3);
    Xaflb=PLB(4);   Xafub=PUB(4);
    clb=PLB(5);     cub=PUB(5);
    klb=PLB(6);     kub=PUB(6);
    
    % Construct the boundary inequality constraint functions
    GInlb=nlb-n;       GInub=n-nub;
    GImlb=mlb-m;       GImub=m-mub;
    GIallb=allb-al;    GIalub=al-alub;
    GIXaflb=Xaflb-Xaf; GIXafub=Xaf-Xafub;
    GIclb=clb-c;       GIcub=c-cub;
    GIklb=klb-k;       GIkub=k-kub;
    % Aggregate all the boundary inequality constraints
    GI=[GInlb; GInub; GImlb; GImub; GIallb; GIalub; GIXaflb; GIXafub; GIclb; GIcub; GIklb; GIkub];
    
    % Generate the values for rho and beta for the penality function:    
    rho=S; % Should be a large number greater than zero
    beta=zeros(size(GI));
    beta(GI>0)=rho*1e20;
    
    % Compute the penalty function for inequality constraints
    PI=sum(beta.*(GI).^2)/numel(GI);
       
    GE1=0; GE2=0; GE3=zeros(2,1); GE4=zeros(2,1); GE5=0; GE6=0; GE7=zeros(2,1);
    %% PRODUCE A PENALTY FUNCTION FOR
    % (1) INITIAL SUSCEPTIBILITY (ANHYSTERETIC)
    %   We can use the above computed anhysteretic/hysteretic equations 
    %   to derive a constraint, namely Xan=dMan/dH at (M=0, H=0) or 
    %   Xan_init.  To do so, we need to rearrange our derivative with the 
    %   understanding that dM/dH=dMan/dH for the initial anhysteretic 
    %   magnetization curve. This allows us to rearrange our earlier 
    %   equation for dMan/dH and determine the initial magnetization 
    %   constraint Xan_init on the anhysteretic curve.
    
    % Grab the M=0,H=0 value of our initial susceptibility 
    % (slope at origin)from our initial (H, Man) curve
    Xan_init=diff(Maninit)./diff(Haninit);
    Xan_init=Xan_init(1);
    
    % Construct the constraint. 
    GE1=Xan_init-(n*m*a - 3*Xaf)/(al*n*m*a + 3*al*Xaf - 3);
    
    %% PRODUCE A PENALTY FUNCTION FOR    
    % (2) CONSTRAINT: INITIAL SUSCEPTIBILITY (HYSTERETIC)
    %   Like the anhysteretic curve, our magnetization curve has an 
    %   initial magnetization slope, X_init=dM/dH for (H=0, M=0). We can 
    %   solve for this by noting that M=0 and Mirr=0 at the origin and 
    %   substituting into the equations in the model. 
    
    % Grab the M=0,H=0 value of our initial susceptibility 
    % (slope at origin) from our initial (H, M) curve
    X_init=diff(Minit)./diff(Hinit);
    X_init=X_init(1);
    
    % Construct the constraint
    GE2=X_init-c.*Xaf + c*n*m*a/3;
        
    %% PRODUCE A PENALTY FUNCTION FOR    
    % (3) CONSTRAINT: LOOP TIP SUSCEPTIBILITY (HYSTERETIC/ANHYSTERETIC)
    %   At the point of saturation, all magnetic domains have been aligned 
    %   and the resulting irreversible field is now fixed in a single 
    %   orientation. Relaxing the applied field now will simply cause 
    %   changes in the reversible field since all the pinned moments are 
    %   now already in the applied field direction. Relaxation will simply 
    %   reverse the direction of the surrounding moments to eventually 
    %   align with the pinned moments. This means that dMirr/dH will not 
    %   change as the field is relaxed from saturation - so on the 
    %   decreasing slope of loop tip and on the increasing slope of the 
    %   loop tip only will the dMirr/dH=0. Note, however, that as soon as 
    %   the H field reverses direction (crosses the zero value) this 
    %   condition is no longer valid ... at this point the irreversible 
    %   components will once again begin to change.
    
    % Define which points form the loop tip on our curve (in this case
    % anything above 2000 [Oe] or below -2000 [Oe]
    Hdectipbnd=((10^3)/(4*pi())).*2000;  Hinctipbnd=-Hdectipbnd;
    Hdectip=Hdec(Hdec>Hdectipbnd);       Hinctip=Hinc((Hinc<Hinctipbnd));
    Mdectip=Mdec(Hdec>Hdectipbnd);       Minctip=Minc((Hinc<Hinctipbnd));
    Handectip=Handec(Handec>Hdectipbnd); Haninctip=Haninc(Haninc<Hinctipbnd);
    Mandectip=Mandec(Handec>Hdectipbnd); Maninctip=Maninc(Haninc<Hinctipbnd);
    
    % if we are traversing a minor loop of the curve, our condition of 2000
    % [Oe] and up being the loop tip may not be true and will lead to an
    % empty matrix for our Hdectip, Mdectip, Handectip, Mandectip values. 
    if (isempty(Hdectip)==0 && isempty(Mdectip)==0 && isempty(Handectip)==0 && isempty(Mandectip)==0)
        % Grab the loop tip of the decreasing curve
        X_dectip=diff(Mdectip)./diff(Hdectip);
        X_dectip=X_dectip(1:2);
        Xan_dectip=diff(Mandectip)./diff(Handectip);
        Xan_dectip=Xan_dectip(1:2);
        
        % Construct the decreasing tip constraint
        GE3=X_dectip-c.*Xan_dectip;
    end

    % if we are traversing a minor loop of the curve, our condition of 2000
    % [Oe] and up being the loop tip may not be true and will lead to an
    % empty matrix for our Hinctip, Minctip, Haninctip, Maninctip values. 
    if (isempty(Hinctip)==0 && isempty(Minctip)==0 && isempty(Haninctip)==0 && isempty(Maninctip)==0)
        % Grab the loop tip of the increasing curve
        X_inctip=diff(Minctip)./diff(Hinctip);
        X_inctip=X_inctip(1:2);
        Xan_inctip=diff(Maninctip)./diff(Haninctip);
        Xan_inctip=Xan_inctip(1:2);

        % Construct the increasing tip constraint
        GE4=X_inctip-c.*Xan_inctip;
    end
    
    %% PRODUCE A PENALTY FUNCTION FOR    
    % (4) CONSTRAINT: COERCIVITY CONSTRAINT
    %   Another constraint can be constructed by taking a closer look at 
    %   H intercept (or coercive) point conditions. As noted earlier the 
    %   that at saturation we have full alignment of all magnetic moments 
    %   in the material. So now when the magnetic field, H, is relaxed to 
    %   zero the only contributors the the magnetization are from the 
    %   irreversible field. As the field direction of H is changed to the 
    %   coercive point our total field, M, goes to zero. This means that 
    %   there is a direct relationship between or irreversible components, 
    %   Mirr and the anhysteretic components, Man. Add to this that our 
    %   susceptibility dMirr/dH falls as the difference between dM/dH and 
    %   dMan/dH on any point on our hysteresis curve we can construct the 
    %   following formulation for a coercive point constraint
    
    % Find two points (H0,M0) & (H1,M1) that are closest to the H axis
    % intercept for the decreasing and increasing curves (for both the
    % hysteretic, M and anhysteretic Man curves)
    KNNdec_c=knnsearch(Mdec,0,'K',2);
    KNNinc_c=knnsearch(Minc,0,'K',2);
    KNNandec_c=knnsearch(Mandec,0,'K',2);
    KNNaninc_c=knnsearch(Maninc,0,'K',2);
    
    % Use linear interpolation to figure out the H-intercept point for
    % the decreasing and increasing loop part for both the hysteretic, M
    % and the anhysteretic, Man curves
    pdec_c=polyfit(Mdec(KNNdec_c),Hdec(KNNdec_c),1);
    Hdec_c=polyval(pdec_c,0);
    pinc_c=polyfit(Minc(KNNinc_c),Hinc(KNNinc_c),1);
    Hinc_c=polyval(pinc_c,0);
 
    pandec_c=polyfit(Mandec(KNNandec_c),Handec(KNNandec_c),1);
    Handec_c=polyval(pandec_c,0);
    paninc_c=polyfit(Maninc(KNNaninc_c),Haninc(KNNaninc_c),1);
    Haninc_c=polyval(paninc_c,0);
    
    % Compute the Man(Hc)value on the anhysteretic curve at the hysteretic
    % curve coercivity point Hc.
    Hdec_Hc=Hdec_c;
    pandec_Hc=polyfit(Handec(KNNandec_c),Mandec(KNNandec_c),1);
    Man_Hc(1)=polyval(pandec_Hc,Hdec_Hc);
    
    Hinc_Hc=Hinc_c;
    paninc_Hc=polyfit(Haninc(KNNaninc_c),Maninc(KNNaninc_c),1);
    Man_Hc(2)=polyval(paninc_Hc,Hinc_Hc);
    
    % Grab the decreasing & increasing part of the loop and compute the 
    % derivative of M and Man.
    X_c(1)=diff(Mdec(KNNdec_c))./diff(Hdec(KNNdec_c));
    X_c(2)=diff(Minc(KNNinc_c))./diff(Hinc(KNNinc_c));
    Xan_c(1)=diff(Mandec(KNNandec_c))./diff(Handec(KNNandec_c));
    Xan_c(2)=diff(Maninc(KNNaninc_c))./diff(Haninc(KNNaninc_c));
    
    % Account for the directionality of the the H field change.
    del_delbar_c=[1/-1, 1/1];
    
    % Construct the constraint
    denominator_c=(1/(1-c)).*X_c - (c/(1-c)).*Xan_c;
    GE5=(k-del_delbar_c.*((Man_Hc./(1-c)).*(1./denominator_c)))';  
        
    %% PRODUCE A PENALTY FUNCTION FOR    
    % (5) CONSTRAINT: REMANANCE CONSTRAINT
    %   Another constraint can be constructed by setting H=0 and finding 
    %   the M axis intercept value. This value will be the remanent
    %   magnetization, Mr. Note that given this condition our equation for
    %   M=cMan + (1-c)Mirr becomes: M(r)=c Man(r) + (1-c) Mirr(r). Likewise
    %   we get the susceptibility: X(r)=c Xan(r) + (1-c) Xirr(r). We can
    %   now solve for M(r) to get our remanent point constraint.
    
    % Find the M axis intercept by finding the nearest neighbouring indices
    % to H=0 for our decreasing and increasing curve. 
    KNNdec_r=knnsearch(Hdec,0,'K',2);
    KNNinc_r=knnsearch(Hinc,0,'K',2);
    KNNandec_r=knnsearch(Handec,0,'K',2);
    KNNaninc_r=knnsearch(Haninc,0,'K',2);
    
    % Use linear interpolation to figure out the M-intercept point for
    % the decreasing and increasing loop part for both the hysteretic, M
    % and the anhysteretic, Man curves
    pdec_r=polyfit(Hdec(KNNdec_r),Mdec(KNNdec_r),1);
    Mdec_r=polyval(pdec_r,0);
    pinc_r=polyfit(Hinc(KNNinc_r),Minc(KNNinc_r),1);
    Minc_r=polyval(pinc_r,0);
 
    pandec_r=polyfit(Handec(KNNandec_r),Mandec(KNNandec_r),1);
    Mandec_r=polyval(pandec_r,0);
    paninc_r=polyfit(Haninc(KNNaninc_r),Maninc(KNNaninc_r),1);
    Maninc_r=polyval(paninc_r,0);
    
    % Grab the Man(Mr) points and the M(Mr) points
    Man_Mr(1)=Mandec_r;
    Man_Mr(2)=Maninc_r;
    M_Mr(1)=Mdec_r;
    M_Mr(2)=Minc_r;
    
    % Grab the decreasing & increasing part of the loop and compute the 
    % derivative of M and Man at the remanence point
    X_r(1)=diff(Mdec(KNNdec_r))./diff(Hdec(KNNdec_r));
    X_r(2)=diff(Minc(KNNinc_r))./diff(Hinc(KNNinc_r));
    Xan_r(1)=diff(Mandec(KNNandec_r))./diff(Handec(KNNandec_r));
    Xan_r(2)=diff(Maninc(KNNaninc_r))./diff(Haninc(KNNaninc_r));
    
    % Account for the directionality of the H field change
    del_delbar_r=[-1/1, 1/1];
    
    % Construct the constraint
    demoninator_r=(1./(X_r-c.*Xan_r) + al./((1-c).*1));
    GE6=(M_Mr-(Man_Mr - del_delbar_r.*(k./demoninator_r)))';  
    
    %% PRODUCE A PENALTY FUNCTION FOR    
    % (6) CONSTRAINT: INITIAL INCREASING TO LOOP TIP CONSTRAINT
    %    We note that when the material is completely demagnetized our
    %    initial magnetization curve will approach saturation near the loop
    %    tips at approximately the same rate as our anhysteretic curve.
    %    That is, X=Xan for the initial curve loop tip to saturation. So
    %    X(it)=Xan(it). Given this if we now look at our equation:
    %    X=cXan+(1-c)Xirr, we can substitute our condition X(it)=Xan(it) to
    %    obtain a final constraint for our parameters. Note that we use the
    %    definitiona Mirr=(M-cMan)/(1-c) to derive:
    %    (Man-Mirr)=(Man-M)/(1-c) and substitute into Xirr=dMirr/dH. 
    
    % Define which points form the loop tip on our initial curve 
    % (in this case anything above 2000 [Oe] or below -2000 [Oe]
    Hinittipbnd=((10^3)/(4*pi())).*2000;  
    Hinittip=Hinit(Hinit>Hinittipbnd);
    Minittip=Minit(Hinit>Hinittipbnd);  
    Haninittip=Haninit(Haninit>Hinittipbnd);
    Maninittip=Maninit(Haninit>Hinittipbnd);
    
    % if we are traversing a minor loop of the curve, our condition of 2000
    % [Oe] and up being the loop tip may not be true and will lead to an
    % empty matrix for our Hinittip, Minittip, Haninittip, Maninittip values. 
    if (isempty(Hinittip)==0 && isempty(Minittip)==0 && isempty(Haninittip)==0 && isempty(Maninittip)==0)    
        % Compute the slope of the initial loop tip toward saturation
        X_inittip=diff(Minittip)./diff(Hinittip);
        X_inittip=X_inittip(1:2);

        Hinittip_=Hinittip(1:2);
        Minittip_=Minittip(1:2);
        Maninittip_=Maninittip(1:2);

        % Construct the decreasing tip constraint
        GE7=Minittip_-(Maninittip_ - (k.*1.*X_inittip.*(1-c))./(1+al.*X_inittip));
    end
        
    %% AGGREGATE THE EQUALITY CONSTRAINTS AND COMPUTE THE PENALTY
    
    % Aggregate our equality constraints.
    GE=[GE1; GE2; GE3; GE4; GE5; GE6; GE7];
    
    % Compute the penalty function for our equality constraints
    PE=sum(rho.*(GE).^2)/numel(GE);
    
    %% OUTPUT VARIABLES
      
    % The power loss per (m^3 . Hz)
    O(1)=W;
    
    % The M-intercept, center, remanence, exchange bias, saturation values
    O(2)=Mdec_r;                % Decreasing curve intercept
    O(3)=Minc_r;                % Increasing curve intercept
    O(4)=(Mdec_r+Minc_r)/2;     % M center position (exchange bias)
    O(5)=Mdec_r-O(4);           % Remanance under hysteresis symmetry assumption
    O(6)=Msp;                   % Saturation Magnetization positive (avg of tip)
    O(7)=Msn;                   % Saturation Magnetization negative (avg of tip)
    O(8)=(abs(Msp)+abs(Msn))/2; % The average saturation magnetization computed from the +/- tip average magnetizations [A/m]
    O(9)=Ms;                    % The average saturation magnetization computed from parameter n*m [A/m]
    
    % The H-intercept, center, coercivity, exhange bias values
    O(10)=Hdec_c;               % Decreasing curve intercept
    O(11)=Hinc_c;               % Increasing curve intercept
    O(12)=(Hdec_c+Hinc_c)/2;    % H center position (exchange bias)
    O(13)=Hinc_c-O(12);         % Coercivity under hysteresis symmetry assumption    
    O(14)=Hsp;                  % Applied field at saturation positive 
    O(15)=Hsn;                  % Applied field at saturation negative (avg of tip)
    O(16)=(abs(Hsp)+abs(Hsn))/2;% The average applied field at saturation computed from the +/- tip average applied fields at saturation [A/m]
    
    %% PRODUCE DEBUGGING PLOTS
    debug=1;
    if(debug)        
        % Plot the known data and the generated curve at those points.
        fighandle=figure(1); clf; hold on; grid on;
        axis([-8e5-1e5 8e5+1e5 -5e5-1e5 5e5+1e5]);

        % TITLE
        title('GENERATED HYSTERESIS CURVE');
        
        p0=plot(Hm,Mm,'k*'); xlabel('H'); ylabel('M');
        p1=plot(Hm,Mt,'kx');
        
        % Plot the generated curve
        p2=plot(H,M,'b-');        % Hysteretic
        p3=plot(H,Man,'r-');      % Anhysteretic
        p4=plot(H,Mirr,'g-');     % Irreversible
        p5=plot(H,Mrev,'c-');     % Reversible
        
        % Plot the rate of change curves
        p6=plot(H(1:end-1),diff(M),'b-.');    % Hysteretic
        p7=plot(H(1:end-1),diff(Man),'r-.');  % Anhysteretic
        p8=plot(H(1:end-1),diff(Mirr),'g-.'); % Irreversible
        p9=plot(H(1:end-1),diff(Mrev),'c-.'); % Reversible
        
        % Plot the X(tip) & Xan(tip) slopes
        p10=plot(Hdectip,Mdectip,'g--');     % Hysteretic
        p10a=plot(Hinctip,Minctip,'g--');    % Hysteretic
        p11=plot(Handectip,Mandectip,'g--'); % Anhysteretic
        p11a=plot(Haninctip,Maninctip,'g--');% Anhysteretic
        set([p10,p10a,p11,p11a],'LineWidth',3)
        
        % Plot the H-intercept Points on the curve
        p12=plot([Hdec_c,Hinc_c],[0,0],'ko');                   % Hysteretic
        p13=plot([Handec_c,Haninc_c],[0,0],'kd');               % Anhysteretic
        p14=plot([Hdec_Hc,Hinc_Hc],[Man_Hc(1),Man_Hc(2)],'ks'); % Anhysteretic
        set([p12,p13,p14],'MarkerSize',10);
        
        % Plot the X(c) and Xan(c) slopes 
        p15=plot(Hdec(KNNdec_c),Mdec(KNNdec_c),'k:');          % Hysteretic
        p15a=plot(Hinc(KNNinc_c),Minc(KNNinc_c),'k:');         % Hysteretic 
        p16=plot(Handec(KNNandec_c),Mandec(KNNandec_c),'k:');  % Anhysteretic
        p16a=plot(Haninc(KNNaninc_c),Maninc(KNNaninc_c),'k:'); % Anhysteretic
        set([p15,p15a,p16,p16a],'LineWidth',3);
        
        % Plot the M-intercept Points on the curve
        p17=plot([0,0],[Mdec_r,Minc_r],'mo');                   % Hysteretic
        p18=plot([0,0],[Mandec_r,Maninc_r],'md');               % Anhysteretic
        p19=plot([0,0],[Man_Mr(1),Man_Mr(2)],'ms');             % Anhysteretic
        set([p17,p18,p19],'MarkerSize',10);
        
        % Plot the X(r) and Xan(r) slopes
        p20=plot(Hdec(KNNdec_r),Mdec(KNNdec_r),'m:');          % Hysteretic
        p20a=plot(Hinc(KNNinc_r),Minc(KNNinc_r),'m:');         % Hysteretic
        p21=plot(Handec(KNNandec_r),Mandec(KNNandec_r),'m:');  % Anhysteretic
        p21a=plot(Haninc(KNNaninc_r),Maninc(KNNaninc_r),'m:'); % Anhysteretic
        set([p20,p20a,p21,p21a],'LineWidth',3);

        % Plot the X(init tip) slopes
        p22=plot(Hinittip,Minittip,'k-.');     % Hysteretic
        set(p22,'LineWidth',3);
        
        % Add legend
        [leghandle,markers]=legend([p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22], ...
                                      'Mm=Measured Data', ...
                                      'Mt=Hysteretic (at Data)', ...
                                      'M=Hysteretic', ...
                                      'Man=Anhysteretic', ...
                                      'Mirr=Irreversible', ...
                                      'Mrev=Reversible', ...
                                      'dMdH', ...
                                      'dManDH', ...
                                      'dMirrDH', ...
                                      'dMrevDH', ...
                                      'X(tip)', ...
                                      'Xan(tip)', ...
                                      'H(c)', ...
                                      'Han(c)', ...
                                      'Man(Hc)', ...
                                      'X(c)', ...
                                      'Xan(c)', ...
                                      'M(r)', ...
                                      'Man(r)', ...
                                      'Man(Mr)', ...
                                      'X(r)', ...
                                      'Xan(r)', ...
                                      'X(init tip)', ...
                                      'Location','SouthEast');
        % Set the legend font size and icon sizesmaller
        set(leghandle,'FontSize',10);
        
        % Remove the menu and toolbar from figure
        set(fighandle, 'MenuBar', 'none');
        set(fighandle, 'ToolBar', 'none');
        
        % Print out the parameter values, error values and power value                          
        lbl='';
        lbl=sprintf('%sn=%0.5g\n',lbl,n);
        lbl=sprintf('%sm=%0.5g\n',lbl,m);
        lbl=sprintf('%sal=%0.5g\n',lbl,al);
        lbl=sprintf('%sXaf=%0.5g\n',lbl,Xaf);
        lbl=sprintf('%sc=%0.5g\n',lbl,c);
        lbl=sprintf('%sk=%0.5g\n\n',lbl,k);
        lbl=sprintf('%sPI=%0.5g\n',lbl,PI);
        lbl=sprintf('%sGE1=%0.5g\n',lbl,GE1);
        lbl=sprintf('%sGE2=%0.5g\n',lbl,GE2);
        lbl=sprintf('%sGE3=%0.5g\n',lbl,GE3(1));
        lbl=sprintf('%sGE3=%0.5g\n',lbl,GE3(2));
        lbl=sprintf('%sGE4=%0.5g\n',lbl,GE4(1));
        lbl=sprintf('%sGE4=%0.5g\n',lbl,GE4(2));
        lbl=sprintf('%sGE5=%0.5g\n',lbl,GE5(1));
        lbl=sprintf('%sGE5=%0.5g\n',lbl,GE5(2));
        lbl=sprintf('%sGE6=%0.5g\n',lbl,GE6(1));
        lbl=sprintf('%sGE6=%0.5g\n',lbl,GE6(2));
        lbl=sprintf('%sGE7=%0.5g\n',lbl,GE7(1));        
        lbl=sprintf('%sGE7=%0.5g\n',lbl,GE7(2));
        lbl=sprintf('%sPE=%0.5g\n',lbl,PE);
        lbl=sprintf('%sS=%0.5g\n\n',lbl,S);
        lbl=sprintf('%sW=%0.5g\n',lbl,O(1));
        lbl=sprintf('%sMdec r=%0.5g\n',lbl,O(2));
        lbl=sprintf('%sMinc r=%0.5g\n',lbl,O(3));
        lbl=sprintf('%sMexr=%0.5g\n',lbl,O(4));
        lbl=sprintf('%sMr sym=%0.5g\n',lbl,O(5));
        lbl=sprintf('%sMsp=%0.5g\n',lbl,O(6));
        lbl=sprintf('%sMsn=%0.5g\n',lbl,O(7));
        lbl=sprintf('%sMs avg=%0.5g\n',lbl,O(8));
        lbl=sprintf('%sMs nm=%0.5g\n',lbl,O(9));
        lbl=sprintf('%sHdec c=%0.5g\n',lbl,O(10));
        lbl=sprintf('%sHinc c=%0.5g\n',lbl,O(11));
        lbl=sprintf('%sHexc=%0.5g\n',lbl,O(12));
        lbl=sprintf('%sHc sym=%0.5g\n',lbl,O(13));
        lbl=sprintf('%sHsp=%0.5g\n',lbl,O(14));
        lbl=sprintf('%sHsn=%0.5g\n',lbl,O(15));
        lbl=sprintf('%sHs avg=%0.5g\n',lbl,O(16));
        lblu='';
        lblu=sprintf('%s[1/m3]\n',lblu);
        lblu=sprintf('%s[J/T]\n',lblu);
        lblu=sprintf('%s[]\n',lblu);
        lblu=sprintf('%s[]\n',lblu);
        lblu=sprintf('%s[]\n',lblu);
        lblu=sprintf('%s[A/m]\n\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n',lblu); 
        lblu=sprintf('%s[(A/m)2]\n',lblu); 
        lblu=sprintf('%s[(A/m)2]\n',lblu); 
        lblu=sprintf('%s[(A/m)2]\n',lblu); 
        lblu=sprintf('%s[(A/m)2]\n',lblu); 
        lblu=sprintf('%s[(A/m)2]\n',lblu); 
        lblu=sprintf('%s[(A/m)2]\n',lblu);
        lblu=sprintf('%s[(A/m)2]\n\n',lblu);
        lblu=sprintf('%s[W/(m3.Hz)]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        lblu=sprintf('%s[A/m]\n',lblu);
        
        text(-8e5,5e5,lbl,'HorizontalAlignment','left','VerticalAlignment','top');
        text(-4.5e5,5e5,lblu,'HorizontalAlignment','left','VerticalAlignment','top');
        
        % Set the position and size of the figure
        scrsize=get(0,'screensize');
        pad=100;
        swidth=scrsize(3); sheight=scrsize(4);
        width=sheight*0.95; height=sheight*0.85;
        set(fighandle, 'Position', [swidth/2-width/2, sheight/2 - height/2, width, height]);
                
        hold off;

        % Write the video frame to file.
        frame=getframe(gcf);
        %writeVideo(VID,frame)
    end    
    
    %% OUTPUT THE SSD ERROR VALUE AFTER ADDING PENALTY FUNCTION VALUES
     % THIS VALUE IS THE ONE THAT WILL BE MINIMIZED
    
    % Compute the SSD 
    S=S + PI + PE;   
end