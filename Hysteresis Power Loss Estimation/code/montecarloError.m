function [Dnew] = montecarloError(D, E, N)   
% Monte Carlo  Monte Carlo dataset creation
%   Create N indices uniformly distributed. Use these N indices (with
%   replacement) to create a new sample set that includes a random 
%   bounded error for each measured value. Use the error bounds in 
%   E to create a normal error distribution for each data point in 
%   D. 
%     N: [1,1]  number of random samples to create out of this data. 
%     M: [1,1]  number of data variables.
%     Q: [1,1]  number of data points used for estimation
%     U: [1,1]  number of data variables, e.g. three if D={H,M,T} and E={eH,eM,eT}
%     D: [Q,M]  the data to use for generating new random data set
%     E: [Q,M]  the error to use for generating a new data set
%     

    % Compute N random integer index values between 1 and Q
    Q=size(D,1);
    indices=randi([1 Q],N,1);
    
    % Sort the indices in ascending order
    indices=sort(indices,1,'ascend');
    
    % Use the unique index values to setup histogram bins. 
    uindices=unique(indices);
    
    % Bin the number of indices that are duplicated and create a count
    % matrix called histo.
    histo= num2cell(hist(indices,uindices)');

    % Generate information on the direction in which the data is going
    % based on the first data variable (assumed to be the H field).
    H=D(:,1);
    dir=sign(diff(H));
    dir(end+1)=dir(end);
    
    % For every unique index in the indices matrix we want to assign a
    % directionality for the sorting process based on the direction of the
    % original data at that index 
    udir=dir(uindices);
    udirltz=find(udir<0);
    udirgtz=find(udir>0);
    sortdir=cell(size(udir));
    sortdir(udirltz)=repmat({'descend'},size(udirltz));
    sortdir(udirgtz)=repmat({'ascend'},size(udirgtz));
    
    % Now generate a sorted array of random numbers matching the counts in
    % the histogram. For instance, if bin 1 contains 3 items then we want
    % to create three random numbers that are sorted relative to each other
    % in descending order --- thus preserving the original decreasing
    % monotonicity of the original data. This monotonic form is important
    % for the differential equation solver
    fun=@(N,sortdir) sort(randn(N,1),sortdir);
    
    % For each data variable extract the uniform rand integer indices
    % computed earlier and perturb them around the error value to produce a
    % new sample that has normally deviated from the mean error for that
    % data point.
    U=size(D,2);
    for i=1:U
        % Generate N normally distributed random numbers sorted
        % monotonically using the @fun anonymous function.
        randomcellarr=cellfun(fun,histo,sortdir,'UniformOutput',false);
        randomcellarr=randomcellarr(:);
        random=cell2mat(randomcellarr);
        
        % find the indices that are the same and arrange the random numbers
        % associated with these indices in descending order. 
        Dnew(:,i)=D(indices,i)+E(indices,i).*random;
    end
    
    % We need to ensure that the data is decreasing or
    % increasing in the same way (locally monotonic depending on location
    % along which we are following the hysteresis curve. A standard sort 
    % across the entire dataset without distinguishing the directionality
    % of the hysteresis curve will not work - we have two curves in one!
    debug=false;
    if (debug)
        Hnew=Dnew(:,1);
        dirnew=sign(diff(Hnew));
        dirnew(end+1)=dirnew(end);   
        test=[Hnew,indices,dirnew]; % monitor this variable 
    end
end