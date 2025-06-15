function [PARAM, SSD, OUT]= paramEstimate (objfun, confun, DATA, PARAM0, ITER, ALGORITHM)    
    % Initialize the first iteration number in i and also initialize the
    % maximum number of iterations to execute.
    i=1; 
    
    % Set the various precisions and tolerances
    precision=1e-6;
   
    % Initialize the count for the number of parameters & data measurement
    % points/samples & the number of data variables.
    %     N:     [1,1]     number of parameters, 
    %     Q:     [1,1]     number of data points used for estimation
    %     U:     [1,1]     number of data variables
    N=numel(PARAM0);
    Q=numel(DATA(:,1));
    U=numel(DATA(1,:));
    
    % Set the parameters for the initial iteration
    PARAM=zeros(ITER,N);
    PARAM(i,:)=PARAM0;

    % Initialize the function SSD VALUE 
    SSD=zeros(ITER,1);

    % This second output variable is simply to store the most recent output
    % of the output function call. 
    currentoutput=0; 
    
    switch(ALGORITHM)
        % Use the Levenberg–Marquardt algorithm (a combination of 
        % gradient descent and newtons method since it is switching 
        % between gradient descent and newtons method)
        % Setup optimization options for the inbuilt MATLAB
        % iteration mechanism.
        % A great explanation of this algorithm is provided by:
        %   Bun - 2009 - Applications of the Levenberg-Marquardt algorithm 
        %                to the inverse problem
        % See the following for a method to implement lagrange for equality
        % and inequality constraints
        %   Dong - 2006 - Methods for Constrained Optimization
        case {'levenberg-marquardt'}
            % See bottom of this main function for a definition of the
            % outputfunction function
            outfun=@(x,o,s) outputfunction(x,o,s);

            % Create the input function - see the bottom of this
            % function file to see a definition for this wrapper
            % function.
            infun=@(x) inputfunction(x);
    
            opt=optimset('Disp','iter', ...
                         'MaxIter',ITER, ...
                         'TolX',precision, ...
                         'OutputFcn',outfun, ...
                         'TolFun',precision, ...
                         'Algorithm','Levenberg-Marquardt', ... 
                         'ScaleProblem','Jacobian', ...
                         'ObjectiveLimit',5, ...
                         'MaxFunEvals',	uint64(1/(precision)^2));
                     
            % Run the matlab inbuilt minimization algorithm
            fminsearch(infun,PARAM0,opt);
            
        % Use a constrained optimization approach. Note that there will be
        % a number of difficulties with this approach for highly non-linear
        % and potentially non-smooth functions. 
        case {'interior-point','sqp'}
                % See bottom of this main function for a definition of the
                % outputfunction function
                outfun=@(x,o,s) outputfunction(x,o,s);
                 
                % Create the input function - see the bottom of this
                % function file to see a definition for this wrapper
                % function.
                infun=@(x) inputfunction(x);
                
                % Create the constraint function - see the bottom of this
                % function file to see a definition for this wrapper
                % function.
                cfun=@(x) constraintfunction(x);
                
                % Setup optimization options for the inbuilt MATLAB
                % iteration mechanism.
                opt=optimset('Disp','iter', ...
                             'MaxIter',ITER, ...
                             'TolX',precision, ...
                             'OutputFcn',outfun, ...
                             'TolFun',precision, ...
                             'Algorithm',ALGORITHM, ...
                             'ScaleProblem','obj-and-constr', ...
                             'AlwaysHonorConstraints','none', ...
                             'MaxFunEvals',	uint64(1/(precision)^2));
                
                % Generate boundary and linear constraints if any.
                [c, ceq, A, B, Aeq, Beq, lb, ub]=confun(PARAM0);
                
                % Run the matlab inbuilt minimization algorithm
                fmincon(infun, PARAM0, A, B, Aeq, Beq, lb, ub, cfun, opt);
                
            case 'genetic'
                % See bottom of this main function for a definition of the
                % outputfunction function
                outfun=@(x,o,s) gaoutputfunction(x,o,s);
                 
                % Create the input function - see the bottom of this
                % function file to see a definition for this wrapper
                % function.
                infun=@(x) -log(inputfunction(x));
                
                % Create the constraint function - see the bottom of this
                % function file to see a definition for this wrapper
                % function.
                cfun=@(x) constraintfunction(x);
                
                opt = gaoptimset('Display','iter', ...
                                 'Generations',ITER, ...
                                 'StallGenLimit', 10, ... 
                                 'PopulationSize',400, ...
                                 'EliteCount',0.2*400, ...
                                 'OutputFcns',outfun);
                             
                % Generate boundary and linear constraints if any.                             
                [c, ceq, A, B, Aeq, Beq, lb, ub]=confun(PARAM0);                             
                
                % Run the matlab inbuilt minimization algorithm
                % Note that when introducing non-linear constraints the
                % time taken increases significantly since for each
                % generation 10's of sub-generational solutions are
                % attempted. Only use non-linear constraints if absolutely
                % necessary. 
                ga(infun,N,A,B,Aeq,Beq,lb,ub,[],opt);
                
        % Throw an error
        otherwise
            error('No known optimization algorithm selected');
    end
    
    % Create an output function to store iteration values for
    % output by our main function if the MATLAB optimization algorithms are
    % being used.
    function stop = outputfunction(param,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
            i=optimvalues.iteration+1;
            PARAM(i,:)=param; 
            SSD(i)=optimvalues.fval;
            OUT(i,:)=currentoutput;
        elseif isequal(state,'done')
            i=optimvalues.iteration+1;
            PARAM=PARAM(1:i,:);
            SSD=SSD(1:i);
            OUT=OUT(1:i,:);
        end
    end

    % Create an input function that wraps the objfun into a form useable by
    % MATLABs built in function call
    function ssd = inputfunction(param)
        [~, ssd, currentoutput]=objfun(DATA,param);
    end

    % Create an output function to store iteration values for
    % output by our main function if the MATLAB optimization algorithms are
    % being used.
    function [state,options,optchanged] = gaoutputfunction(options,state,flag)
        optchanged = false;
        
        if isequal(flag,'init')
            % This will occur when initializing the population
        elseif isequal(flag,'interrupt')
            % This computation will occur multiple times with generations
            % generated if the optimization has non-linear constraints
            % implemented.
        elseif isequal(flag,'iter')
            i=state.Generation+1;
            % Choose the highest fitness member from score 
            % We need to select the fittest member of the population at 
            % the particular generation. This is done by looking at 
            % Best(i) and finding the index in Score. The score index 
            % is then used to retrieve the right population member in 
            % Population    
            [val,ind]=max(state.Score);
            SSD(i)=val;    
            PARAM(i,:)=state.Population(ind,:);    
        elseif isequal(flag,'done')
            PARAM=PARAM(2:i,:);
            SSD=SSD(2:i);
        end
    end

    % Create a constraint function that wraps the confun into a form
    % useable by MATLABs built in function call
    function [c, ceq] = constraintfunction(param)
        [c,ceq,~,~,~,~,~,~]=confun(param);
    end
end