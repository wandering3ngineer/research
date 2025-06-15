function [CON, CONEQ, A, B, Aeq, Beq, lb, ub]=hystereticConstraints(P)
    % Extract the parameters
    n=P(1); m=P(2); al=P(3); Xaf=P(4); c=P(5); k=P(6);

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
    lb=PLB';
    ub=PUB';  
    
    % Package the constraints for output
    CONEQ=[];
    CON=[];
    A=[];
    B=[];
    Aeq=[];
    Beq=[];
end