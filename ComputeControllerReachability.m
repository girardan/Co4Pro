function [V, C, domC] = ComputeControllerReachability(TS, X1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT ARGUMENTS
% TS is an (Nx by Nu by Nv by Np) array. Nx is a total number of symbolic
% states, Nu is the maximal number of control input values available at a
% symbolic state, Nv is the maximal number of successors from a symbolic
% state for a given control value, Np is the number of modes in the
% specification. TS(x, u, w, p) is the successor symbolic 
% state from state x when symbolic control u and simbolic disturbance w are
% applied while the next mode is chosen to be p.
% X1 is a column vector of symbolic terminal states.
% OUTPUT ARGUMENTS
% C is a (Nx by Nu) boolean matrix that represents the reachability
% controller. C(x, u) is true if control input u is enabled at a state x of
% the closed-loop system.
% V is a (Nx by 1) vector that represents the value function. V(x) is the
% maximal number of transitions from state x to reach X1 using controller
% C.
% domC is a (Nx by 1) boolean vector of controllable states. Terminal
% states from X1 are considered controllable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Nx = size(TS, 1);
    Nu = size(TS, 2);
    Nv = size(TS, 3);
    Np = size(TS, 4);
    fixedPointFound = false;
    
    states = (1:Nx)';
    Xu = Nx+1; % Blocking state
    TS(TS == 0) = Xu;
    
    V = Inf(Nx+1, 1);
    V(X1) = 0;
    X = setdiff((1:Nx)', X1);
    
    succ_number = size(TS(1, :, :, :));
    
    while ~fixedPointFound
        Vtmp = V;
        for x_ind = 1:length(X)
            V(X(x_ind)) = 1 + min(max(min(reshape(Vtmp(TS(X(x_ind), :, :, :)), succ_number), [], 4), [], 3), [], 2);
        end
        if Vtmp == V
            fixedPointFound = true;
        end
    end
    
    C = zeros(Nx, Nu, Nv, Np, 'logical');
    C(X1, :, :, :) = false;
    Xreach = (V < Inf) & (V > 0);
    C(Xreach, :, :, :) = 1 + Vtmp(TS(Xreach, :, :, :)) <= repmat(V(Xreach), 1, Nu, Nv, Np);
    
    domC = states(V(1:Nx) < Inf);
end


