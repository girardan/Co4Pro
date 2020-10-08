function [V, C] = ControlProgram(TS, Tasks, R, W0, algorithm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT ARGUMENTS
% TS is a cell array of specification-constrained abstractions. TS{i} is an
% (Nx by Nu by Nv by Np) array. Nz is a total number of symbolic states, Nu is
% the maximal number of control input values available at a symbolic state,
% Nv is the maximal number of successors from a symbolic state for a given
% control value, Np is the number of modes in Tasks{i}. TS(x, u, w, p) is
% the successor symbolic state from state x
% when a symbolic control u and a simbolic disturbance w are applied while
% the next mode is chosen to be p.
% Tasks is a (NT by 1) cell array of control tasks. Each task is a
% structure with fields TaskType (either 'reachability' or 'safety') and X1
% (a column vector of symbolic terminal states). R is a (NT by NT) cell
% array of functions that represents the scheduler. R{i0, i1}(x) is a
% column vector of symbolic states to which the system TS can jump while
% transitioning from task i0 to task i1.
% W0 is a (NT by 1) cell array of column vectors. W0{i} represents the
% subset of Tasks{i}.X1 that are terminal for the control program.
% algorithm is either 'maximal' or 'anytime'.
% OUTPUT ARGUMENTS
% C is a (NT by 1) cell array of task controllers that can be combined into
% program controller. C{i} is (Nx by Nu) boolean matrix
% that represents the respective safety or reachability controller. C(x, u) is
% true if control input u is enabled at a state x of the closed-loop
% system.
% V is a cell array of column vectors. V{i} is a (Nx by 1) vector that
% represents the value function. For a reachability task, V{i}(x) is the
% maximal number of transitions from state x to reach Tasks{i}.X1 using
% controller C{i}. For a safety task, V{i}(x) == 0 if the all the
% trajectories of the closed-loop system are infinite, V{i}(x) == 1 if 
% all the trajectories of the closed-loop system are either infinite or
% reach the terminal set Tasks{i}.X1, V{i}(x) = + Inf otherwise.
% domC is a (Nx by 1) boolean vector of controllable states. Terminal
% states from X1 are considered controllable.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NT = size(TS, 1);
C = cell(NT, 1);
V = cell(NT, 1);
domC = cell(NT, 1);

% First iteration: computing task controllers
Wi = cell(NT, 1);
for i = 1:NT
    if strcmp(algorithm, 'maximal')
        Wi{i} = Tasks{i}.X1;
    else
        Wi{i} = intersect(Tasks{i}.X1, W0{i});
    end
    
    switch Tasks{i}.TaskType
        case 'reachability'
            [V{i}, C{i}, domC{i}] = ComputeControllerReachability(TS{i}, Wi{i});
        case 'safety'
            [V{i}, C{i}, domC{i}] = ComputeControllerSafety(TS{i}, Wi{i});
    end
end

% Fixed point iteration
iteration = 0;
fixedPoint = false;
while ~fixedPoint
    iteration = iteration + 1;
    fprintf('iteration = %d.\n', iteration);
    fixedPoint = true;
    for i = 1:NT
        W = [];
        if strcmp(algorithm, 'maximal')
            states = setdiff(Wi{i}, W0{i});
        else
            states = setdiff(Tasks{i}.X1, union(Wi{i}, W0{i}));
        end
        for k = 1:length(states)
            blockingTerminalState = true;
            
            xk = states(k);
            for iplus = 1:NT
                if ~isempty(intersect(R{i, iplus}(xk), domC{iplus}))
                    blockingTerminalState = false;
                end
            end
            if (blockingTerminalState && strcmp(algorithm, 'maximal')) || ...
                    (~blockingTerminalState && strcmp(algorithm, 'anytime'))
                W = [W; xk]; %#ok<AGROW>
            end
        end
        if ~isempty(W)
            % Updating effective terminal sets
            if strcmp(algorithm, 'maximal')
                Wi{i} = setdiff(Wi{i}, W);
            else
                Wi{i} = union(Wi{i}, W);
            end
            % Calculating new task controllers
            switch Tasks{i}.TaskType
                case 'reachability'
                    [V{i}, C{i}, domC{i}] = ComputeControllerReachability(TS{i}, Wi{i});
                case 'safety'
                    [V{i}, C{i}, domC{i}] = ComputeControllerSafety(TS{i}, Wi{i});
            end
            fixedPoint = false;
        end
    end
end
end

