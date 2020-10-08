%% Toolbox test script #3
% This example is described in the following paper:
% V. Sinyakov and A. Girard. Formal Synthesis from Control Programs.
% Conference on Decision and Control (CDC'20), Jeju Island, Republic of Korea, December 2020.
%% System abstraction
M = 1370;
f_0 = 51.0709;
f_1 = 0.3494;
f_2 = 0.4161;
r = 2;
T0 = 0.25;
d_car = 2;
v0_max = 30;
v_max = 20;
v_takeover = r/T0;

alpha = @(u, v) u - (f_0 + f_1*v + f_2*v^2)/M;
beta = @(u, v, h) v*(u == h | v <= v_takeover) + sqrt((v*T0)^2 - r^2)/T0*(u ~= h & v > v_takeover);
chi = @(v, vmin, vmax) min(max(v, vmin), vmax);

n_x = 4;
Dynamics = @(x, y, u, w) [x(1) + (beta(u(2), x(2), x(4)) - y(3))*T0; ...
    chi(x(2) + alpha(u(1), x(3))*T0*(u(2) == x(4) | x(2) <= v_takeover), 0, v0_max); ...
    chi(y(3) + alpha(w(1), y(3))*T0, 0, v_max); u(2)*(x(2) > v_takeover) + x(4)*(x(2) <= v_takeover)];
F = @(x, y, u, w) Dynamics(x, y, u, w);

Partition = cell(n_x, 1);

tmp = [-Inf, linspace(-100, 40, 71)]';
Partition{1} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{1}(1, 1) = -Inf;

tmp = [0, linspace(0, 30, 31)]';
Partition{2} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{2}(1, 1) = 0;

tmp = linspace(0, v_max, 21)';
Partition{3} = [tmp(1:end-1)+eps(tmp(1:end-1)), tmp(2:end)];
Partition{3}(1, 1) = 0;

Partition{4} = [0.5, 0.5; 1.5, 1.5];

Controls = cell(2, 1);
Controls{1} = [-20, 0, 10];
Controls{2} = [0.5, 1.5];
Disturbance = [-20, 10];

tic;
[TS, TransitionNumber] = AbstractSystemDiscrete(F, Partition, Controls, Disturbance);
toc;

%% Reset system abstraction
Reset = @(x, y, u, w) x + [w(1); 0; -x(3)+w(2); 0]*heaviside(x(1))*heaviside(1-x(4));

Controls_reset = cell(1, 1);
Controls_reset{1} = 0;

Disturbance_reset = [-Inf, -40; 0, v_max];

tic;
[TS_reset, TransitionNumber_reset] = AbstractSystemDiscrete(Reset, Partition, Controls_reset, Disturbance_reset);
toc;

%% Tasks definition and abstraction
delta = 3;
c = 0.1;
l2 = 50;

Spec = cell(3, 1);
Spec{1}.isMonotone = true;
Spec{1}.Inputs = [5, 25];
N_inputs = length(Spec{1}.Inputs);
Np = 2;
Spec{1}.W = cell(Np, Np, N_inputs);
Spec{1}.D = cell(Np, Np, N_inputs);
for i = 1:N_inputs
    Spec{1}.W{1, 1, i} = @(x, x_plus) (-x_plus(2)+x(2)+c)*(Spec{1, 1}.Inputs(i)-x(2)>=delta) + (x_plus(2)-x(2)+c)*(x(2)-Spec{1, 1}.Inputs(i)>=delta);
    Spec{1}.D{1, 1, i} = @(x) max(x(1)+l2, x(4)-1);
    Spec{1}.W{1, 2, i} = @(x, x_plus) (x_plus(2)-x(2)+c)*(x(2)-Spec{1, 1}.Inputs(i)>=delta);
    Spec{1}.D{1, 2, i} = @(x) max([x(1)+d_car, -x(1)-l2, x(4)-1]);
    Spec{1}.W{2, 1, i} = @(x, x_plus) (-x_plus(2)+x(2)+c)*(Spec{1, 1}.Inputs(i)-x(2)>=delta) + (x_plus(2)-x(2)+c)*(x(2)-Spec{1, 1}.Inputs(i)>=delta);
    Spec{1}.D{2, 1, i} = @(x) max(x(1)+l2, x(4)-1);
    Spec{1}.W{2, 2, i} = @(x, x_plus) (x_plus(2)-x(2)+c)*(x(2)-Spec{1, 1}.Inputs(i)>=delta);
    Spec{1}.D{2, 2, i} = @(x) max([x(1)+d_car, -x(1)-l2, x(4)-1]);
end
Spec{1, 1}.TaskType = 'safety';
Spec{1, 1}.X1 = [];
dims1 = [71; 31; 20; 2; 2; 2];
for linearIndexX = 1:prod(dims1)
    i = Ind2Sub(dims1, linearIndexX);
    if Partition{1}(i(1), 2) <= -d_car && Partition{2}(i(2), 2) <= Spec{1}.Inputs(i(6)) && Partition{4}(i(4), 2) <= 1
        Spec{1}.X1 = [Spec{1}.X1; linearIndexX];
    end
end

Spec{2}.isMonotone = true;
Spec{2}.Inputs = 0;
N_inputs = length(Spec{2}.Inputs);
Np = 3;
Spec{2}.W = cell(Np, Np, N_inputs);
Spec{2}.D = cell(Np, Np, N_inputs);
for i = 1:N_inputs
    for p1 = 1:Np
        for p2 = 1:Np
            Spec{2}.W{p1, p2, i} = @(x, x_plus) -1;
            Spec{2}.D{p1, p2, i} = @(x) 1;
            if p1 == 1 && p2 == 1
                Spec{2}.D{p1, p2, i} = @(x) x(1) + d_car;
            elseif p1 == 1 && p2 == 2
                Spec{2}.D{p1, p2, i} = @(x) max(x(1) + d_car, 1 - x(4));
            elseif p1 == 2 && p2 == 2
                Spec{2}.D{p1, p2, i} = @(x) 1 - x(4);
            elseif p1 == 2 && p2 == 3
                Spec{2}.D{p1, p2, i} = @(x) max([-x(1) + d_car, x(1) - 40, 1 - x(4)]);
            elseif p1 == 3 && p2 == 3
                Spec{2}.D{p1, p2, i} = @(x) -x(1) + d_car;
            end
        end
    end
end
Spec{2}.TaskType = 'reachability';
Spec{2}.X1 = [];
dims2 = [71; 31; 20; 2; 3];
for linearIndexX = 1:prod(dims2)
    i = Ind2Sub(dims2, linearIndexX);
	if Partition{1}(i(1), 1) >= d_car && Partition{1}(i(1), 2) <= 40 && Partition{4}(i(4), 2) <= 1 && i(5) == 3
        Spec{2}.X1 = [Spec{2}.X1; linearIndexX];
    end
end

Spec{3}.isMonotone = true;
Spec{3}.Inputs = 0;
N_inputs = length(Spec{3}.Inputs);
Np = 1;
Spec{3}.W = cell(Np, Np, N_inputs);
Spec{3}.D = cell(Np, Np, N_inputs);
Spec{3}.W{1, 1, 1} = @(x, x_plus) -1;
Spec{3}.D{1, 1, 1} = @(x) -1;
Spec{3}.TaskType = 'reachability';
Spec{3}.X1 = [];
dims3 = [71; 31; 20; 2];
for linearIndexX = 1:prod(dims3)
    i = Ind2Sub(dims3, linearIndexX);
    if Partition{1}(i(1), 1) <= -d_car && Partition{4}(i(4), 2) <= 1
        Spec{3}.X1 = [Spec{3}.X1; linearIndexX];
    end
end

TS_with_Spec = cell(3, 1);
tic;
TS_with_Spec{1} = AbstractSpecificationDiscrete(Spec{1}, Partition, TS, TransitionNumber);
TS_with_Spec{2} = AbstractSpecificationDiscrete(Spec{2}, Partition, TS, TransitionNumber);
TS_with_Spec{3} = AbstractSpecificationDiscrete(Spec{3}, Partition, TS_reset, TransitionNumber_reset);
toc;

%% Scheduler and terminal set definition
paren = @(x, varargin) x(varargin{:});

R = cell(3, 3);
R{1, 1} = @(x) prod(dims1)+1;
R{1, 2} = @(x) Sub2Ind(dims2, [paren(Ind2Sub(dims1, x), 1:4); 1]);
R{1, 3} = @(x) prod(dims3)+1;
R{2, 1} = @(x) prod(dims1)+1;
R{2, 2} = @(x) prod(dims2)+1;
R{2, 3} = @(x) Sub2Ind(dims3, paren(Ind2Sub(dims2, x), 1:4));
R{3, 1} = @(x) Sub2Ind(dims1, [paren(Ind2Sub(dims3, x), 1:4); 1; 1]);
R{3, 2} = @(x) prod(dims2)+1;
R{3, 3} = @(x) prod(dims3)+1;

W0 = cell(3, 1);

%% Controller computation using maximal controller
tic;
[Vm, Cm] = ControlProgram(TS_with_Spec, Spec, R, W0, 'maximal');
toc;

%% Controller computation using anytime controller
tic;
[Va, Ca] = ControlProgram(TS_with_Spec, Spec, R, W0, 'anytime');
toc;

%% Let us check if the two algorithms produce the same result
for i = 1:3
    display(sum(Vm{i} ~= Va{i}));
end

%% Controllable sets visualization
[XX, YY, ZZ] = meshgrid(linspace(Partition{2}(2), Partition{2}(end)-1, 30)+1/2, ...
    linspace(Partition{1}(2), Partition{1}(end)-1, 70)+1/2, ...
    linspace(Partition{3}(1), Partition{3}(end)-1, 20)+1/2);

VV = reshape(Vm{2}(1:(end-1)), [71, 31, 20, 2, 3]);
VV = VV(2:end, 2:end, :, 1, 3);
VV(VV == Inf) = 10000;
figure;
h = patch(isosurface(XX, YY, ZZ, VV, 2000));
set(h, 'EdgeColor', 'None');
set(h, 'FaceColor', [1, 0, 0]);
set(h, 'EdgeLighting', 'phong');
set(h, 'FaceLighting', 'phong');
set(h, 'FaceAlpha', 1);
set(h, 'DiffuseStrength', 0.6);
set(h, 'SpecularStrength', 0.9);
VV(1, :, :) = 2000+1e-4;
VV(:, 1, :) = 2000+1e-4;
VV(:, end, :) = 2000+1e-4;
VV(:, :, 1) = 2000+1e-4;
VV(:, :, end) = 2000+1e-4;
h1 = patch(isosurface(XX, YY, ZZ, VV, 2000));
set(h1, 'EdgeColor', 'None');
set(h1, 'FaceColor', [1, 0, 0]);
set(h1, 'EdgeLighting', 'phong');
set(h1, 'FaceLighting', 'phong');
set(h1, 'FaceAlpha', 0.1);
set(h1, 'DiffuseStrength', 0.6);
set(h1, 'SpecularStrength', 0.9);
xlabel('v^0');
ylabel('d^1');
zlabel('v^1');
grid on;
legend('Controllable set');
xlim([Partition{2}(1), Partition{2}(end)]);
ylim([Partition{1}(42), Partition{1}(end)]);
zlim([Partition{3}(1), Partition{3}(end)]);
view([0, 42]);
camlight('right');
view([53, 43]);

%% Now let us compare with the respective reachable set of Task 2 (takeover task)
tic;
[V_task2, C_task2] = ComputeControllerReachability(TS_with_Spec{2}, Spec{2}.X1);
toc;

%% Task 2 reachable set visualization
[XX, YY, ZZ] = meshgrid(linspace(Partition{2}(2), Partition{2}(end)-1, 30)+1/2, ...
    linspace(Partition{1}(2), Partition{1}(end)-1, 70)+1/2, ...
    linspace(Partition{3}(1), Partition{3}(end)-1, 20)+1/2);

VV = reshape(V_task2(1:(end-1)), [71, 31, 20, 2, 3]);
VV = VV(2:end, 2:end, :, 1, 3);
VV(VV == Inf) = 10000;
figure;
h = patch(isosurface(XX, YY, ZZ, VV, 2000));
set(h, 'EdgeColor', 'None');
set(h, 'FaceColor', [1, 0, 0]);
set(h, 'EdgeLighting', 'phong');
set(h, 'FaceLighting', 'phong');
set(h, 'FaceAlpha', 1);
set(h, 'DiffuseStrength', 0.6);
set(h, 'SpecularStrength', 0.9);
VV(1, :, :) = 2000+1e-4;
VV(:, 1, :) = 2000+1e-4;
VV(:, end, :) = 2000+1e-4;
VV(:, :, 1) = 2000+1e-4;
VV(:, :, end) = 2000+1e-4;
h1 = patch(isosurface(XX, YY, ZZ, VV, 2000));
set(h1, 'EdgeColor', 'None');
set(h1, 'FaceColor', [1, 0, 0]);
set(h1, 'EdgeLighting', 'phong');
set(h1, 'FaceLighting', 'phong');
set(h1, 'FaceAlpha', 0.1);
set(h1, 'DiffuseStrength', 0.6);
set(h1, 'SpecularStrength', 0.9);
xlabel('v^0');
ylabel('d^1');
zlabel('v^1');
grid on;
legend('Controllable set');
xlim([Partition{2}(1), Partition{2}(end)]);
ylim([Partition{1}(42), Partition{1}(end)]);
zlim([Partition{3}(1), Partition{3}(end)]);
view([0, 42]);
camlight('right');
view([53, 43]);


%% Trajectory simulation
rng(14042020);
traj_length = 600;
z1 = zeros(6, traj_length+1);
z1(5:6, 1) = [1; 2];

x = zeros(4, traj_length+1);
x(:, 1) = [-60; 20; 15; 0.5];
task = ones(1, traj_length);
u_lin = zeros(1, traj_length);
u = zeros(2, traj_length);

dims = cell(2, 1);
dims{1} = dims1;
dims{2} = dims2;

x2 = zeros(3, traj_length+1);
x2(:, 1) = [-200; 15; 20];
z2 = zeros(3, traj_length+1);
for i = 1:traj_length
    for s = 1:4
        z1(s, i) = find((Partition{s}(:, 1) <= x(s, i)) & (Partition{s}(:, 2) >= x(s, i)), 1, 'first');
    end
    if task(i) == 1
        z_lin1 = Sub2Ind(dims{1}, z1(1:6, i));
    else
        z_lin1 = Sub2Ind(dims{2}, z1(1:5, i));
    end
    
    if ismember(z_lin1, Spec{task(i)}.X1)
        if task(i) == 2
            x(:, i) = Reset(x(:, i), x(:, i), 0, [x2(1, i); x2(3, i)]);
            x2(:, i) = [-150+100*rand; x(3, i); v_max*rand];
            task(i+1) = 1;
            for s = 1:4
                z1(s, i) = find((Partition{s}(:, 1) <= x(s, i)) & (Partition{s}(:, 2) >= x(s, i)), 1, 'first');
            end
            z1(5, i) = 1;
            z_lin1 = Sub2Ind(dims{1}, z1(1:6, i));
        elseif Vm{2}(Sub2Ind(dims{2}, [z1(1:4, i); 1])) <= (-40 - d_car - x2(1, i)) / (v_max*T0)
            task(i+1) = 2;
            
            z1(5, i) = 1;
            z_lin1 = Sub2Ind(dims{2}, z1(1:5, i));
        else
            task(i+1) = task(i);
        end
    else
        task(i+1) = task(i);
    end
    
    u_lin(i) = find(min(max(Cm{task(i+1), 1}(z_lin1, :, :, :), [], 4), [], 3), 1, 'last');
    u_sub = Ind2Sub([3; 2], u_lin(i));
    u(1, i) = Controls{1}(u_sub(1));
    u(2, i) = Controls{2}(u_sub(2));
    
    x(:, i+1) = F(x(:, i), x(:, i), u(:, i), -10+20*rand);
    if task(i+1) == 1
        if x(1, i+1) + l2 <= 0
            z1(5, i+1) = 1;
        else
            z1(5, i+1) = 2;
        end
    else
        if x(1, i+1) >= d_car && x(4, i+1) == 0.5
            z1(5, i+1) = 3;
        elseif x(4, i+1) == 1.5
            z1(5, i+1) = 2;
        else
            z1(5, i+1) = 1;
        end
    end
    if i <= 300
        z1(6, i+1) = 2;
    else
        z1(6, i+1) = 1;
    end
    
    x2(1, i+1) = chi(x2(1, i) + (x2(2, i) - x2(3, i))*T0, -Inf, -d_car);
    x2(2, i+1) = x(3, i+1);
    x2(3, i+1) = chi(x2(3, i) + alpha(-10+20*rand, x2(3, i))*T0, 0, v_max);
end

%% Plot trajectory
tt = (0:traj_length)'*T0;

n_takeover1 = [find(diff(x(4, :)) > 0)+1, traj_length+1];
n_takeover2 = find(diff(x(4, :)) < 0);

figure;
% subplot(3, 1, 1);
plot(tt, -l2*ones(traj_length+1, 1), '-.r');
hold on
plot(tt, d_car*ones(traj_length+1, 1), '-r');
plot(tt(1:n_takeover1(1)), x(1, 1:n_takeover1(1)), 'k');
for j = 1:length(n_takeover1)-1
    plot(tt(n_takeover1(j):n_takeover2(j)), x(1, n_takeover1(j):n_takeover2(j)), 'b');
    plot(tt(n_takeover2(j):n_takeover1(j+1)), x(1, n_takeover2(j):n_takeover1(j+1)), 'k');
end
plot(tt, -d_car*ones(traj_length+1, 1), '-r');
hold off
grid on
xlabel('t')
ylabel('Distance')
legend('Safe distance', 'Car dimensions', 'First lane', 'Second lane', 'Location', 'northeast');
ylim([-200, 20]);

figure;
% subplot(3, 1, 2);
plot(tt, x(2, :), 'b');
hold on;
plot(tt, x(3, :), 'r');
plot(tt, x2(3, :), 'g');
plot(tt, Spec{1}.Inputs(z1(6, :)), 'k');
hold off;
grid on
xlabel('t')
ylabel('v')
ylim([-5, 35]);
legend('v_0', 'v_1', 'v_2', 'v^*', 'Location', 'northeast');

current_mode = z1(5, 1:end-1) + 2*(task(2:end)-1);
current_mode(find(diff(current_mode) == -3)+1) = 5;

figure;
% subplot(3, 1, 3);
plot(tt(1:end-1), current_mode, 'k');
grid on;
xlabel('t');
ylabel('Mode');
ylim([0, 6]);








