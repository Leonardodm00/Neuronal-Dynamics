clear all
close all
clc
%%%
% This is an example of the usage of the PECUZAL algorithm 
% on the Roessler system.

% set parameters for standard Roessler
a = 0.2;
b = 0.2;
c = 5.7;
% Set time series length and sampling
N = 2500;
dt = 0.2;
transient_time = 2500;
tf = (N+transient_time) * dt; 
t = 0:dt:tf;
% Roessler ODE
roe = @(t,x) [-x(2)-x(3); x(1)+a*x(2); b+x(3)*(x(1)-c)];
% solve the ODE
x0 = [1, -2, 0.1]; % initial conditions  
[t, x] = ode45(roe, t, x0);
% remove transients
x(1:transient_time+1,:) = [];
t = t(1:end-transient_time-1);

y = x(:,1);  % bind y-component
mi = mutualinformation(y,50); % compute MI
figure('Units', 'normalized', 'Position', [.3 .3 .7 .5])
subplot(1,2,1)
plot(t(1:1000), y(1:1000))
title('y-component of Roessler time series')
xlabel('time')
grid on
subplot(1,2,2)
plot(mi(:,1), mi(:,2), '.-', 'MarkerSize',15)
xlabel('time delay \tau')
ylabel('mutual information [nats]')
title('Mutual information for y-component of Roessler time series')
grid on



% Now we are ready to go and call the PECUZAL algorithm pecuzal_embedding with a Theiler window determined from the
% first minimum of the mutual information shown in the above Figure, i.e. 'theiler'=7,
% and possible delays ranging in the interval [0:50].  We also use the 'econ' option for an accelerated computation.

[Y_reconstruct, tau_vals, ts_vals, Ls, E] = pecuzal_embedding(y, 0:50, 'theiler', 7, 'econ', true);


subplot(1,2,1)
plot3(Y_reconstruct(:,1), Y_reconstruct(:,2), Y_reconstruct(:,3))
xlabel(sprintf('y(t+%i)', tau_vals(1)))
ylabel(sprintf('y(t+%i)', tau_vals(2)))
zlabel(sprintf('y(t+%i)', tau_vals(3)))
title('PECUZAL reconstructed Roessler system')
grid on
view([-34.6 5.2])
subplot(1,2,2)
plot3(x(:,1), x(:,2), x(:,3))
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')
title('Original Roessler system')
grid on


clf
plot(0:length(E{1})-1, E{1}(:,1), 'LineWidth',3), hold on
scatter(tau_vals(2), E{1}(tau_vals(2)+1,1), 'k', 'filled','HandleVisibility',"off")
plot(0:length(E{2})-1, E{2}(:,1), 'LineWidth',3)
scatter(tau_vals(3), E{2}(tau_vals(3)+1,1), 'k', 'filled','HandleVisibility',"off")
plot(0:length(E{3})-1, E{3}(:,1), 'LineWidth',3)
legend('1st embedding cycle', '2nd embedding cycle', '3rd embedding cycle', 'Location', 'northeastoutside')
title('PECUZAL continuity statistics for Roessler y-component')
xlabel('delay')
ylabel('$\langle \varepsilon^\star \rangle$', 'Interpreter', 'latex')
set(gca,'LineWidth',2)
set(gca,'FontSize',18)
grid on
