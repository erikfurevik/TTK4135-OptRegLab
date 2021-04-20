%% Initialization and model definition
init01                                 % Load helicopter parameters

global nx N

% Continuous system
Ac = [0 1 0 0 0 0;
      0 0 -K_2 0 0 0;
      0 0 0 1 0 0;
      0 0 -K_1*K_pp -K_1*K_pd 0 0;
      0 0 0 0 0 1;
      0 0 0 0 -K_3*K_ep -K_3*K_ed];
Bc = [0 0; 0 0; 0 0; K_1*K_pp 0; 0 0; 0 K_3*K_ep];

% Discrete system (forward euler)
Dt = 0.25;
Ad = Dt*Ac + eye(width(Ac));
Bd = Dt * Bc;

% Initial value
lambda0 = pi;
x0 = [lambda0; 0; 0; 0; 0; 0];

% Number of states and inputs
nx = width(Ad);                         % Number of states
nu = width(Bd);                         % Number of inputs

% Time horizon and initialization
N  = 40;                               % Time horizon for states
% N = 60 % For optional exercise
Nx = N*nx;                           % Total state variables
Nu = N*nu;                           % Total input variables
Nz = N*nx+N*nu;                      % Total number of optim. variables
z0 = ones(Nz,1);                     % Initial guess

alpha = 0.2;
beta = 20;
lambdaT = 2*pi/3;

% Bounds
ul 	    = [-30*pi/180; -inf];            % Lower bound on control
uu 	    = [30*pi/180; inf];              % Upper bound on control

xl      = -Inf*ones(nx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(nx,1);               % Upper bound on states (no bound)
xl(3)   = ul(1);                           % Lower bound on pitch
xu(3)   = uu(1);                           % Upper bound on pitch

% For optional exercise
% maxElevRate = 0.03;
% maxTravRate = 0.4;
% xl(2) = -maxTravRate;
% xu(2) = maxTravRate;
% xl(6) = -maxElevRate;
% xu(6) = maxElevRate;


% Generate constraints on measurements and inputs (the same)
[zlb,zub] = gen_constraints(N,N,xl,xu,ul,uu);
zlb(Nz)  = 0;                    % We want the last input to be zero
zub(Nz)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights) 
Q = zeros(nx);
Q(1,1) = 1;
q = diag([1 1]);                              % Weight on input
G = gen_q(Q,q,N,N);                 % Generate G matrix
c = zeros(Nz,1);                        % There is no linear term.

% Cost function
phi = @(z) z'*G*z;

%% Generate system matrixes for linear model (these take the same shape)
Aeq = gen_aeq(Ad,Bd,N,nx,nu);             % Generate A, hint: gen_aeq
beq = zeros(size(Aeq,1),1);             % Generate b
beq(1:nx) = Ad*x0;

%% Solve QP problem with linear model
opt = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',40000);
z = fmincon(phi,z0,[],[],Aeq,beq,zlb,zub,@nonLinCon,opt);

%% Extract control inputs and states
u1 = [z(Nx+1:nu:Nz);z(Nz-1)];           % Control input 1 from solution
u2 = [z(Nx+2:nu:Nz);z(Nz)];             % Control input 2 from solution
x1 = [x0(1);z(1:nx:Nx)];              % State x1 from solution
x2 = [x0(2);z(2:nx:Nx)];              % State x2 from solution
x3 = [x0(3);z(3:nx:Nx)];              % State x3 from solution
x4 = [x0(4);z(4:nx:Nx)];              % State x4 from solution
x5 = [x0(5);z(5:nx:Nx)];              % State x5 from solution
x6 = [x0(6);z(6:nx:Nx)];              % State x6 from solution

num_variables = 5/Dt;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

% v = [5 sec padding; v; 5 sec padding]
u1  = [zero_padding; u1; zero_padding];
u2  = [zero_padding; u2; zero_padding];
x1  = [lambda0*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

% Create timeseries
t = 0:Dt:Dt*(length(u1)-1);
X_opt = [x1 x2 x3 x4 x5 x6];
U = [t' u1 u2];

%% LQR

Rlqr = 1;
Qlqr = diag([6 2 0.2 0.1 1 1]);

Klqr = dlqr(Ad,Bd,Qlqr,Rlqr);


%% Plotting

figure(3)
subplot(4,2,1),stairs(t,u1),grid,ylabel('u - pitch')
subplot(4,2,2),stairs(t,u2),grid,ylabel('u - elev')
subplot(4,2,3),plot(t,x1','-mo'),grid,ylabel('travel')
subplot(4,2,4),plot(t,x2,'-mo'),grid,ylabel('r')
subplot(4,2,5),plot(t,x3','-mo'),grid,ylabel('pitch')
subplot(4,2,6),plot(t,x4','-mo'),grid,ylabel('pdot')
subplot(4,2,7),plot(t,x5','-mo'),grid,xlabel('tid (s)'),ylabel('elevation')
subplot(4,2,8),plot(t,x6','-mo'),grid,xlabel('tid (s)'),ylabel('elev rate')

% Returns vector of constraints for each z based on non-linear eq.
function [c,ceq] = nonLinCon(z)
    alpha = 0.2;
    %alpha = 0.3; % Optional exercise
    beta = 20;
    lambdaT = 2*pi/3;
    global nx N
    c = zeros(N,1);
    for k=1:N
        lambda_k = z((k-1)*nx + 1); % first element of each state vector
        elev_k = z((k-1)*nx + 5); % fifth element of each state vector
        c(k) = alpha*exp(-beta*(lambda_k-lambdaT)^2) - elev_k;
    end
    ceq = [];
end
