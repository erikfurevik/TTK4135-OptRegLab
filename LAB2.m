%% Initialization
init01                                 % Load helicopter parameters

% Continuous model
Ac = [0 1         0          0;
      0 0      -K_2          0;
      0 0         0          1;
      0 0 -K_1*K_pp -K_1*K_pd];

Bc = [0; 0; 0; K_1*K_pp];

% Problem constants
N  = 100;                               % Time horizon
nx = width(Ac);                         % Number of states
nu = width(Bc);                         % Number of inputs
Nx = N*nx;                              % Total state variables
Nu = N*nu;                              % Total input variables
Nz = N*nx+N*nu;                         % Total number of optim. variables

% Discrete model by forward euler
Dt = 0.25;
Ad = Dt*Ac + eye(nx); 
Bd = Dt*Bc;

% Initialization
lambda0 = pi;
x0 = [lambda0; 0; 0; 0];

%% Generate bounds
xl = -Inf*ones(nx,1);  % Lower bound (assume none)
xu = Inf*ones(nx,1);   % Upper bound (assume none)

pmax = 30*pi/180;      % Maximum allowed pitch either way
ul 	  = -pmax;         % Min pitch ref
uu 	  =  pmax;         % Max pitch ref
xl(3) = -pmax;         % Min pitch
xu(3) =  pmax;         % Max pitch

[zlb,zub] = gen_constraints(N,N,xl,xu,ul,uu);
zlb(Nz)  = 0;                    % We want the last input to be zero
zub(Nz)  = 0;                    % We want the last input to be zero

%% Objective function weights
Q = zeros(nx);
Q(1,1) = 1;                 % Travel term
q = 1;                      % Weight on input
G = 2*gen_q(Q,q,N,N);       % Generate G matrix
c = zeros(Nz,1);            % No linear term in function

%% Equality constraints
Aeq = gen_aeq(Ad,Bd,N,nx,nu);
beq = zeros(size(Aeq,1),1);
beq(1:nx) = Ad*x0;

%% Solve problem using quadprog
opt = optimset('Display','on', 'Diagnostics','off', 'LargeScale','off', 'Algorithm', 'interior-point-convex');
z = quadprog(G,c,[],[],Aeq,beq,zlb,zub,[],opt); 

%% Extract control inputs and states
padLength = 5/Dt;
zeroPad = zeros(padLength,1);
iPad = ones(padLength,1);

% v = [5s padding;   optimal solution; 5s padding]
u   = [zeroPad;      z(Nx+1:Nz); z(Nz); zeroPad];
x1  = [lambda0*iPad; x0(1);z(1:nx:Nx);  zeroPad];
x2  = [zeroPad;      x0(2);z(2:nx:Nx);  zeroPad];
x3  = [zeroPad;      x0(3);z(3:nx:Nx);  zeroPad];
x4  = [zeroPad;      x0(4);z(4:nx:Nx);  zeroPad];
t = 0:Dt:Dt*(length(u)-1);

% For use in simulink:
u = [t' u];

%% Plotting
figure(2)
subplot(511),stairs(t,u(:,2)),grid,ylabel('u')
subplot(512),plot(t,x1,'-mo'),grid,ylabel('lambda')
subplot(513),plot(t,x2,'-mo'),grid,ylabel('r')
subplot(514),plot(t,x3,'-mo'),grid,ylabel('p')
subplot(515),plot(t,x4,'-mo'),grid,ylabel('pdot')
xlabel('tid (s)')