function var_out = fn_lxf_shw()
global nq xq wq

i_exponent = 1;
i_interpolant = 1;
i_scheme =1;
cell_cell_arr_shw = {};
scheme_arr = {'SHW_LXF'};


p = 0;   % Sets the lower domain boundary
q = 100; % Sets the upper domain boundary
g = 9.81; % Defines gravity constant
s = 0.8; % Safety Factor for timestep calculation
runningtime = 10;  % Set time simulation runs for

showplot = 1;
% gauss();
Choice=1;
%========== Sets constants (same for all schemes)==========================

A = 1;
B = 1;      % DO NOT CHANGE THESE VALUES
C = 0;

maxit_arr = [1:4]; % max refinement iterations

for m = 1:length(maxit_arr)
%     dx = 2^(-(maxit_arr(m)+3)); % mesh size
% if m == 1
%     N=501 % Defines the no. of grid points
% else 
%     N = 500*(2*(m-1)) + 1 
% end
%     dx = (q - p)/(N - 1); % Calculates the grid spacing between points
    dx = 100*(2^(-(maxit_arr(m)+7))); % mesh size
    iterations = 0; % Sets initial iterations to 0
    pausetime = 0.1;  % Set time inverval between simulation
    leftboundary = 0;    % Initiates left solid boundary condition if set to 1
    rightboundary = 0;   % Initiates right solid boundary condition if set to 1
    x = p : dx : q ;  % Calculates initial values of x using grid spacing
    N= length(x);
    damwall = 30;
    %==========================================================================
    % for i = 1 : N
    %     h(i) = Hinitial(x(i));        % Sets initial flow depth
    % end
    h= zeros(size(x));
%     h = .1+.1*exp(-.1*(x-50).^2);
    h(x<damwall) =.2;
    h(x>=damwall) = .1;
    
    v= zeros(size(x));
    U1(1:N+2) = zeros;   % Defines array size of future grid points
    U2(1:N+2) = zeros;

    [h, v] = boundary(h,v, leftboundary, rightboundary,N,Choice);  % Sets boundary values for height and velocity
    [un1, un2, F1, F2] = dependent(h,v,g);        % Sets initial dependent variables




    t=0;
    it=1;
    L2L2R = 0; %L2L2 accumulation of ||R||^2
                        L2L2R_arr = zeros(size(x));

    dist_x_pl = dx;
    dist_x_min = dx;
    % arrays used for EOC for bound and error
    bound_arr_eoc = [];%zeros(1,ceil(T/dt));
    error_arr_eoc = [];% zeros(1,ceil(T/dt));
    time_arr = [];% zeros(1,ceil(T/dt));
    dofs_arr = [];% zeros(1,ceil(T/dt));
    EI_index = [];
    dt = .1*dx;
    %=======================  Runs Scheme  ====================================
    while t <= (runningtime-dt/2) % Halts the program when condition is violated
        dt =.1*dx;% timestepinequality(dx,s,v,g,h);       % Calculates timestep dt
        % Calculates time program runs for
        iterations = iterations + 1;               % Adds 1 to the iteration no.
        c = dt/dx;    % Calculates this here to save computation
        
        hold = [h(2:end-1);h(2:end-1).*v(2:end-1)];
        f_h_old = [0.5 * c * (F1(3:end) - F1(1:end-2));0.5 * c *(F2(3:end) - F2(1:end-2))];

        U1(2:end-1) = (A*un1(1:end-2) + B*un1(3:end))/(A+B+C) - 0.5 * c * (F1(3:end) - F1(1:end-2));
        U2(2:end-1)= (A*un2(1:end-2)+ B*un2(3:end))/(A+B+C) - 0.5 * c *(F2(3:end) - F2(1:end-2));
        

        
        if dt < 10^-5  % Terminates loop is scheme becomes unstablej+1
            t = runningtime;
            fprintf('\n\nScheme has become unstable. Simulation terminated.')
        end
        
        [h,v] = redependent(U1,U2,N, leftboundary, rightboundary,Choice);     % Recalculates dependent variables
        

        timestore(iterations) = t; % Stores time at each iteration
        timestep(iterations) = dt; % Stores time step at each iteration
        
        [un1, un2, F1, F2] = dependent(h,v,g);   % Recalculates flux vectors and
        [h,v] = redependent(U1,U2,N, leftboundary, rightboundary,Choice);     % Recalculates dependent variables
        
        hnew = [h(2:end-1); h(2:end-1).*v(2:end-1)];
        f_h_new =  [0.5 * c * (F1(3:end) - F1(1:end-2)); 0.5 * c *(F2(3:end) - F2(1:end-2))];
        
         
        for iq = 1 : nq
            tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
            [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old, IUh] = compute_Rs_vector_temp_3_spatiotemp_3_shw_semidiscrete(x,dist_x_pl,dist_x_min,hold,hnew,tq(iq),t,dt,f_h_old,f_h_new,1,1); %compute R at gauss point
            L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
            L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
           
            
            if (L2L2R<0)
                disp('L2L2R<0')
            end
        end
        
        if showplot ==1
            subplot(2,1,1)
            plot(x, v(2:end-1),'b')  % Plots water height aproximations for each timestep.
            axis([p q min(v) 1.1*max(v)])
            
            %=============== Sets titles and labels depending on choice ===========
            title('Water height for Lax-Friedrichs Scheme','Fontsize',12)
            
            xlabel('x [m]','Fontsize', 12)
            ylabel('h [m]', 'Fontsize', 12)
            
            subplot(2,1,2)
            plot(x, IUh,'b')   % Plots water velocity aproximations for each timestep.
            axis([p q min(IUh) 1.1*max(IUh)])
            %set(gcf, 'Position',  [100, 100, 800, 400])
            
            %=============== Sets titles and labels depending on choice ===========
            title('Recon for Lax-Friedrichs Scheme','Fontsize',12)
            
            xlabel('x [m]','Fontsize', 12)
            ylabel('v [m/s]', 'Fontsize', 12)
            
            
            
            pause(pausetime)   % Shows results for each timestep.
        end
        
        time_arr(it) = it*dt;
        bound_arr_eoc(it) = sqrt(L2L2R);
        % error_arr_eoc(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex));
        error_arr_eoc(it) = 1;%sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,it*dt,ex,c_0_coeff_arr_new));
        dofs_arr(it) = length(x);
        EI_index(it) = bound_arr_eoc(it)./error_arr_eoc(it);
        it = it + 1;
        t = t + dt  ;

        % and matrix of dependent variables
    end
    it
    cell_cell_arr_shw{m}=[time_arr;bound_arr_eoc;error_arr_eoc;EI_index;dofs_arr];

end
    save([scheme_arr{i_scheme},'_cell_arr_file_shw_dam_break.mat'],'cell_cell_arr_shw')

%=========================== Terminates Scheme ============================
pause(5)  % Pauses for 5 seconds before time step graph appears
clf
plot(timestore,timestep,'r')  % Plots time step results
%=============== Sets titles and labels depending on choice ===========
title('Values for dt obtained from heuristic calculation - Lax-Friedrichs','fontsize',14)

xlabel('Time [seconds]','fontsize',18)
ylabel('dt [seconds]','fontsize',18)
axis([0 t 0 max(timestep)])
%====================  Functions of initial conditions ====================
function depth = Hinitial(x)
depth = 0.2;   % <<<<<<<< This is the depth of the Dam (can be changed)
damwall = 30;  % <<< Defines the location of the breached Dam wall
% (do not select a number less than p or larger than q)
if x >= damwall
    depth = 0.1;   % Defines the initial depth of the floor water (do not
end                % use zero).
end
function velocity = Vinitial(x,N)  % Defines initial velocity
for j = 1 : N
    velocity (j) = 0; % <<<< Sets initial velocity to zero (can be changed).
end
end
%======== Function of initial independent and dependent variables =========
function[un1, un2, F1, F2] = dependent(h,v,g)
un1 = h;                 % Function to calculate initial
un2 = h.*v;              % dependent variables.
% un2= v;
F1 = h.*v;
F2 = h.*v.^2 + (g.*h.^2)/2;
% F2 = .5*v.^2 + g*h;

end
%=========== Functions to recalculate height and velocities================
function[h,v] = redependent(U1,U2,N, leftboundary, rightboundary,Choice)
h = U1(2:N+1);             % Calculates new h and v values for Lax-F
v = U2(2:N+1)./U1(2:N+1);  % (Old ghost values are deleted)
% v = U2(2:N+1);  % (Old ghost values are deleted)

[h, v] = boundary(h,v, leftboundary, rightboundary,N,Choice);  % Recalculates boundary values
end
%=========== Function to calculate timestep ==============================
function [dt, max1, max2, maxoverall] = timestepinequality(dx,s,v,g,h)
max1 = max(abs(v + sqrt(g*h)));
max2 = max(abs(v - sqrt(g*h)));
maxoverall = max(max1,max2);      % Calculates maximum wave speed
dt = s * dx / maxoverall;   % Calculates timestep
end
%============= Function to calculate Boundary conditions ==================
function [h v] = boundary(h,v, leftboundary, rightboundary,N,Choice)
% Sets boundary conditions (either zero gradient or reflective)
if leftboundary == 0 && rightboundary == 0
    h = [h(1) h h(N)];
    v = [v(1) v v(N)];
elseif leftboundary ==1 && rightboundary == 0
    h = [h(1) h h(N)];
    v = [-v(3) v v(N)];
elseif leftboundary == 0 && rightboundary == 1
    h = [h(1) h h(N)];
    v = [v(1) v -v(N-2)];
elseif leftboundary == 1 && rightboundary == 1
    h = [h(1) h h(N)];
    v = [-v(3) v -v(N-2)];
end

end

var_out =1;
end