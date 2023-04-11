clear all; close all; clc %#ok<CLALL>
%  !===================================================================
%  ! /* Description: A simple Smoothed-particle hydrodynamics (SPH) code.
%  ! The solver takes into account effects of gravity, pressure, and
%  ! particle-particle gravitational pull.
%  !    --
%  !
%  !    Record of Revisions:
%  !    Date                 Author(s)          Description of Change
%  !    =============         ========           =====================
%  !    April 9, 2023      A. Fakhreddine            Original
%  !===================================================================

% Simulation parameters
N = 100;         % Number of particles
L = 1.0;         % Simulation length
rho0 = 0.0;      % Reference density
h = 0.02;        % Smoothing length
m = 0.1;         % Particle mass
G = 0.0;         % Gravitational constant
dt = 0.0025;     % Time step
tmax = 2000*dt;  % Maximum simulation time
c0 = 1.2;        % Speed of sound

% Particle initialization
x = rand(N, 1)*L; % x-location of particles
y = rand(N, 1)*L; % y-location of particles

% Velocity vector where v(:,1) is u and v(:,2) is v component
v = zeros(N, 2);
v_half = zeros(N, 2);

rho = rho0*ones(N, 1); % Particle density
P = zeros(N, 1);       % Particle pressure contribution

% Media info
save_image_flag = 1; % Choose to save images or not
iter_save = 1;       % Saving frequency

% Main simulation loop
iter_count = 0;
for t = 0:dt:tmax
    % Update particle density and pressure
    for i=1:N
        rho(i) = m*getW(0,h);
        for j=i+1:N
            dx = x(j) - x(i);
            dy = y(j) - y(i);
            rij = [dx, dy];
            rij_norm = sqrt(dx^2+dy^2);
            Wij = getW(rij_norm/h,h);
            rho(i) = rho(i) + m*Wij;
            rho(j) = rho(j) + m*Wij;
        end
        P(i) = c0^2*(rho(i)-rho0);
    end
    x = x + dt*v_half(:, 1);
    y = y + dt*v_half(:, 2);
    
    % Compute acceleration due to gravity
    a_gravity = repmat([0, -G], N, 1);
    
    % Compute acceleration due to pressure and particle mass
    a = zeros(N, 2);
    for i = 1:N
        for j = i+1:N
            dx = x(j) - x(i);
            dy = y(j) - y(i);
            rij = [dx, dy];
            rij_norm = sqrt(dx^2+dy^2);
            dWij = getdW(rij_norm/h,h);
            
            a(i,:) = a(i,:) + m*(P(i)/(rho(i)^2) + P(j)/(rho(j)^2))*dWij*rij/rij_norm;
            a(j,:) = a(j,:) - m*(P(i)/(rho(i)^2) + P(j)/(rho(j)^2))*dWij*rij/rij_norm;
            
            a(i,:) = a(i,:) + 2.5*m*m*rij/(rij_norm^2);
            a(j,:) = a(j,:) - 2.5*m*m*rij/(rij_norm^2);
        end
        a(i,:) = a(i,:) - v(i,:)*10; % Add some damping to stabilize solution
    end
    
    % Update velocity and position using velocity-Verlet integration
    v_half = v + 0.5*dt*(a + a_gravity);
    x = x + dt*v_half(:, 1);
    y = y + dt*v_half(:, 2);
    
    for i = 1:N
        v(i, :) = v_half(i, :) + 0.5*dt*(a(i, :) + a_gravity(i, :));
    end
    
    % Simple wall BC for bottom boundary
    for i=1:N
        if(y(i) < 2*h)
            delta = 2*h - y(i);
            y(i) = y(i) + delta;
            v(i,2) = -v(i,2);
        end
    end
    
    % Plot particles
    clf;
    if(save_image_flag == 1)
        f=figure;
    end
    scatter(x, y, 50, P, 'filled');
    axis equal;
    axis([0 L 0 L]);
    colormap;
    c = colorbar;
    c.FontName = 'Times';
    title(sprintf('$Time = %.3f$', t),'interpreter','latex');
    xlabel('$x$','interpreter','latex','Fontsize',25)
    ylabel('$y$','interpreter','latex','Fontsize',25)
    set(gca,'ticklabelinterpreter','latex','Fontsize',20)
    drawnow;
    
    if(save_image_flag == 1 && mod(iter_count,iter_save) == 0)
        mkdir('Frames')
        if iter_count < 10
            saveas(f,['Frames/frame0000' num2str(iter_count) '.jpg']);
        elseif iter_count < 100
            saveas(f,['Frames/frame000' num2str(iter_count) '.jpg']);
        elseif iter_count < 1000
            saveas(f,['Frames/frame00' num2str(iter_count) '.jpg']);
        end
    end
    iter_count = iter_count + 1;
end


