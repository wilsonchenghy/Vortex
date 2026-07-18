function vortex_sim
% Single-EDF mono-copter -- waypoint flight: takeoff, out, circle, return, land.
%
% Inner loop : Hover LQR  (attitude + altitude + yaw, integral on z, psi)
% Outer loop : Position LQR (xy hold, integral on x, y) -> phi_ref, theta_ref
%
% State (16) : x(1:12) = nonlinear plant state
%              x(13:14)= hover integrators (int_z, int_psi)
%              x(15:16)= position integrators (int_x, int_y)

    clear; clc; close all;

    %% Parameters
    p.m   = 0.8;
    p.g   = 9.81;
    p.Jxx = 0.012; p.Jyy = 0.012; p.Jzz = 0.0015;
    p.Ir  = 8e-5;
    p.l   = 0.15;
    p.r   = 0.08;
    p.kF  = 6.0;
    p.CD  = 0.05;
    p.kf  = 2.0e-6;
    p.kq  = 1.0e-8;
    p.s   = +1;

    %% Equilibrium
    Omega0 = sqrt(p.m*p.g/p.kf);
    delta0 = p.s*p.kq*p.m*p.g/(4*p.r*p.kF*p.kf);
    h0     = p.s*p.Ir*Omega0;
    u_trim = [+delta0; -delta0; -delta0; +delta0; Omega0];

    fprintf('Omega_0 = %.1f rad/s (%.0f RPM)\n', Omega0, Omega0*60/(2*pi));
    fprintf('delta_0 = %.2f deg\n', rad2deg(delta0));
    fprintf('h_0     = %.4f kg m^2/s\n', h0);
    fprintf('Nutation freq = %.2f Hz\n\n', h0/sqrt(p.Jxx*p.Jyy)/(2*pi));

    %% Hover LQR (inner loop, 8 states + 2 integrators)
    A = zeros(8); B = zeros(8,5);
    A(1,4)=1; A(2,5)=1; A(3,6)=1; A(7,8)=1;
    A(4,5) = -h0/p.Jxx;
    A(5,4) = +h0/p.Jyy;
    B(4,[1 3]) =  p.l*p.kF/p.Jxx;
    B(5,[2 4]) = -p.l*p.kF/p.Jyy;
    B(6,1) =  p.r*p.kF/p.Jzz;  B(6,2) = -p.r*p.kF/p.Jzz;
    B(6,3) = -p.r*p.kF/p.Jzz;  B(6,4) =  p.r*p.kF/p.Jzz;
    B(6,5) = -2*p.s*p.kq*Omega0/p.Jzz;
    B(8,5) =  2*p.kf*Omega0*(1-p.CD)/p.m;

    G_hov = [0 0 0 0 0 0 1 0; 0 0 1 0 0 0 0 0];
    Aa_h  = [A zeros(8,2); G_hov zeros(2,2)];
    Ba_h  = [B; zeros(2,5)];

    Q_h = diag([ 1/0.1^2  1/0.1^2  1/1^2   ...
                 1/1^2    1/1^2    1/2^2   ...
                 1/0.25^2 1/1^2            ...
                 1/0.15^2 1/0.5^2 ]);
    R_h = diag([(1/deg2rad(10))^2*ones(1,4)  (1/(1000*2*pi/60))^2]);
    K_hov = lqr(Aa_h, Ba_h, Q_h, R_h);

    %% Position LQR (outer loop, 4 states + 2 integrators)
    A_p = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    B_p = [0 0; 0 0; 0 p.g; -p.g 0];        % [phi_ref, theta_ref]
    G_p = [1 0 0 0; 0 1 0 0];
    Aa_p = [A_p zeros(4,2); G_p zeros(2,2)];
    Ba_p = [B_p; zeros(2,2)];

    Q_p = diag([ 1/0.4^2  1/0.4^2  ...      % position error
                 1/0.8^2  1/0.8^2  ...      % velocity
                 1/1.5^2  1/1.5^2 ]);       % integral
    R_p = diag([(1/deg2rad(8))^2  (1/deg2rad(8))^2]);
    K_pos = lqr(Aa_p, Ba_p, Q_p, R_p);

    %% Initial condition: on the pad, near-zero altitude
    x0 = zeros(12,1);
    x0(6) = 0.05;            % 5 cm above ground
    xfull0 = [x0; 0; 0; 0; 0];

    %% Simulate
    T_end = 22;
    tspan = linspace(0, T_end, T_end*30);
    opts  = odeset('RelTol',1e-6,'AbsTol',1e-8);
    [t, X] = ode45(@(t,x) closed_loop(t,x,p,K_hov,K_pos,u_trim), tspan, xfull0, opts);

    %% Build reference trace for plotting
    R_ref = arrayfun(@(tt) get_reference(tt), t);
    x_ref = arrayfun(@(s) s.x, R_ref);
    y_ref = arrayfun(@(s) s.y, R_ref);
    z_ref = arrayfun(@(s) s.z, R_ref);

    %% Plot
    plot_states(t, X, x_ref, y_ref, z_ref);

    %% Animate
    animate(t, X, x_ref, y_ref, z_ref);
end

%% ============================================================
%  Waypoint sequence  --  smooth cosine ramps between phases
%   0-2 s     takeoff to 1.5 m at origin
%   2-3 s     hover
%   3-5.5 s   fly out to circle entry (+R, 0)
%   5.5-6.5 s hover on circle entry
%   6.5-13.5 s  horizontal circle (R = 1 m, centered at origin, CCW)
%   13.5-16.5 s return to origin along +x
%   16.5-17.5 s hover
%   17.5-22 s land
%% ============================================================
function ref = get_reference(t)
    z_hover = 1.5;  z_pad = 0.05;  R = 1.0;

    if     t < 2
        ref.x = 0;  ref.y = 0;
        ref.z = z_pad + (z_hover - z_pad)*sramp(t, 0, 2);
    elseif t < 3
        ref.x = 0;  ref.y = 0;  ref.z = z_hover;
    elseif t < 5.5
        s = sramp(t, 3, 5.5);
        ref.x = R*s;  ref.y = 0;  ref.z = z_hover;
    elseif t < 6.5
        ref.x = R;  ref.y = 0;  ref.z = z_hover;
    elseif t < 13.5
        theta = 2*pi*sramp(t, 6.5, 13.5);
        ref.x = R*cos(theta);
        ref.y = R*sin(theta);
        ref.z = z_hover;
    elseif t < 16.5
        s = sramp(t, 13.5, 16.5);
        ref.x = R*(1 - s);  ref.y = 0;  ref.z = z_hover;
    elseif t < 17.5
        ref.x = 0;  ref.y = 0;  ref.z = z_hover;
    else
        ref.x = 0;  ref.y = 0;
        ref.z = z_pad + (z_hover - z_pad)*(1 - sramp(t, 17.5, 22));
    end
    ref.psi = 0;
end

function s = sramp(t, t0, t1)
    if     t <= t0, s = 0;
    elseif t >= t1, s = 1;
    else,           s = 0.5*(1 - cos(pi*(t - t0)/(t1 - t0)));
    end
end

%% ============================================================
function dx = closed_loop(t, x, p, K_hov, K_pos, u_trim)
    ref = get_reference(t);

    % World-frame velocity for the position controller
    R = R_zyx(x(1), x(2), x(3));
    v_w = R * x(10:12);

    % --- Position controller (outer) -> phi_ref, theta_ref --------------
    pos_state = [x(4) - ref.x; x(5) - ref.y; v_w(1); v_w(2); x(15); x(16)];
    pos_u = -K_pos * pos_state;          % [phi_ref; theta_ref]
    phi_ref   = pos_u(1);
    theta_ref = pos_u(2);

    att_lim = deg2rad(15);
    phi_ref   = max(min(phi_ref,   att_lim), -att_lim);
    theta_ref = max(min(theta_ref, att_lim), -att_lim);

    % --- Hover controller (inner) ----------------------------------------
    xhov = [x(1); x(2); x(3); x(7); x(8); x(9); x(6); x(12); x(13); x(14)];
    xref_hov = [phi_ref; theta_ref; ref.psi; 0; 0; 0; ref.z; 0; 0; 0];
    xerr = xhov - xref_hov;
    xerr(3) = wrap_pi(xerr(3));

    u = u_trim + (-K_hov * xerr);
    u(1:4) = max(min(u(1:4), deg2rad(15)), -deg2rad(15));
    u(5)   = max(u(5), 0.4*u_trim(5));

    % --- Integrator dynamics ---------------------------------------------
    dInt_hov = [x(6) - ref.z;  wrap_pi(x(3) - ref.psi)];
    dInt_pos = [x(4) - ref.x;  x(5) - ref.y];

    % --- Plant -----------------------------------------------------------
    dx12 = monocopter_eom(x(1:12), u, p);

    % Soft ground constraint: don't allow z to go below zero with vz < 0
    if x(6) < 0 && dx12(6) < 0
        dx12(6)  = 0;
        dx12(12) = max(dx12(12), 0);
    end

    dx = [dx12; dInt_hov; dInt_pos];
end

%% ============================================================
function dx = monocopter_eom(x, u, p)
    phi=x(1); th=x(2); psi=x(3);
    wx=x(7); wy=x(8); wz=x(9);
    vx=x(10); vy=x(11); vz=x(12);
    d1=u(1); d2=u(2); d3=u(3); d4=u(4); Om=u(5);

    Ft = p.kf*Om^2; Fd = Ft*p.CD;
    F1=p.kF*d1; F2=p.kF*d2; F3=p.kF*d3; F4=p.kF*d4;
    h  = p.s*p.Ir*Om;

    tx = (F1+F3)*p.l;
    ty = -(F2+F4)*p.l;
    tz = (F1-F2-F3+F4)*p.r - p.s*p.kq*Om^2;

    wxdot = (tx + (p.Jyy-p.Jzz)*wy*wz - h*wy)/p.Jxx;
    wydot = (ty + (p.Jzz-p.Jxx)*wx*wz + h*wx)/p.Jyy;
    wzdot = (tz + (p.Jxx-p.Jyy)*wx*wy)/p.Jzz;

    cphi=cos(phi); sphi=sin(phi); cth=cos(th); sth=sin(th); tth=tan(th);
    phidot = wx + sphi*tth*wy + cphi*tth*wz;
    thdot  = cphi*wy - sphi*wz;
    psidot = (sphi*wy + cphi*wz)/cth;

    vxdot = (F2+F4)/p.m + p.g*sth;
    vydot = (F1+F3)/p.m - p.g*cth*sphi;
    vzdot = (Ft - Fd)/p.m - p.g*cth*cphi;

    R = R_zyx(phi, th, psi);
    pdot = R * [vx; vy; vz];

    dx = [phidot; thdot; psidot; pdot;
          wxdot; wydot; wzdot;
          vxdot; vydot; vzdot];
end

%% ============================================================
function R = R_zyx(phi, th, psi)
    cph=cos(phi); sph=sin(phi); cth=cos(th); sth=sin(th);
    cps=cos(psi); sps=sin(psi);
    R = [cps*cth, cps*sth*sph - sps*cph, cps*sth*cph + sps*sph;
         sps*cth, sps*sth*sph + cps*cph, sps*sth*cph - cps*sph;
        -sth,     cth*sph,               cth*cph];
end

function y = wrap_pi(a), y = mod(a+pi, 2*pi) - pi; end

%% ============================================================
function plot_states(t, X, x_ref, y_ref, z_ref)
    figure('Color','w','Position',[100 100 950 700]);
    subplot(3,1,1);
    plot(t, rad2deg(X(:,1:3)), 'LineWidth',1.2); grid on;
    ylabel('Attitude [deg]'); legend('\phi','\theta','\psi');
    title('Vortex flight: takeoff -> out -> circle -> return -> land');

    subplot(3,1,2);
    plot(t, X(:,4), 'LineWidth',1.4); hold on;
    plot(t, x_ref, '--', 'LineWidth',1);
    plot(t, X(:,5), 'LineWidth',1.2);
    plot(t, y_ref, '--', 'LineWidth',1);
    grid on; ylabel('Horizontal [m]'); legend('x','x_{ref}','y','y_{ref}');

    subplot(3,1,3);
    plot(t, X(:,6), 'LineWidth',1.4); hold on;
    plot(t, z_ref, '--', 'LineWidth',1);
    grid on; ylabel('Altitude [m]'); xlabel('Time [s]'); legend('z','z_{ref}');
end

%% ============================================================
function animate(t, X, x_ref, y_ref, z_ref)
    fig = figure('Color','w','Position',[150 80 960 720], ...
                 'MenuBar','none','ToolBar','none','Resize','off');
    ax = axes('Parent',fig,'Position',[0.08 0.08 0.88 0.86]);
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    xlim(ax,[-1.3 1.8]); ylim(ax,[-1.3 1.3]); zlim(ax,[0 2.2]);
    view(ax, 40, 18);
    xlabel(ax,'x_w [m]'); ylabel(ax,'y_w [m]'); zlabel(ax,'z_w [m]');
    title(ax,'Vortex -- Takeoff / Out / Circle / Return / Land');

    % Ground + pad
    [Xg,Yg] = meshgrid(linspace(-1.5,1.8,9), linspace(-1.5,1.5,9));
    surf(ax,Xg,Yg,zeros(size(Xg)),'FaceColor',[0.92 0.95 0.92], ...
         'EdgeColor',[0.7 0.7 0.7],'FaceAlpha',0.3);
    plot3(ax, 0, 0, 0.001, 'ks', 'MarkerSize',12, 'MarkerFaceColor',[0.6 0.9 0.6]);

    % Reference trace (full path)
    plot3(ax, x_ref, y_ref, z_ref, ':', ...
          'Color',[0.6 0.6 0.6],'LineWidth',1);

    % Vehicle geometry
    [xc,yc,zc] = cylinder(0.05, 30); zc = zc*0.30 - 0.15;
    hBody = surf(ax,xc,yc,zc,'FaceColor',[0.7 0.7 0.85], ...
                 'EdgeColor','none','FaceAlpha',0.75);
    vanes = gobjects(4,1);
    for k = 1:4
        a = (k-1)*pi/2; c=cos(a); s=sin(a);
        V = [ 0.02*c           0.02*s          -0.15;
              (0.02+0.08)*c   (0.02+0.08)*s    -0.15;
              (0.02+0.08)*c   (0.02+0.08)*s    -0.18;
              0.02*c           0.02*s          -0.18];
        vanes(k) = patch(ax,'Vertices',V,'Faces',[1 2 3 4], ...
                         'FaceColor',[0.3 0.5 0.8],'EdgeColor','k');
    end
    hThrust = plot3(ax,[0 0],[0 0],[0 -0.45],'r-','LineWidth',3);
    hAngMom = plot3(ax,[0 0],[0 0],[0  0.35],'g-','LineWidth',3);
    hXax    = plot3(ax,[0 0.15],[0 0],[0 0],'r:','LineWidth',1);
    hYax    = plot3(ax,[0 0],[0 0.15],[0 0],'b:','LineWidth',1);

    tBody = hgtransform('Parent', ax);
    set([hBody; vanes; hThrust; hAngMom; hXax; hYax], 'Parent', tBody);

    hTrail = animatedline(ax,'Color',[0 0.4 0.75],'LineWidth',1.5);
    hTxt = uicontrol('Parent',fig,'Style','text', ...
        'Units','normalized','Position',[0.01 0.82 0.22 0.16], ...
        'BackgroundColor',[1 1 1],'HorizontalAlignment','left', ...
        'FontName','Consolas','FontSize',10);

    drawnow; pause(0.3);

    out_path = fullfile(pwd, 'vortex_flight.mp4');
    vw = VideoWriter(out_path,'MPEG-4');
    vw.FrameRate = 30; vw.Quality = 95;
    open(vw);

    img1 = print(fig,'-RGBImage','-r100');
    H = size(img1,1) - mod(size(img1,1),2);
    W = size(img1,2) - mod(size(img1,2),2);

    for k = 1:length(t)
        phi=X(k,1); th=X(k,2); psi=X(k,3); pw=X(k,4:6);
        T = [R_zyx(phi,th,psi), pw(:); 0 0 0 1];
        set(tBody,'Matrix',T);
        addpoints(hTrail, pw(1), pw(2), pw(3));

        phase = flight_phase(t(k));
        set(hTxt,'String', sprintf( ...
            't = %5.2f s\nphase: %s\nz   = %5.2f m\nx   = %+5.2f m\ny   = %+5.2f m\nphi = %+5.1f deg\ntheta = %+5.1f deg', ...
            t(k), phase, pw(3), pw(1), pw(2), rad2deg(phi), rad2deg(th)));
        drawnow;

        img = print(fig,'-RGBImage','-r100');
        img = img(1:H,1:W,:);
        writeVideo(vw, img);
    end
    close(vw);
    fprintf('Saved video to: %s\n', out_path);
end

function s = flight_phase(t)
    if     t < 2,    s = 'TAKEOFF';
    elseif t < 3,    s = 'HOVER';
    elseif t < 5.5,  s = 'OUTBOUND';
    elseif t < 6.5,  s = 'HOVER';
    elseif t < 13.5, s = 'CIRCLE';
    elseif t < 16.5, s = 'RETURN';
    elseif t < 17.5, s = 'HOVER';
    else,            s = 'LAND';
    end
end