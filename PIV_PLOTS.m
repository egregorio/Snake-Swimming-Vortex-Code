%% Vector field data processing 
clear all 
clc  
set(0,'DefaultAxesFontSize',26);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextarrowshapeInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultPolaraxesTickLabelInterpreter','latex');
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
%% Select trial and load formatted data
trial = 'fwd2';   % <-- change to 'turn1', 'fwd2', or 'fwd6'

ddptv = load('/path/to/the/file/ddptv_data.mat'); % <-- CHANGE TO MATCH YOUR COMPUTER
vf = ddptv.data.(trial).vectorfield;

x = vf.x; X = vf.X; y = vf.y; Y = vf.Y; z = vf.z; Z = vf.Z;
u = vf.u; U = vf.U; v = vf.v; V = vf.V; w = vf.w; W = vf.W;
RawImages = ddptv.data.(trial).images;
nbimages  = size(RawImages, ndims(RawImages));

sigma = 0.6;
[X,Y,Z,z,Vel_norm,Vort_X,Vort_Y,Vort_Z,Vort_norm,Qcrit,front,top,side,ImageX,ImageY,ImageZ] = piv_functions.HydroParam(x,X,y,z,U,V,W,sigma);

%% Initialization

iso4 = ddptv.data.(trial).iso(1);
iso  = ddptv.data.(trial).iso(2);
iso6 = ddptv.data.(trial).iso(3);

rho = 1000;

% Space/Time domain
time = [1:nbimages]/14*1000;
x1 = 1; % boundaries
x2 = 74;
y1 = 1;
y2 = 74;
z1 = 1;
z2 = 43;

% Coordinates
XX = X(y1:y2,x1:x2,z1:z2);
YY = Y(y1:y2,x1:x2,z1:z2);
ZZ = Z(y1:y2,x1:x2,z1:z2);
Coord = [XX(:) YY(:) ZZ(:)]*0.001;
dx=(x(2)-x(1))*0.001; %m 
dy=(y(2)-y(1))*0.001;
dz=(z(2)-z(1))*0.001;

%% Impulse/KE computation

mi=1;
mf=nbimages;
Imp=[];
jj=0;
for m = mi:mf
    jj=jj+1;
    Wx = Vort_X(y1:y2,x1:x2,z1:z2,m); % Vorticity in selected volume
    Wy = Vort_Y(y1:y2,x1:x2,z1:z2,m);
    Wz = Vort_Z(y1:y2,x1:x2,z1:z2,m);
    
%%% Vortex isolation for Q = 5
    Vortfilt = Qcrit(y1:y2,x1:x2,z1:z2,m) > iso;
    Un = Vel_norm(y1:y2,x1:x2,z1:z2,m);
    Unfilt = Un.*Vortfilt;
    Vort = [Wx(:) Wy(:) Wz(:)];
    Vort2 = Vort.*[Vortfilt(:) Vortfilt(:) Vortfilt(:)]; % Vorticité filtré
    
%%% Vortex isolation for Q = 4
    Vortfilt4 = Qcrit(y1:y2,x1:x2,z1:z2,m) > iso4;
    Vort4 = Vort.*[Vortfilt4(:) Vortfilt4(:) Vortfilt4(:)]; % Vorticité filtré
    
%%% Vortex isolation for Q = 6
    Vortfilt6 = Qcrit(y1:y2,x1:x2,z1:z2,m) > iso6;
    Vort6 = Vort.*[Vortfilt6(:) Vortfilt6(:) Vortfilt6(:)]; % Vorticité filtré
    
%%% Vortex center position
    PosVortX(jj) = mean(nonzeros(Vortfilt(:).*XX(:)))*0.001;
    PosVortY(jj) = mean(nonzeros(Vortfilt(:).*YY(:)))*0.001;
    PosVortZ(jj) = mean(nonzeros(Vortfilt(:).*ZZ(:)))*0.001;
    Coord2 = [Coord(:,1)-PosVortX(jj) Coord(:,2)-PosVortY(jj) Coord(:,3)-PosVortZ(jj)];
    
%%% Calcul Impulsion et Ec
    Cr2 = cross(Coord2,Vort2);
    Cr4 = cross(Coord2,Vort4);
    Cr6 = cross(Coord2,Vort6);
    Impvort(:,jj) = rho/2*sum(Cr2*dx*dy*dz,'omitnan');
    Impvort4(:,jj) = rho/2*sum(Cr4*dx*dy*dz,'omitnan');
    Impvort6(:,jj) = rho/2*sum(Cr6*dx*dy*dz,'omitnan');
    KEvort(jj) = rho/2*sum(Unfilt(:).^2*dx*dy*dz,'omitnan');
end
implim = max(max([Impvort(1,:) Impvort(2,:) Impvort(3,:)]),abs(min([Impvort(1,:) Impvort(2,:) Impvort(3,:)]))); 

%% Impulse plot for full recording

if strcmp(trial, 'turn1')
    v_time1 = ddptv.data.(trial).v_time(1);
    v_time2 = ddptv.data.(trial).v_time(2);
    v_time3 = ddptv.data.(trial).v_time(3);
end


% for error bars
Posx = abs(Impvort(1,:) - Impvort6(1,:));
Negx = abs(Impvort(1,:) - Impvort4(1,:));
Posy = abs(Impvort(2,:) - Impvort4(2,:));
Negy = abs(Impvort(2,:) - Impvort6(2,:));
Posz = abs(Impvort(3,:) - Impvort4(3,:));
Negz = abs(Impvort(3,:) - Impvort6(3,:));

fullfig('color','w')
subplot(3,1,1)
hold on
if strcmp(trial, 'turn1')
    xline([v_time1, v_time2, v_time3],':','linewidth',2,'Color',[0 0 0]+0.5)
end
errorbar(time/1000,Impvort(1,:),Negx,Posx,'d-','linewidth',3,'markersize',7,'color',"#005A6A")
box on 
% grid on 
ylabel('$J_x$ (N/s)')
ylim([-0.0006 0.0006])
xlim([0.0 4.0])
set(gca,'linewidth',2)

subplot(3,1,2)
hold on
if strcmp(trial, 'turn1')
    xline([v_time1, v_time2, v_time3],':','linewidth',2,'Color',[0 0 0]+0.5)
end
errorbar(time/1000,Impvort(2,:),Negy,Posy,'d-','linewidth',3,'markersize',7,'color',"#0098B3")
box on 
% grid on 
ylabel('$J_y$ (N/s)')
ylim([-0.0006 0.0006])
xlim([0.0 4.0])
set(gca,'linewidth',2)

subplot(3,1,3)
hold on
if strcmp(trial, 'turn1')
    xline([v_time1, v_time2, v_time3],':','linewidth',2,'Color',[0 0 0]+0.5)
end
errorbar(time/1000,Impvort(3,:),Negz,Posz,'d-','linewidth',3,'markersize',7,'color',"#00C4E7")
box on 
% grid on 
xlabel('$t$ (s)')
ylabel('$J_z$ (N/s)')
ylim([-0.0006 0.0006])
xlim([0.0 4.0])
set(gca,'linewidth',2)

%% Vortex Force Interpolation Plot

vr = ddptv.data.(trial).vortex_range;
n_events = numel(vr);

event_colors = {"#DE0049", "#F59300", "#118C4D", "#5A3E96", "#00A19A", "#C77800"};

time_events = cell(1, n_events);
fx = cell(1, n_events); fy = cell(1, n_events); fz = cell(1, n_events);
gx = cell(1, n_events); gy = cell(1, n_events); gz = cell(1, n_events);

for e = 1:n_events
    t = vr{e}/14;
    time_events{e} = t;

    ppx = polyfit(t, Impvort(1, vr{e}), 2);
    ppy = polyfit(t, Impvort(2, vr{e}), 2);
    ppz = polyfit(t, Impvort(3, vr{e}), 2);

    fx{e} = polyval(ppx, t); fy{e} = polyval(ppy, t); fz{e} = polyval(ppz, t);

    px = polyder(ppx); py = polyder(ppy); pz = polyder(ppz);
    gx{e} = polyval(px, t); gy{e} = polyval(py, t); gz{e} = polyval(pz, t);
end

fullfig('color','w')
subplot(3,1,1)
hold on
for e = 1:n_events
    plot(time_events{e}, fx{e}, '-', 'linewidth',3, 'color', event_colors{e})
    % plot(time_events{e}, gx{e}, '--', 'linewidth',3, 'color', event_colors{e})
end
plot(time/1000, Impvort(1,:), 'd', 'linewidth',3, 'markersize',7, 'color', "#005A6A")
plot(time/1000, Impvort(1,:), 'd', 'linewidth',3, 'markersize',7, 'color', "#005A6A")
if strcmp(trial, 'turn1')
    xline([v_time1, v_time2, v_time3],':','linewidth',2,'Color',[0 0 0]+0.5)
end
box on
ylabel('$J_x$ (N/s)')
ylim([-0.0006 0.0006])
xlim(ddptv.data.(trial).x_lim)
set(gca,'linewidth',2)

subplot(3,1,2)
hold on
for e = 1:n_events
    plot(time_events{e}, fy{e}, '-', 'linewidth',3, 'color', event_colors{e})
    % plot(time_events{e}, gy{e}, '--', 'linewidth',3, 'color', event_colors{e})
end
plot(time/1000, Impvort(2,:), 'd', 'linewidth',3, 'markersize',7, 'color', "#0098B3")
if strcmp(trial, 'turn1')
    xline([v_time1, v_time2, v_time3],':','linewidth',2,'Color',[0 0 0]+0.5)
end
box on
ylabel('$J_y$ (N/s)')
ylim([-0.0006 0.0006])
xlim(ddptv.data.(trial).x_lim)
set(gca,'linewidth',2)

subplot(3,1,3)
hold on
for e = 1:n_events
    plot(time_events{e}, fz{e}, '-', 'linewidth',3, 'color', event_colors{e})
    % plot(time_events{e}, gz{e}, '--', 'linewidth',3, 'color', event_colors{e})
end
plot(time/1000, Impvort(3,:), 'd', 'linewidth',3, 'markersize',7, 'color', "#00C4E7")
if strcmp(trial, 'turn1')
    xline([v_time1, v_time2, v_time3],':','linewidth',2,'Color',[0 0 0]+0.5)
end
box on
xlabel('$t$ (s)')
ylabel('$J_z$ (N/s)')
ylim([-0.0006 0.0006])
xlim(ddptv.data.(trial).x_lim)
set(gca,'linewidth',2)


%% Images of vortex with snake -- m = 14:25 in Figure 2


% fil = fullfile(ExpName+"_vid_Imp.avi");
fig = fullfig('color','white'); 


for jj=14
    clf
    axes()

    hold on
    surf(ImageX,ImageY,ImageZ,'CData',4095-RawImages(:,:,jj),'FaceColor','texturemap'), alpha 1, colormap gray;

    box on
    grid on
    xlabel('$X$'), ylabel('$Y$'),zlabel('$Z$'), axis equal
    set(gca, 'color', 'none');
    axis equal
    xlim([x(x1) x(x2)]),ylim([y(y1) y(y2)]),zlim([z(z1) z(z2)])

    ax2 = axes()%'Position',[0.06 0.3 0.35 0.35]);
    [RR,GG,BB] = piv_functions.isosurfacecolor(Vort_X,jj);
    t = patch(isosurface(X,Y,Z,Qcrit(:,:,:,jj),iso));
    isonormals(X,Y,Z,Qcrit(:,:,:,jj),t)
    isocolors(X,Y,Z,RR(:,:,:,jj),GG(:,:,:,jj),BB(:,:,:,jj),t)
    t.FaceColor = 'interp';
    t.EdgeColor = 'none';

    caxis([-20 20])
    t = patch(isosurface(XX,YY,ZZ,Qcrit(y1:y2,x1:x2,z1:z2,jj),iso));
    isonormals(XX,YY,ZZ,Qcrit(y1:y2,x1:x2,z1:z2,jj),t)
    t.FaceColor = 'interp';
    t.EdgeColor = 'none';
    camlight('headlight')
    lighting gouraud
    box on
    grid on
    xlabel('$X^*$'), ylabel('$Y^*$'),zlabel('$Z^*$'), axis equal
    set(gca, 'color', 'none');
    axis equal
    xlim([x(x1) x(x2)]),ylim([y(y1) y(y2)]),zlim([z(z1) z(z2)])
    drawnow
    frame = getframe(fig);
    im = frame2im(frame); 
end


