%% Set times
Tints1 = irf.tint('2023-04-24T02:29:11.85Z/2023-04-24T02:29:11.85Z');
Tints3 = irf.tint('2023-04-24T04:02:03.384Z/2023-04-24T04:02:04.384Z');
Tints4 = irf.tint('2023-04-24T04:15:25.0Z/2023-04-24T04:15:25.0Z');
Tints5 = irf.tint('2023-04-24T04:21:00.00Z/2023-04-24T04:21:00.00Z');

Vsn1 = 50;
Vsn3 = 60;
Vsn4 = -50;
Vsn5 = 40;

ndistframp = [-500 2000];
Tint1 = Tints1+[ndistframp(1)/Vsn1 ndistframp(2)/Vsn1];
Tint3 = Tints3+[ndistframp(1)/Vsn3 ndistframp(2)/Vsn3];
Tint4 = Tints4+[ndistframp(2)/Vsn4 ndistframp(1)/Vsn4];
Tint5 = Tints5+[ndistframp(1)/Vsn5 ndistframp(2)/Vsn5];

%% Load data
ic = 1;
iPDist1 = mms.get_data('PDi_fpi_brst_l2',Tint1,ic);
iPDist3 = mms.get_data('PDi_fpi_brst_l2',Tint3,ic);
iPDist4 = mms.get_data('PDi_fpi_brst_l2',Tint4,ic);
%iPDist5 = mms.get_data('PDi_fpi_brst_l2',Tint5,ic);

Bxyz1 = mms.get_data('B_dmpa_brst_l2',Tint1,ic);
Bxyz3 = mms.get_data('B_dmpa_brst_l2',Tint3,ic);
Bxyz4 = mms.get_data('B_dmpa_brst_l2',Tint4,ic);
%Bxyz5 = mms.get_data('B_dmpa_brst_l2',Tint5,ic);

%% Define Coordinate systems
ndsl = [0.9686    0.0056   -0.2487];
t1dsl = [-0.2263   -0.4382   -0.8700];
t2dsl = [0.2695    0.8301   -0.4882];
ndsl1 = ndsl/norm(ndsl);
t1dsl1 = t1dsl/norm(t1dsl);
t2dsl1 = t2dsl/norm(t2dsl);

ndsl = [0.8871   -0.4602   -0.0356];
t1dsl = [-0.3382   -0.5956   -0.7286];
t2dsl = [0.3142    0.6584   -0.6840];
ndsl3 = ndsl/norm(ndsl);
t1dsl3 = t1dsl/norm(t1dsl);
t2dsl3 = t2dsl/norm(t2dsl);

ndsl = [0.9489   -0.2611   -0.1774];
t1dsl = [-0.3127   -0.7001   -0.6419];
t2dsl = [0.0434    0.6646   -0.7460];
ndsl4 = ndsl/norm(ndsl);
t1dsl4 = t1dsl/norm(t1dsl);
t2dsl4 = t2dsl/norm(t2dsl);

ndsl = [0.8027   -0.5795    0.1412];
t1dsl = [-0.3078   -0.6052   -0.7341];
t2dsl = [0.5108    0.5458   -0.6642];
ndsl5 = ndsl/norm(ndsl);
t1dsl5 = t1dsl/norm(t1dsl);
t2dsl5 = t2dsl/norm(t2dsl);

%% Rotate B field
Bntt1 = irf_newxyz(Bxyz1,ndsl1,t1dsl1,t2dsl1);
Bntt3 = irf_newxyz(Bxyz3,ndsl3,t1dsl3,t2dsl3);
Bntt4 = irf_newxyz(Bxyz4,ndsl4,t1dsl4,t2dsl4);

%% Compute 1D reduced ion distributions
iPDist1.data(:,1:11,:,:) = 0;
iPDist3.data(:,1:11,:,:) = 0;
iPDist4.data(:,1:11,:,:) = 0;
%iPDist5.data(:,1:11,:,:) = 0;

vlimn = [-1800,1400];
vlimt2 = [-1600,1600];

nMC = 5e2;
vg1Dn = linspace(vlimn(1),vlimn(2),300);
vg1Dt2 = linspace(vlimt2(1),vlimt2(2),300);

nvec1 = irf.ts_vec_xyz(iPDist1.time,repmat(ndsl1,length(iPDist1.time),1));
t2vec1 = irf.ts_vec_xyz(iPDist1.time,repmat(t2dsl1,length(iPDist1.time),1));

nvec3 = irf.ts_vec_xyz(iPDist3.time,repmat(ndsl3,length(iPDist3.time),1));
t2vec3 = irf.ts_vec_xyz(iPDist3.time,repmat(t2dsl3,length(iPDist3.time),1));

nvec4 = irf.ts_vec_xyz(iPDist4.time,repmat(ndsl4,length(iPDist4.time),1));
t2vec4 = irf.ts_vec_xyz(iPDist4.time,repmat(t2dsl4,length(iPDist4.time),1));

%nvec5 = irf.ts_vec_xyz(iPDist5.time,repmat(ndsl5,length(iPDist5.time),1));
%t2vec5 = irf.ts_vec_xyz(iPDist5.time,repmat(t2dsl5,length(iPDist5.time),1));

f1Dn1 = iPDist1.reduce('1D',nvec1,'vg',vg1Dn,'nMC',nMC);
f1Dt21 = iPDist1.reduce('1D',t2vec1,'vg',vg1Dt2,'nMC',nMC);

f1Dn3 = iPDist3.reduce('1D',nvec3,'vg',vg1Dn,'nMC',nMC);
f1Dt23 = iPDist3.reduce('1D',t2vec3,'vg',vg1Dt2,'nMC',nMC);

f1Dn4 = iPDist4.reduce('1D',nvec4,'vg',vg1Dn,'nMC',nMC);
f1Dt24 = iPDist4.reduce('1D',t2vec4,'vg',vg1Dt2,'nMC',nMC);

%f1Dn5 = iPDist5.reduce('1D',nvec5,'vg',vg1Dn,'nMC',nMC);
%f1Dt25 = iPDist5.reduce('1D',t2vec5,'vg',vg1Dt2,'nMC',nMC);

%% Shift distributions to NIF
f1Dn1.depend{1,1} = f1Dn1.depend{1,1}-50;
f1Dn1.ancillary.v_edges = f1Dn1.ancillary.v_edges-50;

f1Dt21.depend{1,1} = f1Dt21.depend{1,1}+70;
f1Dt21.ancillary.v_edges = f1Dt21.ancillary.v_edges+70;

f1Dn3.depend{1,1} = f1Dn3.depend{1,1}-60;
f1Dn3.ancillary.v_edges = f1Dn3.ancillary.v_edges-60;

f1Dt23.depend{1,1} = f1Dt23.depend{1,1}+80;
f1Dt23.ancillary.v_edges = f1Dt23.ancillary.v_edges+80;

f1Dn4.depend{1,1} = f1Dn4.depend{1,1}+50;
f1Dn4.ancillary.v_edges = f1Dn4.ancillary.v_edges+50;

f1Dt24.depend{1,1} = f1Dt24.depend{1,1}-80;
f1Dt24.ancillary.v_edges = f1Dt24.ancillary.v_edges-80;


%% Load Model data 
load('histstruct_shock1.mat');
modeldists1 = histstruct;
moms1 = calculatemodelmoments(modeldists1);

load('histstruct_shock3.mat');
modeldists3 = histstruct;
moms3 = calculatemodelmoments(modeldists3);

load('histstruct_shock4.mat');
modeldists4 = histstruct;
moms4 = calculatemodelmoments(modeldists4);

%load('histstruct_shock5.mat');
%modeldists5 = histstruct;

%% calculate bulk proton velocites from peaks
% Observations 
[~,idxv] = max(f1Dn1.data,[],2);
vnvec = f1Dn1.depend{1,1}(1,:);
vn1 = vnvec(idxv);
vn1 = irf.ts_scalar(f1Dn1.time,vn1);

[~,idxv] = max(f1Dt21.data,[],2);
vt2vec = f1Dt21.depend{1,1}(1,:);
vt21 = vt2vec(idxv);
vt21 = irf.ts_scalar(f1Dt21.time,vt21);

[~,idxv] = max(f1Dn3.data,[],2);
vnvec = f1Dn3.depend{1,1}(1,:);
vn3 = vnvec(idxv);
vn3 = irf.ts_scalar(f1Dn3.time,vn3);

[~,idxv] = max(f1Dt23.data,[],2);
vt2vec = f1Dt23.depend{1,1}(1,:);
vt23 = vt2vec(idxv);
vt23 = irf.ts_scalar(f1Dt23.time,vt23);

[~,idxv] = max(f1Dn4.data,[],2);
vnvec = f1Dn4.depend{1,1}(1,:);
vn4 = vnvec(idxv);
vn4 = irf.ts_scalar(f1Dn4.time,vn4);

[~,idxv] = max(f1Dt24.data,[],2);
vt2vec = f1Dt24.depend{1,1}(1,:);
vt24 = vt2vec(idxv);
vt24 = irf.ts_scalar(f1Dt24.time,vt24);

% Model 
if 0,
[~,idxv] = max(modeldists1.n1Dvx,[],2);
vn1m = modeldists1.vxpositions(idxv);

[~,idxv] = max(modeldists1.n1Dvz,[],2);
vt21m = modeldists1.vzpositions(idxv);

[~,idxv] = max(modeldists3.n1Dvx,[],2);
vn3m = modeldists1.vxpositions(idxv);

[~,idxv] = max(modeldists3.n1Dvz,[],2);
vt23m = modeldists3.vzpositions(idxv);

[~,idxv] = max(modeldists4.n1Dvx,[],2);
vn4m = modeldists4.vxpositions(idxv);

[~,idxv] = max(modeldists4.n1Dvz,[],2);
vt24m = modeldists4.vzpositions(idxv);
end

vn1m = moms1.vx;
vt21m = moms1.vz;

vn3m = moms3.vx;
vt23m = moms3.vz;

vn4m = moms4.vx;
vt24m = moms4.vz;

%% Color map
c = [55,137,187;...
  106,193,165;...
  172,220,166;...
  230,244,157;...
  255,254,194;...
  253,223,144;...
  251,173,104;...
  242,109,074;...
  211,064,082]/255;
cmap = interp1(linspace(1,64,size(c,1)),c,1:64);

%% Plot Figure

h=irf_plot(18,'newfigure');
%h=irf_figure(540+ic,8);
xSize=600; ySize=800;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.40;
ywidth = 0.09;

set(h(1),'position',[0.1 0.975-ywidth xwidth ywidth]);
set(h(2),'position',[0.1 0.975-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.1 0.975-3*ywidth xwidth ywidth]);

set(h(4),'position',[0.1 0.915-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.1 0.915-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.1 0.915-6*ywidth xwidth ywidth]);

set(h(7),'position',[0.1 0.855-7*ywidth xwidth ywidth]);
set(h(8),'position',[0.1 0.855-8*ywidth xwidth ywidth]);
set(h(9),'position',[0.1 0.855-9*ywidth xwidth ywidth]);

set(h(10),'position',[0.515 0.975-ywidth xwidth ywidth]);
set(h(11),'position',[0.515 0.975-2*ywidth xwidth ywidth]);
set(h(12),'position',[0.515 0.975-3*ywidth xwidth ywidth]);

set(h(13),'position',[0.515 0.915-4*ywidth xwidth ywidth]);
set(h(14),'position',[0.515 0.915-5*ywidth xwidth ywidth]);
set(h(15),'position',[0.515 0.915-6*ywidth xwidth ywidth]);

set(h(16),'position',[0.515 0.855-7*ywidth xwidth ywidth]);
set(h(17),'position',[0.515 0.855-8*ywidth xwidth ywidth]);
set(h(18),'position',[0.515 0.855-9*ywidth xwidth ywidth]);

h(1)=irf_panel('Bntt1');
irf_plot(h(1),Bntt1);
ylabel(h(1),{'B (nT)'},'Interpreter','tex','fontsize',12);
irf_legend(h(1),{'B_{n}','B_{t1}','B_{t2}'},[0.70 0.7],'fontsize',12)
irf_legend(h(1),'(a)',[0.02 0.98],'color','k','fontsize',12)
irf_zoom(h(1),'y',[-15 80])
title(h(1),'Shock 1')

specrecn = f1Dn1.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(2)=irf_panel('f1Dn1');
irf_spectrogram(h(2),specrecn,'donotshowcolorbar');
hold(h(2),'on')
irf_plot(h(2),vn1,'k')
hold(h(2),'off')
ylabel(h(2),{'v_n (km s^{-1})'},'interpreter','tex','fontsize',12)
colormap(h(2),cmap)
irf_zoom(h(2),'y',vlimn)
caxis(h(2),[-3 2])
irf_legend(h(2),'(b)',[0.02 0.98],'color','k','fontsize',12)

specrect2 = f1Dt21.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2}','(s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(3)=irf_panel('f1Dt21');
irf_spectrogram(h(3),specrect2,'donotshowcolorbar');
hold(h(3),'on')
irf_plot(h(3),vt21,'k')
hold(h(3),'off')
ylabel(h(3),{'v_{t2} (km s^{-1})'},'interpreter','tex','fontsize',12)
colormap(h(3),cmap)
irf_zoom(h(3),'y',vlimt2)
caxis(h(3),[-3 2])
irf_legend(h(3),'(c)',[0.02 0.98],'color','k','fontsize',12)

irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',Tint1);
irf_timeaxis(h(1:3),'nodate');


h(4)=irf_panel('Bntt3');
irf_plot(h(4),Bntt3);
ylabel(h(4),{'B (nT)'},'Interpreter','tex','fontsize',12);
%irf_legend(h(4),{'B_{n}','B_{t1}','B_{t2}'},[0.70 0.7],'fontsize',12)
irf_legend(h(4),'(g)',[0.02 0.98],'color','k','fontsize',12)
irf_zoom(h(4),'y',[-15 80])
title(h(4),'Shock 3')

specrecn = f1Dn3.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(5)=irf_panel('f1Dn3');
irf_spectrogram(h(5),specrecn,'donotshowcolorbar');
hold(h(5),'on')
irf_plot(h(5),vn3,'k')
hold(h(5),'off')
ylabel(h(5),{'v_n (km s^{-1})'},'interpreter','tex','fontsize',12)
colormap(h(5),cmap)
irf_zoom(h(5),'y',vlimn)
caxis(h(5),[-3 2])
irf_legend(h(5),'(h)',[0.02 0.98],'color','k','fontsize',12)

specrect2 = f1Dt23.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2}','(s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(6)=irf_panel('f1Dt23');
irf_spectrogram(h(6),specrect2,'donotshowcolorbar');
hold(h(6),'on')
irf_plot(h(6),vt23,'k')
hold(h(6),'off')
ylabel(h(6),{'v_{t2} (km s^{-1})'},'interpreter','tex','fontsize',12)
colormap(h(6),cmap)
irf_zoom(h(6),'y',vlimt2)
caxis(h(6),[-3 2])
irf_legend(h(6),'(i)',[0.02 0.98],'color','k','fontsize',12)

irf_plot_axis_align(h(4:6));
irf_zoom(h(4:6),'x',Tint3);
irf_timeaxis(h(4:6),'nodate');


h(7)=irf_panel('Bntt4');
irf_plot(h(7),Bntt4);
ylabel(h(7),{'B (nT)'},'Interpreter','tex','fontsize',12);
%irf_legend(h(7),{'B_{n}','B_{t1}','B_{t2}'},[0.70 0.7],'fontsize',12)
irf_legend(h(7),'(m)',[0.99 0.98],'color','k','fontsize',12)
irf_zoom(h(7),'y',[-15 80])
title(h(7),'Shock 4')

specrecn = f1Dn4.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(8)=irf_panel('f1Dn4');
irf_spectrogram(h(8),specrecn,'donotshowcolorbar');
hold(h(8),'on')
irf_plot(h(8),vn4,'k')
hold(h(8),'off')
ylabel(h(8),{'v_n (km s^{-1})'},'interpreter','tex','fontsize',12)
colormap(h(8),cmap)
irf_zoom(h(8),'y',vlimn)
caxis(h(8),[-3 2])
irf_legend(h(8),'(n)',[0.99 0.98],'color','k','fontsize',12)

specrect2 = f1Dt24.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2}','(s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(9)=irf_panel('f1Dt24');
irf_spectrogram(h(9),specrect2,'donotshowcolorbar');
hold(h(9),'on')
irf_plot(h(9),vt24,'k')
hold(h(9),'off')
ylabel(h(9),{'v_{t2} (km s^{-1})'},'interpreter','tex','fontsize',12)
colormap(h(9),cmap)
irf_zoom(h(9),'y',vlimt2)
caxis(h(9),[-3 2])
irf_legend(h(9),'(o)',[0.99 0.98],'color','k','fontsize',12)

irf_plot_axis_align(h(7:9));
irf_zoom(h(7:9),'x',Tint4);
irf_timeaxis(h(7:9),'nodate');

set(h([2:3 5:6 8:9]),'Color',0.75*[1 1 1]);


Units = irf_units;
mu0 = Units.mu0;
e = Units.e;
mi = Units.mp;

n0 = modeldists1.n0; 
n1 = modeldists1.n1;
B0 = modeldists1.B0;
B1 = modeldists1.B1;
P0 = modeldists1.P0;
P1 = modeldists1.P1;
Bn = modeldists1.Bn;
xposp = modeldists1.xpositions;
vxvec = modeldists1.vxpositions;
vzvec = modeldists1.vzpositions;

xvec = [-20e2:1:5e2]*1e3;

l = modeldists1.l;
By = -B0*tanh(xvec/l)+B1;
ni = -n0*tanh(xvec/l)+n1;
Jz = -B0*sech(xvec/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni)+P0*sech(xvec/l).^2./(e*ni*l);

plot(h(10),xvec/1e3,By*1e9,'b')
hold(h(10),'on')
plot(h(10),xvec/1e3,Bn*ones(size(xvec))*1e9,'k')
plot(h(10),xvec/1e3,zeros(size(xvec))*1e9,'r')
hold(h(10),'off')
irf_legend(h(10),{'B_{n}','B_{t1}','B_{t2}'},[0.70 0.7],'fontsize',12)
grid(h(10),'on')
axis(h(10),[-2000 500 -15 80])
set(h(10),'xticklabel',[])
set(h(10),'yticklabel',[])
irf_legend(h(10),'(d)',[0.02 0.98],'color','k','fontsize',12)

ftemp = modeldists1.n1Dvx;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(11),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(11),'flat')
hold(h(11),'on')
plot(h(11),modeldists1.xpositions/1e3,vn1m/1e3,'color','k')
hold(h(11),'off')
axis(h(11),[-1999 500 -1800 1400])
ylabel(h(11),' ','FontSize',12)
caxis(h(11),[-3 2]);
c = colorbar(h(11),'position',[0.920 0.985-3*ywidth 0.013 0.1600]);
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(11),'XTickLabel',[])
colormap(h(11),cmap)
set(h(11),'yticklabel',[])
irf_legend(h(11),'(e)',[0.02 0.98],'color','k','fontsize',12)

ftemp = modeldists1.n1Dvz;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(12),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(12),'flat')
hold(h(12),'on')
plot(h(12),modeldists1.xpositions/1e3,vt21m/1e3,'color','k')
hold(h(12),'off')
axis(h(12),[-1999 500 -1600 1600])
ylabel(h(12),' ','FontSize',12)
caxis(h(12),[-3 2]);
xlabel(h(12),'n (km)')
colormap(h(12),cmap)
set(h(12),'yticklabel',[])
irf_legend(h(12),'(f)',[0.02 0.98],'color','k','fontsize',12)

irf_plot_axis_align(h(10:12));
set(h(10:12),'xdir','reverse')



n0 = modeldists3.n0; 
n1 = modeldists3.n1;
B0 = modeldists3.B0;
B1 = modeldists3.B1;
P0 = modeldists3.P0;
P1 = modeldists3.P1;
Bn = modeldists3.Bn;
xposp = modeldists3.xpositions;
vxvec = modeldists3.vxpositions;
vzvec = modeldists3.vzpositions;

xvec = [-20e2:1:5e2]*1e3;

l = modeldists3.l;
By = -B0*tanh(xvec/l)+B1;
ni = -n0*tanh(xvec/l)+n1;
Jz = -B0*sech(xvec/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni)+P0*sech(xvec/l).^2./(e*ni*l);

plot(h(13),xvec/1e3,By*1e9,'b')
hold(h(13),'on')
plot(h(13),xvec/1e3,Bn*ones(size(xvec))*1e9,'k')
plot(h(13),xvec/1e3,zeros(size(xvec))*1e9,'r')
hold(h(13),'off')
%ylabel(h(13),{'B (nT)'},'fontsize',14)
grid(h(13),'on')
axis(h(13),[-1999 500 -15 80])
set(h(13),'xticklabel',[])
set(h(13),'yticklabel',[])
irf_legend(h(13),'(j)',[0.02 0.98],'color','k','fontsize',12)

ftemp = modeldists3.n1Dvx;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(14),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(14),'flat')
hold(h(14),'on')
plot(h(14),modeldists3.xpositions/1e3,vn3m/1e3,'color','k')
hold(h(14),'off')
axis(h(14),[-1999 500 -1800 1400])
ylabel(h(14),' ','FontSize',12)
caxis(h(14),[-3 2]);
c = colorbar(h(14),'position',[0.920 0.925-6*ywidth 0.013 0.1600]);
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(14),'XTickLabel',[])
colormap(h(14),cmap)
set(h(14),'yticklabel',[])
irf_legend(h(14),'(k)',[0.02 0.98],'color','k','fontsize',12)

ftemp = modeldists3.n1Dvz;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(15),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(15),'flat')
hold(h(15),'on')
plot(h(15),modeldists3.xpositions/1e3,vt23m/1e3,'color','k')
hold(h(15),'off')
axis(h(15),[-1999 500 -1600 1600])
ylabel(h(15),' ','FontSize',12)
caxis(h(15),[-3 2]);
xlabel(h(15),'n (km)')
colormap(h(15),cmap)
set(h(15),'yticklabel',[])
irf_legend(h(15),'(l)',[0.02 0.98],'color','k','fontsize',12)

irf_plot_axis_align(h(13:15));
set(h(13:15),'xdir','reverse')



n0 = modeldists4.n0; 
n1 = modeldists4.n1;
B0 = modeldists4.B0;
B1 = modeldists4.B1;
P0 = modeldists4.P0;
P1 = modeldists4.P1;
Bn = modeldists4.Bn;
xposp = modeldists4.xpositions;
vxvec = modeldists4.vxpositions;
vzvec = modeldists4.vzpositions;

xvec = [-20e2:1:5e2]*1e3;

l = modeldists4.l;
By = -B0*tanh(xvec/l)+B1;
ni = -n0*tanh(xvec/l)+n1;
Jz = -B0*sech(xvec/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni)+P0*sech(xvec/l).^2./(e*ni*l);

plot(h(16),xvec/1e3,By*1e9,'b')
hold(h(16),'on')
plot(h(16),xvec/1e3,Bn*ones(size(xvec))*1e9,'k')
plot(h(16),xvec/1e3,zeros(size(xvec))*1e9,'r')
hold(h(16),'off')
%ylabel(h(16),{'B (nT)'},'fontsize',14)
grid(h(16),'on')
axis(h(16),[-2000 500 -15 80])
set(h(16),'xticklabel',[])
set(h(16),'yticklabel',[])
irf_legend(h(16),'(p)',[0.99 0.98],'color','k','fontsize',12)

ftemp = modeldists4.n1Dvx;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(17),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(17),'flat')
hold(h(17),'on')
plot(h(17),modeldists4.xpositions/1e3,vn4m/1e3,'color','k')
hold(h(17),'off')
axis(h(17),[-1999 500 -1800 1400])
ylabel(h(17),' ','FontSize',12)
caxis(h(17),[-3 2]);
c = colorbar(h(17),'position',[0.920 0.865-9*ywidth 0.013 0.1600]);
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(17),'XTickLabel',[])
colormap(h(17),cmap)
set(h(17),'yticklabel',[])
irf_legend(h(17),'(q)',[0.99 0.98],'color','k','fontsize',12)

ftemp = modeldists4.n1Dvz;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(18),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(18),'flat')
hold(h(18),'on')
plot(h(18),modeldists4.xpositions/1e3,vt24m/1e3,'color','k')
hold(h(18),'off')
axis(h(18),[-1999 500 -1600 1600])
ylabel(h(18),' ','FontSize',12)
caxis(h(18),[-3 2]);
xlabel(h(18),'n (km)')
colormap(h(18),cmap)
set(h(18),'yticklabel',[])
irf_legend(h(18),'(r)',[0.99 0.98],'color','k','fontsize',12)

irf_plot_axis_align(h(16:18));



set(h([11:12 14:15 17:18]),'Color',0.75*[1 1 1]);

set(h(1:18),'fontsize',12);

