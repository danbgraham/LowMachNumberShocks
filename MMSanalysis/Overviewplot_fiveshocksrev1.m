Tint = irf.tint('2023-04-24T02:00:00.000Z/2023-04-24T04:40:00.000Z');

ic = 1;
iPDist = mms.get_data('PDi_fpi_fast_l2',Tint,ic);
Bxyz=mms.get_data('B_gse_srvy_l2',Tint,ic);
ne = mms.get_data('Ne_fpi_fast_l2',Tint,ic);
Vi = mms.get_data('Vi_dbcs_fpi_fast_l2',Tint,ic);

iPDistomni = iPDist.omni.deflux;

%% Constants
Units = irf_units;
e = Units.e;
me = Units.me;
mi = Units.mp;

Vabs = Vi.abs;
EVi = 0.5*mi*Vi.abs2*1e6/e;

%% Load burst mode data

Tint1 = irf.tint('2023-04-24T02:28:50.000Z/2023-04-24T02:30:10.000Z');
Tint2 = irf.tint('2023-04-24T03:49:33.000Z/2023-04-24T03:50:30.000Z');
Tint3 = irf.tint('2023-04-24T04:01:40.000Z/2023-04-24T04:03:00.000Z');
Tint4 = irf.tint('2023-04-24T04:14:31.000Z/2023-04-24T04:15:51.000Z');
Tint5 = irf.tint('2023-04-24T04:20:35.000Z/2023-04-24T04:21:55.000Z');

iPDist1 = mms.get_data('PDi_fpi_brst_l2',Tint1,ic);
iPDist2 = mms.get_data('PDi_fpi_brst_l2',Tint2,ic);
iPDist3 = mms.get_data('PDi_fpi_brst_l2',Tint3,ic);
iPDist4 = mms.get_data('PDi_fpi_brst_l2',Tint4,ic);
iPDist5 = mms.get_data('PDi_fpi_brst_l2',Tint5,ic);

Bxyz1 = mms.get_data('B_dmpa_brst_l2',Tint1,ic);
Bxyz2 = mms.get_data('B_dmpa_brst_l2',Tint2,ic);
Bxyz3 = mms.get_data('B_dmpa_brst_l2',Tint3,ic);
Bxyz4 = mms.get_data('B_dmpa_brst_l2',Tint4,ic);
Bxyz5 = mms.get_data('B_dmpa_brst_l2',Tint5,ic);

%% Coordinate systems
ndsl = [0.9203   -0.3886   -0.0445];
t1dsl = [-0.2271   -0.4383   -0.8696];
t2dsl = [0.3184    0.8105   -0.4917];
ndsl1 = ndsl/norm(ndsl);
t1dsl1 = t1dsl/norm(t1dsl);
t2dsl1 = t2dsl/norm(t2dsl);

ndsl = [0.9015    0.1799   -0.3935];
t1dsl = [-0.1538   -0.7168   -0.6801];
t2dsl = [-0.4045    0.6736   -0.6186];
ndsl2 = ndsl/norm(ndsl);
t1dsl2 = t1dsl/norm(t1dsl);
t2dsl2 = t2dsl/norm(t2dsl);

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
Bntt2 = irf_newxyz(Bxyz2,ndsl2,t1dsl2,t2dsl2);
Bntt3 = irf_newxyz(Bxyz3,ndsl3,t1dsl3,t2dsl3);
Bntt4 = irf_newxyz(Bxyz4,ndsl4,t1dsl4,t2dsl4);
Bntt5 = irf_newxyz(Bxyz5,ndsl5,t1dsl5,t2dsl5);

%% 1D reduced ion distributions
iPDist1.data(:,1:11,:,:) = 0;
iPDist2.data(:,1:11,:,:) = 0;
iPDist3.data(:,1:11,:,:) = 0;
iPDist4.data(:,1:11,:,:) = 0;
iPDist5.data(:,1:11,:,:) = 0;

vlimn = [-1800,1400];
vlimt2 = [-1600,1600];

nMC = 5e2;
vg1Dn = linspace(vlimn(1),vlimn(2),300);
vg1Dt2 = linspace(vlimt2(1),vlimt2(2),300);

nvec1 = irf.ts_vec_xyz(iPDist1.time,repmat(ndsl1,length(iPDist1.time),1));
t2vec1 = irf.ts_vec_xyz(iPDist1.time,repmat(t2dsl1,length(iPDist1.time),1));

nvec2 = irf.ts_vec_xyz(iPDist2.time,repmat(ndsl2,length(iPDist2.time),1));
t2vec2 = irf.ts_vec_xyz(iPDist2.time,repmat(t2dsl2,length(iPDist2.time),1));

nvec3 = irf.ts_vec_xyz(iPDist3.time,repmat(ndsl3,length(iPDist3.time),1));
t2vec3 = irf.ts_vec_xyz(iPDist3.time,repmat(t2dsl3,length(iPDist3.time),1));

nvec4 = irf.ts_vec_xyz(iPDist4.time,repmat(ndsl4,length(iPDist4.time),1));
t2vec4 = irf.ts_vec_xyz(iPDist4.time,repmat(t2dsl4,length(iPDist4.time),1));

nvec5 = irf.ts_vec_xyz(iPDist5.time,repmat(ndsl5,length(iPDist5.time),1));
t2vec5 = irf.ts_vec_xyz(iPDist5.time,repmat(t2dsl5,length(iPDist5.time),1));

f1Dn1 = iPDist1.reduce('1D',nvec1,'vg',vg1Dn,'nMC',nMC);
f1Dt21 = iPDist1.reduce('1D',t2vec1,'vg',vg1Dt2,'nMC',nMC);

f1Dn2 = iPDist2.reduce('1D',nvec2,'vg',vg1Dn,'nMC',nMC);
f1Dt22 = iPDist2.reduce('1D',t2vec2,'vg',vg1Dt2,'nMC',nMC);

f1Dn3 = iPDist3.reduce('1D',nvec3,'vg',vg1Dn,'nMC',nMC);
f1Dt23 = iPDist3.reduce('1D',t2vec3,'vg',vg1Dt2,'nMC',nMC);

f1Dn4 = iPDist4.reduce('1D',nvec4,'vg',vg1Dn,'nMC',nMC);
f1Dt24 = iPDist4.reduce('1D',t2vec4,'vg',vg1Dt2,'nMC',nMC);

f1Dn5 = iPDist5.reduce('1D',nvec5,'vg',vg1Dn,'nMC',nMC);
f1Dt25 = iPDist5.reduce('1D',t2vec5,'vg',vg1Dt2,'nMC',nMC);

%%
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

rr = interp1([1 64 128 192 256],[0 90 190 250 255]/255,1:256,'pchip');
gg = interp1([1 64 128 192 256],[0 15 55 140 255]/255,1:256,'pchip');
bb = interp1([1 64 128 192 256],[0 110 80 10 0]/255,1:256,'pchip');
cmapinf = [rr' gg' bb'];

%% Plot Figure

h=irf_plot(18,'newfigure');
%h=irf_figure(540+ic,8);
xSize=900; ySize=1000;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.92;
ywidth = 0.09;
ywidth2 = 0.095;
set(h(1),'position',[0.07 0.98-ywidth xwidth ywidth]);
set(h(2),'position',[0.07 0.98-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.07 0.98-3*ywidth xwidth ywidth]);

set(h(4),'position',[0.08 0.935-4*ywidth2 0.28 ywidth2]);
set(h(5),'position',[0.08 0.935-5*ywidth2 0.28 ywidth2]);
set(h(6),'position',[0.08 0.935-6*ywidth2 0.28 ywidth2]);

set(h(7),'position',[0.37 0.935-4*ywidth2 0.28 ywidth2]);
set(h(8),'position',[0.37 0.935-5*ywidth2 0.28 ywidth2]);
set(h(9),'position',[0.37 0.935-6*ywidth2 0.28 ywidth2]);

set(h(10),'position',[0.66 0.935-4*ywidth2 0.32 ywidth2]);
set(h(11),'position',[0.66 0.935-5*ywidth2 0.32 ywidth2]);
set(h(12),'position',[0.66 0.935-6*ywidth2 0.32 ywidth2]);

set(h(13),'position',[0.225 0.03+2*ywidth2 0.28 ywidth2]);
set(h(14),'position',[0.225 0.03+ywidth2 0.28 ywidth2]);
set(h(15),'position',[0.225 0.03 0.28 ywidth2]);

set(h(16),'position',[0.515 0.03+2*ywidth2 0.32 ywidth2]);
set(h(17),'position',[0.515 0.03+ywidth2 0.32 ywidth2]);
set(h(18),'position',[0.515 0.03 0.32 ywidth2]);

h(1)=irf_panel('idist');
irf_spectrogram(h(1),iPDistomni.specrec,'log');
colormap(h(1),cmapinf);
hold(h(1),'on')
irf_plot(h(1),EVi,'k')
hold(h(1),'off')
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',14);
set(h(1),'yscale','log');
set(h(1),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(1),'E_{i} (eV)','fontsize',14,'Interpreter','tex');
grid(h(1),'off')
title(h(1),'MMS1')

Bxyzmag = irf.ts_scalar(Bxyz.time,[Bxyz.data Bxyz.abs.data]);

h(2)=irf_panel('Bxyz');
irf_plot(h(2),Bxyzmag);
ylabel(h(2),{'B (nT)'},'Interpreter','tex');
irf_legend(h(2),{'B_{x}','B_{y}','B_{z}','|B|'},[0.65 0.75],'fontsize',14)
irf_legend(h(2),'(b)',[0.99 0.8],'color','k','fontsize',14)
irf_zoom(h(2),'y',[-70 80])

h(3)=irf_panel('ne');
irf_plot(h(3),ne);
ylabel(h(3),{'n_e (cm^{-3})'},'Interpreter','tex');
irf_legend(h(3),'(c)',[0.99 0.98],'color','k','fontsize',14)

irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',Tint);

% Shock 1
h(4)=irf_panel('Bntt1');
irf_plot(h(4),Bntt1);
ylabel(h(4),{'B (nT)'},'Interpreter','tex');
irf_legend(h(4),{'B_{n}','B_{t1}','B_{t2}'},[0.70 0.7],'fontsize',14)
irf_legend(h(4),'(d)',[0.99 0.7],'color','k','fontsize',14)
irf_zoom(h(4),'y',[-15 80])
title(h(4),'Shock 1')

specrecn = f1Dn1.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(5)=irf_panel('f1Dn1');
irf_spectrogram(h(5),specrecn,'donotshowcolorbar');
ylabel(h(5),{'v_n (km s^{-1})'},'interpreter','tex')
colormap(h(5),cmap)
irf_zoom(h(5),'y',vlimn)
caxis(h(5),[-3 2])
irf_legend(h(5),'(e)',[0.99 0.98],'color','k','fontsize',14)

specrect2 = f1Dt21.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2}','(s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(6)=irf_panel('f1Dt21');
irf_spectrogram(h(6),specrect2,'donotshowcolorbar');
ylabel(h(6),{'v_{t2} (km s^{-1})'},'interpreter','tex')
colormap(h(6),cmap)
irf_zoom(h(6),'y',vlimt2)
caxis(h(6),[-3 2])
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',14)

irf_plot_axis_align(h(4:6));
irf_zoom(h(4:6),'x',Tint1);
irf_timeaxis(h(4:6),'nodate');

% Shock 3 
h(7)=irf_panel('Bntt3');
irf_plot(h(7),Bntt3);
ylabel(h(7),' ');
set(h(7),'yticklabel',[])
irf_legend(h(7),'(g)',[0.99 0.7],'color','k','fontsize',14)
irf_zoom(h(7),'y',[-15 80])
title(h(7),'Shock 3')

specrecn = f1Dn3.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(8)=irf_panel('f1Dn3');
irf_spectrogram(h(8),specrecn,'donotshowcolorbar');
ylabel(h(8),{' '})
set(h(8),'yticklabel',[])
colormap(h(8),cmap)
irf_zoom(h(8),'y',vlimn)
caxis(h(8),[-3 2])
irf_legend(h(8),'(h)',[0.99 0.98],'color','k','fontsize',14)

specrect2 = f1Dt23.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2}','(s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(9)=irf_panel('f1Dt23');
irf_spectrogram(h(9),specrect2,'donotshowcolorbar');
ylabel(h(9),{' '})
colormap(h(9),cmap)
irf_zoom(h(9),'y',vlimt2)
set(h(9),'yticklabel',[])
caxis(h(9),[-3 2])
irf_legend(h(9),'(i)',[0.99 0.98],'color','k','fontsize',14)

irf_plot_axis_align(h(7:9));
irf_zoom(h(7:9),'x',Tint3);
irf_timeaxis(h(7:9),'nodate');

% Shock 5
h(10)=irf_panel('Bntt5');
irf_plot(h(10),Bntt5);
ylabel(h(10),' ');
set(h(10),'yticklabel',[])
irf_legend(h(10),'(j)',[0.99 0.7],'color','k','fontsize',14)
irf_zoom(h(10),'y',[-15 80])
title(h(10),'Shock 5')

specrecn = f1Dn5.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n} (s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(11)=irf_panel('f1Dn5');
irf_spectrogram(h(11),specrecn);
ylabel(h(11),{' '})
set(h(11),'yticklabel',[])
colormap(h(11),cmap)
irf_zoom(h(11),'y',vlimn)
caxis(h(11),[-3 2])
irf_legend(h(11),'(k)',[0.99 0.98],'color','k','fontsize',14)

specrect2 = f1Dt25.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2} (s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(12)=irf_panel('f1Dt25');
irf_spectrogram(h(12),specrect2);
ylabel(h(12),{' '})
colormap(h(12),cmap)
irf_zoom(h(12),'y',vlimt2)
set(h(12),'yticklabel',[])
caxis(h(12),[-3 2])
irf_legend(h(12),'(l)',[0.99 0.98],'color','k','fontsize',14)

irf_plot_axis_align(h(10:12));
irf_zoom(h(10:12),'x',Tint5);
irf_timeaxis(h(10:12),'nodate');

% Shock 2
h(13)=irf_panel('Bntt2');
irf_plot(h(13),Bntt2);
ylabel(h(13),{'B (nT)'},'Interpreter','tex');
irf_legend(h(13),'(m)',[0.99 0.98],'color','k','fontsize',14)
irf_zoom(h(13),'y',[-15 80])
title(h(13),'Shock 2')

specrecn = f1Dn2.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(14)=irf_panel('f1Dn2');
irf_spectrogram(h(14),specrecn,'donotshowcolorbar');
ylabel(h(14),{'v_n (km s^{-1})'},'interpreter','tex')
colormap(h(14),cmap)
irf_zoom(h(14),'y',vlimn)
caxis(h(14),[-3 2])
irf_legend(h(14),'(n)',[0.99 0.98],'color','k','fontsize',14)

specrect2 = f1Dt22.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2}','(s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(15)=irf_panel('f1Dt22');
irf_spectrogram(h(15),specrect2,'donotshowcolorbar');
ylabel(h(15),{'v_{t2} (km s^{-1})'},'interpreter','tex')
colormap(h(15),cmap)
irf_zoom(h(15),'y',vlimt2)
caxis(h(15),[-3 2])
irf_legend(h(15),'(o)',[0.99 0.98],'color','k','fontsize',14)

irf_plot_axis_align(h(13:15));
irf_zoom(h(13:15),'x',Tint2);
irf_timeaxis(h(13:15),'nodate');

% Shock 4
h(16)=irf_panel('Bntt4');
irf_plot(h(16),Bntt4);
ylabel(h(16),' ');
set(h(16),'yticklabel',[])
irf_legend(h(16),'(p)',[0.99 0.98],'color','k','fontsize',14)
irf_zoom(h(16),'y',[-15 80])
title(h(16),'Shock 4')

specrecn = f1Dn4.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n} (s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(17)=irf_panel('f1Dn4');
irf_spectrogram(h(17),specrecn);
ylabel(h(17),{' '})
set(h(17),'yticklabel',[])
colormap(h(17),cmap)
irf_zoom(h(17),'y',vlimn)
caxis(h(17),[-3 2])
irf_legend(h(17),'(q)',[0.99 0.98],'color','k','fontsize',14)

specrect2 = f1Dt24.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2} (s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(18)=irf_panel('f1Dt24');
irf_spectrogram(h(18),specrect2);
ylabel(h(18),{' '})
colormap(h(18),cmap)
irf_zoom(h(18),'y',vlimt2)
set(h(18),'yticklabel',[])
caxis(h(18),[-3 2])
irf_legend(h(18),'(r)',[0.99 0.98],'color','k','fontsize',14)

irf_plot_axis_align(h(16:18));
irf_zoom(h(16:18),'x',Tint4);
irf_timeaxis(h(16:18),'nodate');


set(h([5:6 8:9 11:12 14:15 17:18]),'Color',0.75*[1 1 1]);

irf_plot_zoomin_lines_between_panels(h(3),h(4));
c_eval('irf_pl_mark(h(?),irf_time(Tint1,''epochtt>epoch'')'',[230 200 0]/255)',1:3);

irf_plot_zoomin_lines_between_panels(h(3),h(7));
c_eval('irf_pl_mark(h(?),irf_time(Tint3,''epochtt>epoch'')'',[230 200 0]/255)',1:3);

irf_plot_zoomin_lines_between_panels(h(3),h(10));
c_eval('irf_pl_mark(h(?),irf_time(Tint5,''epochtt>epoch'')'',[230 200 0]/255)',1:3);

c_eval('irf_pl_mark(h(?),irf_time(Tint2,''epochtt>epoch'')'',[230 200 0]/255)',1:3);
c_eval('irf_pl_mark(h(?),irf_time(Tint4,''epochtt>epoch'')'',[230 200 0]/255)',1:3);

set(h(1:18),'fontsize',14);
xtickangle(h(1:18),0)