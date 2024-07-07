ic = 1;

Tint1 = irf.tint('2023-04-24T02:28:50.000Z/2023-04-24T02:30:50.000Z');
Tint2 = irf.tint('2023-04-24T04:01:40.000Z/2023-04-24T04:03:30.000Z');
Tint3 = irf.tint('2023-04-24T04:20:35.000Z/2023-04-24T04:22:20.000Z');

iPDist1 = mms.get_data('PDi_fpi_brst_l2',Tint1,ic);
iPDist2 = mms.get_data('PDi_fpi_brst_l2',Tint2,ic);
iPDist3 = mms.get_data('PDi_fpi_brst_l2',Tint3,ic);

Bxyz1 = mms.get_data('B_dmpa_srvy_l2',Tint1,ic);
Bxyz2 = mms.get_data('B_dmpa_srvy_l2',Tint2,ic);
Bxyz3 = mms.get_data('B_dmpa_srvy_l2',Tint3,ic);

Bxyzlf1 = Bxyz1.filt(0,0.25,16,5);
Bxyzlf2 = Bxyz2.filt(0,0.25,16,5);
Bxyzlf3 = Bxyz3.filt(0,0.25,16,5);

Bxyz1mag = Bxyz1.abs;
Bxyz2mag = Bxyz2.abs;
Bxyz3mag = Bxyz3.abs;

Bxyz1maglf = Bxyzlf1.abs;
Bxyz2maglf = Bxyzlf2.abs;
Bxyz3maglf = Bxyzlf3.abs;

%%
ndsl = [0.9203   -0.3886   -0.0445];
t1dsl = [-0.2271   -0.4383   -0.8696];
t2dsl = [0.3184    0.8105   -0.4917];
ndsl1 = ndsl/norm(ndsl);
t1dsl1 = t1dsl/norm(t1dsl);
t2dsl1 = t2dsl/norm(t2dsl);

ndsl = [0.8871   -0.4602   -0.0356];
t1dsl = [-0.3382   -0.5956   -0.7286];
t2dsl = [0.3142    0.6584   -0.6840];
ndsl2 = ndsl/norm(ndsl);
t1dsl2 = t1dsl/norm(t1dsl);
t2dsl2 = t2dsl/norm(t2dsl);

ndsl = [0.8027   -0.5795    0.1412];
t1dsl = [-0.3078   -0.6052   -0.7341];
t2dsl = [0.5108    0.5458   -0.6642];
ndsl3 = ndsl/norm(ndsl);
t1dsl3 = t1dsl/norm(t1dsl);
t2dsl3 = t2dsl/norm(t2dsl);
%%

iPDist1.data(:,1:11,:,:) = 0;
iPDist2.data(:,1:11,:,:) = 0;
iPDist3.data(:,1:11,:,:) = 0;

vlimn = [-1000,500];

nMC = 5e2;
vg1Dn = linspace(vlimn(1),vlimn(2),300);

nvec1 = irf.ts_vec_xyz(iPDist1.time,repmat(ndsl1,length(iPDist1.time),1));
nvec2 = irf.ts_vec_xyz(iPDist2.time,repmat(ndsl2,length(iPDist2.time),1));
nvec3 = irf.ts_vec_xyz(iPDist3.time,repmat(ndsl3,length(iPDist3.time),1));

f1Dn1 = iPDist1.reduce('1D',nvec1,'vg',vg1Dn,'nMC',nMC);

f1Dn2 = iPDist2.reduce('1D',nvec2,'vg',vg1Dn,'nMC',nMC);

f1Dn3 = iPDist3.reduce('1D',nvec3,'vg',vg1Dn,'nMC',nMC);

%% calculate bulk proton velocites from peaks
% Observations 
[~,idxv] = max(f1Dn1.data,[],2);
vnvec = f1Dn1.depend{1,1}(1,:);
vn1 = vnvec(idxv);
vn1 = irf.ts_scalar(f1Dn1.time,vn1);

[~,idxv] = max(f1Dn2.data,[],2);
vnvec = f1Dn2.depend{1,1}(1,:);
vn2 = vnvec(idxv);
vn2 = irf.ts_scalar(f1Dn2.time,vn2);

[~,idxv] = max(f1Dn3.data,[],2);
vnvec = f1Dn3.depend{1,1}(1,:);
vn3 = vnvec(idxv);
vn3 = irf.ts_scalar(f1Dn3.time,vn3);

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

%% 

h=irf_plot(6,'newfigure');
%h=irf_figure(540+ic,8);
xSize=700; ySize=700;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.86;
ywidth = 0.123;
set(h(1),'position',[0.10 0.97-ywidth xwidth ywidth]);
set(h(2),'position',[0.10 0.97-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.10 0.89-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.10 0.89-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.10 0.81-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.10 0.81-6*ywidth xwidth ywidth]);

% B minima times
tint1B1 = irf.tint('2023-04-24T02:29:17.5Z/2023-04-24T02:29:17.5Z');
tint1B2 = irf.tint('2023-04-24T02:29:25.5Z/2023-04-24T02:29:25.5Z');
tint1B3 = irf.tint('2023-04-24T02:29:32.0Z/2023-04-24T02:29:32.0Z');
tint1B4 = irf.tint('2023-04-24T02:29:38.0Z/2023-04-24T02:29:38.0Z');
tint1B5 = irf.tint('2023-04-24T02:29:45.0Z/2023-04-24T02:29:45.0Z');
tint1B6 = irf.tint('2023-04-24T02:29:50.2Z/2023-04-24T02:29:50.2Z');
tint1B7 = irf.tint('2023-04-24T02:29:57.0Z/2023-04-24T02:29:57.0Z');
tint1B8 = irf.tint('2023-04-24T02:30:03.4Z/2023-04-24T02:30:03.4Z');
tint1B9 = irf.tint('2023-04-24T02:30:07.5Z/2023-04-24T02:30:07.5Z');
tint1B10 = irf.tint('2023-04-24T02:30:12.2Z/2023-04-24T02:30:12.2Z');
tint1B11 = irf.tint('2023-04-24T02:30:18.0Z/2023-04-24T02:30:18.0Z');


h(1)=irf_panel('Bxyz1');
irf_plot(h(1),Bxyz1mag);
hold(h(1),'on')
irf_plot(h(1),Bxyz1maglf,'r')
c_eval('irf_plot(h(1),[tint1B?.epochUnix [-100; 100]],''b--'')',1:11)
hold(h(1),'off')
ylabel(h(1),{'|B| (nT)'},'Interpreter','tex','fontsize',14);
irf_zoom(h(1),'y',[0 80]);
irf_legend(h(1),'(a)',[0.99 0.7],'color','k','fontsize',14)
title(h(1),'Shock 1');

specrecn = f1Dn1.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(2)=irf_panel('f1Dn1');
irf_spectrogram(h(2),specrecn,'donotshowcolorbar');
colormap(h(2),cmap)
irf_zoom(h(2),'y',vlimn)
caxis(h(2),[-3 2])
hold(h(2),'on')
irf_plot(h(2),vn1,'k')
c_eval('irf_plot(h(2),[tint1B?.epochUnix [-1000; 1000]],''b--'')',1:11)
hold(h(2),'off')
ylabel(h(2),{'v_n (km s^{-1})'},'interpreter','tex')
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',14)

irf_zoom(h(1:2),'x',Tint1);
irf_timeaxis(h(1:2),'nodate');


% B minima times
tint2B1 = irf.tint('2023-04-24T04:02:10.3Z/2023-04-24T04:02:10.3Z');
tint2B2 = irf.tint('2023-04-24T04:02:19.0Z/2023-04-24T04:02:19.0Z');
tint2B3 = irf.tint('2023-04-24T04:02:27.3Z/2023-04-24T04:02:27.3Z');
tint2B4 = irf.tint('2023-04-24T04:02:36.2Z/2023-04-24T04:02:36.2Z');
tint2B5 = irf.tint('2023-04-24T04:02:47.0Z/2023-04-24T04:02:47.0Z');
tint2B6 = irf.tint('2023-04-24T04:03:00.0Z/2023-04-24T04:03:00.0Z');
tint2B7 = irf.tint('2023-04-24T04:03:25.0Z/2023-04-24T04:03:25.0Z');

h(3)=irf_panel('Bxyz2');
irf_plot(h(3),Bxyz2mag);
irf_zoom(h(3),'y',[0 80]);
hold(h(3),'on')
irf_plot(h(3),Bxyz2maglf,'r')
c_eval('irf_plot(h(3),[tint2B?.epochUnix [-100; 100]],''b--'')',1:7)
hold(h(3),'off')
ylabel(h(3),{'|B| (nT)'},'Interpreter','tex','fontsize',14);
irf_legend(h(3),'(c)',[0.99 0.7],'color','k','fontsize',14)
title(h(3),'Shock 3');

specrecn = f1Dn2.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(4)=irf_panel('f1Dn2');
irf_spectrogram(h(4),specrecn,'donotshowcolorbar');
colormap(h(4),cmap)
irf_zoom(h(4),'y',vlimn)
caxis(h(4),[-3 2])
hold(h(4),'on')
irf_plot(h(4),vn2,'k')
c_eval('irf_plot(h(4),[tint2B?.epochUnix [-1000; 1000]],''b--'')',1:7)
hold(h(4),'off')
ylabel(h(4),{'v_n (km s^{-1})'},'interpreter','tex')
irf_legend(h(4),'(d)',[0.99 0.98],'color','k','fontsize',14)

irf_zoom(h(3:4),'x',Tint2);
irf_timeaxis(h(3:4),'nodate');


% B minima times
tint3B1 = irf.tint('2023-04-24T04:21:13.5Z/2023-04-24T04:21:13.5Z');
tint3B2 = irf.tint('2023-04-24T04:21:25.0Z/2023-04-24T04:21:25.0Z');
tint3B3 = irf.tint('2023-04-24T04:21:38.0Z/2023-04-24T04:21:38.0Z');

h(5)=irf_panel('Bxyz');
irf_plot(h(5),Bxyz3mag);
irf_zoom(h(5),'y',[0 80]);
hold(h(5),'on')
irf_plot(h(5),Bxyz3maglf,'r')
c_eval('irf_plot(h(5),[tint3B?.epochUnix [-100; 100]],''b--'')',1:3)
hold(h(5),'off')
ylabel(h(5),{'|B| (nT)'},'Interpreter','tex','fontsize',14);
irf_legend(h(5),'(e)',[0.99 0.7],'color','k','fontsize',14)
title(h(5),'Shock 5');

specrecn = f1Dn3.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(6)=irf_panel('f1Dn3');
irf_spectrogram(h(6),specrecn,'donotshowcolorbar');
colormap(h(6),cmap)
irf_zoom(h(6),'y',vlimn)
caxis(h(6),[-3 2])
hold(h(6),'on')
irf_plot(h(6),vn3,'k')
c_eval('irf_plot(h(6),[tint3B?.epochUnix [-1000; 1000]],''b--'')',1:3)
hold(h(6),'off')
ylabel(h(6),{'v_n (km s^{-1})'},'interpreter','tex')
irf_legend(h(6),'(f)',[0.99 0.98],'color','k','fontsize',14)

irf_zoom(h(5:6),'x',Tint3);
irf_plot_axis_align(h(1:6));
set(h(1:6),'fontsize',14);