%% Shock times
Tint1 = irf.tint('2023-04-24T02:28:13.00Z/2023-04-24T02:31:03.00Z');
Tint2 = irf.tint('2023-04-24T03:49:33.00Z/2023-04-24T03:51:03.00Z');
Tint3 = irf.tint('2023-04-24T04:00:03.00Z/2023-04-24T04:03:33.00Z');
Tint4 = irf.tint('2023-04-24T04:14:03.00Z/2023-04-24T04:16:33.00Z');
Tint5 = irf.tint('2023-04-24T04:20:03.00Z/2023-04-24T04:22:23.00Z');

Tints1 = irf.tint('2023-04-24T02:29:11.85Z/2023-04-24T02:29:11.85Z');
Tints2 = irf.tint('2023-04-24T03:50:12.3Z/2023-04-24T03:50:12.3Z');
Tints3 = irf.tint('2023-04-24T04:02:03.384Z/2023-04-24T04:02:04.384Z');
Tints4 = irf.tint('2023-04-24T04:15:25.0Z/2023-04-24T04:15:25.0Z');
Tints5 = irf.tint('2023-04-24T04:21:00.00Z/2023-04-24T04:21:00.00Z');

Vsn1 = 50;
Vsn2 = -140;
Vsn3 = 60;
Vsn4 = -50;
Vsn5 = 40;

%% Shock coordinate systems
ndsl = [0.9686    0.0056   -0.2487];
t1dsl = [-0.2263   -0.4382   -0.8700];
t2dsl = [0.2695    0.8301   -0.4882];
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

%% Load data
Bxyz1 = mms.get_data('B_dmpa_brst_l2',Tint1,1);
Bxyz2 = mms.get_data('B_dmpa_brst_l2',Tint2,1);
Bxyz3 = mms.get_data('B_dmpa_brst_l2',Tint3,1);
Bxyz4 = mms.get_data('B_dmpa_brst_l2',Tint4,1);
Bxyz5 = mms.get_data('B_dmpa_brst_l2',Tint5,1);

Exyz1 = mms.get_data('E_dsl_edp_brst_l2',Tint1,1);
Exyz2 = mms.get_data('E_dsl_edp_brst_l2',Tint2,1);
Exyz3 = mms.get_data('E_dsl_edp_brst_l2',Tint3,1);
Exyz4 = mms.get_data('E_dsl_edp_brst_l2',Tint4,1);
Exyz5 = mms.get_data('E_dsl_edp_brst_l2',Tint5,1);

Ve1 = mms.get_data('Ve_dbcs_fpi_brst_l2',Tint1,1);
Ve2 = mms.get_data('Ve_dbcs_fpi_brst_l2',Tint2,1);
Ve3 = mms.get_data('Ve_dbcs_fpi_brst_l2',Tint3,1);
Ve4 = mms.get_data('Ve_dbcs_fpi_brst_l2',Tint4,1);
Ve5 = mms.get_data('Ve_dbcs_fpi_brst_l2',Tint5,1);


%% Rotate of vectors

Bntt1 = irf_newxyz(Bxyz1,ndsl1,t1dsl1,t2dsl1);
Bntt2 = irf_newxyz(Bxyz2,ndsl2,t1dsl2,t2dsl2);
Bntt3 = irf_newxyz(Bxyz3,ndsl3,t1dsl3,t2dsl3);
Bntt4 = irf_newxyz(Bxyz4,ndsl4,t1dsl4,t2dsl4);
Bntt5 = irf_newxyz(Bxyz5,ndsl5,t1dsl5,t2dsl5);

Entt1 = irf_newxyz(Exyz1,ndsl1,t1dsl1,t2dsl1);
Entt2 = irf_newxyz(Exyz2,ndsl2,t1dsl2,t2dsl2);
Entt3 = irf_newxyz(Exyz3,ndsl3,t1dsl3,t2dsl3);
Entt4 = irf_newxyz(Exyz4,ndsl4,t1dsl4,t2dsl4);
Entt5 = irf_newxyz(Exyz5,ndsl5,t1dsl5,t2dsl5);

Ventt1 = irf_newxyz(Ve1,ndsl1,t1dsl1,t2dsl1);
Ventt2 = irf_newxyz(Ve2,ndsl2,t1dsl2,t2dsl2);
Ventt3 = irf_newxyz(Ve3,ndsl3,t1dsl3,t2dsl3);
Ventt4 = irf_newxyz(Ve4,ndsl4,t1dsl4,t2dsl4);
Ventt5 = irf_newxyz(Ve5,ndsl5,t1dsl5,t2dsl5);

Bntt1 = Bntt1.resample(Entt1);
Bntt2 = Bntt2.resample(Entt2);
Bntt3 = Bntt3.resample(Entt3);
Bntt4 = Bntt4.resample(Entt4);
Bntt5 = Bntt5.resample(Entt5);

Ventt1 = Ventt1.resample(Entt1);
Ventt2 = Ventt2.resample(Entt2);
Ventt3 = Ventt3.resample(Entt3);
Ventt4 = Ventt4.resample(Entt4);
Ventt5 = Ventt5.resample(Entt5);

Entt1 = Entt1.filt(0,50,8192,5);
Entt2 = Entt2.filt(0,100,8192,5);
Entt3 = Entt3.filt(0,10,8192,5);
Entt4 = Entt4.filt(0,5,8192,5);
Entt5 = Entt5.filt(0,5,8192,5);

%% Convert position to time
Times1 = Bntt1.time.epochUnix - Tints1(1).epochUnix;
Times2 = Bntt2.time.epochUnix - Tints2(1).epochUnix;
Times3 = Bntt3.time.epochUnix - Tints3(1).epochUnix;
Times4 = Bntt4.time.epochUnix - Tints4(1).epochUnix;
Times5 = Bntt5.time.epochUnix - Tints5(1).epochUnix;

npos1 = -Times1*Vsn1;
npos2 = -Times2*Vsn2;
npos3 = -Times3*Vsn3;
npos4 = -Times4*Vsn4;
npos5 = -Times5*Vsn5;

%% Figure 

h=irf_plot(5,'newfigure');
%h=irf_figure(540+ic,8);
xSize=700; ySize=600;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.88;
ywidth = 0.18;
set(h(1),'position',[0.095 0.985-ywidth xwidth ywidth]);
set(h(2),'position',[0.095 0.980-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.095 0.975-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.095 0.970-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.095 0.965-5*ywidth xwidth ywidth]);


h(1)=irf_panel('Bn');
plot(h(1),npos1,Bntt1.x.data,'k','linewidth',1.5)
hold(h(1),'on')
plot(h(1),npos2,Bntt2.x.data,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(1),npos3,Bntt3.x.data,'r','linewidth',1.5)
plot(h(1),npos4,Bntt4.x.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(1),npos5,Bntt5.x.data,'m','linewidth',1.5)
hold(h(1),'off')
axis(h(1),[-600 600 -8 30])
ylabel(h(1),'B_{n} (nT)')
set(h(1),'xticklabel',[])
irf_legend(h(1),'(a)',[0.99 0.97])

h(2)=irf_panel('Bt1');
plot(h(2),npos1,Bntt1.y.data,'k','linewidth',1.5)
hold(h(2),'on')
plot(h(2),npos2,Bntt2.y.data,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(2),npos3,Bntt3.y.data,'r','linewidth',1.5)
plot(h(2),npos4,Bntt4.y.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(2),npos5,Bntt5.y.data,'m','linewidth',1.5)
hold(h(2),'off')
axis(h(2),[-600 600 0 80])
ylabel(h(2),'B_{t1} (nT)')
set(h(2),'xticklabel',[])
irf_legend(h(2),'(b)',[0.99 0.97])

h(3)=irf_panel('Bt2');
plot(h(3),npos1,Bntt1.z.data,'k','linewidth',1.5)
hold(h(3),'on')
plot(h(3),npos2,Bntt2.z.data,'color',[0.3 0.3 1],'linewidth',1.5)
plot(h(3),npos3,Bntt3.z.data,'r','linewidth',1.5)
plot(h(3),npos4,Bntt4.z.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(3),npos5,Bntt5.z.data,'m','linewidth',1.5)
hold(h(3),'off')
axis(h(3),[-600 600 -14 14])
ylabel(h(3),'B_{t2} (nT)')
set(h(3),'xticklabel',[])
irf_legend(h(3),'(c)',[0.99 0.97])

h(4)=irf_panel('En');
plot(h(4),npos1,Entt1.x.data,'k','linewidth',1.5)
hold(h(4),'on')
plot(h(4),npos2,Entt2.x.data,'b','linewidth',1.5)
plot(h(4),npos3,Entt3.x.data,'r','linewidth',1.5)
plot(h(4),npos4,Entt4.x.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(4),npos5,Entt5.x.data,'m','linewidth',1.5)
hold(h(4),'off')
axis(h(4),[-600 600 -50 220])
ylabel(h(4),'E_n (mV m^{-1})')
legend(h(4),{'Shock 1','Shock 2','Shock 3','Shock 4','Shock 5'},'location','northwest')
set(h(4),'xticklabel',[])
irf_legend(h(4),'(d)',[0.99 0.97])

h(5)=irf_panel('Vet2');
plot(h(5),npos1,Ventt1.z.data,'k','linewidth',1.5)
hold(h(5),'on')
plot(h(5),npos2,Ventt2.z.data,'b','linewidth',1.5)
plot(h(5),npos3,Ventt3.z.data,'r','linewidth',1.5)
plot(h(5),npos4,Ventt4.z.data,'color',[0 0.6 0],'linewidth',1.5)
plot(h(5),npos5,Ventt5.z.data,'m','linewidth',1.5)
hold(h(5),'off')
axis(h(5),[-600 600 -500 2400])
ylabel(h(5),'V_{e,t2} (km s^{-1})')
xlabel(h(5),'n (km)')
irf_legend(h(5),'(e)',[0.99 0.97])

set(h(1:5),'fontsize',14);
