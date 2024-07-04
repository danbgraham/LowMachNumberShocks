%% load data
load('histstruct_shock1_l5_Bn0.mat')
modeldists1 = histstruct;
load('histstruct_shock1_l10_Bn0.mat')
modeldists2 = histstruct;
load('histstruct_shock1_l15_Bn0.mat')
modeldists3 = histstruct;
load('histstruct_shock1_l20_Bn0.mat')
modeldists4 = histstruct;
load('histstruct_shock1_l25_Bn0.mat')
modeldists5 = histstruct;
load('histstruct_shock1_l30_Bn0.mat')
modeldists6 = histstruct;
load('histstruct_shock1_l50_Bn0.mat')
modeldists7 = histstruct;
%% Calculate moemnts
moments1 = calculatemodelmoments(modeldists1);
moments2 = calculatemodelmoments(modeldists2);
moments3 = calculatemodelmoments(modeldists3);
moments4 = calculatemodelmoments(modeldists4);
moments5 = calculatemodelmoments(modeldists5);
moments6 = calculatemodelmoments(modeldists6);
moments7 = calculatemodelmoments(modeldists7);

Ts1 = moments1.Ts;
Ts2 = moments2.Ts;
Ts3 = moments3.Ts;
Ts4 = moments4.Ts;
Ts5 = moments5.Ts;
Ts6 = moments6.Ts;
Ts7 = moments7.Ts;

%% Calculate magnetic field
B01 = modeldists1.B0;
B11 = modeldists1.B1;
l1 = modeldists1.l;

B03 = modeldists3.B0;
B13 = modeldists3.B1;
l3 = modeldists3.l;

B06 = modeldists6.B0;
B16 = modeldists6.B1;
l6 = modeldists6.l;

xvec = [-20e2:1:5e2]*1e3;

By1 = -B01*tanh(xvec/l1)+B11;
By3 = -B03*tanh(xvec/l3)+B13;
By6 = -B06*tanh(xvec/l6)+B16;

xposp = modeldists1.xpositions;
vxvec = modeldists1.vxpositions;
vzvec = modeldists1.vzpositions;

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

h=irf_plot(8,'newfigure');
%h=irf_figure(540+ic,8);
xSize=700; ySize=800;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.80;
ywidth = 0.11;

set(h(1),'position',[0.1 0.990-ywidth xwidth ywidth]);
set(h(2),'position',[0.1 0.990-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.1 0.990-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.1 0.990-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.1 0.990-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.1 0.990-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.1 0.990-7*ywidth xwidth ywidth]);
set(h(8),'position',[0.1 0.940-8*ywidth xwidth 0.16]);

plot(h(1),xvec/1e3,By1*1e9,'color',[0 0.4470 0.7410])
hold(h(1),'on')
plot(h(1),xvec/1e3,By3*1e9,'color',[0.9290 0.6940 0.1250])
plot(h(1),xvec/1e3,By6*1e9,'color',[0.3010 0.7450 0.9330])
hold(h(1),'off')
%ylabel(h(1),{'B (nT)'},'fontsize',14)
grid(h(1),'on')
axis(h(1),[-2000 500 -5 80])
set(h(1),'xticklabel',[])
ylabel(h(1),'B_{t1} (nT)','fontsize',14)
irf_legend(h(1),'(a)',[0.01 0.98],'color','k','fontsize',14)
legend(h(1),'l = 5 km','l = 15 km','l = 30 km','location','southeast')

ftemp = modeldists1.n1Dvx;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(2),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(2),'flat')
axis(h(2),[-2000 500 -1800 1400])
ylabel(h(2),' ','FontSize',14)
caxis(h(2),[-3 2]);
colormap(h(2),cmap)
ylabel(h(2),'v_n (km s^{-1})','fontsize',14)
set(h(2),'xticklabel',[])
irf_legend(h(2),'(b)',[0.01 0.98],'color','k','fontsize',14)
irf_legend(h(2),'l = 5 km',[0.05 0.8],'color','k','fontsize',14)

ftemp = modeldists1.n1Dvz;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(3),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(3),'flat')
axis(h(3),[-2000 500 -1600 1600])
ylabel(h(3),' ','FontSize',14)
caxis(h(3),[-3 2]);
colormap(h(3),cmap)
ylabel(h(3),'v_{t2} (km s^{-1})','fontsize',14)
set(h(3),'xticklabel',[])
irf_legend(h(3),'(c)',[0.01 0.98],'color','k','fontsize',14)

ftemp = modeldists3.n1Dvx;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(4),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(4),'flat')
axis(h(4),[-2000 500 -1800 1400])
ylabel(h(4),' ','FontSize',14)
caxis(h(4),[-3 2]);
colormap(h(4),cmap)
ylabel(h(4),'v_n (km s^{-1})','fontsize',14)
set(h(4),'xticklabel',[])
irf_legend(h(4),'(d)',[0.01 0.98],'color','k','fontsize',14)
irf_legend(h(4),'l = 15 km',[0.05 0.8],'color','k','fontsize',14)

ftemp = modeldists3.n1Dvz;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(5),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(5),'flat')
axis(h(5),[-2000 500 -1600 1600])
ylabel(h(5),' ','FontSize',14)
caxis(h(5),[-3 2]);
colormap(h(5),cmap)
ylabel(h(5),'v_{t2} (km s^{-1})','fontsize',14)
set(h(5),'xticklabel',[])
irf_legend(h(5),'(e)',[0.01 0.98],'color','k','fontsize',14)

ftemp = modeldists6.n1Dvx;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(6),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(6),'flat')
axis(h(6),[-2000 500 -1800 1400])
ylabel(h(6),' ','FontSize',14)
caxis(h(6),[-3 2]);
colormap(h(6),cmap)
ylabel(h(6),'v_n (km s^{-1})','fontsize',14)
set(h(6),'xticklabel',[])
irf_legend(h(6),'(f)',[0.01 0.98],'color','k','fontsize',14)
irf_legend(h(6),'l = 30 km',[0.05 0.8],'color','k','fontsize',14)

ftemp = modeldists6.n1Dvz;
ftemp(ftemp < 1e-4) = NaN;
pcolor(h(7),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(7),'flat')
axis(h(7),[-2000 500 -1600 1600])
ylabel(h(7),' ','FontSize',14)
caxis(h(7),[-3 2]);
c = colorbar(h(7),'position',[0.910 0.990-7*ywidth+0.01 0.013 0.11*6-0.02]);
c.Label.String = 'log_{10}f_i (s m^{-4})';
colormap(h(7),cmap)
ylabel(h(7),'v_{t2} (km s^{-1})','fontsize',14)
set(h(7),'xticklabel',[])
irf_legend(h(7),'(g)',[0.01 0.98],'color','k','fontsize',14)

set(h([2:7]),'Color',0.75*[1 1 1]);

plot(h(8),xposp/1e3,Ts1,'color',[0 0.4470 0.7410])
hold(h(8),'on')
plot(h(8),xposp/1e3,Ts2,'color',[0.8500 0.3250 0.0980])
plot(h(8),xposp/1e3,Ts3,'color',[0.9290 0.6940 0.1250])
plot(h(8),xposp/1e3,Ts4,'color',[0.4940 0.1840 0.5560])
plot(h(8),xposp/1e3,Ts5,'color',[0.4660 0.6740 0.1880])
plot(h(8),xposp/1e3,Ts6,'color',[0.3010 0.7450 0.9330])
plot(h(8),xposp/1e3,Ts7,'color',[0.6350 0.0780 0.1840])
hold(h(8),'off')
%ylabel(h(8),{'B (nT)'},'fontsize',14)
grid(h(8),'on')
axis(h(8),[-2000 500 0 40])
%set(h(8),'xticklabel',[])
xlabel(h(8),'n (km)')
ylabel(h(8),'T_{p} (eV)','fontsize',14)
irf_legend(h(8),'(h)',[0.01 0.98],'color','k','fontsize',14)
%legend(h(8),'l = 5 km','l = 10 km','l = 15 km','l = 20 km','l = 25 km','l = 30 km','l = 50 km','location','northeast')
irf_legend(h(8),'l = 5 km',[1.01 0.99],'color',[0 0.4470 0.7410],'fontsize',12)
irf_legend(h(8),'l = 10 km',[1.01 0.84],'color',[0.8500 0.3250 0.0980],'fontsize',12)
irf_legend(h(8),'l = 15 km',[1.01 0.69],'color',[0.9290 0.6940 0.1250],'fontsize',12)
irf_legend(h(8),'l = 20 nT',[1.01 0.54],'color',[0.4940 0.1840 0.5560],'fontsize',12)
irf_legend(h(8),'l = 25 nT',[1.01 0.30],'color',[0.4660 0.6740 0.1880],'fontsize',12)
irf_legend(h(8),'l = 30 nT',[1.01 0.15],'color',[0.3010 0.7450 0.9330],'fontsize',12)
irf_legend(h(8),'l = 50 nT',[1.01 0.00],'color',[0.6350 0.0780 0.1840],'fontsize',12)

irf_plot_axis_align(h(1:8));

set(h(1:8),'xdir','reverse')