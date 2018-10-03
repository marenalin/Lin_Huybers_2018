addpath('./TRMM_interp_3am/');
addpath('./Functions');
addpath('./Data');
addpath('./Interp_IMD');
load IMD_rf;%loads variables rf25,stn25,y25,m25,d25,lonr,latr, which give
            %the 0.25x0.25 daily data from 1901-2016
load IMD_rf_1x1;%loads variables ind_rf_cut,
                %ind_stn_cut,lonm_cut,latm_cut, mm,yy,dd for the
                %1x1 Rajeevan data, which only use stations with
                %90% of the stations reporting
load TRMM_mic;
ym=unique(TRMM_mic(1).yy);
load stn25_interp_dist_all; %loads all of the interpolation
                            %distances daily for each gridbox


%%check temporal changes to station counts for the 0.25 degree and
%%1 degree data over Central India
ax=find(lonr>=74.5 & lonr<=86.5);
ay=find(latr>=16.5 & latr<=26.5);
sk=intersect(ax,ay); %sk is then the gridboxes in 0.25x0.25 data
                     %that are in central India
ax=find(lonm_cut>=74.5 & lonm_cut<=86.5);
ay=find(latm_cut>=16.5 & latm_cut<=26.5);
tk=intersect(ax,ay); %tk is then the gridboxes in 1x1 data
                     %that are in central India

s25=nansum(stn25(sk,:),1);
s1=nansum(ind_stn_cut(:,tk),2);
yx=1951:2016;
mt=find(m25>5 & m25<10);
jt=find(mm>5 & mm<10);
for ct1=1:length(yx)
    z1=find(yx(ct1)==y25);
    z1=intersect(z1,mt);
    san(ct1,1)=nanmean(s25(z1));
    
    z1=find(yx(ct1)==yy);
    z1=intersect(z1,jt);
    san(ct1,2)=nanmean(s1(z1));
end

gg=figure;
gg.Units='inches';
gg.Position=[2 2 14 8]; 

plot(yx,san,'linewidth',2);
axis tight;
hold on;
xlabel('Monsoon season','fontsize',14);
ylabel('Average daily station total over Central India','fontsize', ...
       14);
set(gca,'fontsize',16);
legend('0.25^\circ x0.25^\circ','1^\circ x1^\circ');
figpos = get(gcf, 'Position');
paperWidth = figpos(3);
paperHeight = figpos(4);
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print('-dpdf','Station_totals_IMD');

%first regrid TRMM onto the IMD grid using linear
%interpolation, integrates 3 AM UTC to 3 AM UTC, assuming 12 AM 
TRMM_rg=struct();
nstrl=[];
for ct1=1:length(TRMM_mic(1).yy)
    nstrl(ct1)=datenum(TRMM_mic(1).yy(ct1),TRMM_mic(1).mm(ct1), ...
                       TRMM_mic(1).dd(ct1),TRMM_mic(1).hh(ct1),0,0);
end

[AA BB]=sort(nstrl);

%TRMM_mic(1).yy(BB) is now ordered chronologically

vk=find(m25>5 & m25<10);
for ctm=1:length(ym)
    yk=find(y25==ym(ctm));     
    yk=intersect(yk,vk);        
    rk=find(ym(ctm)==TRMM_mic(1).yy);
    for ct1=1:length(yk)
        ex=find(TRMM_mic(1).mm(rk)==m25(yk(ct1)));
        dx=find(TRMM_mic(1).dd(rk)==d25(yk(ct1)));
        ex=intersect(ex,dx);
        dx=find(TRMM_mic(1).hh(rk)==0); %midnight UTC on the
                                        %collection date, so
        ex=intersect(ex,dx);
        ex=rk(ex);
        ex=find(BB==ex); %finds the location of this index in
                         %the reordering
        
        mx=BB((ex-6):ex); %middle indices for hours midnight of
                          %collection date, and preceding 5
                          %3-hr windows
        ns=nansum(3*TRMM_mic(1).vals(:,mx),2);
        %bookends, hours 12am-1:30 am of the collection date,
        %and 1:30 am-3 am of the previous date
        ns=nansum([ns 1.5*TRMM_mic(1).vals(:,BB(ex-7)) 1.5*TRMM_mic(1).vals(:,BB(ex+1))],2);
        rf25_TRMM(ct1,1:length(lon_TRMM))=ns;
    end


    rf25_TRMM_rg=nan(length(yk),length(lonr));
    for ct1=1:length(yk)
        temp=rf25_TRMM(ct1,:);
        kl=find(isnan(temp)==0);
        rf25_TRMM_rg(ct1,1:length(lonr))= griddata(lon_TRMM(kl), ...
                                                   lat_TRMM(kl), ...
                                                   double(temp(kl)),lonr,latr);
        %linear interpolation of the TRMM data onto the IMD grid
    end
    TRMM_rg(ctm).rf=rf25_TRMM_rg;
    ctm
end

    

%load all produced data and look at statistics of missing
%events and using the 99.6th percentiles
rfal=struct();
for ct1=1:length(ym)
    pz=num2str(ym(ct1));
    load(['TRMM_interp_3am_' pz '.mat']);
    rfal(ct1).grd=rf25_interp;
end


l1=find(m25>5 & m25<10);
l2=find(y25>1950);
l1=intersect(l1,l2);

yx=1951:2016;
ax=find(lonr>=74.5 & lonr<=86.5);
ay=find(latr>=16.5 & latr<=26.5);
ak=intersect(ax,ay);
for ct1=1:length(ym)
    temp=rfal(ct1).grd(:,ak);
    temp=temp(:);
    lm=prctile(temp,[99.6]); %find the 99.6th percentile of each
                             %set of data

    for ct2=1:length(yx)
        lr=find(yx(ct2)==y25(l1));
        temp=rfal(ct1).grd(lr,ak);
        temp=temp(:);
        vs(ct1,ct2)=length(find(temp>lm));
        ds(ct1,ct2)=length(find(temp>100));
    end
end

%anomaly of extreme event counts >99.6th percentile rainfall for
%each year yx
mva=vs-repmat(mean(vs,2),1,length(yx));


mvt=mean(mva);

bpps=2:(length(yx)-1);
for ct1=1:length(bpps);
    m1=mean(mvt(1:bpps(ct1)));
    m2=mean(mvt((bpps(ct1)+1):end));
    res=[mvt(1:bpps(ct1))-m1 mvt((bpps(ct1)+1):end)-m2];
    srr(ct1)=sum(res.^2);
end


%check the breakpoint of each of the time-series separately
for ctx=1:size(mva,1)
    mvk=mva(ctx,:);
    bpps=2:(length(yx)-1);
    srr=[];
    for ct1=1:length(bpps);
        m1=mean(mvk(1:bpps(ct1)));
        m2=mean(mvk((bpps(ct1)+1):end));
        res=[mvk(1:bpps(ct1))-m1 mvk((bpps(ct1)+1):end)-m2];
        srr(ct1)=sum(res.^2);
    end
    bpx(ctx)=find(min(srr)==srr);
    ma1(ctx)=mean(mvk(1:25));
    ma2(ctx)=mean(mvk(26:end));;
end

%find the index with the largest offset, across 1975
find((ma2-ma1)==max(ma2-ma1))



%find the number of stations 
mf=find(m25>5 & m25<10);
yf=find(y25>1950);
yf=intersect(yf,mf);
tstn=stn25(ak,yf);

%fraction of gridbox days that have no stations
length(find(tstn(:))==0)/length(tstn(:))

addpath('./Interp_IMD/');
%load the re-interpolated IMD data
ax=find(lonr>=74.5 & lonr<=86.5);
ay=find(latr>=16.5 & latr<=26.5);
ak=intersect(ax,ay);
yb=zeros(length(yx),1);
cvals=nan(length(yx),length(yx),2);
cval_anom=nan(length(yx),length(yx),2);
for ctn=1:length(yx)
    ry=0;
    try
        load(['IMD_interp_rep_n_' num2str(yx(ctn)) '.mat']);
        ry=1;
    end
    if ry==1

        for ct1=1:length(yx)
            temp=rfint(ct1,:,ak);
            temp=temp(:);
            cvals(ctn,ct1,1)=length(find(temp>100));
            cvals(ctn,ct1,2)=length(find(temp>150));
        end
        cval_anom(ctn,:,1)=cvals(ctn,:,1)-mean(cvals(ctn,:,1));
        cval_anom(ctn,:,2)=cvals(ctn,:,2)-mean(cvals(ctn,:,2));
        yb(ctn)=1;
    end
end



%check count of boxes with interpolation distances over 0.5^\circ and find the breakpoint
for ct1=1:length(ym)
    rx=find(ym(ct1)==y25);
    rx=intersect(yf,rx);
    temp=stn25_avd(ak,rx);
    temp=temp(:);
    mdis(ct1)=mean(temp);
    frav(ct1)=length(find(temp>0.5))/length(temp);
end


bx=2:length(frav);
X=[ones(length(ym),1) (ym-mean(ym))'];
for ct1=1:length(bx)
    v1=frav(1:bx(ct1));
    v2=frav((bx(ct1)+1):end);
    Yp=[mean(v1)*ones(length(v1),1); mean(v2)*ones(length(v2), ...
                                                   1)];
    res=frav'-Yp;
    ssx(ct1)=sum(res.^2);
end

%%FIG 3
MX=prctile(squeeze(cval_anom(:,:,1)),[2.5 97.5]);
MA=max(squeeze(cval_anom(:,:,1)));
MI=min(squeeze(cval_anom(:,:,1)));
gg=figure;
gg.Units='inches';
gg.Position=[2 2 14 8]; 
px=gca;%subplot(2,2,[1 2]);
kv=get(px,'Position');
set(px,'Position',[kv(1) kv(2)*1.1 kv(3:4)]);
Ygg=[MI fliplr(MA)];
Xgg=[yx fliplr(yx)];
gp=patch(Xgg,Ygg,[0.7 0.7 0.7]);
set(gp,'EdgeColor','None');
Ygg=[MX(1,:) fliplr(MX(2,:))];
Xgg=[yx fliplr(yx)];
gp=patch(Xgg,Ygg,[0.8 0.8 0.8]);
set(gp,'EdgeColor','None');


hold on;
cl=brewermap(length(1998:2015),'Spectral');
for ct1=1:length(ym)
plot(yx,mva(ct1,:),'linewidth',1.5,'Color',cl(ct1,:))
end
box on; axis tight;
plot(yx,mvt,'linewidth',4,'Color',[0 0 0]);

m2=mean(mvt(26:end));
m1=mean(mvt(1:25));

plot([yx(1) yx(25)],[m1 m1],'k--','linewidth',4,'Color',[0 0 0]);
plot([yx(25) yx(end)],[m2 m2],'k--','linewidth',4,'Color',[0 0 0]); 

set(gca,'linewidth',2,'fontsize',16);



h=colorbar;
colormap(gca,cl);
caxis([1997.5 2015.5]);
set(h,'YTick',1998:2015,'YTickLabel',1998:2015);
ylabel(h,'TRMM monsoon season');
set(h,'Linewidth',2,'fontsize',16);
xlabel('Year of IMD station data','fontsize',18);
ylabel('Extreme event count (anomaly)','fontsize',18);
%text(1946,160,'(a)','fontsize',18);

figpos = get(gcf, 'Position');
paperWidth = figpos(3);
paperHeight = figpos(4);
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print('-dpdf','Fig_3');


mimd=mean(squeeze(cval_anom(:,:,1)),1);
bx=2:length(mimd);
X=[ones(length(ym),1) (ym-mean(ym))'];
for ct1=1:length(bx)
    v1=mimd(1:bx(ct1));
    v2=mimd((bx(ct1)+1):end);
    Yp=[mean(v1)*ones(length(v1),1); mean(v2)*ones(length(v2), ...
                                                   1)];
    res=mimd'-Yp;
    ssx(ct1)=sum(res.^2);
end


%figure out what year gives you the highest offset in
%the IMD data, using 1973 as the breakpoint
os1=mean(squeeze(cval_anom(:,1:23,1)),2); %using 1973 as a breakpoint
os2=mean(squeeze(cval_anom(:,24:end,1)),2);
find(max(os2-os1)==(os2-os1)) %count of 91



mk=find(m25>5& m25<10);
lk=find(y25>1950);
lk=intersect(mk,lk);
for ct1=1:length(yx);
    oj=find(yx(ct1)==y25);
    oj=intersect(mk,oj);
    pj=stn25_avd(ak,oj);
    pj=pj(:);
    vprc(ct1,1:3)=prctile(pj,[50 70 90]);
    vmean(ct1)=mean(pj);
end
avr=nanmean(rf25(:,lk),2); %average daily monsoon rainfall





%check variance of the original time-series
for ct1=1:length(ym)
    temp=TRMM_rg(ct1).rf(:,ak);  
    temp=temp(:);
    gva(ct1)=nanvar(temp);
end


%changes to mean rainfall in TRMM 
for ctk=1:length(ym)
    for ct1=1:length(yx)
        rb=find(yx(ct1)==y25(lk));
        temp=rfal(ctk).grd(rb,ak);
        temp=temp(:);
        anvr(ctk,ct1)=nanmean(temp);
    end
end

X=[ones(size(yx(:))) yx(:)-mean(yx)];
for ctx=1:length(ym)
     sl=multilin_lsq(X,anvr(ctx,:)');
     IMDmn(ctx)=nanmean(nanmean(rf25(ak,l1(y25(l1)==ym(ctx)))));
     TRMMmn(ctx)=nanmean(nanmean(TRMM_rg(ctx).rf(:,ak))); 
     scl(ctx)=sl(2)*50*(IMDmn(ctx)/TRMMmn(ctx)); %scale for
                                                 %proportionally
                                                 %lower rainfall
                                                 %values in TRMM
                                                 %versus in IMD
end

%check whether the jump/offset also occurs in the average
%interpolation distance
bs=2:length(yx);
for ct1=1:length(bs);
    i1=1:bs(ct1);
    i2=(bs(ct1)+1):length(yx);
    v1=mean(vmean(i1));
    v2=mean(vmean(i2));
    prv=[v1*ones(1,length(i1)) v2*ones(1,length(i2))];
    rrsq(ct1)=sum((vmean-prv).^2);
end
find(rrsq==min(rrsq))
%corresponds to year 1975

%count of >100 mm events for each gridpoint in central India
temp=rf25(:,lk);
for ct1=1:length(lonr)
    rb=find(temp(ct1,:)>100);
    hl(ct1)=length(rb);
end

p1=find(m25<6);
p2=find(m25>9);
p1=union(p1,p2);
p1=intersect(p1,find(y25>1950));
p1=find(y25>1950);

%proportion of total annual rainfall that falls during the
%monsoon
nansum(nansum(rf25(ak,lk)))/nansum(nansum(rf25(ak,p1)))

for ct1=1:length(lk)
    %count the 150mm 'events' on every day
    temp=rf25(ak,lk(ct1));
    ct(ct1)=length(find(temp>150));
end

dd=find(y25(lk)==2006);
cumsum(fliplr(sort(ct(dd))))
[AA BB]=sort(ct(dd));
AA=fliplr(AA);
BB=fliplr(BB);


dd=find(y25(lk)==2006);
ff=find(m25(lk)>6 & m25(lk)<9);
ff=intersect(dd,ff);
gg=find(d25(lk)>1& d25(lk)<6);
gg=intersect(gg,ff);
sum(ct(gg)); %this gives the count of E_y that can be attributed
%to two monsoon depressions, July 2-5 and Aug 2-5, 2006





b2=find(y25>1950 & y25<1976);
b2=intersect(yf,b2);
b3=find(y25>1975);
b3=intersect(yf,b3);

av1=nanmean(stn25_avd(:,b2),2);
av2=nanmean(stn25_avd(:,b3),2);
gg=figure;
gg.Units='inches';
gg.Position=[2 2 14 20]; 

pl1=subplot(3,2,5);
hold on;
nax=gca;
lc=get(nax,'Position');
ww=lc(3);
hh=lc(4);
vals=av1;
s1=[0 0.6];
tv=linspace(s1(1),s1(2),100);
blp=brewermap(100,'Spectral');
col=[];
for ct2=1:length(lonr)
    fg=find(vals(ct2)<=tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end
yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);


yy=[25 30 30 25 25];
xx=[70 70 90 90 70];
pd=patch(xx,yy,[1 1 1]);
set(pd,'EdgeColor','none');
colw=[1 1 1];
colw=repmat(colw,length(lonr),1,1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,colw,pl1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);    

cY=[16.5 16.5 26.5 26.5 16.5];
cX=[86.5 74.5 74.5 86.5 86.5];

plot(cX,cY,'linewidth',3,'Color',[0 0 0]);

grid off;
box on;
axis([72 87 15 27]);

set(gca,'linewidth',4,'fontsize',16);
df=get(nax,'Xtick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);
text(70,27.5,'(d) 1951-1975 average','fontsize',18);
set(gca,'Position',[lc(1) lc(2)*0.5 ww hh],'linewidth',4);

pl1=subplot(3,2,6);
hold on;
fc=get(gca,'Position');
nax=gca;
set(gca,'Position',[fc(1) fc(2) hh ww],'linewidth',4);
vals=av2;
s1=[0 0.6];
tv=linspace(s1(1),s1(2),100);
blp=brewermap(100,'Spectral');
col=[];
for ct2=1:length(lonr)
    fg=find(vals(ct2)<=tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end

yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);

yy=[25 30 30 25 25];
xx=[70 70 90 90 70];
pd=patch(xx,yy,[1 1 1]);
set(pd,'EdgeColor','none');

colw=[1 1 1];
colw=repmat(colw,length(lonr),1,1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,colw,pl1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);

plot(cX,cY,'linewidth',3,'Color',[0 0 0]);
set(gca,'linewidth',2,'fontsize',16);
box on; grid off;
colormap(brewermap(100,'Spectral'));
caxis([s1(1) s1(2)]);
axis([72 87 15 27]);
grid on; box on;

set(gca,'linewidth',2,'fontsize',16);
df=get(gca,'Xtick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);


pc=colorbar('Location','Northoutside');
set(gca,'Position',[fc(1) fc(2)*0.5 ww hh],'linewidth',4);


jc=get(pc,'Position');
set(pc,'Position',[lc(1) jc(2)*1.05 ww*2.32 0.02 ]);
set(pc,'linewidth',2,'Xtick',0:0.1:0.6,'XtickLabel',{'0^\circ' '0.1^\circ' '0.2^\circ' ...
'0.3^\circ' '0.4^\circ' '0.5^\circ' ...
'0.6^\circ'},'fontsize',14,'linewidth',2);

text(70,27.5,'(e) 1976-2016 average','fontsize',18);

xlabel(pc,'Average interpolation distance');

px=subplot(3,3,4:6);
hold on;
kv=get(px,'Position');
set(px,'Position',[kv(1) kv(2) kv(3) kv(4)*1.1]); 
box on; 
plot(yx,vmean,'linewidth',4,'Color',[0 0 0]);
axis tight;
i1=1:25;
i2=26:length(yx);
v1=mean(vmean(i1));
v2=mean(vmean(i2));
prv=[v1*ones(1,length(i1)) v2*ones(1,length(i2))];

plot(yx(1:25),prv(1:25),'k--','linewidth',4);
plot(yx(26:end),prv(26:end),'k--','linewidth',4);
set(gca,'linewidth',2,'fontsize',16);
xlabel('Year','fontsize',18);
ylabel('Average interpolation distance','fontsize',18);
yt=get(gca,'Ytick');
Ylab=[];
for ct1=1:length(yt)
    Ylab{ct1}=[num2str(yt(ct1)) '^\circ'];
end
set(gca,'Yticklabel',Ylab);
text(1944,0.35,'(c)','fontsize',18);

pl1=subplot(3,2,1);
hold on;
nax=gca;
lc=get(gca,'Position');
set(gca,'Position',[lc(1) lc(2) ww hh ]);
vals=avr;
s1=[4 12];
tv=linspace(s1(1),s1(2),100);
blp=brewermap(100,'Spectral');
col=[];
for ct2=1:length(lonr)
    fg=find(vals(ct2)<=tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end

yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);

yy=[25 30 30 25 25];
xx=[70 70 90 90 70];
pd=patch(xx,yy,[1 1 1]);
set(pd,'EdgeColor','none');

colw=[1 1 1];
colw=repmat(colw,length(lonr),1,1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,colw,pl1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);    
cY=[16.5 16.5 26.5 26.5 16.5];
cX=[86.5 74.5 74.5 86.5 86.5];

plot(cX,cY,'linewidth',3,'Color',[0 0 0]);

grid off;
box on;
axis([72 87 15 27]);

set(gca,'linewidth',4,'fontsize',16);
df=get(nax,'Xtick');
dfs=[];
for ct2=1:length(df)
    dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
    dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);
text(69,28,'(a)','fontsize',18);
h=colorbar('location','NorthOutside');
xlabel(h,'Average daily rainfall (mm)','fontsize',14);
caxis([s1(1) s1(2)]);
colormap(brewermap(100,'Spectral'));
set(h,'linewidth',2,'fontsize',16);
cpos=get(h,'Position');
set(h,'Position',[cpos(1) cpos(2:3) 0.02]);
set(gca,'Position',[lc(1) lc(2)*0.95 ww hh ]);

pl1=subplot(3,2,2);
hold on;
fc=get(gca,'Position');
nax=gca;
set(gca,'Position',[fc(1) fc(2) ww hh],'linewidth',4);
vals=hl;
s1=[0 80];
tv=linspace(s1(1),s1(2),100);
blp=brewermap(100,'Spectral');
col=[];
for ct2=1:length(lonr)
    fg=find(vals(ct2)<=tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end

yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);

yy=[25 30 30 25 25];
xx=[70 70 90 90 70];
pd=patch(xx,yy,[1 1 1]);
set(pd,'EdgeColor','none');



colw=[1 1 1];
colw=repmat(colw,length(lonr),1,1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,colw,pl1);
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);

plot(cX,cY,'linewidth',3,'Color',[0 0 0]);


set(gca,'linewidth',2,'fontsize',16);
box on; grid off;
colormap(brewermap(100,'Spectral'));
axis([72 87 15 27]);
grid on; box on;

df=get(gca,'Xtick');
dfs=[];
for ct2=1:length(df)
    dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
    dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);



text(69,28,'(b)','fontsize',18);
h=colorbar('location','NorthOutside');
xlabel(h,'Total count of extreme events (>100 mm day^{-1})','fontsize',13);
caxis([s1(1) s1(2)]);
set(h,'linewidth',2,'fontsize',16);
cpos=get(h,'Position');
set(h,'Position',[cpos(1) cpos(2:3) 0.02]);
set(gca,'Position',[fc(1) fc(2)*0.95 ww hh],'linewidth',4);

figpos = get(gcf, 'Position');
paperWidth = figpos(3);
paperHeight = figpos(4);
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print('-dpdf','Fig_1');

%%Fig 4

rf25(:,ak);
bb=find(m25>5 & m25<10);
for ct1=1:length(yx)
    aa=find(y25==yx(ct1));
    aa=intersect(aa,bb);
    ms=rf25(ak,aa);
    ctex(ct1,1)=length(find(ms(:)>100));
    ctex(ct1,2)=length(find(ms(:)>150));
end


%an offset of 130 would lead to an insignificant trend
nn=10000;
bp=1974;
d1=find(yx>bp);
dp=find(yx<=bp);
Cte=[ctex(dp,1); ctex(d1,1)-130];
X=[ones(size(yx(:))) yx(:)-mean(yx(:))];
Xt=X;

sl1=multilin_lsq(Xt,Cte);
Yd1=X(d1,:)*sl1;

slt=nan(length(nn),1);
for ct1=1:nn
    Y1=randsample(Cte,length(Cte),'true');
    sl=multilin_lsq(X,Y1);
    slt(ct1)=sl(2);
end

rcr(1:2)=prctile(slt,[2.5 97.5]);




   

nn=10000;
%plot up the shuffling time-series
d1=find(yx>=1973); 
X=[ones(size(yx(:))) yx(:)-mean(yx(:))];
Xt=X;
Xt(d1,2)=X(d1,2)-mean(X(d1,2));
sl1=multilin_lsq(X(d1,:),ctex(d1,1));
Yd1=X(d1,:)*sl1;
sl2=multilin_lsq(X(d1,:),ctex(d1,2));
Yd2=X(d1,:)*sl2;


for ct1=1:nn
    Y1=randsample(ctex(d1,1),length(d1),'true');
    Y2=randsample(ctex(d1,2),length(d1),'true');
    sl=multilin_lsq(X(d1,:),Y1);
    slt(ct1,1)=sl(2);
    sl=multilin_lsq(X(d1,:),Y2);
    slt(ct1,2)=sl(2);
end

length(find(slt(:,1)>sl1(2)))*2/nn  
length(find(slt(:,2)>sl2(2)))*2/nn  

rc(1,1:2)=prctile(slt(:,1),[5 95]);
rc(2,1:2)=prctile(slt(:,2),[5 95]);

rc(1,1:2)=prctile(slt(:,1),[2.5 97.5]);
rc(2,1:2)=prctile(slt(:,2),[2.5 97.5]);


gg=figure;
gg.Units='inches';
gg.Position=[2 2 10 5]; 
hold on;

mv1=mean(ctex(d1,1));
fp(1:2,1)=mv1+rc(1,1:2)*Xt(d1(end),2);
fp(1:2,2)=mv1+rc(1,1:2)*Xt(d1(1),2);
Xtr=[yx(d1(1)) yx(d1(1)) mean(yx(d1))];
Ytr=[fp(1,1) fp(2,1) mean(ctex(d1,1))];
hg=patch(Xtr,Ytr,[0.8 0.8 0.8]);
set(hg,'EdgeColor','None');

Xtr=[mean(yx(d1)) yx(d1(end)) yx(d1(end)) ];
Ytr=[mean(ctex(d1,1)) fp(1,2) fp(2,2) ];
hg=patch(Xtr,Ytr,[0.8 0.8 0.8]);
set(hg,'EdgeColor','None');
a1=plot(yx,ctex(:,1),'linewidth',2);
plot(yx(d1),Yd1,'k--','linewidth',2);

%742 events attributed to July2-5, Aug2-5
plot([2006 2006], [max(ctex(:,1)) max(ctex(:,1))-742],'k', ...
     'linewidth',2);

scatter([2006],[max(ctex(:,1))],50,[0 0 ...
                    0],'filled','v');

scatter(2006,max(ctex(:,1))-742,50,[0 0 0 ],'^','filled');
text(2002,1300,'Jul 2-5, 2006','fontsize',16,'Rotation',90);
text(2004, 1300, 'Aug 2-5, 2006','fontsize',16,'Rotation',90);


plot([2006 2006], [max(ctex(:,2)) max(ctex(:,2))-383],'k', ...
     'linewidth',2);

scatter([2006],[max(ctex(:,2))],50,[0 0 ...
                    0],'filled','v');

scatter(2006,max(ctex(:,2))-383,50,[0 0 0 ],'^','filled');



mv2=mean(ctex(d1,2));
fp(1:2,1)=mv2+rc(2,1:2)*Xt(d1(end),2);
fp(1:2,2)=mv2+rc(2,1:2)*Xt(d1(1),2);



Xtr=[yx(d1(1)) yx(d1(1)) mean(yx(d1))];
Ytr=[fp(1,1) fp(2,1) mean(ctex(d1,2))];
hg=patch(Xtr,Ytr,[0.8 0.8 0.8]);
set(hg,'EdgeColor','None');

Xtr=[mean(yx(d1)) yx(d1(end)) yx(d1(end)) ];
Ytr=[mean(ctex(d1,2)) fp(1,2) fp(2,2) ];
hg=patch(Xtr,Ytr,[0.8 0.8 0.8]);
set(hg,'EdgeColor','None');

a2=plot(yx,ctex(:,2),'linewidth',2);
set(gca,'linewidth',2,'fontsize',16);
plot(yx(d1),Yd2,'k--','linewidth',2);
box on;
xlabel('Year','fontsize',16);
ylabel('Extreme event count','fontsize',16);
py=legend([a1 a2],'>100 mm day^{-1}','>150 mm day^{-1}');
set(py,'fontsize',16,'linewidth',2,'Location','Best');
figpos = get(gcf, 'Position');
paperWidth = figpos(3);
paperHeight = figpos(4);
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print('-dpdf',['Fig_4']);




%%FIG 2

load TRMM_interp_3am_2007;
tk25=ak;
[yt mt dt]=datevec(datenum(2000,6,1):datenum(2000,9,30));
gg=figure;
gg.Units='inches';
gg.Position=[2 2 14 8]; 
pl1=subplot(2,3,1);
dm=8;
rb=0
dn=7;
hx=find(mt==dm);
kx=find(dt==dn);
hx=intersect(hx,kx);

vals=TRMM_rg(10).rf(hx,:);
kf=find(vals(tk25)>100);
cmm(1)=length(find(vals(tk25)>100)) %60
cc=find(isnan(vals)==0);
rlim=120;
s1=[rb rlim];
csch='YlGnBu';
nax=gca;
tv=linspace(s1(1),s1(2),100);
blp=brewermap(100,csch);
col=[];
cc=find(isnan(vals)==0);
for ct2=1:length(lonr)
    fg=find(vals(ct2)<=tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end
yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);
hold on;
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);   
hold on;
plot(lonr(tk25(kf)), latr(tk25(kf)),'rx','Markersize',5,'MarkerFaceColor',[1 1 1]);
grid off;
box on;
axis([72 88 16 30]);
set(gca,'fontsize',15,'linewidth',2);
posg=get(gca,'Position');
h=colorbar('Location','NorthOutside');
set(h,'linewidth',2,'fontsize',16);
xlabel(h,['August ' num2str(dn) ' rainfall (mm day^{-1})']);
colormap(gca,blp);
text(70,32,'(a)','fontsize',18);
text(82,29,'TRMM 2007','fontsize',18);
text(82,28,[num2str(cmm(1)) ' extremes'],'fontsize',18);
caxis(s1);
set(gca,'Position',[posg(1)*0.35 posg(2)*0.1 posg(3:4)*1.25]);
cpos=get(h,'Position');
set(h,'Position',[cpos(1) cpos(2)*1.2 cpos(3:4)]);

df=get(nax,'Xtick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);


pl1=subplot(2,3,2);
nax=gca;
yp=find(y25==1971); dp=find(d25==dn);
cp=find(m25==dm);cp=intersect(cp,dp); cp=intersect(cp,yp);
vals=stn25(:,cp);
cc=find(isnan(vals)==0);
s1=[0 4];
tv=1:4;
blp=brewermap(4,'OrRd');
col=[];
cc=find(isnan(vals)==0);
for ct2=1:length(lonr)
    fg=find(vals(ct2)==tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end
col(vals==0,:)=0;
yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);
hold on;
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);   
dpos=get(gca,'Position');
hold on;
grid off;
box on;
text(67.5,30,'(b)','fontsize',18);
text(84,28.5,'1971','fontsize',18);
axis([72 88 16 30]);
set(gca,'fontsize',15,'linewidth',2);
dpos=get(gca,'Position');
h=colorbar('Location','SouthOutside');
set(h,'linewidth',2,'fontsize',16);
set(h,'Position',[cpos(1) cpos(2)*1.5 cpos(3:4)]); 
xlabel(h,['IMD August ' num2str(dn) ' station count']);
set(gca,'Position',[dpos(1)*0.95 dpos(2)*0.98 dpos(3:4)*1.2]);
blp=[0 0 0 ; blp];
colormap(gca,blp);
caxis([-0.5 s1(2)+0.5]);

df=get(nax,'Xtick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);

%station panel 1, 
yp=find(y25==1993); dp=find(d25==dn);
cp=find(m25==dm);cp=intersect(cp,dp); cp=intersect(cp,yp);
pl1=subplot(2,3,3);
vals=stn25(:,cp);
cc=find(isnan(vals)==0);
s1=[0 4];
nax=gca;
tv=1:4;
blp=brewermap(4,'OrRd');
col=[];
cc=find(isnan(vals)==0);
for ct2=1:length(lonr)
    fg=find(vals(ct2)==tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end
col(vals==0,:)=0;
yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);
hold on;

[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);   
hold on;
grid off;
box on;
text(67.5,30,'(c)','fontsize',18);
text(84,28.5,'1993','fontsize',18);
axis([72 88 16 30]);
set(gca,'fontsize',15,'linewidth',2);
rpos=get(gca,'Position');
set(gca,'Position',[rpos(1)*1.05 rpos(2)*0.98 dpos(3:4)*1.2])

df=get(nax,'Xtick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);

aa=find(y25>1950);
bb=find(m25>5 & m25<10);
aa=intersect(aa,bb);


a1=find(y25(aa)==1971);
a2=find(m25(aa)==dm);
a1=intersect(a1,a2);
a3=find(d25(aa)==dn);
a1=intersect(a1,a3);
vals=rf25_interp(a1,:);
kf=find(vals(tk25)>100);
cmm(2)=length(find(vals(tk25)>100)) %==33
a1
pl1=subplot(2,3,5);
cc=find(isnan(vals)==0);
s1=[rb rlim];
nax=gca;

tv=linspace(s1(1),s1(2),100);
blp=brewermap(100,csch);
col=[];
cc=find(isnan(vals)==0);
for ct2=1:length(lonr)
    fg=find(vals(ct2)<=tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end
yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);
hold on;
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);   
hold on;
plot(lonr(tk25(kf)), latr(tk25(kf)),'rx','Markersize',5,'MarkerFaceColor',[1 1 1]);
grid off;
text(81.5,28.5,[ num2str(cmm(2)) ' extremes'],'fontsize',18);
rpos=get(gca,'Position');
set(gca,'Position',[rpos(1)*0.95 rpos(2)*0.4 dpos(3:4)*1.2])
box on;
axis([72 88 16 30]);
set(gca,'fontsize',15,'linewidth',2);
posg=get(gca,'Position');
df=get(nax,'Xtick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);
text(67.5,30,'(d)','fontsize',18);

df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);

text(78,31.5,'\downarrow','fontsize',50,'fontweight','bold');

a1=find(y25(aa)==1993);
a2=find(m25(aa)==dm);
a1=intersect(a1,a2);
a3=find(d25(aa)==dn);
a1=intersect(a1,a3);
vals=rf25_interp(a1,:);
a1
kf=find(vals(tk25)>100);
cmm(3)=length(find(vals(tk25)>100)) %==52
pl1=subplot(2,3,6);
cc=find(isnan(vals)==0);
s1=[rb rlim];
nax=gca;
tv=linspace(s1(1),s1(2),100);
blp=brewermap(100,csch);
col=[];
cc=find(isnan(vals)==0);
for ct2=1:length(lonr)
    fg=find(vals(ct2)<=tv);
    if length(fg)>0
        col(ct2,1:3)=blp(fg(1),:);
    else;
        col(ct2,1:3)=blp(end,:);
    end
end
yy=[5 24 24 5 5];
xx=[70 70 90 90 70];
patch(xx,yy,[0.3 0.3 0.3]);
hold on;
[pl1] = gridded_map_India(lonr,latr,0.25/2,col,pl1);   
hold on;
plot(lonr(tk25(kf)), latr(tk25(kf)),'rx','Markersize',5, ...
     'MarkerFaceColor',[1 1 1]);
text(67.5,30,'(e)','fontsize',18);
grid off;
text(81.5,28.5,[ num2str(cmm(3)) ' extremes'],'fontsize',18);
box on;
axis([72 88 16 30]);
set(gca,'fontsize',15,'linewidth',2);
posg=get(gca,'Position');
set(gca,'Position',[posg(1)*1.05 posg(2)*0.4 dpos(3:4)*1.2]);
df=get(nax,'Xtick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ E'];
end
set(gca,'XTickLabel',dfs);


df=get(nax,'Ytick');
dfs=[];
for ct2=1:length(df)
dfs{ct2}=[num2str(df(ct2)) '^\circ N'];
end
set(gca,'YTickLabel',dfs);

text(78,31.5,'\downarrow','fontsize',50,'fontweight','bold');

figpos = get(gcf, 'Position');
paperWidth = figpos(3);
paperHeight = figpos(4);
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print('-dpdf',['Fig_2']);

