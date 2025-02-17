addpath('../Functions/');
addpath('../Data');

% pull up the 1 degree data
load ind_rain_all_1951_2015_rev.mat;
ax=find(lonm_cut>=74.5 & lonm_cut<=86.5);
ay=find(latm_cut>=16.5 & latm_cut<=26.5);
tk=intersect(ax,ay); %tk is then the gridboxes in 1x1 data
                     %that are in central India


rr=find(mm>5 & mm<10);
%get distribution of station count

%for each gridbox, I need to determine the number of times the
%station count falls below 4
for ct1=1:length(latm_cut)
    cc=find(ind_stn_cut(rr,ct1)<4);
    rbd(ct1)=length(cc);
end

load IMD_rf;
%relate each 1 deg gridbox to the 0.25 degree grid 
bd=0.5;
figure;
hold on;
stdv=struct();
for ct1=1:length(lonm_cut);
    xx=lonm_cut(ct1);
    xx=[xx-bd xx-bd xx-bd xx xx+bd xx+bd xx+bd xx xx-bd];
    yy=latm_cut(ct1);
    yy=[yy-bd yy yy+bd yy+bd yy+bd yy yy-bd yy-bd yy-bd];
    dd=inpolygon(lonr,latr,xx,yy);
    %plot(xx,yy);
    stdv(ct1).stx=find(dd==1); %indices of the IMD 0.25 grid that
                               %are within the corresponding 1 deg gridbox;
end

%find all the distances between each 1 deg gridpoint and all 0.25
%stations
st_dist_all=struct();
for ct1=1:length(lonm_cut)
    xx=lonm_cut(ct1);
    yy=latm_cut(ct1);
    zz=sqrt((xx-lonr).^2+(yy-latr).^2);
    st_dist_all(ct1).dist=zz;
end

for ct1=1:length(lonm_cut)
    tic,
    xx=lonm_cut(ct1);
    yy=latm_cut(ct1);
    dist=sqrt((lonm_cut-xx).^2+(latm_cut-yy).^2);
    xx=find(dist<=2);
    [AA BB]=sort(dist(xx));
    for ct2=1:length(rr) %go through each day to determine where
                         %these stations could be
        stc=ind_stn_cut(rr(ct2),xx(BB)); %station counts of each of
                                         %the locations xx(BB)
        
        %then you need to do a cumsum of these, find <=4;
        pz=cumsum(stc);
        pz=find(pz<=4);
        axc=find(stc>0);
        axc=intersect(pz,axc);
        if length(axc)==0;%only admit members of the cumsum<=4 that
                          %actually have stations, such that you
                          %end at a gridbox with stations         
        else
            pz=axc(end);
        end
        if length(pz)>0;
            if stc(pz(end))==0 %if the gridbox we are on has a station
                          %count of 0, and then cumsum indicates
                          %that the next gridbox has a sufficient
                          %number of stations, then set this index
                          %to the next one (so long as that next
                          %index exists
                if (pz(end)+1)<=length(stc)
                    pz=pz(end)+1;
                end
            else
                pz=pz(end);
            end
          
            iz=find(stc(1:pz)>0);
            stlo=xx(BB(iz)); %indices of station locations in the
                             %1x1 grid
            stc=ind_stn_cut(rr(ct2),xx(BB(iz)));
            if sum(stc)<5
                stnx=[];
                for ct3=1:length(stc)
                    if stc(ct3)<=length(stdv(stlo(ct3)).stx)
                        stnx=[stnx; randsample(stdv(stlo(ct3)).stx, ...
                                               stc(ct3))];
                    else
                        stnx=[stnx; stdv(stlo(ct3)).stx(:)];
                        
                    end
                end
            else
                if length(stdv(stlo).stx)>3
                    stnx=randsample(stdv(stlo).stx,4);
                else
                    stnx=randsample(stdv(stlo).stx,length(stdv(stlo).stx));
                end
            end
            station_network(ct1).tstep(ct2).stnx=stnx;
            station_network(ct1).tstep(ct2).dist= ...
                st_dist_all(ct1).dist(stnx);
        else
            ct3=ct1;
            if length(stdv(ct3).stx)>3
                if stc(1)>=4
                    stnx=randsample(stdv(ct3).stx,4);
                else
                    stnx=randsample(stdv(ct3).stx,stc(1));
                end
            else 
                stnx=stdv(ct3).stx;
            end
            station_network(ct1).tstep(ct2).stnx=stnx;
            station_network(ct1).tstep(ct2).dist= ...
                st_dist_all(ct1).dist(stnx);
        end
    
    end
    toc
end

tout=nan(length(lonm_cut),length(rr));
for ct1=1:length(lonm_cut)
    for ct2=1:length(station_network(ct1).tstep)
        ds=station_network(ct1).tstep(ct2).dist;
        if length(ds(ds>2))>0
            xx=length(ds(ds>2));
            tout(ct1,ct2)=1;
        end
    end
    gct(ct1)=length(find(tout(ct1,:)>0))/length(rr);
end 


tout=nan(length(lonm_cut),length(rr));
for ct1=1:length(lonm_cut)
    gct(ct1)=length(find(tout(ct1,:)>0))/length(rr);
end 

%each day, we can also compute the cos weight that goes into the
%interpolation
for ct1=1:length(lonm_cut)
    tic,
    for ct2=1:length(station_network(ct1).tstep)
        ds=station_network(ct1).tstep(ct2).dist;
        ds(ds==0)=0.01; %if we happened to have selected an exactly
                        %identical location, put it 1 km
                        %away... just so the inverse distance
                        %weighting doesn't just give infinity
        mx=2;
        if max(ds)>2
            mx=max(ds); %if the maximum distance is greater than 2
                        %deg, the interpolation radius, then set it
                        %to the max interpolation distance.
        end
        if length(ds)>1
            dws=ds; %these are the distances in degrees, and
                    %the next two lines of code modify the
                    %distances with slight changes to the
                    %weighting according to Rajeevan 2003... it's not that 
                    %different from just using the distance. 0.5
                    %degrees is a threshold 
            mdd=mx/3;
            dws(ds<=mdd)=1./ds(ds<=mdd);
            dws(ds>mdd)=(1/(4*mx))*(27*(ds(ds>mdd)/mx-1).^2);
            comt=zeros(size(ds));
            
            %compute the cosine between each of these points, so
            %that the weights of these stations are balanced
            t=nan(length(ds),length(ds));
            XX=lonm_cut(ct1);
            YY=latm_cut(ct1);
            xii=lonr(station_network(ct1).tstep(ct2).stnx);
            yii=latr(station_network(ct1).tstep(ct2).stnx);
            for ctr=1:length(ds)
                for cts=2:length(ds)
                    t(ctr,cts)=((XX-xii(ctr))*(XX-xii(cts)) + ...
                                (YY-yii(ctr))*(YY-yii(cts)))/(ds(ctr)*ds(cts));
                    t(cts,ctr)=t(ctr,cts);
                end
            end
            comt=zeros(size(ds));
            for ctr=1:length(comt)
                t(ctr,ctr)=1;
                ri=t(ctr,:);
                fs=dws.*(1-ri');
                if nansum(dws(fs>0))>0               
                    comt(ctr)=nansum(fs)/nansum(dws(fs>0));
                end
            end

            station_network(ct1).tstep(ct2).shep_weights=(dws.^2).*(1+comt);
        end
        
    end
    toc
end


%check that there are stations to interpolate from for every
%time-step
for ct1=1:length(lonm_cut)
    for ct2=1:length(station_network(ct1).tstep)
        if length(station_network(ct1).tstep(ct2).stnx)>0
            lvs(ct1,ct2)=1;
        end
    end
end

%rr are the indices for the dates that ar interpolated
%station_network
%save('rf1x1_stnnetwork','station_network','rr');

for ct1=1:length(ym)
    Cos_interp_1x1(ym(ct1))
end


%now analyze what the extremes look like in each of these... 99.6th
%percentile
for ctg=1:18
    load(['TRMM_interp_cos_1x1_' num2str(ym(ctg))]);
    temp=TRMM_rg(ctg).rf(:,tk25);
    vlim=prctile(temp(:),99.6);
    ry=1951:2015;
    for ct1=1:length(ry)
        pz=find(yy(rr)==ry(ct1));
        xx=rf1_interp(pz,tk);
        evc(ctg,ct1)=length(find(xx(:)>=vlim));
    end
end

%pull out the trends for each of these timeseries
X=[ones(size(ry(:))) ry(:)-mean(ry)];
for ct1=1:length(ym)
    yT=evc(ct1,:);
    sl=multilin_lsq(X,yT(:));
    ux(ct1)=sl(2)*10; %events per decade
end

for ct1=1:length(ry)
    pz=find(yy(rr)==ry(ct1));
    alts=ind_rf_cut(rr(pz),tk);
    accss(ct1)=length(find(alts>100));
end

sl=multilin_lsq(X,accss(:));


%find the aggregate count when each gridpoint has a different
%threshold
for ctg=1:length(ym)
    load(['TRMM_interp_cos_1x1_' num2str(ym(ctg))]);
    temp=TRMM_rg(ctg).rf;
    cts=[];
    for cty=1:length(ry)
        pz=find(yy(rr)==ry(cty));
        xx=rf1_interp(pz,tk);
        for ct1=1:length(tk)
            temp2=temp(:,stdv(tk(ct1)).stx);
            vlim(ct1)=prctile(temp2(:),99.6);
            cts(cty,ct1)=length(find(xx(:,ct1)>vlim(ct1)));
        end
    end
    evclim(1:length(ry),ctg)=sum(cts');
end


for ctg=1:length(ym)
    load(['TRMM_interp_cos_1x1_' num2str(ym(ctg))]);
    temp=TRMM_rg(ctg).rf;
    cts=[];
    for cty=1:length(ry)
        pz=find(yy(rr)==ry(cty));
        xx=rf1_interp(pz,tk);
        anmean(ctg,cty)=nanmean(xx(:));
    end
    anmean_anom(ctg,1:length(ry))=anmean(ctg,:)-nanmean(anmean(ctg,:));
end

for ctg=1:18
    evc_anom(ctg,1:length(ry))=evc(ctg,:)-nanmean(evc(ctg,:));
end

figure;
plot(ry, evc_anom,'linewidth',2);
set(gca,'fontsize',16,'linewidth',2);

%find the same location in the interpolated TRMM data
ax=find(lonr>=74.5 & lonr<=86.5);
ay=find(latr>=16.5 & latr<=26.5);
tk25=intersect(ax,ay); %tk is then the gridboxes in 1x1 data
                     %that are in central India


%plot average annual station count over Central India
for ct1=1:length(ry)
    bpp=find(ry(ct1)==yy(rr));
    tct(ct1)=sum(sum(ind_stn_cut(rr(bpp),tk)))/122;
end



%ok look at each season of TRMM and find the gridpoint counts
%>99.5%
for ctg=1:length(TRMM_rg)
    load(['TRMM_interp_cos_1x1_' num2str(ym(ctg))]);
    temp=rf1_interp;
    tempk=TRMM_rg1(ctg).rf(tk,:);
    vlim=prctile(tempk(:),99.6);
    for ct1=1:length(rr);
        dc(ctg,ct1)=length(find(rf1_interp(ct1,tk)>vlim));
        stct(ctg,ct1)=sum(sum(ind_stn_cut(rr(ct1),tk)));
    end
    
end
stct=stct(1,:);
anomd=(dc-repmat(nanmean(dc),18,1));

%find the fraction of each year that is an underestimate
ax=actct(:);
aa=find(ax==0);
for ctg=1:length(ry)
    xx=find(yy(rr)==ry(ctg));
    az=anomd(:,xx);
    dm=az(:);
    bb=find(dm>0);
    bb=intersect(aa,bb);
    overdays(ctg)=length(bb)/length(aa); %fraction of all true
                                         %0-count days that have
                                         %extreme events
    fracover(ctg)=length(find(az(:)>0))/length(az(:));
    stan(ctg)=mean(stct(xx));
end


%griddata onto the 1x1 grid using the nearest neighbor
%interpolation for all years of TRMM data
for ctg=1:length(ym)
    tic,
    temp=TRMM_rg(ctg).rf;
    for ct1=1:size(temp,1)
        Vq = griddata(lonr,latr,temp(ct1,:),lonm_cut,latm_cut,'nearest');
        TRMM_rg1(ctg).rf(1:length(lonm_cut),ct1)=Vq;
    end
    toc
end

for ctg=1:length(ym)    
    temp=TRMM_rg1(ctg).rf(tk,:);
    vlim=prctile(temp(:),99.6);
    for ct1=1:size(TRMM_rg1(ctg).rf,2)
        temp2=TRMM_rg1(ctg).rf(tk,ct1);
        actct(ctg,ct1)=length(find(temp2>vlim));
    end
end

actct_rep=repmat(actct,1,65);
danom=dc-actct_rep;

%find the fraction of each year that is an underestimate
for ctg=1:length(ry)
    xx=find(yy(rr)==ry(ctg));
    az=danom(:,xx);
    fracover(ctg)=length(find(az(:)>0))/length(az(:));
    fracunder(ctg)=length(find(az(:)<0))/length(az(:));
    omitted(ctg)=sum(az(find(az(:)<0)))/18;
    stan(ctg)=mean(stct(xx));
end



gg=figure;
gg.Units='inches';
gg.Position=[2 2 15 20]; 
letn={'(a)' '(b)' '(c)' '(d)' '(e)' '(f)' '(g)' '(h)' '(i)' '(j)' ...
      '(k)' '(l)' '(m)' '(n)' '(o)' '(p)' '(q)' '(r)' '(s)' '(t)' ...
      '(u)' '(v)' '(w)'};
for ct1=1:length(ym)
    subplot(6,3,ct1);
    plot(ry,evc(ct1,:),'linewidth',2);
    axis tight;
    set(gca,'fontsize',16,'linewidth',2);
    %title([num2str(ym(ct1))],'fontsize',16);
    ylabel('# extremes','fontsize',14);
    xlabel('Year','fontsize',14);
    uy=get(gca,'ylim');
    text(1935,uy(2),letn{ct1},'fontsize',18);
end
figpos = get(gcf, 'Position');
paperWidth = figpos(3);
paperHeight = figpos(4);
set(gcf, 'paperunits', 'inches');
set(gcf, 'papersize', [paperWidth paperHeight]);
set(gcf, 'PaperPosition', [0    0   paperWidth paperHeight]);
print('-dpng',['All_TRMM_timeseries_1x1']);

%fraction of all extreme events omitted, if there were stations at
%the center of every 1 degree gridbox
mean(omitted)/(122*length(tk)*0.004);


for ct1=1:length(ry)
    yz=find(ry(ct1)==yy(rr));
    rrc=ind_rf_cut(rr(yz),tk);
    imd_cts(ct1)=length(find(rrc(:)>100));
end



%check rf1_interp
ry=1951:2015;
for ctg=1:length(ym)
    load(['TRMM_interp_cos_1x1_' num2str(ym(ctg))]);
    for ctk=1:length(ry)
        cs=find(ry(ctk)==yy(rr));
        anmn(ctg,ctk,1:length(tk))=nanmean(rf1_interp(cs,tk));
    end
    %TRMMn(ctg).rf=rf1_interp;
end

llx=squeeze(nanmean(anmn,3));
X=[ones(size(ry(:))) ry(:)-mean(ry)];

for ctm=1:length(ym)
    sl=multilin_lsq(X,llx(ctm,:)');
    mdec(ctm)=sl(2)*10; %mm per decade
end

acrf=nan(size(ry));
for ct1=1:length(ry)
    jz=find(yy(rr)==ry(ct1));
    acrf(ct1)=mean(mean(ind_rf_cut(rr(jz),tk)'));
end

slc=multilin_lsq(X,acrf(:));

%-0.0054
%-0.0544 per decade
%gotta scale this for the mean diff with TRMM
%TRMM_rg1 is the nearest neighbor interp of TRMM_rg 
%at mean(mean(ind_rf_cut(rr,tk))) is the mean for the IMD data
%7.1850 mm per day
nmt=nan(size(ym));
for ct1=1:length(ym)
    nmt(ct1)=nanmean(nanmean((TRMM_rg1(ct1).rf(tk,:))))
end

%7.185/mean(nmt) = 1.4146
%scaling this trend by something similar 

%this is the scaling of trends to the IMD data
1.4146*mdec./(slc(2)*10)