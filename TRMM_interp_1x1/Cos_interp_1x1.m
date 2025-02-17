function [] = Cos_interp_1x1(yy1);
%function [] = Cos_interp_1x1(yy1);
%yy1 is the year of TRMM data to be interpolated for all years of
%IMD station networks. This function saves interpolated data as
%TRMM_interp_cos_XXXX, where XXXX is the TRMM year.
    addpath('../Data');
    addpath('../Functions');
    load TRMM_mic;
    ym=unique(TRMM_mic(1).yy); %individual years of available TRMM data
    load IMD_rf; %all of the IMD rainfall data, rf25 are the
                 %gridded rainfall values
    load ind_rain_all_1951_2015_rev.mat;
    load rf1x1_stnnetwork;
    dd1_g=dd;
    mm1_g=mm;
    yy1_g=yy;
    %first regrid TRMM onto the IMD grid using linear
    %interpolation, integrates 3 AM UTC to 3 AM UTC, assuming 12 AM 
    TRMM_rg=struct();
    nstrl=[];
    for ct1=1:length(TRMM_mic(1).yy)
        nstrl(ct1)=datenum(TRMM_mic(1).yy(ct1),TRMM_mic(1).mm(ct1), ...
                           TRMM_mic(1).dd(ct1),TRMM_mic(1).hh(ct1),0,0);
    end
    
    [AA BB]=sort(nstrl);
    
    %TRMM_mic(1).yy(BB) is now ordered exactly
    
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
            ns=nansum([ns 1.5*TRMM_mic(1).vals(:,BB(ex-7))+1.5*TRMM_mic(1).vals(:,BB(ex+1))],2);
            rf25_TRMM(ct1,1:length(lon_TRMM))=ns;
        end


        rf25_TRMM_rg=nan(length(yk),length(lonr));
        for ct1=1:length(yk)
            temp=rf25_TRMM(ct1,:);
            kl=find(isnan(temp)==0);
            rf25_TRMM_rg(ct1,1:length(lonr))= griddata(lon_TRMM(kl), ...
                                                       lat_TRMM(kl),double(temp(kl)),lonr,latr);
        end
        TRMM_rg(ctm).rf=rf25_TRMM_rg;
        ctm
    end
    
    
    
   
    l1=rr;

    [ytm mtm dtm]=datevec(datenum(2000,6,1):datenum(2000,9,30));
    %for lg=1:length(ym)
    lg=find(ym==yy1);
        %cth=find(ym==yy1); 
        cth=ym(lg);
        rtem=TRMM_rg(lg).rf;
        rf1_interp=nan(length(l1),length(lonm_cut));
        for ctp=1:length(lonm_cut)
            for ctn=1:length(station_network(ctp).tstep)
                if length(station_network(ctp).tstep(ctn).dist)>0
                    a2=l1(ctn);
                    k1=find(mm1_g(a2)==mtm);
                    k2=find(dd1_g(a2)==dtm);
                    k1=intersect(k1,k2);
                    rft=rtem(k1,:)';
                    iv=station_network(ctp).tstep(ctn).stnx;
                    if length(station_network(ctp).tstep(ctn).stnx)>1
                        dws=station_network(ctp).tstep(ctn).shep_weights; 
                        rf1_interp(ctn,ctp)=sum(rft(iv).*(dws))/ ...
                            sum(dws);
                    else
                        rf1_interp(ctn,ctp)=rft(iv);
                    end
                end
            end
            ctp
        end

        save(['TRMM_interp_cos_1x1_' num2str(ym(lg))],'rf1_interp',['-' ...
                            'v7.3']);
        %end
end

