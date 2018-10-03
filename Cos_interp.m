function [] = Cos_interp(yy1);
%function [] = Cos_interp(yy1);
%yy1 is the year of TRMM data to be interpolated for all years of
%IMD station networks. This function saves interpolated data as
%TRMM_interp_cos_XXXX, where XXXX is the TRMM year.
    addpath('./Data');
    addpath('./Functions');
    load TRMM_mic;
    ym=unique(TRMM_mic(1).yy); %individual years of available TRMM data
    load IMD_rf; %all of the IMD rainfall data, rf25 are the
                 %gridded rainfall values
    

    
    %search the radius for every cell first
    ix=struct();
    for ct1=1:length(lonr)
        XX=lonr(ct1);
        YY=latr(ct1);
        rx=sqrt((lonr-XX).^2+(latr-YY).^2);
        xl=find(rx<=1.5);
        dist=rx(xl);
        %rank the 'stations' from nearest to farthest
        [AA BB]=sort(dist);
        ix(ct1).dist=AA;
        ix(ct1).ix=xl(BB); %get the direction factor out here
        
        xii=lonr(ix(ct1).ix);
        yii=latr(ix(ct1).ix);
        for ctr=1:length(ix(ct1).dist)
            %this computes the cosine between all two possible
            %distance vectors into matrix t
            for cts=2:length(ix(ct1).dist)
                t(ctr,cts)=((XX-xii(ctr))*(XX-xii(cts)) + ...
                            (YY-yii(ctr))*(YY-yii(cts)))/(ix(ct1).dist(ctr)*ix(ct1).dist(cts));
                t(cts,ctr)=t(ctr,cts);
            end
        end
        ix(ct1).cos_t=t;
    end

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
    
    
    
    
    l1=find(m25>5 & m25<10);
    l2=find(y25>1950);
    l1=intersect(l1,l2);

    [ytm mtm dtm]=datevec(datenum(2000,6,1):datenum(2000,9,30));
    
    cth=find(ym==yy1); 

    rtem=TRMM_rg(cth).rf;
    rf25_interp=nan(length(l1),length(lonr));
    for ctn=1:length(l1)
        a2=l1(ctn);
        k1=find(m25(a2)==mtm);
        k2=find(d25(a2)==dtm);
        k1=intersect(k1,k2);
        rft=rtem(k1,:)';
        stt=stn25(:,a2); %station count in each gridbox
        rfn=nan(size(rft));
        rfn(stt>0)=rft(stt>0);
        qp=find(isnan(rfn)==1); %these are where there are no
                                %stations, and the locations that need
                                %to be 'interpolated' based on the
                                %gridpoints that do have stations
        for ct1=1:length(qp)    
            iv=ix(qp(ct1)).ix;
            kz=find(stt(iv)>0);
            if length(kz)>0
                stk=stt(iv(kz));
                stx=cumsum(stk); %this takes the cumulative sum of the
                                 %station counts of gridpoints that
                                 %have already been sorted by distance
                                 %from the gridpoint
                rb=find(stx<=4);
                if length(rb)==0
                    rb=1; %just take the first index, as this indicates
                          %that the nearest 'station' is more than four stations
                end
                sw=kz(1:rb(end));
                iv=iv(sw);
                ds=ix(qp(ct1)).dist(kz(1:rb(end)));
                if length(sw)==1
                    rfn(qp(ct1))=rfn(iv); %fixed on June 19th
                elseif length(sw)>1
                    dws=ds; %these are the distances in degrees, and
                            %the next two lines of code modify the
                            %distances with slight changes to the
                            %weighting according to Rajeevan 2003... it's not that 
                            %different from just using the distance. 0.5 degrees is a threshold 
                    dws(ds<=0.5)=1./ds(ds<=0.5);
                    dws(ds>0.5)=(1/(4*1.5))*(27*(ds(ds>0.5)/1.5-1).^2);
                    comt=zeros(size(ds));
                    cos_t=ix(qp(ct1)).cos_t;
                    for ctr=1:length(comt)
                        ri=cos_t(sw(ctr),sw);
                        fs=dws.*(1-ri');
                        if nansum(dws(fs>0))>0               
                            comt(ctr)=sum(fs)/nansum(dws(fs>0));
                        end
                    end
                    
                    dws=(dws.^2).*(1+comt);
                    rfn(qp(ct1))=sum(stk(rb).*rfn(iv).*(dws))/sum(stk(rb).*dws);                        
                end
                
            end
        end
        ctn, cth
        rf25_interp(ctn,1:length(lonr))=rfn;
    end

    save(['TRMM_interp_cos_' num2str(yy1)],'rf25_interp','-v7.3');
end

