function [ out ] = gasex( datapath, startyear, startmon, endmon )
    
    ngasex=30;
    nmon=length(startmon:endmon);
    updateint=3;
    
    [T, S, H, P, FICE, XKW, TRACERS, LAT, LON] = loaddata(datapath);
    
    nlat=length(LAT);
    nlon=length(LON);
    ntrac=length(TRACERS);
    
    mldepth=squeeze(sum(H,2));

%     out.time=zeros(nmon,1);
    time=startyear;
    fill1D=single(zeros(nmon,1));
    fill3D=single(zeros(nmon,nlat,nlon));
    
    for i=1:length(TRACERS)
        out(i).name=TRACERS(i).name;
        out(i).lat=LAT;
        out(i).lon=LON;
        out(i).mldepth=mldepth;
        out(i).time=fill1D;
        out(i).conc=fill3D;
        out(i).satval=fill3D;
        out(i).relsat=fill3D;
    end
    
    tmpconc=zeros(nlat,nlon);
    
    for mon=startmon:endmon        
        midx=mod(mon,12)+1; % Index advanced by 1?
        if midx==0
            midx=12;
        end
        dt=eomday(1,midx)*86400;
        dt_gasex=dt/ngasex;
        timedt=dt_gasex/31536000;
%         disp(dt_gasex)
%         pause
        surfT=squeeze(T(midx,1,:,:));
        surfS=squeeze(S(midx,1,:,:));
        surfMLH=squeeze(mldepth(midx,:,:));
        surfP=squeeze(P(midx,:,:));
        surfXKW=squeeze(XKW(midx,:,:));
        surfFICE=squeeze(FICE(midx,:,:));
%         disp([surfT(100,100) surfS(100,100)])
        if mod(mon,updateint)==0
%             disp(sprintf('\nTime: %f',time))
%             disp(sprintf('Tracer\t\tRelSat\t\tConcMixedLayer'));
        end
        for trac=2:2
            Sc=calc_schmidt( TRACERS(trac).sch_coeffs,surfT);
            sol=calc_sol( TRACERS(trac).sol_coeffs,surfS,surfT);
            Sk=(1-surfFICE).*surfXKW.*(Sc./660).^(-.5)./surfMLH*dt_gasex;
            if mon>1
                tmpconc=squeeze(out(trac).conc(mon-1,:,:));
            else
                tmpconc=zeros(nlat,nlon);
            end
            for t=1:ngasex
%                 disp(sprintf('%g',time+timedt*t))
                atmcon=calc_atmcon(TRACERS(trac),time+timedt*t,LAT);
                atmcon=reshape(atmcon,length(LAT),1);
                atmcon=repmat(atmcon,1,length(LON));
                satval=atmcon.*surfP.*sol;
                flux=Sk.*(satval-tmpconc);
                tmpconc=tmpconc+flux;
                
            end
            out(trac).satval(mon,:,:)=single(satval);
            out(trac).conc(mon,:,:)=single(tmpconc);
            out(trac).relsat(mon,:,:)=single(tmpconc./satval);
            out(trac).time(mon)=single(time+timedt*ngasex);
            
            disp(sprintf('Lat/Lon %f %f',LAT(100),LON(100)));
            disp(sprintf('\nTime: %f',time+timedt*ngasex))
            disp(sprintf('Atmospheric Concentration=%6.4f',atmcon(100,100)))
            disp(sprintf('T=%g S=%g P=%g',surfT(100,100),surfS(100,100),surfP(100,100)))
            disp(sprintf('Solubility=%g',sol(100,100)))
            disp(sprintf('Mixed Layer Concentration=%g',satval(100,100)))
            
            if mod(mon,updateint)==0
%                 disp(sprintf('%s\t\t%g\t%g',out(trac).name,out(trac).relsat(mon,100,100),...
%                     out(trac).conc(mon,100,100)));
            end
            
            %             disp(Sk(100,100));
        end
        time=time+(timedt*ngasex);
    pause
    end
    
    
end
% BEGIN SUBROUTINES
function [T, S, H, P, FICE, XKW, TRACERS,LAT,LON] = loaddata(datapath)
    start4D=[0 0 0 0];
    count4D=[12 4 -1 -1];
    start3D=[0 0 0];
    count3D=[12 -1 -1];
    
    T=nc_varget([datapath filesep 'ts.nc'],'TEMPCLIM',start4D,count4D);
    S=nc_varget([datapath filesep 'ts.nc'],'SALTCLIM',start4D,count4D);
    H=nc_varget([datapath filesep 'H-clim.nc'],'HCLIM',start4D,count4D);
    LAT=nc_varget([datapath filesep 'H-clim.nc'],'yh');
    LON=nc_varget([datapath filesep 'H-clim.nc'],'xh');
    P=nc_varget([datapath filesep 'gasx_ocmip2_himgrid.nc'],'OCMIP_ATMP', ...
        start3D,count3D);
    FICE=nc_varget([datapath filesep 'gasx_ocmip2_himgrid.nc'],'OCMIP_FICE',...
        start3D,count3D);
    XKW=nc_varget([datapath filesep 'gasx_ocmip2_himgrid.nc'],'OCMIP_XKW',...
        start3D,count3D);
    load([datapath 'cfc11.mat']);
    TRACERS(1)=cfc11;
    load([datapath 'cfc12.mat']);
    TRACERS(2)=cfc12;
    load([datapath 'sf6.mat']);
    TRACERS(3)=sf6;
end


function [ Sc ] = calc_schmidt( C,T )
    
    Sc=C(1)-C(2)*T+C(3)*T.^2-C(4)*T.^3;
    
end

function [ trac_sol ] = calc_sol(sol_coeffs,S,T)
    TempK=T+273.15;
    trac_sol = sol_coeffs(1) + sol_coeffs(2)*(100 ./ TempK)+ ...
    sol_coeffs(3) * log(TempK ./ 100) + sol_coeffs(4) * (TempK ./ 100).^2 ...
    + S .* (sol_coeffs(5) + sol_coeffs(6) * (TempK ./ 100) ...
    + sol_coeffs(7) * (TempK ./ 100).^2);
trac_sol = exp(trac_sol);
end

function [ atmcon ] = calc_atmcon( tracer, time, lat )
    
    nlat=10.2;
    slat=-10.2;
    nidx=find(lat>nlat);
    sidx=find(lat<slat);
    eidx=find(lat<=nlat & lat>=slat);
    
    if time<tracer.year(length(tracer.year))
        
        nval=interp1(tracer.year,tracer.Nval,time);
        sval=interp1(tracer.year,tracer.Sval,time);
        
        atmcon=zeros(size(lat));
    else
        nval=tracer.Nval(length(tracer.year));
        sval=tracer.Sval(length(tracer.year));
    end
    atmcon(nidx)=nval;
    atmcon(sidx)=sval;
    atmcon(eidx)=interp1([nlat slat],[nval sval],lat(eidx));
    
end

