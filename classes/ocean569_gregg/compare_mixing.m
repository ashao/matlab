% Compares the output from GOLD runs where the mixing parameterizations
% were altered
%% Load model output
prop.path = '/ltraid2/ashao/uw-apl/models/gold/mixing_exp/';
prop.expnames = {'control','noGM','noMLRESTRAT','weakmix'};
prop.nexp = length(prop.expnames);
prop.fields = {'temp','uh','vh','uhgm','vhgm','uhml','vhml','h','salt'};
prop.nfields = length(prop.fields);

% Extract model output
for ifield = 9:prop.nfields
    for iexp=1:prop.nexp
        fprintf('Experiment %s Field %s\n', ...
            prop.expnames{iexp},prop.fields{ifield})
        
        % Ecception  here because disabling MLRESTRAT makes GOLD output nan for
        % uhML, vhML
        if sum(ifield==6) + sum(ifield==7) > 0 & .... 
                (strcmp(prop.expnames{iexp},'noMLRESTRAT'))
            
            output.(prop.expnames{iexp}).(prop.fields{ifield})= ...
                zeros(size(output.control.temp));
        else
            filename = [prop.path 'ocean_year.' prop.expnames{iexp} '.nc'];
            output.(prop.expnames{iexp}).(prop.fields{ifield})= ...
                nc_varget(filename,prop.fields{ifield});
        end
        
    end
end
%% Calculate total transport (UH + UHGM + UHML)
[ntime nz nlat nlon]=size(output.control.temp);
for iexp = 1:prop.nexp
    if iexp==2 % Exception to deal with noGM case       
        output.(prop.expnames{iexp}).uhgm=zeros([ntime,nz,nlat,nlon]);
        output.(prop.expnames{iexp}).vhgm=zeros([ntime,nz,nlat,nlon]);
    end
    
    uhtot.(prop.expnames{iexp})=...
        output.(prop.expnames{iexp}).uh + ...
        output.(prop.expnames{iexp}).uhgm + ...
        output.(prop.expnames{iexp}).uhml;
    vhtot.(prop.expnames{iexp})= ...
        output.(prop.expnames{iexp}).vh + ...
        output.(prop.expnames{iexp}).vhgm + ...
        output.(prop.expnames{iexp}).vhml;
end

%% Calculate anomalies of transport in drake passage

Integrate in depth
for iexp=1:prop.nexp
    temp_array=uhtot.(prop.expnames{iexp});
    drakepassage.(prop.expnames{iexp}).uhtot_depth= ...
        sum(temp_array(:,:,mask.drakepassage.idx),2);
    temp_array=vhtot.(prop.expnames{iexp});
    drakepassage.(prop.expnames{iexp}).vhtot_depth= ...
        sum(temp_array(:,:,mask.drakepassage.idx),2);
end

mask.drakepassage.spatial2d = metrics.geolon.data == -66.5 & ...
    metrics.geolat.data <= -55.5;
mask.drakepassage.idx = find(mask.drakepassage.spatial2d);
colors = {'b','k','g'};

%% Plot yearly anomaly from control run
figure; hold on;
plot((drakepassage.control.uhtot_depth-drakepassage.noGM.uhtot_depth)./1e6, ...
    colors{1})
plot((drakepassage.control.uhtot_depth-drakepassage.noMLRESTRAT.uhtot_depth)./1e6, ...
    colors{2})
plot((drakepassage.control.uhtot_depth-drakepassage.weakmix.uhtot_depth)./1e6, ...
    colors{3})
grid on; axis square;
xlabel('Model Year')
ylabel('(Sv)')
box on
title('Drake Passage Transport Anomaly (120Sv)')
legend('noGM','noMLRESTRAT','WeakVMix')
% print -depsc -painters drakepassagetransport.eps


%% Mixed layer depth
m_proj('Mercator','lon',[-240 -160],'lat',[20 50])
for iexp=1:4
    subplot(2,2,iexp)
    m_pcolor(metrics.geolon.data,metrics.geolat.data, ... 
        mean(sum(output.(prop.expnames{iexp}).h(9:10,1:2,:,:),2)));
    shading interp;    
    m_grid;
    caxis([0 100])
    m_coast('patch',[0 0 0]);
    colorbar SouthOutside
    colormap(othercolor('BuDRd_12'))
    title(prop.expnames{iexp})
end
% print -depsc -painters MixedLayerDepth.eps
%% Meridional Heat Transpor
plotyear = 10;
for iexp = 1:prop.nexp    
    disp(iexp)
    temp_array=vhtot.(prop.expnames{iexp});
    cp = sw_cp(output.(prop.expnames{iexp}).salt(plotyear,:,:,:), ...
        output.(prop.expnames{iexp}).temp(plotyear,:,:,:), ...
        2000*ones(1,63,nlat,nlon));
    cp=reshape(cp,[63,210,360]);
    cp=cp.*squeeze(temp_array(plotyear,:,:,:)).* ...
        squeeze(output.(prop.expnames{iexp}).temp(plotyear,:,:,:)+273.15);
    cp= sum(...
        cp.*reshape(sw_pden(output.(prop.expnames{iexp}).salt(plotyear,:,:,:), ...
        output.(prop.expnames{iexp}).temp(plotyear,:,:,:), ...
        2000*ones(1,63,nlat,nlon),0),[63 210 360]),1);    
    mht.(prop.expnames{iexp})=sum(squeeze(cp),2);
end

%% Plot MHT
clf
hold on

plot(metrics.lath.data,(mht.noGM-mht.control)./1e12,'b')
plot(metrics.lath.data,(mht.noMLRESTRAT-mht.control)./1e12,'g')
plot(metrics.lath.data,(mht.weakmix-mht.control)./1e12,'k')
ylim([-100 200])
xlim([-90 90])
xlabel('Latitude')
ylabel('Heat Flux Anomaly (TW)')
grid on
axis square
legend('noGM','noMLRESTRAT','WeakMix')
box on
% saveas(gcf,'MHT.eps','epsc')