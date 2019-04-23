function [efflux, sitesC, sitesCg, sitesC2] = EquilibriumConcentrationFlux2(iniConcentrationGas, waterVolume,gasVolume,sitesC,sitesC2,sitesCg,HenryDomain,InvasedIsland,AlphaM,sitesNEquil,changeN)

efflux = zeros(1,5);
ListGas = [1 2 6 7 8];
%oxygen, CO2, NH3
for iG =1:5
    iN = ListGas(iG);
    if iN == 7 %HONO separately due to the index problem
        
        totalN = sitesC2{6}.*waterVolume + sitesCg{iN}.*gasVolume;
        sitesCtemp = sitesNEquil{iN}.*InvasedIsland + (InvasedIsland==0).*HenryDomain{iN}.*totalN./AlphaM{iN};
        sitesCgtemp = iniConcentrationGas(iN).*InvasedIsland + (InvasedIsland==0).*totalN./AlphaM{iN};
        newTotN = sitesCtemp.*waterVolume + sitesCgtemp.*gasVolume;
        deltaN = -1*(newTotN - totalN).*InvasedIsland;
        efflux(iG) = sum(deltaN(:)) + sum(sum(changeN(:,:,iN)));%./HenryDomain{iN}));
        %deltaN = (totalN - AlphaM{iN}.*iniConcentrationGas(iN)).*InvasedIsland./HenryDomain{iN};
        %efflux(iG) = sum(deltaN(:)) + sum(sum(changeN(:,:,iN)./HenryDomain{iN}));
        sitesC2{6} = sitesCtemp;
        sitesCg{iN,1} = sitesCgtemp;
        
    else if iN == 6 % ammonia
            
            totalN = sitesC2{3}.*waterVolume + sitesCg{iN}.*gasVolume;
            sitesCtemp = sitesNEquil{iN}.*InvasedIsland + (InvasedIsland==0).*HenryDomain{iN}.*totalN./AlphaM{iN};
            sitesCgtemp = iniConcentrationGas(iN).*InvasedIsland + (InvasedIsland==0).*totalN./AlphaM{iN};
            newTotN = sitesCtemp.*waterVolume + sitesCgtemp.*gasVolume;
            deltaN = -1*(newTotN - totalN).*InvasedIsland;
            efflux(iG) = sum(deltaN(:)) + sum(sum(changeN(:,:,iN)));%./HenryDomain{iN}));
            %deltaN = (totalN - AlphaM{iN}.*iniConcentrationGas(iN)).*InvasedIsland./HenryDomain{iN};
            %efflux(iG) = sum(deltaN(:)) + sum(sum(changeN(:,:,iN)./HenryDomain{iN}));
            sitesC2{3} = sitesCtemp;
            sitesCg{iN,1} = sitesCgtemp;
        else
            totalN = sitesC{iN}.*waterVolume + sitesCg{iN}.*gasVolume;
            sitesCtemp = sitesNEquil{iN}.*InvasedIsland + (InvasedIsland==0).*HenryDomain{iN}.*totalN./AlphaM{iN};
            sitesCgtemp = iniConcentrationGas(iN).*InvasedIsland + (InvasedIsland==0).*totalN./AlphaM{iN};
            newTotN = sitesCtemp.*waterVolume + sitesCgtemp.*gasVolume;
            deltaN = -1*(newTotN - totalN).*InvasedIsland;
            efflux(iG) = sum(deltaN(:)) + sum(sum(changeN(:,:,iN)));%./HenryDomain{iN}));
            %deltaN = (totalN - AlphaM{iN}.*iniConcentrationGas(iN)).*InvasedIsland./HenryDomain{iN};
            %efflux(iG) = sum(deltaN(:)) + sum(sum(changeN(:,:,iN)./HenryDomain{iN}));
            sitesC{iN,1} = sitesCtemp;
            sitesCg{iN,1} = sitesCgtemp;
        end
    end
end

end
