function [ cfc11, cfc12 ] = off_piecewise_atmospheric ( cfc11, cfc12 )
    
    cfc11.Nval=interp_to_1970(cfc11.year,cfc11.Nval);
    cfc11.Sval=interp_to_1970(cfc11.year,cfc11.Sval);
    cfc12.Nval=interp_to_1970(cfc12.year,cfc12.Nval);
    cfc12.Sval=interp_to_1970(cfc12.year,cfc12.Sval);

end

function atmval = interp_to_1970(year,atmval)
    sidx=min(find(atmval>0));
    eidx=min(find(year>1970));
    ridx=sidx:eidx;
    
    idx_endpoints=[sidx eidx];
    endpts_year=year(idx_endpoints);
    endpts_val=atmval(idx_endpoints);
    atmval(ridx)=interp1(endpts_year,endpts_val,year(ridx));
end