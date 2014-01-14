function [ ] = off_lintracval(cfc11,cfc12,outpath)
%     Performs a least squares fit to a polynomial of degree 1 for the
%     given tracer

%%CFC-11    
cfc11=calc_linvals( cfc11 );
write_trac_tsv(cfc11,[outpath filesep 'cfc11_lin.tsv']);

%%CFC-12
cfc11=calc_linvals( cfc12 );
write_trac_tsv(cfc11,[outpath filesep 'cfc12_lin.tsv']);

end

function array = calc_linvals( array )
% Get the linear least squares fit of the atmospheric increase from
% 1970-1990
    idx=find(array.year>1970 & array.year<1990);
    Nlinfit=polyfit(array.year(idx),array.Nval(idx),1);
    Slinfit=polyfit(array.year(idx),array.Sval(idx),1);
% Replace the values after 1990
    replaceidx=find(array.year>1990);
    Nlinvals=polyval(Nlinfit,array.year(replaceidx));
    Slinvals=polyval(Slinfit,array.year(replaceidx));
    array.Nval(replaceidx)=Nlinvals;
    array.Sval(replaceidx)=Slinvals;
end