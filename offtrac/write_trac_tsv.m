function [ ] = write_trac_tsv( array, outfile )
    
    dlmwrite(outfile,length(array.year),'delimiter','\t');
    dlmwrite(outfile,array.sch_coeffs,'delimiter','\t','-append')
    dlmwrite(outfile,array.sol_coeffs,'delimiter','\t','-append')
    data = [array.year array.Nval array.Sval];
    dlmwrite(outfile,data,'delimiter','\t','-append')
    
end