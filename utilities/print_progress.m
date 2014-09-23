function print_progress(iter, totiter, outfreq)

if (mod(iter,outfreq)==0) || (iter == 1)
    
    fprintf('Iteration %d/%d\n',iter,totiter);
    
end