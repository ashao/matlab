clear

mphealth = [150 225 326 457 639 895 1253 1755 2457 3439];
basehealth = mphealth(1);
ntrials = 100;

for plevel = 1:length(mphealth)
        
    probdrop = 0.1*plevel;
    runtime = zeros(ntrials,1);
    counter = zeros(ntrials,1);
    for trial = 1:ntrials;
        drop = false;
        while ~drop            
            counter(trial) = counter(trial) +1;
%             disp(counter(trial))
            runtime(trial) = runtime(trial)+mphealth(plevel)/basehealth;
            if rand(1) < probdrop
                drop = true;
            end
        end
        
    end
    
    time(plevel)=mean(runtime);
    totruns(plevel)=mean(counter);
    
    
end

[null optlevel]=min(time);
fprintf('Optimal MP Level: %d Number of runs: %d\n',optlevel,totruns(optlevel))