load him_wet.mat
count=0;
nwet=length(find(wet));

fill2d=zeros(nwet,2);
p_ptemp=fill2d;
p_mldepth=fill2d;
p_ssh=fill2d;

for i=1:210
    disp(sprintf('Latitude %d/%d',i,210))
    for j=1:360
%     randy=randi(210,1);
%     randx=randi(360,1);
    if wet(i,j)
       count=count+1;
        p_ptemp(count,:)=polyfit(array.year,squeeze(array.ptemp(:,2,i,j)'),1);
        p_mldepth(count,:)=polyfit(array.year,squeeze(array.mldepth(:,i,j))',1);
        p_ssh(count,:)=polyfit(array.year,squeeze(array.ssh(:,i,j))',1);
    end
    end
end