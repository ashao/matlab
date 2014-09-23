inpath = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/';
files = dir([inpath filesep 'ocean_month.*.nc']);
nfiles = length(files);
outpath = '/ltraid3/ashao/uw-apl/models/HIM/hindcast/';
p = 2000*ones([12 49 210 360]);
pr = zeros([12 49 210 360]);
years = 1948:2007;

load metrics
wet = logical(metrics.wet.data);

for i=1:nfiles
    
    infile = [inpath files(i).name];
    fprintf(infile);
    tic;
    temp = nc_varget(infile,'temp');
    salt = nc_varget(infile,'salt');
    pden = nan([12 49 210 360]);
    pden(:,:,wet) = single(sw_pden(salt(:,:,wet),temp(:,:,wet), ...
        p(:,:,wet),pr(:,:,wet)));
    outfile = [outpath sprintf('sigmatheta.%d.mat',years(i))];
    save(outfile,'pden');
    fprintf(': %f secs\n',toc);
    
end