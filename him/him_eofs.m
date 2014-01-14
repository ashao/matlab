function [ pcs, eofs, pervar, data_proc ] = him_eofs( time, data_proc, neofs )
%% Calculates eofs on HIM gridded data
% Applies the landmask, weights by area, and filters out annual harmonics
% Input:
%   time: timestamp in days
%   data: physical fields [time,lat,lon]
%   neofs: number of modes to compute (default: 20)
% Output:
%   pcs: Principal Components
%   eofs: Modes [time lat lon]
%   pervar: Normalized eigenvalues


% Set default number of modes
if nargin<3
    neofs=20;
end

% Load HIM metrics file
load metrics.mat

disp('Preprocessing Data')
% Make a landmask
dry=logical(~metrics.wet.data);

[ntime nlat nlon]=size(data_proc);
Ah=metrics.Ah.data;

% Reshape matrices so that they are just [ntime nlat*nlon]
Ah=reshape(Ah,[numel(Ah), 1]);
data_proc=reshape(data_proc,[ntime numel(Ah)]);

% Delete all data on land so as not to bias results
Ah(dry)=[];
data_proc(:,dry)=[];

% Weight the data by area
total_area=sum(sum(Ah));
Ahgrid=repmat(Ah',[ntime 1]);
data_proc=data_proc.*Ahgrid./total_area;
[ntime ncells]=size(data_proc);
meandata(1,:)=mean(data_proc);
data_proc=data_proc-repmat(meandata,[ntime 1]);
disp('Removing seasonal cycle')

% Remove annual cycle and harmonics
parfor datapt=1:ncells
    
%     data_proc(:,datapt)=data_proc(:,datapt)-mean(data_proc(:,datapt));    
    if mod(datapt,1000)==0
        disp(sprintf('Gridpoint %d/%d',datapt,ncells));
    end
    [null null null null data_proc(:,datapt)]=annual_harms( ... 
        data_proc(:,datapt),time,5,5,365);
    
    
end

% Perform SVD decomp to estimate EOFs
disp('Calculating PCs and EOFs')
[U S V]=svds(data_proc,neofs);
pcs=U*S;
pervar=diag(S).^2/sum(diag(S).^2);
eofs=zeros([neofs nlat nlon]);
Vtemp=zeros(nlat,nlon);

disp('Reshaping modes')
for mode=1:neofs
    
    Vtemp(~dry)=V(:,mode);
    eofs(mode,:,:)=reshape(Vtemp,[nlat,nlon]);
end

end