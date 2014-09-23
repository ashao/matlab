parfor tracknum=1:254

	datapath = '/ltraid4/aviso/alongtrack/sla/vxxc/'
	outpath = '/ltraid4/aviso/matlab/vxxc/';
	track = extract_monomission_sla(datapath,outpath,tracknum);
end
