function array = read_psd(url)

tempfile = 'temp.psd';
urlwrite(url,tempfile);

delete(tempfile)

end