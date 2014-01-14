% function [ fronts ] = orsi_fronts( );
inpath = 'C:\Users\ashao\data\orsi_fronts\';
data = csvread([inpath filesep 'saf.csv']);
lon = data(:,1)';
lat = data(:,2)';
flipidx = lon<0;
fronts.saf.lon = [lon(~flipidx) lon(flipidx)+360];
fronts.saf.lat = [lat(~flipidx) lat(flipidx)];