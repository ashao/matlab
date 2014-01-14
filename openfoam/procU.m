function [ data ] = procU(filename)
array=load(filename);
data.x=array(:,1);
data.y=array(:,2);
data.z=array(:,3);
data.u=array(:,4);
data.v=array(:,5);
data.w=array(:,6);
data.vx=array(:,7);
data.vy=array(:,8);
data.vz=array(:,9);
end