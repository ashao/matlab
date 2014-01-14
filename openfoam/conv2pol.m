function [ vt vr vz ] = conv2pol( ux, uy, uz, r )

    vr=sqrt(ux.^2+uy.^2);
    vt=atan(uy./ux).*abs(r);
    vz=uz;
    
end