function[x,y,z]= mid2corners(x_mid,y_mid,z_mid,W,L,t)
    x=x_mid-(0.5*W);
    y=y_mid-(0.5*L);
    z=z_mid-(0.5*t);
end