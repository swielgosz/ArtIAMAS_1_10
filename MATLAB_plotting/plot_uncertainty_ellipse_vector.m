function plot_vecs = plot_uncertainty_ellipse_vector(FIM, target)
    Cov_Mat = inv(FIM);
    sx2 = Cov_Mat(1,1); sy2 = Cov_Mat(2,2); sxy = Cov_Mat(1,2);
    sx = sx2^(1/2); sy = sy2^(1/2);
    
    a2 = (sx^2+sy^2)/2 + (((sx^2-sy^2)^2 )/ 4 + sxy^2)^(1/2);
    b2 = (sx^2+sy^2)/2 - (((sx^2-sy^2)^2 )/ 4 + sxy^2)^(1/2);
    
    a = a2^(1/2); b = b2^(1/2);
    theta = (1/2)*atan2(2*sxy, sx^2-sy^2);
    u = target(1); v = target(2);
    
    t = 0:0.01:2*pi;
    x_plot = u+a*cos(t)*cos(theta) - b*sin(theta)*sin(t);
    y_plot = v+a*cos(t)*sin(theta) + b*cos(theta)*sin(t);

    plot_vecs = [x_plot; y_plot];
end