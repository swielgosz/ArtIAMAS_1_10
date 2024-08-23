function plot_vecs = plot_wsn_radii(sensor, radii)
    u = sensor(1); v = sensor(2);
    t = 0:0.01:2*pi;
    x_plot = u+radii*cos(t);
    y_plot = v+radii*sin(t);
    plot_vecs = [x_plot; y_plot];
end