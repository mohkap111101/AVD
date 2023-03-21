function geometry = createGeometry(t_s, t, h, b,wing)

    geometry = wing;
    geometry.t_s = t_s;
    geometry.t = t;
    geometry.h = h;
    geometry.b = b;
    geometry.L = 0.5;

end