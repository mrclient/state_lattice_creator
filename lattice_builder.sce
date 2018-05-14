// Abbreviations used:
// res = result

clear ;
path = get_absolute_file_path('lattice_builder.sce');
exec(path+'utils.sci');

k0 = 0.5;
xf = 3;
yf = 5;
theta_f = %pi/2;
kf = 0.0;

[a, b, c, d, sf, status] = find_curve(k0, xf, yf, theta_f, kf);

// plotting
[x,y] = points(a, b, c, d, sf);
plot2d(x,y);
