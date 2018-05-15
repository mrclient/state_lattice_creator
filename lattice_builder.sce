// Abbreviations used:
// res = result

clear ;
path = get_absolute_file_path('lattice_builder.sce');
exec(path + 'utils.sci');

// INPUT DATA
N = 7;
step = 0.25;
theta_des = [2,0; 2,2; 0,2];
k = [0, -4, 4];

sz = size(theta_des);
theta_num = sz(1);

counter = 1
kolor = 1;
for yf = 0:floor(N/2)
    for xf = 0:floor(N/2)
        if xf == 0 & yf == 0 then
            continue;
        end
        for i = 1:theta_num
            theta_0 = atan(theta_des(i,2), theta_des(i,1));
            new_xy = rotate([xf*step;yf*step], -theta_0);
            new_xf = new_xy(1);
            new_yf = new_xy(2);
            for j = 1:theta_num
                theta_f = atan(theta_des(j,2), theta_des(j,1));

                // the ban of S-shaped turns
                if clean([cos(theta_0), sin(theta_0)] * [-yf;xf] * [cos(theta_f), sin(theta_f)] * [-yf;xf]) > 0 then
                    continue;
                end

                theta_f = theta_f - theta_0;
                for k0 = k
                    for kf = k
                        counter = counter + 1;
                        [a, b, c, d, sf, status] = find_curve(k0, new_xf, new_yf, theta_f, kf);
                        // plotting
                        if status == 1 then
                            [x, y] = points(a, b, c, d, sf);
                            xy = rotate([x;y], theta_0)
                            plot2d(xy(1,:), xy(2,:), kolor);
                        else
                            printf("Bad status for kolor = %d\n", kolor);
                            printf("k0 = %f, xf = %f, yf = %f, theta_f = %f, kf = %f\n", k0, xf, yf, theta_f+theta_0, kf);
                            printf("k0 = %f, xf = %f, yf = %f, theta_f = %f, kf = %f\n", k0, new_xf, new_yf, theta_f, kf);
                            printf("----------------\n");
                        end
                    end
                end
            end
            kolor = kolor + 1;
            if kolor == 8 then
                kolor = kolor + 1
            end
        end
    end
end

