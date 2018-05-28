// Abbreviations used:
// res = result

clear ;
path = get_absolute_file_path('lattice_builder.sce');
exec(path + 'utils.sci');

// INPUT DATA
N = 7;
step = 0.25;
k = [0, -4, 4];

// to make an array of all needed angle
i = 0;
for x = 0:floor(N/2)
    for y = 0:floor(N/2)
        if x == 0 & y == 0 then
            continue;
        end
        new_angle = atan(y, x);
        if i==0 | find(theta_des==new_angle) == [] then
            i = i + 1;
            theta_des(:,i) = new_angle;
            if new_angle >= %pi/4 & new_angle < %pi/2 then
                i = i + 1;
                theta_des(:,i) = %pi - new_angle;
            end
            if new_angle <= %pi/4 & new_angle > 0 then
                i = i + 1;
                theta_des(:,i) = -new_angle;
            end
        end
    end
end
theta_des = gsort(theta_des, 'c', 'i');
theta_num = ceil(i/2);

total = 1;
success = 1;
kolor = 1;
v_file = mopen(path + "vertices.txt", "wt");
p_file = mopen(path + "points.txt", "wt");
for i = 1:theta_num
    theta_0 = theta_des(i + floor(theta_num/2));
    for k0 = k
        for yf = 0:floor(N/2)
            for xf = 0:floor(N/2)
                if xf == 0 & yf == 0 then
                    continue;
                end
                new_xy = rotate([xf*step;yf*step], -theta_0);
                new_xf = new_xy(1);
                new_yf = new_xy(2);

                // NoT = number of theta
                NoT = find(theta_des==atan(yf, xf));
                NoT = NoT - floor(theta_num/2);
                inf = NoT;
                sup = NoT + theta_num - 1;
                if NoT == 1 then
                    sup = sup - floor(theta_num/2);
                elseif NoT == theta_num
                    inf = inf + floor(theta_num/2);
                end
                for j = inf:sup
                    theta_f = theta_des(j);
                    // the ban of S-shaped turns
                    if clean([cos(theta_0), sin(theta_0)] * [-yf;xf] * [cos(theta_f), sin(theta_f)] * [-yf;xf]) > 0 then
                        continue;
                    end
                    theta_f = theta_f - theta_0;
                    for kf = k
                        total = total + 1;
                        [a, b, c, d, sf, status] = find_curve(k0, new_xf, new_yf, theta_f, kf);

                        // checking for results; saving; plotting
                        if status == 1 then
                            [x, y, theta] = poses(a, b, c, d, sf);
                            xy = rotate([x;y], theta_0)
                            theta = theta + theta_0;
                            write_vertices(v_file, theta_0, k0, theta_f+theta_0, kf, xf*step, yf*step, sf, success);
                            write_points(p_file, success, xy, theta)
                            plot2d(xy(1,:), xy(2,:), kolor);
                            success = success + 1;
                        else
                            printf("Bad status for kolor = %d\n", kolor);
                            printf("k0 = %f, xf = %f, yf = %f, theta_f = %f, kf = %f\n", k0, xf, yf, theta_f+theta_0, kf);
                            printf("k0 = %f, xf = %f, yf = %f, theta_f = %f, kf = %f\n", k0, new_xf, new_yf, theta_f, kf);
                            printf("----------------\n");
                        end
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
mclose(p_file);
mclose(v_file);
printf("total = %d\nsuccess = %d\npercentage = %.2f\n", total, success, 100*success/total);
