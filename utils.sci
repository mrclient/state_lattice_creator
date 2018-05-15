//////////////////////////////////////////////////////////
//// The file contains some functions for calculating ////
////////// generalized Fresnel integrals and /////////////
///////////////// some other things //////////////////////
//////////////////////////////////////////////////////////

function out = cosTheta(n, a, b, c, d, s)
    out = s^n * cos(a*s + b*s^2/2 + c*s^3/3 + d*s^4/4);
endfunction

function out = sinTheta(n, a, b, c, d, s)
    out = s^n * sin(a*s + b*s^2/2 + c*s^3/3 + d*s^4/4);
endfunction

function out = C(n, a, b, c, d, sf)
    deff('[out]=aux(s)','out=cosTheta(n, a, b, c, d, s)')
    [out,  err, ierr] = intg(0, sf, aux);
endfunction

function out = S(n, a, b, c, d, sf)
    deff('[out]=aux(s)','out=sinTheta(n, a, b, c, d, s)')
    [out, err, ierr] = intg(0, sf, aux);
endfunction

// function of constraints that must be equal to zero
function [y] = g(x, k0, xf, yf, theta_f, kf)
    a = k0;
    b = x(1);
    c = x(2);
    d = x(3);
    sf = x(4);
    y(1) = a + b*sf + c*sf^2 + d*sf^3 - kf;
    y(2) = a*sf + b*sf^2/2 + c*sf^3/3 + d*sf^4/4 - theta_f;
    y(3) = C(0, a, b, c, d, sf) - xf;
    y(4) = S(0, a, b, c, d, sf) - yf;
endfunction

// the function gradient
function [y] = grad_g(x, k0, xf, yf, theta_f, kf)
    a = k0;
    b = x(1);
    c = x(2);
    d = x(3);
    sf = x(4);
    y(1,:) = [sf, sf^2, sf^3, b + 2*c*sf + 3*d*sf^2];
    y(2,:) = [sf^2/2, sf^3/3, sf^4/4, a + b*sf + c*sf^2 + d*sf^3];
    y(3,:) = [-S(1, a, b, c, d, sf), -S(2, a, b, c, d, sf)/2, -S(3, a, b, c, d, sf)/3, cosTheta(0, a, b, c, d, sf)];
    y(4,:) = [C(1, a, b, c, d, sf), C(2, a, b, c, d, sf)/2, C(3, a, b, c, d, sf)/3, sinTheta(0, a, b, c, d, sf)];
endfunction

function [a, b, c, d, sf, status] = find_curve(k0, xf, yf, theta_f, kf)
    res_0 = [0; 0; 0; sqrt(xf^2+yf^2)];
    [res, g_val, status] = fsolve(res_0, list(g, k0, xf, yf, theta_f, kf), list(grad_g, k0, xf, yf, theta_f, kf));
//    disp(res,'coeffs from fsolve:');
//    disp(g_val,'value of g from fsolve:');
//    disp(status,'info from fsolve:');
    a = k0;
    b = res(1);
    c = res(2);
    d = res(3);
    sf = res(4);
endfunction


function [x, y, theta] = poses(a, b, c, d, sf)
    i = 1;
    for s = 0:(0.01*sf):sf
        x(:,i) = C(0, a, b, c, d, s);
        y(:,i) = S(0, a, b, c, d, s);
        theta(:,i) = a*s + b*s^2/2 + c*s^3/3 + d*s^4/4;
        i = i + 1;
    end
endfunction


function write_points(des, num, xy, theta)

    mfprintf(des, "%d\n", num);

    for i = 1:max(size(xy))
        mfprintf(des, "%f ", xy(1, i));
    end
    mfprintf(des, "\n");

    for i = 1:max(size(xy))
        mfprintf(des, "%f ", xy(2, i));
    end
    mfprintf(des, "\n");

    for i = 1:max(size(xy))
        mfprintf(des, "%f ", theta(i));
    end
    mfprintf(des, "\n");
endfunction


function write_vertices(des, theta_0, k0, theta_f, kf, dx, dy, sf, num)
    mfprintf(des, "%f %f %f %f %f %f %f %d %d\n", theta_0, k0, theta_f, kf, dx, dy, sf, num, 0)
    mfprintf(des, "%f %f %f %f %f %f %f %d %d\n", %pi-theta_0, -k0, %pi-theta_f, -kf, -dx, dy, sf, num, 1);
    mfprintf(des, "%f %f %f %f %f %f %f %d %d\n", -theta_0, -k0, -theta_f, -kf, dx, -dy, sf, num, 2);
    mfprintf(des, "%f %f %f %f %f %f %f %d %d\n", theta_0-%pi, k0, theta_f-%pi, kf, -dx, -dy, sf, num, 3);
endfunction
