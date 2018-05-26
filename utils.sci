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
    d = x(1);
    sf = x(2);
    b = ((-4*kf - 8*k0)*sf + d*sf^4 + 12*theta_f) / (2 * sf^2);
    c = -((-6*kf-6*k0)*sf + 3*d*sf^4 + 12*theta_f) / (2*sf^3);
    y(1) = C(0, a, b, c, d, sf) - xf;
    y(2) = S(0, a, b, c, d, sf) - yf;
endfunction


function [a, b, c, d, sf, status] = find_curve(k0, xf, yf, theta_f, kf)
    res_0 = [0; (theta_f^2/5 + 1)*sqrt(xf^2+yf^2)];
    [res, g_val, status] = fsolve(res_0, list(g, k0, xf, yf, theta_f, kf));
//    disp(res,'coeffs from fsolve:');
//    disp(g_val,'value of g from fsolve:');
//    disp(status,'info from fsolve:');
    a = k0;
    d = res(1);
    sf = res(2);
    b = ((-4*kf - 8*k0)*sf + d*sf^4 + 12*theta_f) / (2 * sf^2);
    c = -((-6*kf-6*k0)*sf + 3*d*sf^4 + 12*theta_f) / (2*sf^3);
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
