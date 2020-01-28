#include "node.h"

icy::Node::Node()
{
    id = altId = globalNodeId = -1;
    x0 = y0 = z0 = 0;
    cx = cy = cz = 0;
    ux = uy = uz = 0;
    vx = vy = vz = 0;
    ax = ay = az = 0;
    fx = fy = fz = 0;
    unx = uny = unz = 0;
    vnx = vny = vnz = 0;
    anx = any = anz = 0;
    tx = ty = tz = 0;
    dux = duy = duz = 0;
    anchored = isSurface = false;
}

void icy::Node::Initialize(double x, double y, double z, int _id)
{
    id=_id;
    cx = x0 = x;
    cy = y0 = y;
    cz = z0 = z;
}



void icy::Node::AcceptTentativeValues(double h)
{
    vnx = (unx - ux) / h;
    vny = (uny - uy) / h;
    vnz = (unz - uz) / h;

    ax = (vnx - vx) / h;
    ay = (vny - vy) / h;
    az = (vnz - vz) / h;

    vx = vnx; vy = vny; vz = vnz;

    ux = unx; uy = uny; uz = unz;
    cx = x0 + ux; cy = y0 + uy; cz = z0 + uz;
    dux = duy = duz = 0;
}


void icy::Node::InferTentativeValues(double h, double beta, double gamma)
{
    auto InferTentativeUVA = [](double _du, double _u, double _v, double _a, double _h,
            double &_un, double &_vn, double &_an, double _beta, double _gamma)
    {
        _un = _u + _du;
        _an = _a * (1.0 - 1.0 / (2*_beta)) + _du / (_h*_h*_beta) + _v * (-1.0/(_h*_beta));
        _vn = _v + _h * ((1.0 - _gamma) * _a + _gamma * _an);
    };

    InferTentativeUVA(dux, ux, vx, ax, h, unx, vnx, anx, beta, gamma);
    InferTentativeUVA(duy, uy, vy, ay, h, uny, vny, any, beta, gamma);
    InferTentativeUVA(duz, uz, vz, az, h, unz, vnz, anz, beta, gamma);
    tx = x0 + unx;
    ty = y0 + uny;
    tz = z0 + unz;
}
