#include "numbercrunching.h"
#include <stdexcept>
#include <cfloat>
#include <cmath>


const double icy::NumberCrunching::EPS = 1e-10;

std::vector<int> icy::NumberCrunching::resultingList;
std::unordered_set<long long> icy::NumberCrunching::NL2set;
std::vector<long long> icy::NumberCrunching::NL2vector;
std::vector<icy::CPResult> icy::NumberCrunching::cprList;
SymmetricEigensolver3x3<double> icy::NumberCrunching::solver;


void icy::NumberCrunching::NarrowPhase(std::vector<Element*> &broadList, MeshCollection &mc)
{
    int nTetra = broadList.size();
    if(nTetra % 2 != 0) throw std::runtime_error("broarList.size is odd");
    int nPairs = nTetra / 2;

    resultingList.reserve(nPairs+1000);
    resultingList.resize(nPairs);

#pragma omp parallel for
    for(int i=0;i<nPairs;i++) resultingList[i] = NarrowPhaseTwoElems(broadList[i * 2], broadList[i * 2 + 1]);

    NL2set.clear();
    NL2set.reserve(nPairs);

    // remove duplications from resultingList
    for (int i = 0; i < nPairs; i++) {
        int bits = resultingList[i];
        if(bits == 0) continue;
        if ((bits & 1) != 0) NL2set.insert((long long) (broadList[i * 2 + 1]->vrts[0]->globalNodeId) << 32 | (broadList[i * 2]->globalElementId));
        if ((bits & 2) != 0) NL2set.insert((long long) (broadList[i * 2 + 1]->vrts[1]->globalNodeId) << 32 | (broadList[i * 2]->globalElementId));
        if ((bits & 4) != 0) NL2set.insert((long long) (broadList[i * 2 + 1]->vrts[2]->globalNodeId) << 32 | (broadList[i * 2]->globalElementId));
        if ((bits & 8) != 0) NL2set.insert((long long) (broadList[i * 2 + 1]->vrts[3]->globalNodeId) << 32 | (broadList[i * 2]->globalElementId));
        if ((bits & 16) != 0) NL2set.insert((long long) (broadList[i * 2]->vrts[0]->globalNodeId) << 32 | (broadList[i * 2 + 1]->globalElementId));
        if ((bits & 32) != 0) NL2set.insert((long long) (broadList[i * 2]->vrts[1]->globalNodeId) << 32 | (broadList[i * 2 + 1]->globalElementId));
        if ((bits & 64) != 0) NL2set.insert((long long) (broadList[i * 2]->vrts[2]->globalNodeId) << 32 | (broadList[i * 2 + 1]->globalElementId));
        if ((bits & 128) != 0) NL2set.insert((long long) (broadList[i * 2]->vrts[3]->globalNodeId) << 32 | (broadList[i * 2 + 1]->globalElementId));
    }

    nPairs = NL2set.size();
    cprList.resize(nPairs);
    if (nPairs == 0) return;

    NL2vector.clear();
    NL2vector.reserve(nPairs);
    NL2vector.insert(NL2vector.end(), NL2set.begin(),NL2set.end());

#pragma omp parallel for
    for(int i=0;i<nPairs;i++) {
        long long value = NL2vector[i];
        int nodeIdx = (int)(value >> 32);
        int elemIdx = (int)(value & 0xffffffff);
        Node *nd = mc.allNodes[nodeIdx];
        cprList[i].nd = nd;
        cprList[i].fc = FindClosestFace(nd, mc.surfaceElements[elemIdx]);
    }

}

double icy::NumberCrunching::dtn(
            double f1x, double f1y, double f1z,
            double f2x, double f2y, double f2z,
            double f3x, double f3y, double f3z,
        double ndx, double ndy, double ndz)
{

    double t0[3], t1[3], t2[3];
    double edge0[3], edge1[3];
    double v0[3], sourcePosition[3];

    t0[0] = f1x; t0[1] = f1y; t0[2] = f1z;
    t1[0] = f2x; t1[1] = f2y; t1[2] = f2z;
    t2[0] = f3x; t2[1] = f3y; t2[2] = f3z;
    sourcePosition[0] = ndx; sourcePosition[1] = ndy; sourcePosition[2] = ndz;

    SUB(edge0, t1, t0);
    SUB(edge1, t2, t0);
    SUB(v0, t0, sourcePosition);

    double a = DOT(edge0, edge0);
    double b = DOT(edge0, edge1);
    double c = DOT(edge1, edge1);
    double d = DOT(edge0, v0);
    double e = DOT(edge1, v0);

    double det = a * c - b * b;
    double s = b * e - c * d;
    double t = b * d - a * e;

    if (s + t < det)
    {
        if (s < 0)
        {
            if (t < 0)
            {
                if (d < 0)
                {
                    s = clamp(-d / a);
                    t = 0;
                }
                else
                {
                    s = 0;
                    t = clamp(-e / c);
                }
            }
            else
            {
                s = 0;
                t = clamp(-e / c);
            }
        }
        else if (t < 0)
        {
            s = clamp(-d / a);
            t = 0;
        }
        else
        {
            double invDet = 1.0 / det;
            s *= invDet;
            t *= invDet;
        }
    }
    else
    {
        if (s < 0)
        {
            double tmp0 = b + d;
            double tmp1 = c + e;
            if (tmp1 > tmp0)
            {
                double numer = tmp1 - tmp0;
                double denom = a - 2 * b + c;
                s = clamp(numer / denom);
                t = 1 - s;
            }
            else
            {
                t = clamp(-e / c);
                s = 0;
            }
        }
        else if (t < 0)
        {
            if (a + d > b + e)
            {
                double numer = c + e - b - d;
                double denom = a - 2 * b + c;
                s = clamp(numer / denom);
                t = 1 - s;
            }
            else
            {
                s = clamp(-e / c);
                t = 0;
            }
        }
        else
        {
            double numer = c + e - b - d;
            double denom = a - 2 * b + c;
            s = clamp(numer / denom);
            t = 1 - s;
        }
    }

    double d1[3];

    d1[0] = t0[0] + s * edge0[0] + t * edge1[0] - sourcePosition[0];
    d1[1] = t0[1] + s * edge0[1] + t * edge1[1] - sourcePosition[1];
    d1[2] = t0[2] + s * edge0[2] + t * edge1[2] - sourcePosition[2];

    double sqdist = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];

    return sqdist; // squared
}

icy::Face* icy::NumberCrunching::FindClosestFace(Node *nd, Element *elem)
{
    double smallestDistance = DBL_MAX;
    Face *closestFace = nullptr;
    if(elem->adjFaces.size() == 0) throw std::runtime_error("adjacent faces of teh element not initialized");

    for(auto const &f : elem->adjFaces) {
        double dist = dtn(f->vrts[0]->tx, f->vrts[0]->ty, f->vrts[0]->tz,
        f->vrts[1]->tx, f->vrts[1]->ty, f->vrts[1]->tz,
            f->vrts[2]->tx, f->vrts[2]->ty, f->vrts[2]->tz,
        nd->tx, nd->ty, nd->tz);
        if (smallestDistance >= dist)
        {
            smallestDistance = dist;
            closestFace = f;
        }
    }
    if(closestFace == nullptr) throw std::runtime_error("closest face not found");
    return closestFace;
}








int icy::NumberCrunching::NarrowPhaseTwoElems(Element *tetra1, Element *tetra2)
  {
      double nds[24];
      int nd_idxs[8];
      for (int i = 0; i < 4; i++)
      {
          nd_idxs[i] = tetra1->vrts[i]->globalNodeId;
          nd_idxs[i + 4] = tetra2->vrts[i]->globalNodeId;
          nds[i * 3 + 0] = tetra1->vrts[i]->tx;
          nds[i * 3 + 1] = tetra1->vrts[i]->ty;
          nds[i * 3 + 2] = tetra1->vrts[i]->tz;

          nds[(i + 4) * 3 + 0] = tetra2->vrts[i]->tx;
          nds[(i + 4) * 3 + 1] = tetra2->vrts[i]->ty;
          nds[(i + 4) * 3 + 2] = tetra2->vrts[i]->tz;
      }

      // verify that elements are non-adjacent
      for (int i = 0; i < 4; i++)
          for (int j = 4; j < 8; j++)
              if (nd_idxs[i] == nd_idxs[j]) return 0;

      // b-values for elements
      double bv0[9];
      double bv1[9];

      Bvalues(nds[0], nds[1], nds[2],
      nds[3], nds[4], nds[5],
      nds[6], nds[7], nds[8],
      nds[9], nds[10], nds[11],
      bv0[0], bv0[1], bv0[2],
      bv0[3], bv0[4], bv0[5],
      bv0[6], bv0[7], bv0[8]);

      Bvalues(nds[12], nds[13], nds[14],
      nds[15], nds[16], nds[17],
      nds[18], nds[19], nds[20],
      nds[21], nds[22], nds[23],
      bv1[0], bv1[1], bv1[2],
      bv1[3], bv1[4], bv1[5],
      bv1[6], bv1[7], bv1[8]);

      int result = 0;
      // perform tests
      bool bres;
      double x0 = nds[0];
      double y0 = nds[1];
      double z0 = nds[2];

      // test if 1 element nodes are inside element 0
      bres = ctest(bv0[0], bv0[1], bv0[2],
      bv0[3], bv0[4], bv0[5],
      bv0[6], bv0[7], bv0[8],
      nds[12] - x0, nds[13] - y0, nds[14] - z0);
      if (bres) result |= (1);

      bres = ctest(bv0[0], bv0[1], bv0[2],
      bv0[3], bv0[4], bv0[5],
      bv0[6], bv0[7], bv0[8],
      nds[15] - x0, nds[16] - y0, nds[17] - z0);
      if (bres) result |= (2);

      bres = ctest(bv0[0], bv0[1], bv0[2],
      bv0[3], bv0[4], bv0[5],
      bv0[6], bv0[7], bv0[8],
      nds[18] - x0, nds[19] - y0, nds[20] - z0);
      if (bres) result |= (4);

      bres = ctest(bv0[0], bv0[1], bv0[2],
      bv0[3], bv0[4], bv0[5],
      bv0[6], bv0[7], bv0[8],
      nds[21] - x0, nds[22] - y0, nds[23] - z0);
      if (bres) result |= (8);

      // test if 0 element nodes are inside element 1
      x0 = nds[12];
      y0 = nds[13];
      z0 = nds[14];

      bres = ctest(bv1[0], bv1[1], bv1[2],
      bv1[3], bv1[4], bv1[5],
      bv1[6], bv1[7], bv1[8],
      nds[0] - x0, nds[1] - y0, nds[2] - z0);
      if (bres) result |= (16);

      bres = ctest(bv1[0], bv1[1], bv1[2],
      bv1[3], bv1[4], bv1[5],
      bv1[6], bv1[7], bv1[8],
      nds[3] - x0, nds[4] - y0, nds[5] - z0);
      if (bres) result |= (32);

      bres = ctest(bv1[0], bv1[1], bv1[2],
      bv1[3], bv1[4], bv1[5],
      bv1[6], bv1[7], bv1[8],
      nds[6] - x0, nds[7] - y0, nds[8] - z0);
      if (bres) result |= (64);

      bres = ctest(bv1[0], bv1[1], bv1[2],
      bv1[3], bv1[4], bv1[5],
      bv1[6], bv1[7], bv1[8],
      nds[9] - x0, nds[10] - y0, nds[11] - z0);
      if (bres) result |= (128);

      return result;
  }

bool icy::NumberCrunching::ctest(double b11, double b12, double b13,
    double b21, double b22, double b23,
    double b31, double b32, double b33,
    double x, double y, double z)
{

    double y1 = x * b11 + y * b12 + z * b13;
    double y2 = x * b21 + y * b22 + z * b23;
    double y3 = x * b31 + y * b32 + z * b33;
    return (y1 > EPS && y2 > EPS && y3 > EPS && (y1 + y2 + y3) < (1 - EPS));
}


void icy::NumberCrunching::Bvalues(double x0, double y0, double z0,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3,
    double &b11, double &b12, double &b13,
    double &b21, double &b22, double &b23,
    double &b31, double &b32, double &b33)
{
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;
    a11 = x1 - x0;
    a12 = x2 - x0;
    a13 = x3 - x0;
    a21 = y1 - y0;
    a22 = y2 - y0;
    a23 = y3 - y0;
    a31 = z1 - z0;
    a32 = z2 - z0;
    a33 = z3 - z0;

    // inverse
    double det = a31 * (-a13 * a22 + a12 * a23) + a32 * (a13 * a21 - a11 * a23) + a33 * (-a12 * a21 + a11 * a22);
    b11 = (-a23 * a32 + a22 * a33) / det;
    b12 = (a13 * a32 - a12 * a33) / det;
    b13 = (-a13 * a22 + a12 * a23) / det;
    b21 = (a23 * a31 - a21 * a33) / det;
    b22 = (-a13 * a31 + a11 * a33) / det;
    b23 = (a13 * a21 - a11 * a23) / det;
    b31 = (-a22 * a31 + a21 * a32) / det;
    b32 = (a12 * a31 - a11 * a32) / det;
    b33 = (-a12 * a21 + a11 * a22) / det;
}


// point-triangle distance and derivatives
double pt(double (&x)[12], double (&fd)[12], double (&sd)[12][12], double &zeta2, double &zeta3);



// =========================================== collision response
void icy::NumberCrunching::OneCollision(double distanceEpsilonSqared, double k, CPResult &res)
{
    Node *nd = res.nd;
    Face *fc = res.fc;
    res.Clear();

    double x[12];
    x[0] = nd->tx;
    x[1] = nd->ty;
    x[2] = nd->tz;
    x[3] = fc->vrts[0]->tx;
    x[4] = fc->vrts[0]->ty;
    x[5] = fc->vrts[0]->tz;
    x[6] = fc->vrts[1]->tx;
    x[7] = fc->vrts[1]->ty;
    x[8] = fc->vrts[1]->tz;
    x[9] = fc->vrts[2]->tx;
    x[10] = fc->vrts[2]->ty;
    x[11] = fc->vrts[2]->tz;

    res.idxs[0] = nd->altId;
    res.idxs[1] = fc->vrts[0]->altId;
    res.idxs[2] = fc->vrts[1]->altId;
    res.idxs[3] = fc->vrts[2]->altId;

    res.nds[0] = nd;
    res.nds[1] = fc->vrts[0];
    res.nds[2] = fc->vrts[1];
    res.nds[3] = fc->vrts[2];

    // exclude colliding rigid surfaces
    if (res.idxs[0] < 0 && res.idxs[1] < 0 && res.idxs[2] < 0 && res.idxs[3] < 0) return;

    double fd[12];
    double sd[12][12]; //12x12
    double w[3];
    double dsq = pt(x, fd, sd, w[1], w[2]);

    if (dsq < distanceEpsilonSqared) return;

    w[0] = 1 - (w[1] + w[2]);

    double fx, fy, fz;
    fx = k * 0.5 * fd[0];
    fy = k * 0.5 * fd[1];
    fz = k * 0.5 * fd[2];

    res.fi[0] = -fx;
    res.fi[1] = -fy;
    res.fi[2] = -fz;
    res.fi[3] = w[0] * fx;
    res.fi[4] = w[0] * fy;
    res.fi[5] = w[0] * fz;
    res.fi[6] = w[1] * fx;
    res.fi[7] = w[1] * fy;
    res.fi[8] = w[1] * fz;
    res.fi[9] = w[2] * fx;
    res.fi[10] = w[2] * fy;
    res.fi[11] = w[2] * fz;

    for (int i = 0; i < 12; i++)
        for (int j = i; j < 12; j++)
            res.dfi[i][j] = res.dfi[j][i] = k * sd[i][j] / 2;
}

void icy::NumberCrunching::CollisionResponse(LinearSystem &ls, double DistanceEpsilon, double k)
{
    double distanceEpsilonSqared = DistanceEpsilon*DistanceEpsilon;
    int N = (int)cprList.size();

//#pragma omp parallel for
    for(int i=0;i<N;i++) OneCollision(distanceEpsilonSqared, k, cprList[i]);

    // distribute the values into linear system
    for(auto const &res : cprList)
    {
        for(int r=0;r<4;r++) {
            int ni = res.idxs[r];
            double fx = res.fi[r * 3 + 0];
            double fy = res.fi[r * 3 + 1];
            double fz = res.fi[r * 3 + 2];

            ls.AddToRHS(ni, fx, fy, fz);

            // add to node's force
            res.nds[r]->fx += fx;
            res.nds[r]->fy += fy;
            res.nds[r]->fz += fz;

            for (int c = 0; c < 4; c++)
            {
                int nj = res.idxs[c];
                ls.AddToLHS_Symmetric(ni, nj,
                res.dfi[r * 3 + 0][c * 3 + 0], res.dfi[r * 3 + 0][c * 3 + 1], res.dfi[r * 3 + 0][c * 3 + 2],
                res.dfi[r * 3 + 1][c * 3 + 0], res.dfi[r * 3 + 1][c * 3 + 1], res.dfi[r * 3 + 1][c * 3 + 2],
                res.dfi[r * 3 + 2][c * 3 + 0], res.dfi[r * 3 + 2][c * 3 + 1], res.dfi[r * 3 + 2][c * 3 + 2]);
            }
        }
    }
}


//======================= LINEAR TETRAHEDRON
void icy::NumberCrunching::AssembleElems(LinearSystem &ls, std::vector<Element*> &elasticElements, ModelPrms &prms)
{
}

void icy::NumberCrunching::ElementElasticity(Element *elem)
{

}


void icy::NumberCrunching::fastRotationMatrix(
        double p0x, double p0y, double p0z,
        double p1x, double p1y, double p1z,
        double p2x, double p2y, double p2z,
        double &r11, double &r12, double &r13,
        double &r21, double &r22, double &r23,
        double &r31, double &r32, double &r33)
{
    double d10x = p1x - p0x;
    double d10y = p1y - p0y;
    double d10z = p1z - p0z;

    double mag = sqrt(d10x*d10x + d10y*d10y + d10z*d10z);
    r11 = d10x / mag;
    r21 = d10y / mag;
    r31 = d10z / mag;

    // p2-p0
    double wx = p2x - p0x;
    double wy = p2y - p0y;
    double wz = p2z - p0z;

    // cross product
    double cx = -d10z * wy + d10y * wz;
    double cy = d10z * wx - d10x * wz;
    double cz = -d10y * wx + d10x * wy;

    mag = sqrt(cx*cx + cy*cy + cz*cz);
    r12 = cx / mag;
    r22 = cy / mag;
    r32 = cz / mag;

    r13 = r22 * r31 - r21 * r32;
    r23 = -r12 * r31 + r11 * r32;
    r33 = r12 * r21 - r11 * r22;
    mag = sqrt(r13*r13 + r23*r23 + r33*r33);
    r13 /= mag;
    r23 /= mag;
    r33 /= mag;
}

void icy::NumberCrunching::multABd(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33,
        double b11, double b12, double b13,
        double b21, double b22, double b23,
        double b31, double b32, double b33,
        double &m11, double &m12, double &m13,
        double &m21, double &m22, double &m23,
        double &m31, double &m32, double &m33)
{
    m11 = a11*b11 + a12*b21 + a13*b31; m12 = a11*b12 + a12*b22 + a13*b32; m13 = a11*b13 + a12*b23 + a13*b33;
    m21 = a21*b11 + a22*b21 + a23*b31; m22 = a21*b12 + a22*b22 + a23*b32; m23 = a21*b13 + a22*b23 + a23*b33;
    m31 = a31*b11 + a32*b21 + a33*b31; m32 = a31*b12 + a32*b22 + a33*b32; m33 = a31*b13 + a32*b23 + a33*b33;
}

void icy::NumberCrunching::multAX(
        double a11, double a12, double a13,
        double a21, double a22, double a23,
        double a31, double a32, double a33,
        double x1, double x2, double x3,
        double &y1, double &y2, double &y3)
{
    y1 = x1 * a11 + x2 * a12 + x3 * a13;
    y2 = x1 * a21 + x2 * a22 + x3 * a23;
    y3 = x1 * a31 + x2 * a32 + x3 * a33;
}


void icy::NumberCrunching::F_and_Df_Corotational(
    const double(&x0)[12], const double(&xc)[12],
    double(&f)[12], double(&Df)[12][12], double &V,
double(&sigma)[6], double (&principal_str)[3], const double (&E)[6][6])
{
    // Colorational formulation:
    // f = RK(Rt xc - x0)
    // Df = R K Rt

    // calculate K
    double x12, x13, x14, x23, x24, x34;
    double y12, y13, y14, y23, y24, y34;
    double z12, z13, z14, z23, z24, z34;
    double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4;
    double Jdet;

    x12 = x0[0] - x0[3]; x13 = x0[0] - x0[6]; x14 = x0[0] - x0[9]; x23 = x0[3] - x0[6]; x24 = x0[3] - x0[9]; x34 = x0[6] - x0[9];
    y12 = x0[1] - x0[4]; y13 = x0[1] - x0[7]; y14 = x0[1] - x0[10]; y23 = x0[4] - x0[7]; y24 = x0[4] - x0[10]; y34 = x0[7] - x0[10];
    z12 = x0[2] - x0[5]; z13 = x0[2] - x0[8]; z14 = x0[2] - x0[11]; z23 = x0[5] - x0[8]; z24 = x0[5] - x0[11]; z34 = x0[8] - x0[11];

    Jdet = -(x12 * (y23 * z34 - y34 * z23) + x23 * (y34 * z12 - y12 * z34) + x34 * (y12 * z23 - y23 * z12));
    V = Jdet / 6.;

    a1 = y24 * z23 - y23 * z24; b1 = x23 * z24 - x24 * z23; c1 = x24 * y23 - x23 * y24;
    a2 = y13 * z34 - y34 * z13; b2 = x34 * z13 - x13 * z34; c2 = x13 * y34 - x34 * y13;
    a3 = y24 * z14 - y14 * z24; b3 = x14 * z24 - x24 * z14; c3 = x24 * y14 - x14 * y24;
    a4 = -y13 * z12 + y12 * z13; b4 = -x12 * z13 + x13 * z12; c4 = -x13 * y12 + x12 * y13;

    a1 /= Jdet; a2 /= Jdet; a3 /= Jdet; a4 /= Jdet;
    b1 /= Jdet; b2 /= Jdet; b3 /= Jdet; b4 /= Jdet;
    c1 /= Jdet; c2 /= Jdet; c3 /= Jdet; c4 /= Jdet;

    double B[6][12] = {
        { a1, 0, 0, a2, 0, 0, a3, 0, 0, a4, 0, 0 },
        { 0, b1, 0, 0, b2, 0, 0, b3, 0, 0, b4, 0 },
        { 0, 0, c1, 0, 0, c2, 0, 0, c3, 0, 0, c4 },
        { b1, a1, 0, b2, a2, 0, b3, a3, 0, b4, a4, 0 },
        { 0, c1, b1, 0, c2, b2, 0, c3, b3, 0, c4, b4 },
        { c1, 0, a1, c2, 0, a2, c3, 0, a3, c4, 0, a4 } };

    double BtE[12][6] = {}; // result of multiplication (Bt x E)
    for (int r = 0; r < 12; r++)
        for (int c = 0; c < 6; c++)
            for (int i = 0; i < 6; i++) BtE[r][c] += B[i][r] * E[i][c];

    // K = Bt x E x B x V
    double K[12][12] = {};
    for (int r = 0; r < 12; r++)
        for (int c = 0; c < 12; c++)
            for (int i = 0; i < 6; i++) K[r][c] += BtE[r][i] * B[i][c] * V;

    double R0[3][3], R1[3][3], R[3][3];
    fastRotationMatrix(
        x0[0], x0[1], x0[2],
        x0[3], x0[4], x0[5],
        x0[6], x0[7], x0[8],
        R0[0][0], R0[0][1], R0[0][2],
        R0[1][0], R0[1][1], R0[1][2],
        R0[2][0], R0[2][1], R0[2][2]);

    fastRotationMatrix(
        xc[0], xc[1], xc[2],
        xc[3], xc[4], xc[5],
        xc[6], xc[7], xc[8],
        R1[0][0], R1[0][1], R1[0][2],
        R1[1][0], R1[1][1], R1[1][2],
        R1[2][0], R1[2][1], R1[2][2]);

    multABd(
        R1[0][0], R1[0][1], R1[0][2],
        R1[1][0], R1[1][1], R1[1][2],
        R1[2][0], R1[2][1], R1[2][2],
        R0[0][0], R0[1][0], R0[2][0],
        R0[0][1], R0[1][1], R0[2][1],
        R0[0][2], R0[1][2], R0[2][2],
        R[0][0], R[0][1], R[0][2],
        R[1][0], R[1][1], R[1][2],
        R[2][0], R[2][1], R[2][2]);

    double RK[12][12] = {};
    double RKRt[12][12] = {};

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            // RK = R * K
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++)
                        RK[3 * i + k][3 * j + l] += R[k][m] * K[3 * i + m][3 * j + l];
                }

            // RKRT = RK * R^T
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++)
                        RKRt[3 * i + k][3 * j + l] += RK[3 * i + k][3 * j + m] * R[l][m];
                }
        }

    // xr = Rt xc
    double xr[12] = {};
    multAX(R[0][0], R[1][0], R[2][0],
        R[0][1], R[1][1], R[2][1],
        R[0][2], R[1][2], R[2][2],
        xc[0], xc[1], xc[2],
        xr[0], xr[1], xr[2]);
    multAX(R[0][0], R[1][0], R[2][0],
        R[0][1], R[1][1], R[2][1],
        R[0][2], R[1][2], R[2][2],
        xc[3], xc[4], xc[5],
        xr[3], xr[4], xr[5]);
    multAX(R[0][0], R[1][0], R[2][0],
        R[0][1], R[1][1], R[2][1],
        R[0][2], R[1][2], R[2][2],
        xc[6], xc[7], xc[8],
        xr[6], xr[7], xr[8]);
    multAX(R[0][0], R[1][0], R[2][0],
        R[0][1], R[1][1], R[2][1],
        R[0][2], R[1][2], R[2][2],
        xc[9], xc[10], xc[11],
        xr[9], xr[10], xr[11]);

    for (int i = 0; i < 12; i++) xr[i] -= x0[i];

    // f = RK(Rt pm - mx)
    // Df = RKRt
    for (int i = 0; i < 12; i++)
        for (int j = 0; j < 12; j++) {
            f[i] += RK[i][j] * xr[j];
            Df[i][j] = RKRt[i][j];
        }

    // calculation of strain (rotation excluded) e = B.xr
    // B[6][12]
    double e[6] = {};
    for (int j = 0; j < 6; j++)
        for (int i = 0; i < 12; i++)
            e[j] += B[j][i] * xr[i];

    // calculation of stress (rotation excluded) s = E.e
    // E[6][6]
    sigma[0]=sigma[1]=sigma[2]=sigma[3]=sigma[4]=sigma[5]=0;
    for (int j = 0; j < 6; j++)
        for (int i = 0; i < 6; i++)
            sigma[j] += E[i][j] * e[i];

    // compute principal stresses
    Eigenvalues(sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5], principal_str);
}

void icy::NumberCrunching::Eigenvalues(
        double xx, double yy, double zz,
        double xy, double yz, double zx,
        double* eigenvalues)
{
    std::array<double, 3> eval;
    std::array<std::array<double, 3>, 3> evec;
    solver(xx, xy, zx, yy, yz, zz, false, -1, eval, evec);
    eigenvalues[0] = eval[0];
    eigenvalues[1] = eval[1];
    eigenvalues[2] = eval[2];
}



