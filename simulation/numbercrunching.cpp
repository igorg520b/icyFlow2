#include "numbercrunching.h"
#include <stdexcept>
#include <cfloat>

const double icy::NumberCrunching::EPS = 1e-10;

std::vector<int> icy::NumberCrunching::resultingList;
std::unordered_set<long long> icy::NumberCrunching::NL2set;
std::vector<long long> icy::NumberCrunching::NL2vector;
std::vector<icy::CPResult> icy::NumberCrunching::cprList;

icy::NumberCrunching::NumberCrunching()
{

}

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
        long nodeIdx = (long)(value >> 32);
        long elemIdx = (long)value;
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

#pragma omp parallel for
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




