#ifndef KDOP24_H
#define KDOP24_H

namespace icy { class kDOP24; }

class icy::kDOP24
{
    /*
      k=18 (Aabb + 12 diagonal planes that "cut off" some space of the edges):
  (-1,0,0) and (1,0,0)  -> indices 0 and 9
  (0,-1,0) and (0,1,0)  -> indices 1 and 10
  (0,0,-1) and (0,0,1)  -> indices 2 and 11
  (-1,-1,0) and (1,1,0) -> indices 3 and 12
  (-1,0,-1) and (1,0,1) -> indices 4 and 13
  (0,-1,-1) and (0,1,1) -> indices 5 and 14
  (-1,1,0) and (1,-1,0) -> indices 6 and 15
  (-1,0,1) and (1,0,-1) -> indices 7 and 16
  (0,-1,1) and (0,1,-1) -> indices 8 and 17
*/
public:
    double d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11; // lower boundaries
    double d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23; // higher boundaries
    kDOP24();
    void Reset();
    bool Overlaps(kDOP24 &b);
    void Expand(double x, double y, double z);
    void Expand(kDOP24 &b);
    void Dimensions(double &dx, double &dy, double &dz);
    double centerX() { return (d0 + d12) / 2; }
    double centerY() { return (d1 + d13) / 2; }
    double centerZ() { return (d2 + d14) / 2; }

private:
    inline void MinMax(double p, double &mi, double &ma);
    inline void MinMax(double a, double b, double &mi, double &ma);
    double min(double a, double b) { return a < b ? a : b; }
    double max(double a, double b) { return a > b ? a : b; }
};

#endif // KDOP24_H
