#ifndef MARCHINGTILES_H
#define MARCHINGTILES_H

#include "array3.h"
#include "hashtable.h"
#include "vec.h"

struct MarchingTiles
{
   std::vector<Vec3i> tri;
   std::vector<Vec3f> x;
   Vec3f origin;
   float dx;
   Array3f phi;

   explicit MarchingTiles(float dx_=1)
      : origin(0), dx(dx_)
   {}

   explicit MarchingTiles(const Vec3f &origin_, float dx_=1)
      : origin(origin_), dx(dx_)
   {}

   void contour_grid(void);

   private:
   HashTable<Vec6i,unsigned int> edge_cross; // stores vertices that have been created already at given edge crossings
   void contour_tile(int i, int j, int k); // add triangles for contour in the given tile (starting at grid point (4*i,4*j,4*k))
   void contour_tet(const Vec3i& x0, const Vec3i& x1, const Vec3i& x2, const Vec3i& x3, float p0, float p1, float p2, float p3);
   int find_edge_cross(const Vec3i& x0, const Vec3i& x1, float p0, float p1);
};

#endif
