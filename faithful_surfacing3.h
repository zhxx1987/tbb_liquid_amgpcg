#ifndef FAITHFUL_SURFACING3_H
#define FAITHFUL_SURFACING3_H

#include "marching_tiles.h"
#include "sparse_matrix.h"

struct FaithfulSurfacing3
{
   // defining state
   std::vector<Vec3f> particle_x;
   float inner_radius, outer_radius;
   float grid_dx;
   float gs_damping;
   int shrink_steps;
   int smooth_steps;
   int free_smooth_steps;

   // useful stuff, including tri and x (the output mesh)
   Vec3f origin;
   int nx, ny, nz;
   MarchingTiles march;
   std::vector<Vec3i> &tri;
   std::vector<Vec3f> &x;
   Array3<std::vector<unsigned int> > particle_bucket;
   SparseMatrixf A, B;

   FaithfulSurfacing3(void);
   void run_surfacing(void);

   private:
   void init(void);
   void build_phi_and_particle_buckets(void);
   void build_A(void);
   void shrink_mesh(bool enforce_fidelity);
   void build_B(void);
   void smooth_mesh(bool enforce_fidelity);
   void make_point_faithful(Vec3f& xp);
};

#endif
