#include "marching_tiles.h"

// definition of acute tile
static const int num_nodes=40;
static Vec3i node[num_nodes]={
   Vec3i(1,0,0),
   Vec3i(2,2,0),
   Vec3i(1,4,0),
   Vec3i(3,4,0),
   Vec3i(2,6,0),
   Vec3i(1,0,4),
   Vec3i(3,0,4),
   Vec3i(2,1,2),
   Vec3i(0,2,1),
   Vec3i(0,2,3),
   Vec3i(2,2,4),   // 10
   Vec3i(1,4,4),
   Vec3i(3,4,4),
   Vec3i(0,4,2),
   Vec3i(2,3,2),
   Vec3i(2,5,2),
   Vec3i(4,2,1),
   Vec3i(4,2,3),
   Vec3i(5,4,4),
   Vec3i(4,4,2),
   Vec3i(0,6,1),   // 20
   Vec3i(0,6,3),
   Vec3i(2,6,4),
   Vec3i(2,7,2),
   Vec3i(4,6,1),
   Vec3i(4,6,3),
   Vec3i(0,2,5),
   Vec3i(2,3,6),
   Vec3i(2,5,6),
   Vec3i(4,2,5),
   Vec3i(4,4,6),   // 30
   Vec3i(4,6,5),
   Vec3i(2,1,6),
   Vec3i(0,0,2),
   Vec3i(5,0,4),
   Vec3i(4,0,2),
   Vec3i(3,0,0),
   Vec3i(4,0,6),
   Vec3i(5,0,0),
   Vec3i(5,4,0)
};
static const int num_tets=46;
static Vec4i tet[num_tets]={
   Vec4i(2,3,15,14),
   Vec4i(2,15,13,14),
   Vec4i(6,29,17,10),
   Vec4i(6,17,34,35),
   Vec4i(12,14,17,10),
   Vec4i(0,36,1,7),
   Vec4i(29,12,18,17),
   Vec4i(3,14,16,19),
   Vec4i(3,15,14,19),
   Vec4i(14,16,1,3),
   Vec4i(0,7,8,33),
   Vec4i(7,16,36,1),
   Vec4i(9,7,5,33),
   Vec4i(8,7,9,33),
   Vec4i(14,12,17,19),
   Vec4i(12,29,10,17),
   Vec4i(14,9,13,8),
   Vec4i(8,2,14,1),
   Vec4i(14,17,16,19),
   Vec4i(17,14,7,10),
   Vec4i(16,7,36,35),
   Vec4i(17,6,7,35),
   Vec4i(11,15,14,13),
   Vec4i(3,2,1,14),
   Vec4i(9,14,7,8),
   Vec4i(9,14,11,10),
   Vec4i(1,8,0,7),
   Vec4i(14,8,1,7),
   Vec4i(12,15,19,14),
   Vec4i(19,39,16,3),
   Vec4i(17,12,18,19),
   Vec4i(14,9,11,13),
   Vec4i(14,12,11,10),
   Vec4i(2,8,14,13),
   Vec4i(6,17,7,10),
   Vec4i(5,9,26,10),
   Vec4i(14,9,7,10),
   Vec4i(11,9,10,26),
   Vec4i(7,9,5,10),
   Vec4i(17,6,34,29),
   Vec4i(6,7,5,10),
   Vec4i(15,12,11,14),
   Vec4i(16,14,1,7),
   Vec4i(7,17,35,16),
   Vec4i(38,35,16,36),
   Vec4i(14,16,17,7)
};

void MarchingTiles::
contour_grid(void)
{
   tri.resize(0);
   x.resize(0);
   edge_cross.clear();
   for(int k=0; 4*k+6<phi.nk; ++k) for(int j=0; 4*j+7<phi.nj; ++j) for(int i=0; 4*i+5<phi.ni; ++i)
      contour_tile(i,j,k);
}

void MarchingTiles::
contour_tile(int i, int j, int k)
{
   Vec3i cell(4*i,4*j,4*k);
   for(int t=0; t<num_tets; ++t){
      int a, b, c, d; assign(tet[t], a, b, c, d);
      contour_tet(cell+node[a], cell+node[b], cell+node[c], cell+node[d],
                  phi(cell[0]+node[a][0], cell[1]+node[a][1], cell[2]+node[a][2]),
                  phi(cell[0]+node[b][0], cell[1]+node[b][1], cell[2]+node[b][2]),
                  phi(cell[0]+node[c][0], cell[1]+node[c][1], cell[2]+node[c][2]),
                  phi(cell[0]+node[d][0], cell[1]+node[d][1], cell[2]+node[d][2]));
   }
}

// contour the tet with given grid point vertices and function values
// --- corners arranged so that 0-1-2 uses right-hand-rule to get to 3
void MarchingTiles::
contour_tet(const Vec3i& x0, const Vec3i& x1, const Vec3i& x2, const Vec3i& x3, float p0, float p1, float p2, float p3)
{
   // guard against topological degeneracies
   if(p0==0) p0=1e-30f;
   if(p1==0) p1=1e-30f;
   if(p2==0) p2=1e-30f;
   if(p3==0) p3=1e-30f;

   if(p0<0){
      if(p1<0){
         if(p2<0){
            if(p3<0){
               return; // no contour here
            }else // p3>=0
               tri.push_back(Vec3i(find_edge_cross(x0,x3,p0,p3),
                                   find_edge_cross(x1,x3,p1,p3),
                                   find_edge_cross(x2,x3,p2,p3)));
         }else{ // p2>=0
            if(p3<0)
               tri.push_back(Vec3i(find_edge_cross(x0,x2,p0,p2),
                                   find_edge_cross(x3,x2,p3,p2),
                                   find_edge_cross(x1,x2,p1,p2)));
            else{ // p3>=0
               tri.push_back(Vec3i(find_edge_cross(x0,x3,p0,p3),
                                   find_edge_cross(x1,x3,p1,p3),
                                   find_edge_cross(x0,x2,p0,p2)));
               tri.push_back(Vec3i(find_edge_cross(x1,x3,p1,p3),
                                   find_edge_cross(x1,x2,p1,p2),
                                   find_edge_cross(x0,x2,p0,p2)));
            }
         }
      }else{ // p1>=0
         if(p2<0){
            if(p3<0)
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x2,x1,p2,p1),
                                   find_edge_cross(x3,x1,p3,p1)));
            else{ // p3>=0
               tri.push_back(Vec3i(find_edge_cross(x0,x3,p0,p3),
                                   find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x2,x3,p2,p3)));
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x2,x1,p2,p1),
                                   find_edge_cross(x2,x3,p2,p3)));
            }
         }else{ // p2>=0
            if(p3<0){
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x0,x2,p0,p2),
                                   find_edge_cross(x3,x2,p3,p2)));
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x3,x2,p3,p2),
                                   find_edge_cross(x3,x1,p3,p1)));
            }else // p3>=_0
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x0,x2,p0,p2),
                                   find_edge_cross(x0,x3,p0,p3)));
         }
      }
   }else{ // p0>=0
      if(p1<0){
         if(p2<0){
            if(p3<0)
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x0,x3,p0,p3),
                                   find_edge_cross(x0,x2,p0,p2)));
            else{ // p3>=0
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x3,x1,p3,p1),
                                   find_edge_cross(x3,x2,p3,p2)));
               tri.push_back(Vec3i(find_edge_cross(x3,x2,p3,p2),
                                   find_edge_cross(x0,x2,p0,p2),
                                   find_edge_cross(x0,x1,p0,p1)));
            }
         }else{ // p2>=0
            if(p3<0){
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x0,x3,p0,p3),
                                   find_edge_cross(x3,x2,p3,p2)));
               tri.push_back(Vec3i(find_edge_cross(x0,x1,p0,p1),
                                   find_edge_cross(x3,x2,p3,p2),
                                   find_edge_cross(x2,x1,p2,p1)));
            }else // p3>=0
               tri.push_back(Vec3i(find_edge_cross(x1,x0,p1,p0),
                                   find_edge_cross(x1,x3,p1,p3),
                                   find_edge_cross(x1,x2,p1,p2)));
         }
      }else{ // p1>=0
         if(p2<0){
            if(p3<0){
               tri.push_back(Vec3i(find_edge_cross(x1,x3,p1,p3),
                                   find_edge_cross(x0,x3,p0,p3),
                                   find_edge_cross(x0,x2,p0,p2)));
               tri.push_back(Vec3i(find_edge_cross(x1,x3,p1,p3),
                                   find_edge_cross(x0,x2,p0,p2),
                                   find_edge_cross(x1,x2,p1,p2)));
            }else // p3>=0
               tri.push_back(Vec3i(find_edge_cross(x0,x2,p0,p2),
                                   find_edge_cross(x1,x2,p1,p2),
                                   find_edge_cross(x3,x2,p3,p2)));
         }else{ // p2>=0
            if(p3<0)
               tri.push_back(Vec3i(find_edge_cross(x0,x3,p0,p3),
                                   find_edge_cross(x2,x3,p2,p3),
                                   find_edge_cross(x1,x3,p1,p3)));
            else{ // p3>=0
               return; // assume no degenerate cases (where some of the p's are zero)
            }
         }
      }
   }
}

// return the vertex of the edge crossing (create it if necessary) between given grid points and function values
int MarchingTiles::
find_edge_cross(const Vec3i& x0, const Vec3i& x1, float p0, float p1)
{
   unsigned int vertex_index;
   if(edge_cross.get_entry(Vec6i(x0.v[0], x0.v[1], x0.v[2], x1.v[0], x1.v[1], x1.v[2]), vertex_index)){
      return vertex_index;
   }else if(edge_cross.get_entry(Vec6i(x1.v[0], x1.v[1], x1.v[2], x0.v[0], x0.v[1], x0.v[2]), vertex_index)){
      return vertex_index;
   }else{
      float a=p1/(p1-p0), b=1-a;
      vertex_index=(int)x.size();
      x.push_back(Vec3f(origin[0]+dx*(a*x0[0]+b*x1[0]),
                        origin[1]+dx*(a*x0[1]+b*x1[1]),
                        origin[2]+dx*(a*x0[2]+b*x1[2])));
      edge_cross.add(Vec6i(x0.v[0], x0.v[1], x0.v[2], x1.v[0], x1.v[1], x1.v[2]), vertex_index);
      return vertex_index;
   }
}

