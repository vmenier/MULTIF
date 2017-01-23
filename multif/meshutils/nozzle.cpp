#include <Python.h>

#if PY_MAJOR_VERSION >= 3
  #define MOD_ERROR_VAL NULL
  #define MOD_SUCCESS_VAL(val) val
  #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
  #define MOD_ERROR_VAL
  #define MOD_SUCCESS_VAL(val)
  #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
  #define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

static char module_docstring[] =
    "This module provides an interface for generating a rocket nozzle mesh in AERO-S format using Gmsh.";
static char generate_docstring[] =
    "Generate thermal and structural meshes for given parameters defined in files NOZZLE.txt and BOUNDARY.txt.";
static char convert_docstring[] =
    "Convert temperature output file TEMP.1 from thermal analysis to input file TEMPERATURES.txt for structural analysis.";

static PyObject *nozzle_generate(PyObject *self, PyObject *args);
static PyObject *nozzle_convert(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"generate", nozzle_generate, METH_VARARGS, generate_docstring},
    {"convert", nozzle_convert, METH_VARARGS, convert_docstring},
    {NULL, NULL, 0, NULL}
};

MOD_INIT(_nozzle_module)
{
  PyObject *m;

  MOD_DEF(m, "_nozzle_module", module_docstring, module_methods)

  if(m == NULL)
    return MOD_ERROR_VAL;

  return MOD_SUCCESS_VAL(m);
}

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <numeric>
#include <set>
#include <string>
#include <vector>
#include "OS.h"
#include "MVertex.h"
#include "MElement.h"
#include "Gmsh.h"
#include "GModel.h"
#include "GEntity.h"
//#include "Options.h"
#include "option.h"
#include "Context.h"
#include "gmshLevelset.h"

struct PointData { std::vector<double> xyz; double dydx; std::vector<double> t, tt; double ts, ws; };
struct VertexData { int p, mb; double wb; int nb; double tb; };
struct SegmentData { std::vector<int> m; int ns, ms, nn, mn, sn; std::vector<int> mt, tn; int ln; };
struct MaterialData
{
  enum { ISOTROPIC=0, ANISOTROPIC } type;
  union { double E; double E1; };
  union { double nu; double nu12; };
  double rho;
  union { double w; double w1; };
  double k;
  double h;
  double E2, G12, mu1, mu2, w2, w12;
};
struct BoundaryData { double x, P, T, Ta; };

inline bool cmp(const BoundaryData &a, const double &b) { return a.x < b; }

struct Cmp {
  const std::vector<PointData> &points;
  double tol;
  Cmp(const std::vector<PointData> &p, double t = std::numeric_limits<double>::epsilon()) : points(p), tol(t) {}
  bool operator()(const VertexData &a, const double &b) { return points[a.p].xyz[0] < (b+tol); }
};

struct Cmp2 {
  double tol;
  Cmp2(double t = std::numeric_limits<double>::epsilon()) : tol(t) {}
  bool operator()(const PointData &p, const double &b) { return p.xyz[0] < (b+tol); }
};

void writeNode(MVertex *m, FILE *fp, double scalingFactor)
{
  if(m->getIndex() < 0) return;

  double x1 = m->x() * scalingFactor;
  double y1 = m->y() * scalingFactor;
  double z1 = m->z() * scalingFactor;

  fprintf(fp, "%-8d %-16.9G %-16.9G %-16.9G\n", m->getIndex(), x1, y1, z1);
}

void writeElem(MElement *m, FILE *fp, int type)
{
  int n = m->getNumVertices();

  fprintf(fp, "%-8d %-8d", m->getNum(), type);
  for(int i = 0; i < n; i++) {
    fprintf(fp, " %-8d", m->getVertexBDF(i)->getIndex());
  }
  fprintf(fp, "\n");
}

void writeFace(MElement *m, FILE *fp, int type)
{
  int n = m->getNumVertices();

  fprintf(fp, "%-8d %-8d", m->getNum(), type);
  for(int i = 0; i < n; i++) {
    fprintf(fp, " %-8d", m->getVertexBDF(i)->getIndex());
  }
  fprintf(fp, "\n");
}

void writeAttr(MElement *m, FILE *fp, double scalingFactor, int physicalTag, const std::vector<PointData> &points,
               const std::vector<VertexData> &vertices, const std::vector<SegmentData> &segments) 
{
  // find x-coordinate of element centroid
  double x = 0;
  for(int i = 0; i < m->getNumVertices(); i++) {
    x += m->getVertexBDF(i)->x() * scalingFactor;
  }
  x /= m->getNumVertices();

  // find the material id
  Cmp cmp(points);
  std::vector<VertexData>::const_iterator it = std::lower_bound(vertices.begin(), vertices.end(), x, cmp);
  int materialId;
  switch(physicalTag) {
    case 0 : { // thermal inner layer
      materialId = 1+(segments.begin()+std::distance(vertices.begin(),it-1))->mt[0];
    } break;
    case 999 : { // thermal outer layer
      materialId = 1+(segments.begin()+std::distance(vertices.begin(),it-1))->mt[1];
    } break;
    case 1 : case 2 :  case 3 : { // load layer
      materialId = 1+(segments.begin()+std::distance(vertices.begin(),it-1))->m[0];
    } break;
  } 

  fprintf(fp, "%-8d %-8d\n", m->getNum(), materialId);
}

void writeMate(MElement *m, FILE *fp, double scalingFactor, const std::vector<BoundaryData> &boundaries,
               int physicalTag, const std::vector<MaterialData> &materials, const std::vector<PointData> &points,
               const std::vector<VertexData> &vertices, const std::vector<SegmentData> &segments)
{
  // find x-coordinate of triangle centroid
  double x = 0;
  for(int i = 0; i < 3; i++) {
    x += m->getVertexBDF(i)->x() * scalingFactor;
  }
  x /= 3;

  // interpolate ambient temperature at centroid from boundary data
  std::vector<BoundaryData>::const_iterator it = std::lower_bound(boundaries.begin(), boundaries.end(), x, cmp);
  double Taval = (it == boundaries.begin()) ? it->Ta : ((it-1)->Ta + (it->Ta - (it-1)->Ta)*(x - (it-1)->x)/(it->x - (it-1)->x));

  switch(physicalTag) {

    case 0 : { // thermal layer

    } break;

    case 1 : case 2 : case 3 : { // load layer (layered composite, material properties defined elsewhere)
      fprintf(fp, "%d 0 %g %g %g 0 0 %g 0 %g 0 %g 0 0 0\n",
              m->getNum(), 0., 0., 0., 0., Taval, 0.);
    } break;

    case 4 : { // stringers
      double tval; 
      MaterialData mat;
      Cmp cmp(points);
      std::vector<VertexData>::const_iterator it = std::lower_bound(vertices.begin(), vertices.end(), x, cmp);
      Cmp2 cmp2;
      std::vector<PointData>::const_iterator it2 = std::lower_bound(points.begin(), points.end(), x, cmp2);

      // interpolate stringer thickness from point data
      tval = (it2-1)->ts + (it2->ts - (it2-1)->ts)*(x - (it2-1)->xyz[0])/(it2->xyz[0] - (it2-1)->xyz[0]);
      mat = materials[(segments.begin()+std::distance(vertices.begin(),it-1))->ms];

      fprintf(fp, "%d 0 %g %g %g 0 0 %g 0 %g 0 %g 0 0 0\n",
              m->getNum(), mat.E, mat.nu, mat.rho, tval, Taval, mat.w);
    } break;

    default : { // baffles

      double tval;
      MaterialData mat;
      Cmp cmp(points);
      std::vector<VertexData>::const_iterator it = std::lower_bound(vertices.begin(), vertices.end(), x, cmp); // returns it such that it->p.xyz[0] >= x

      if(it == vertices.begin()) {
        tval = it->tb;
        mat = materials[it->mb]; 
      }
      else {
        tval = (it-1)->tb;
        mat = materials[(it-1)->mb];
      }

      fprintf(fp, "%d 0 %g %g %g 0 0 %g 0 %g 0 %g 0 0 0\n",
              m->getNum(), mat.E, mat.nu, mat.rho, tval, Taval, mat.w);
    } break;
  }
}

void writeComp(MElement *m, FILE *fp, double scalingFactor, const std::vector<PointData> &points,
               const std::vector<VertexData> &vertices, const std::vector<SegmentData> &segments,
               const std::vector<MaterialData> &materials)
{
  // find x-coordinate of triangle centroid
  double x = 0;
  for(int i = 0; i < 3; i++) {
    x += m->getVertexBDF(i)->x() * scalingFactor;
  } 
  x /= 3;

  fprintf(fp, "LAYC %d\n", m->getNum());
  // interpolate each layer thickness at centroid from vertex data
  double tval;
  Cmp cmp(points);
  std::vector<VertexData>::const_iterator it = std::lower_bound(vertices.begin(), vertices.end(), x, cmp);
  int segment_id = std::min<long>(segments.size()-1,std::distance(vertices.begin(), it));
  Cmp2 cmp2;
  std::vector<PointData>::const_iterator it2 = std::lower_bound(points.begin(), points.end(), x, cmp2);
  for(unsigned int k=0; k<it2->t.size(); ++k) {
    tval = (it2 == points.begin()) ? it2->t[k] 
         : ((it2-1)->t[k] + (it2->t[k] - (it2-1)->t[k])*(x - (it2-1)->xyz[0])/(it2->xyz[0] - (it2-1)->xyz[0]));
    const MaterialData &mat = materials[segments[segment_id].m[k]];
    if(mat.type == MaterialData::ISOTROPIC) {
      double G = mat.E/(2*(1+mat.nu));
      fprintf(fp, "%d %g %g %g %g %g %g %g %g %g %g %g %g\n",
              k+1, mat.E, mat.E, mat.nu, G, 0., 0., mat.rho, tval, 0., mat.w, mat.w, 0.);
    }
    else {
      fprintf(fp, "%d %g %g %g %g %g %g %g %g %g %g %g %g\n",
              k+1, mat.E1, mat.E2, mat.nu12, mat.G12, mat.mu1, mat.mu2, mat.rho, tval, 0., mat.w1, mat.w2, mat.w12);
    }
  }
}

void writeCfra(MElement *m, FILE *fp)
{
  fprintf(fp, "%d 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0\n", m->getNum());
}

void writeDisp(MVertex *m, FILE *fp, double scalingFactor)
{
  if(m->getIndex() < 0) return;

  fprintf(fp, "%-8d 1 0.0\n", m->getIndex());
  fprintf(fp, "%-8d 2 0.0\n", m->getIndex());
  fprintf(fp, "%-8d 3 0.0\n", m->getIndex());
  fprintf(fp, "%-8d 4 0.0\n", m->getIndex());
  fprintf(fp, "%-8d 5 0.0\n", m->getIndex());
  fprintf(fp, "%-8d 6 0.0\n", m->getIndex());
}

void writePres(MElement *m, FILE *fp, double scalingFactor, const std::vector<BoundaryData> &boundaries)
{
  // find x-coordinate of triangle centroid
  double x = 0;
  for(int i = 0; i < 3; i++) {
    x += m->getVertexBDF(i)->x() * scalingFactor;
  }
  x /= 3;

  // interpolate pressure at centroid from boundary data
  std::vector<BoundaryData>::const_iterator it = std::lower_bound(boundaries.begin(), boundaries.end(), x, cmp);
  double pval = (it == boundaries.begin()) ? it->P : ((it-1)->P + (it->P - (it-1)->P)*(x - (it-1)->x)/(it->x - (it-1)->x));

  fprintf(fp, "%-8d %-16.9G\n", m->getNum(), pval);
}

void writeTemp(MVertex *m, FILE *fp, double scalingFactor, const std::vector<BoundaryData> &boundaries)
{
  if(m->getIndex() < 0) return;

  double x = m->x() * scalingFactor;

  // interpolate temperature at node from boundary data
  std::vector<BoundaryData>::const_iterator it = std::lower_bound(boundaries.begin(), boundaries.end(), x, cmp);
  double tval = (it == boundaries.begin()) ? it->T : ((it-1)->T + (it->T - (it-1)->T)*(x - (it-1)->x)/(it->x - (it-1)->x));

  fprintf(fp, "%-8d %-16.9G\n", m->getIndex(), tval);
}

void writeConvec(MElement *m, FILE *fp, double scalingFactor, double h, const std::vector<BoundaryData> &boundaries)
{
  // compute the nodal area
  int pOrder = 2;
  double A[4], s[4]; 

  int npts; IntPt *gp;
  m->getIntegrationPoints(pOrder, &npts, &gp);
  for (int i = 0; i < 4; ++i) A[i] = 0;
  for (int i = 0; i < npts; i++) {
    double u = gp[i].pt[0];
    double v = gp[i].pt[1];
    double w = gp[i].pt[2];
    double weight = gp[i].weight;
    double detuvw = m->getJacobianDeterminant(u, v, w);
    m->getShapeFunctions(u, v, w, s);
    for (int j = 0; j < 4; ++j) A[j] += weight*detuvw*s[j];
  }

  for(int i = 0; i < 4; ++i) {
    // interpolate ambient temperature at node from boundary data
    double x = m->getVertexBDF(i)->x();
    std::vector<BoundaryData>::const_iterator it = std::lower_bound(boundaries.begin(), boundaries.end(), x, cmp);
    double Ta = (it == boundaries.begin()) ? it->Ta : ((it-1)->Ta + (it->Ta - (it-1)->Ta)*(x - (it-1)->x)/(it->x - (it-1)->x));

    fprintf(fp, "%-8d %-16.9G %-16.9G %-16.9G\n", m->getVertexBDF(i)->getIndex(), h, A[i], Ta);
  }
}

int writeAEROS(GModel *g,
               const std::vector<PointData> &points,
               const std::vector<VertexData> &vertices,
               const std::vector<SegmentData> &segments,
               const std::vector<MaterialData> &materials,
               const std::vector<BoundaryData> &boundaries,
               const std::vector<std::pair<GEntity::GeomType,int> > &surfaceTags,
               const std::vector<std::pair<GEntity::GeomType,int> > &boundaryTags,
               int lf, const std::string &name,
               int elementTagType=1, bool saveAll=false, double scalingFactor=1.0, bool temp=false)
{
  FILE *fp = Fopen(name.c_str(), "w");
  if(!fp){
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  if(g->noPhysicalGroups()) saveAll = true;

  g->indexMeshVertices(saveAll);

  std::vector<GEntity*> entities;
  g->getEntities(entities);

  // nodes (vertices of all the triangles in the model)
  std::set<int> s;
  FILE *fp1 = Fopen("GEOMETRY.txt", "w");
  fprintf(fp1, "NODES\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      MElement *m = entities[i]->getMeshElement(j);
      const char *str = m->getStringForBDF();
      if(str && std::strcmp(str,"CTRIA3") == 0) {
        for(int k = 0; k < m->getNumVertices(); ++k) {
          int index = m->getVertexBDF(k)->getIndex();
          if(s.find(index) == s.end()) {
            writeNode(m->getVertexBDF(k), fp1, scalingFactor);
            s.insert(index);
          }
        }
      }
    }
  fclose(fp1);

  // elements (all the triangles in the model)
  FILE *fp2 = Fopen("TOPOLOGY.txt", "w");
  fprintf(fp2, "TOPOLOGY\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
      if(str && std::strcmp(str,"CTRIA3") == 0) {
        writeElem(entities[i]->getMeshElement(j), fp2, 15);
      }
    }
  fclose(fp2);

  // attributes (each element has its own material to enable the thickness and ambient temperature to vary)
  FILE *fp3 = Fopen("ATTRIBUTES.txt", "w");
  fprintf(fp3, "ATTRIBUTES\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    bool b = std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 2) != entities[i]->physicals.end();
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      MElement *m = entities[i]->getMeshElement(j);
      const char *str = m->getStringForBDF();
      if(str && std::strcmp(str,"CTRIA3") == 0) {
        if(b) fprintf(fp3, "%-8d %-8d %-8d %-8d\n", m->getNum(), m->getNum(), m->getNum(), m->getNum());
        else fprintf(fp3, "%-8d %-8d\n", m->getNum(), m->getNum());
      }
    }
  }
  fclose(fp3);

  // material properties
  FILE *fp8 = Fopen("MATERIAL.txt", "w");
  fprintf(fp8, "MATERIALS\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    int physicalTag = entities[i]->physicals.empty() ? -1 : entities[i]->physicals[0];
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      MElement *m = entities[i]->getMeshElement(j);
      const char *str = m->getStringForBDF();
      if(str && std::strcmp(str,"CTRIA3") == 0) {
        writeMate(m, fp8, scalingFactor, boundaries, physicalTag, materials, points, vertices, segments);
      }
    }
  }
  fclose(fp8);

  // composite material properties
  FILE *fp9 = Fopen("COMPOSITE.txt", "w");
  fprintf(fp9, "COMPOSITE\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 2) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        MElement *m = entities[i]->getMeshElement(j);
        const char *str = m->getStringForBDF();
        if(str && std::strcmp(str,"CTRIA3") == 0) {
          writeComp(m, fp9, scalingFactor, points, vertices, segments, materials);
        }
      }
    }
  }
  fclose(fp9);

  // composite frames
  FILE *fp10 = Fopen("CFRAMES.txt", "w");
  fprintf(fp10, "CFRAMES\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 2) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        MElement *m = entities[i]->getMeshElement(j);
        const char *str = m->getStringForBDF();
        if(str && std::strcmp(str,"CTRIA3") == 0) {
          writeCfra(m, fp10);
        }
      }
    }
  }
  fclose(fp10);

  // displacement (all nodes on the boundary edges)
  FILE *fp4 = Fopen("DISPLACEMENTS.txt", "w");
  fprintf(fp4, "DISPLACEMENTS\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(boundaryTags.begin(), boundaryTags.end(), std::make_pair(entities[i]->geomType(),entities[i]->tag())) != boundaryTags.end()) {
      for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
        writeDisp(entities[i]->mesh_vertices[j], fp4, scalingFactor);
    }
  }
  fclose(fp4);

  // pressure (all the triangles on the mid-surface of the load layer)
  FILE *fp5 = Fopen("PRESSURES.txt", "w");
  fprintf(fp5, "PRESSURE\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 2) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
        if(str && std::strcmp(str,"CTRIA3") == 0)
          writePres(entities[i]->getMeshElement(j), fp5, scalingFactor, boundaries);
      }
    }
  }
  fclose(fp5);

  // temperatures (all the nodes on the mid-surface of the load layer)
  if(temp) {
    FILE *fp6 = Fopen("TEMPERATURES.txt", "w");
    fprintf(fp6, "TEMPERATURE\n");
    for(unsigned int i = 0; i < entities.size(); i++) {
      if(std::find(surfaceTags.begin(), surfaceTags.end(), std::make_pair(entities[i]->geomType(),entities[i]->tag())) != surfaceTags.end()) {
        for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
          writeTemp(entities[i]->mesh_vertices[j], fp6, scalingFactor, boundaries);
      }
    }
    fclose(fp6);
  }

  // surface topology (all the triangles on the mid-surface of the load layer)
  FILE *fp7 = Fopen("SURFACETOPO.txt", "w");
  fprintf(fp7, "SURFACETOPO 2\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 2) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
        if(str && std::strcmp(str,"CTRIA3") == 0)
          writeFace(entities[i]->getMeshElement(j), fp7, 3);
      }
    }
  }
  // surface topology (all the triangles on the mid-surface of the stringers)
  int stringerCount = 0;
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 4) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
        if(str && std::strcmp(str,"CTRIA3") == 0) {
          if(stringerCount == 0) fprintf(fp7, "SURFACETOPO 4\n");
          writeFace(entities[i]->getMeshElement(j), fp7, 3);
          stringerCount++;
        }
      }
    }
  }
  // surface topology (all the triangles on the mid-surface of the baffles)
  int baffleCount = 0;
  for(std::vector<VertexData>::const_iterator it = vertices.begin(); it != vertices.end(); ++it) if(it->wb > 0) baffleCount++;
  for(int k = 0; k < baffleCount; ++k) {
    fprintf(fp7, "SURFACETOPO %d\n", 5+k);
    for(unsigned int i = 0; i < entities.size(); i++) {
      if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 5+k) != entities[i]->physicals.end()) {
        for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
          const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
          if(str && std::strcmp(str,"CTRIA3") == 0)
            writeFace(entities[i]->getMeshElement(j), fp7, 3);
        }
      }
    }
  }
  fclose(fp7);

  // main file
  fprintf(fp, "STATICS\n");
  fprintf(fp, "sparse\n");
  fprintf(fp, "*\n");
  if(lf == 0) {
    fprintf(fp, "NONLINEAR\n");
    fprintf(fp, "maxit 15\n");
    fprintf(fp, "*\n");
  }
  fprintf(fp, "OUTPUT\n");
  fprintf(fp, "gdisplac 14 7 \"DISP\" 1\n");

  // von mises stress
  fprintf(fp, "stressvm 14 7 \"STRESS\" 1 lower\n");
  fprintf(fp, "stressvm 14 7 \"STRESS.1\" 1 NG 2 lower nodalpartialgroup\n");
  fprintf(fp, "stressvm 14 7 \"STRESS.2\" 1 NG 2 median nodalpartialgroup\n");
  fprintf(fp, "stressvm 14 7 \"STRESS.3\" 1 NG 2 upper nodalpartialgroup\n");
  if(stringerCount > 0)
    fprintf(fp, "stressvm 14 7 \"STRESS.4\" 1 NG 4 median nodalpartialgroup\n");
  for(int k = 0; k < baffleCount; ++k) 
    fprintf(fp, "stressvm 14 7 \"STRESS.%d\" 1 NG %d median nodalpartialgroup\n", 5+k, 5+k);
  if(lf == 1) {
    fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS\" 1 lower thermal\n");
    fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS.1\" 1 NG 2 lower nodalpartialgroup thermal\n");
    fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS.2\" 1 NG 2 median nodalpartialgroup thermal\n");
    fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS.3\" 1 NG 2 upper nodalpartialgroup thermal\n");
    if(stringerCount > 0)
      fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS.4\" 1 NG 4 median nodalpartialgroup thermal\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS.%d\" 1 NG %d median nodalpartialgroup thermal\n", 5+k, 5+k);
    fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS\" 1 lower mechanical\n");
    fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS.1\" 1 NG 2 lower nodalpartialgroup mechanical\n");
    fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS.2\" 1 NG 2 median nodalpartialgroup mechanical\n");
    fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS.3\" 1 NG 2 upper nodalpartialgroup mechanical\n");
    if(stringerCount > 0)
      fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS.4\" 1 NG 4 median nodalpartialgroup mechanical\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS.%d\" 1 NG %d median nodalpartialgroup mechanical\n", 5+k, 5+k);
  }
  // 1st principal stress
  fprintf(fp, "stressp1 14 7 \"STRESSP1\" 1 lower\n");
  fprintf(fp, "stressp1 14 7 \"STRESSP1.1\" 1 NG 2 lower nodalpartialgroup\n");
  fprintf(fp, "stressp1 14 7 \"STRESSP1.2\" 1 NG 2 median nodalpartialgroup\n");
  fprintf(fp, "stressp1 14 7 \"STRESSP1.3\" 1 NG 2 upper nodalpartialgroup\n");
  if(stringerCount > 0)
    fprintf(fp, "stressp1 14 7 \"STRESSP1.4\" 1 NG 4 median nodalpartialgroup\n");
  for(int k = 0; k < baffleCount; ++k)
    fprintf(fp, "stressp1 14 7 \"STRESSP1.%d\" 1 NG %d median nodalpartialgroup\n", 5+k, 5+k);
  if(lf == 1) {
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1\" 1 lower thermal\n");
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.1\" 1 NG 2 lower nodalpartialgroup thermal\n");
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.2\" 1 NG 2 median nodalpartialgroup thermal\n");
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.3\" 1 NG 2 upper nodalpartialgroup thermal\n");
    if(stringerCount > 0)
      fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.4\" 1 NG 4 median nodalpartialgroup thermal\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.%d\" 1 NG %d median nodalpartialgroup thermal\n", 5+k, 5+k);
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1\" 1 lower mechanical\n");
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.1\" 1 NG 2 lower nodalpartialgroup mechanical\n");
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.2\" 1 NG 2 median nodalpartialgroup mechanical\n");
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.3\" 1 NG 2 upper nodalpartialgroup mechanical\n");
    if(stringerCount > 0)
      fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.4\" 1 NG 4 median nodalpartialgroup mechanical\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.%d\" 1 NG %d median nodalpartialgroup mechanical\n", 5+k, 5+k);
  }
  // 2nd principal stress
  fprintf(fp, "stressp2 14 7 \"STRESSP2\" 1 lower\n");
  fprintf(fp, "stressp2 14 7 \"STRESSP2.1\" 1 NG 2 lower nodalpartialgroup\n");
  fprintf(fp, "stressp2 14 7 \"STRESSP2.2\" 1 NG 2 median nodalpartialgroup\n");
  fprintf(fp, "stressp2 14 7 \"STRESSP2.3\" 1 NG 2 upper nodalpartialgroup\n");
  if(stringerCount > 0)
    fprintf(fp, "stressp2 14 7 \"STRESSP2.4\" 1 NG 4 median nodalpartialgroup\n");
  for(int k = 0; k < baffleCount; ++k)
    fprintf(fp, "stressp2 14 7 \"STRESSP2.%d\" 1 NG %d median nodalpartialgroup\n", 5+k, 5+k);
  if(lf == 1) {
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2\" 1 lower thermal\n");
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.1\" 1 NG 2 lower nodalpartialgroup thermal\n");
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.2\" 1 NG 2 median nodalpartialgroup thermal\n");
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.3\" 1 NG 2 upper nodalpartialgroup thermal\n");
    if(stringerCount > 0)
      fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.4\" 1 NG 4 median nodalpartialgroup thermal\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.%d\" 1 NG %d median nodalpartialgroup thermal\n", 5+k, 5+k);
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2\" 1 lower mechanical\n");
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.1\" 1 NG 2 lower nodalpartialgroup mechanical\n");
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.2\" 1 NG 2 median nodalpartialgroup mechanical\n");
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.3\" 1 NG 2 upper nodalpartialgroup mechanical\n");
    if(stringerCount > 0)
      fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.4\" 1 NG 4 median nodalpartialgroup mechanical\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.%d\" 1 NG %d median nodalpartialgroup mechanical\n", 5+k, 5+k);
  }
  // 3rd principal stress
  fprintf(fp, "stressp3 14 7 \"STRESSP3\" 1 lower\n");
  fprintf(fp, "stressp3 14 7 \"STRESSP3.1\" 1 NG 2 lower nodalpartialgroup\n");
  fprintf(fp, "stressp3 14 7 \"STRESSP3.2\" 1 NG 2 median nodalpartialgroup\n");
  fprintf(fp, "stressp3 14 7 \"STRESSP3.3\" 1 NG 2 upper nodalpartialgroup\n");
  if(stringerCount > 0)
    fprintf(fp, "stressp3 14 7 \"STRESSP3.4\" 1 NG 4 median nodalpartialgroup\n");
  for(int k = 0; k < baffleCount; ++k)
    fprintf(fp, "stressp3 14 7 \"STRESSP3.%d\" 1 NG %d median nodalpartialgroup\n", 5+k, 5+k);
  if(lf == 1) {
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3\" 1 lower thermal\n");
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.1\" 1 NG 2 lower nodalpartialgroup thermal\n");
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.2\" 1 NG 2 median nodalpartialgroup thermal\n");
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.3\" 1 NG 2 upper nodalpartialgroup thermal\n");
    if(stringerCount > 0)
      fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.4\" 1 NG 4 median nodalpartialgroup thermal\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.%d\" 1 NG %d median nodalpartialgroup thermal\n", 5+k, 5+k);
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3\" 1 lower mechanical\n");
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.1\" 1 NG 2 lower nodalpartialgroup mechanical\n");
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.2\" 1 NG 2 median nodalpartialgroup mechanical\n");
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.3\" 1 NG 2 upper nodalpartialgroup mechanical\n");
    if(stringerCount > 0)
      fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.4\" 1 NG 4 median nodalpartialgroup mechanical\n");
    for(int k = 0; k < baffleCount; ++k)
      fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.%d\" 1 NG %d median nodalpartialgroup mechanical\n", 5+k, 5+k);
  }
  // stresses and strains in the x and y directions of the material coordinate frame (layers 1 and 3 only)
  fprintf(fp, "stressxx 14 7 \"STRESSXX.1\" 1 NG 2 lower nodalpartialgroup matfrm\n");
  fprintf(fp, "stressxx 14 7 \"STRESSXX.3\" 1 NG 2 upper nodalpartialgroup matfrm\n");
  fprintf(fp, "stressyy 14 7 \"STRESSYY.1\" 1 NG 2 lower nodalpartialgroup matfrm\n");
  fprintf(fp, "stressyy 14 7 \"STRESSYY.3\" 1 NG 2 upper nodalpartialgroup matfrm\n");
  fprintf(fp, "strainxx 14 7 \"STRAINXX.1\" 1 NG 2 lower nodalpartialgroup matfrm\n");
  fprintf(fp, "strainxx 14 7 \"STRAINXX.3\" 1 NG 2 upper nodalpartialgroup matfrm\n");
  fprintf(fp, "strainyy 14 7 \"STRAINYY.1\" 1 NG 2 lower nodalpartialgroup matfrm\n");
  fprintf(fp, "strainyy 14 7 \"STRAINYY.3\" 1 NG 2 upper nodalpartialgroup matfrm\n");
  if(stringerCount > 0) {
    fprintf(fp, "stressxx 14 7 \"STRESSXX.4\" 1 NG 4 median nodalpartialgroup matfrm\n");
    fprintf(fp, "stressyy 14 7 \"STRESSYY.4\" 1 NG 4 median nodalpartialgroup matfrm\n");
    fprintf(fp, "strainxx 14 7 \"STRAINXX.4\" 1 NG 4 median nodalpartialgroup matfrm\n");
    fprintf(fp, "strainyy 14 7 \"STRAINYY.4\" 1 NG 4 median nodalpartialgroup matfrm\n");
  }

  fprintf(fp, "*\n");
  fprintf(fp, "GROUPS\n");
  fprintf(fp, "N surface 2 2\n");
  if(stringerCount > 0)
    fprintf(fp, "N surface 4 4\n");
  for(int k = 0; k < baffleCount; ++k)
    fprintf(fp, "N surface %d %d\n", 5+k, 5+k);
  fprintf(fp, "*\n");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "GEOMETRY.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TOPOLOGY.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "ATTRIBUTES.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "MATERIAL.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "COMPOSITE.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "CFRAMES.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "DISPLACEMENTS.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "PRESSURES.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TEMPERATURES.txt"); // this file is generated using thermal analysis if temp is false
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "SURFACETOPO.txt");
  fprintf(fp, "END\n");

  fclose(fp);

  return 1;
}

int writeAEROH(GModel *g,
               const std::vector<PointData> &points,
               const std::vector<VertexData> &vertices,
               const std::vector<SegmentData> &segments,
               const std::vector<MaterialData> &materials,
               const std::vector<BoundaryData> &boundaries,
               const std::vector<std::pair<GEntity::GeomType,int> > &interiorBoundaryTags,
               const std::vector<std::pair<GEntity::GeomType,int> > &exteriorBoundaryTags,
               const std::string &name,
               int elementTagType=1, bool saveAll=false, double scalingFactor=1.0)
{
  FILE *fp = Fopen(name.c_str(), "w");
  if(!fp){
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  if(g->noPhysicalGroups()) saveAll = true;

  g->indexMeshVertices(saveAll);

  std::vector<GEntity*> entities;
  g->getEntities(entities);

  // nodes
  FILE *fp1 = Fopen("GEOMETRY.txt.thermal", "w");
  fprintf(fp1, "NODES\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
      writeNode(entities[i]->mesh_vertices[j], fp1, scalingFactor);
  fclose(fp1);

  // elements (all hexas and exterior boundary quads)
  FILE *fp2 = Fopen("TOPOLOGY.txt.thermal", "w");
  fprintf(fp2, "TOPOLOGY\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    bool b = (std::find(exteriorBoundaryTags.begin(), exteriorBoundaryTags.end(),
              std::make_pair(entities[i]->geomType(),entities[i]->tag())) != exteriorBoundaryTags.end());
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
      if(str && std::strcmp(str,"CHEXA") == 0) 
        writeElem(entities[i]->getMeshElement(j), fp2, 51);
      else if(str && b && std::strcmp(str,"CQUAD4") == 0)
        writeElem(entities[i]->getMeshElement(j), fp2, 48); 
    }
  }
  fclose(fp2);

  // attributes
  FILE *fp3 = Fopen("ATTRIBUTES.txt.thermal", "w");
  fprintf(fp3, "ATTRIBUTES\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    bool b = (std::find(exteriorBoundaryTags.begin(), exteriorBoundaryTags.end(),
              std::make_pair(entities[i]->geomType(),entities[i]->tag())) != exteriorBoundaryTags.end());
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
      if(str && (std::strcmp(str,"CHEXA") == 0 || (b && std::strcmp(str,"CQUAD4") == 0))) {
        writeAttr(entities[i]->getMeshElement(j), fp3, scalingFactor, entities[i]->physicals[0], points, vertices,
                  segments);
      }
    }
  }
  fclose(fp3);

  // material properties
  std::ofstream fout("MATERIAL.txt.thermal");
  fout << "MATERIALS\n";
  for(std::vector<MaterialData>::const_iterator it = materials.begin(); it != materials.end(); ++it) {
    fout << std::distance(materials.begin(),it)+1 << " 0 0 0 " << it->rho
         << " " << it->h << " " << it->k << " 0 0 0 0 0 0 0 0\n";
  }
  fout.close();

  // temperatures (all the nodes on the interior boundary)
  FILE *fp4 = Fopen("TEMPERATURES.txt.thermal", "w");
  fprintf(fp4, "TEMPERATURE\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(interiorBoundaryTags.begin(), interiorBoundaryTags.end(),
                 std::make_pair(entities[i]->geomType(),entities[i]->tag())) != interiorBoundaryTags.end()) {
      for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
        writeTemp(entities[i]->mesh_vertices[j], fp4, scalingFactor, boundaries);
    }
  }
  fclose(fp4);

  // convection (all the nodes on the exterior boundary)
  FILE *fp6 = Fopen("CONVECTION.txt.thermal", "w");
  fprintf(fp6, "CONVECTION\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(exteriorBoundaryTags.begin(), exteriorBoundaryTags.end(),
                 std::make_pair(entities[i]->geomType(),entities[i]->tag())) != exteriorBoundaryTags.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
        if(str && std::strcmp(str,"CQUAD4") == 0) {
          int m = entities[i]->physicals[0];
          double h = materials[m-1].h;
          writeConvec(entities[i]->getMeshElement(j), fp6, scalingFactor, h, boundaries);
        }
      }
    }
  }
  fclose(fp6);

  // surface topology (all the triangles on the mid-surface of the load layer, used to ouput temperatures for structural model
  //                   also, all the quads on the upper and lower surfaces of the load layer)
  FILE *fp7 = Fopen("SURFACETOPO.txt.thermal", "w");

  fprintf(fp7, "SURFACETOPO 1\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 1) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
        if(str && std::strcmp(str,"CQUAD4") == 0)
          writeFace(entities[i]->getMeshElement(j), fp7, 1);
      }
    }
  }
  fprintf(fp7, "SURFACETOPO 2\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 2) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
        if(str && std::strcmp(str,"CTRIA3") == 0)
          writeFace(entities[i]->getMeshElement(j), fp7, 3);
      }
    }
  }
  fprintf(fp7, "SURFACETOPO 3\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(entities[i]->physicals.begin(), entities[i]->physicals.end(), 3) != entities[i]->physicals.end()) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
        if(str && std::strcmp(str,"CQUAD4") == 0)
          writeFace(entities[i]->getMeshElement(j), fp7, 1);
      }
    }
  }

  fclose(fp7);

  // set of nodes used to define group 4
  std::set<int> s;
  for(unsigned int i = 0; i < entities.size(); i++) {
    int physicalTag = entities[i]->physicals.empty() ? -1 : entities[i]->physicals[0];
    if(physicalTag == 0) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        MElement *m = entities[i]->getMeshElement(j);
        const char *str = m->getStringForBDF();
        if(str && std::strcmp(str,"CHEXA") == 0) {
          for(int k = 0; k < m->getNumVertices(); ++k) {
            int index = m->getVertexBDF(k)->getIndex();
            if(s.find(index) == s.end()) {
              s.insert(index);
            }
          }
        }
      }
    }
  }

  // main file
  fprintf(fp, "STATICS\n");
  fprintf(fp, "sparse\n");
  fprintf(fp, "*\n");
  fprintf(fp, "OUTPUT\n");
  fprintf(fp, "gtempera 22 15 \"TEMP\" 1\n");
  fprintf(fp, "gtempera 22 15 \"TEMP.0\" 1 NG 4\n");
  fprintf(fp, "gtempera 22 15 \"TEMP.1\" 1 NG 1\n");
  fprintf(fp, "gtempera 22 15 \"TEMP.2\" 1 NG 2\n");
  fprintf(fp, "gtempera 22 15 \"TEMP.3\" 1 NG 3\n");
  fprintf(fp, "*\n");
  fprintf(fp, "GROUPS\n");
  fprintf(fp, "N surface 1 1\n");
  fprintf(fp, "N surface 2 2\n");
  fprintf(fp, "N surface 3 3\n");
  for(std::set<int>::iterator it = s.begin(); it != s.end(); ++it) fprintf(fp, "N %d 4\n", *it);
  fprintf(fp, "*\n");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "GEOMETRY.txt.thermal");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TOPOLOGY.txt.thermal");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "ATTRIBUTES.txt.thermal");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "MATERIAL.txt.thermal");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TEMPERATURES.txt.thermal");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "CONVECTION.txt.thermal");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "SURFACETOPO.txt.thermal");
  fprintf(fp, "END\n");

  fclose(fp);

  return 1;
}

int writeAEROS2(GModel *g,
                const std::vector<PointData> &points,
                const std::vector<VertexData> &vertices,
                const std::vector<SegmentData> &segments,
                const std::vector<MaterialData> &materials,
                const std::vector<BoundaryData> &boundaries,
                const std::vector<std::pair<GEntity::GeomType,int> > &boundaryTags,
                int lf, const std::string &name,
                int elementTagType=1, bool saveAll=false, double scalingFactor=1.0)
{
  FILE *fp = Fopen(name.c_str(), "w");
  if(!fp){
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  if(g->noPhysicalGroups()) saveAll = true;

  g->indexMeshVertices(saveAll);

  std::vector<GEntity*> entities;
  g->getEntities(entities);

  // nodes
  std::set<int> s;
  FILE *fp1 = Fopen("GEOMETRY.txt.cmc", "w");
  fprintf(fp1, "NODES\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    int physicalTag = entities[i]->physicals.empty() ? -1 : entities[i]->physicals[0];
    if(physicalTag == 0) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        MElement *m = entities[i]->getMeshElement(j);
        const char *str = m->getStringForBDF();
        if(str && std::strcmp(str,"CHEXA") == 0) {
          for(int k = 0; k < m->getNumVertices(); ++k) {
            int index = m->getVertexBDF(k)->getIndex();
            if(s.find(index) == s.end()) {
              writeNode(m->getVertexBDF(k), fp1, scalingFactor);
              s.insert(index);
            }
          }
        }
      }
    }
  }
  fclose(fp1);

  // elements (all cmc hexas)
  FILE *fp2 = Fopen("TOPOLOGY.txt.cmc", "w");
  fprintf(fp2, "TOPOLOGY\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    int physicalTag = entities[i]->physicals.empty() ? -1 : entities[i]->physicals[0];
    if(physicalTag == 0) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        MElement *m = entities[i]->getMeshElement(j);
        const char *str = m->getStringForBDF();
        if(str && std::strcmp(str,"CHEXA") == 0) {
          writeElem(m, fp2, 17);
        }
      }
    }
  }
  fclose(fp2);

  // attributes
  FILE *fp3 = Fopen("ATTRIBUTES.txt.cmc", "w");
  fprintf(fp3, "ATTRIBUTES\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    int physicalTag = entities[i]->physicals.empty() ? -1 : entities[i]->physicals[0];
    if(physicalTag == 0) {
      for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
        MElement *m = entities[i]->getMeshElement(j);
        const char *str = m->getStringForBDF();
        if(str && std::strcmp(str,"CHEXA") == 0) {
          writeAttr(m, fp3, scalingFactor, physicalTag, points, vertices,
                    segments);
        }
      }
    }
  }
  fclose(fp3);

  // material properties
  std::ofstream fout("MATERIAL.txt.cmc");
  fout << "MATERIALS\n";
  for(std::vector<MaterialData>::const_iterator it = materials.begin(); it != materials.end(); ++it) {
    fout << std::distance(materials.begin(),it)+1 << " 0 " << it->E << " " << it->nu << " " << it->rho
         << " 0 0 0 0 " << boundaries.begin()->Ta << " 0 " << it->w << " 0 0 0\n"; // XXX
  }
  fout.close();

  // displacement (all nodes on the boundary edges)
  FILE *fp4 = Fopen("DISPLACEMENTS.txt.cmc", "w");
  fprintf(fp4, "DISPLACEMENTS\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    if(std::find(boundaryTags.begin(), boundaryTags.end(), std::make_pair(entities[i]->geomType(),entities[i]->tag())) != boundaryTags.end()) {
      for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
        writeDisp(entities[i]->mesh_vertices[j], fp4, scalingFactor);
    }
  }
  fclose(fp4);

  // main file
  fprintf(fp, "STATICS\n");
  fprintf(fp, "sparse\n");
  fprintf(fp, "*\n");
  if(lf == 0) {
    fprintf(fp, "NONLINEAR\n");
    fprintf(fp, "maxit 15\n");
    fprintf(fp, "*\n");
  }
  fprintf(fp, "OUTPUT\n");
  fprintf(fp, "gdisplac 14 7 \"DISP.cmc\" 1\n");

  // von mises stress
  fprintf(fp, "stressvm 14 7 \"STRESS.cmc\" 1\n");
  fprintf(fp, "stressvm 14 7 \"STRESS.0\" 1 NG 1\n");
  if(lf == 1) {
    fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS.cmc\" 1 thermal\n");
    fprintf(fp, "stressvm 14 7 \"THERMAL_STRESS.0\" 1 NG 1 thermal\n");
    fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS.cmc\" 1 mechanical\n");
    fprintf(fp, "stressvm 14 7 \"MECHANICAL_STRESS.0\" 1 NG 1 mechanical\n");
  }
  // 1st principal stress
  fprintf(fp, "stressp1 14 7 \"STRESSP1.cmc\" 1 lower\n");
  fprintf(fp, "stressp1 14 7 \"STRESSP1.0\" 1 NG 1 lower\n");
  if(lf == 1) {
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.cmc\" 1 thermal\n");
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.0\" 1 NG 1 thermal\n");
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.cmc\" 1 mechanical\n");
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.0\" 1 NG 1 mechanical\n");
  }
  // 2nd principal stress
  fprintf(fp, "stressp2 14 7 \"STRESSP2.cmc\" 1 lower\n");
  fprintf(fp, "stressp2 14 7 \"STRESSP2.0\" 1 NG 1 lower\n");
  if(lf == 1) {
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.cmc\" 1 thermal\n");
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.0\" 1 NG 1 thermal\n");
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.cmc\" 1 mechanical\n");
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.0\" 1 NG 1 mechanical\n");
  }
  // 3rd principal stress
  fprintf(fp, "stressp3 14 7 \"STRESSP3.cmc\" 1 lower\n");
  fprintf(fp, "stressp3 14 7 \"STRESSP3.0\" 1 NG 1 lower\n");
  if(lf == 1) {
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.cmc\" 1 thermal\n");
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.0\" 1 NG 1 thermal\n");
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.cmc\" 1 mechanical\n");
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.0\" 1 NG 1 mechanical\n");
  }

  fprintf(fp, "*\n");
  fprintf(fp, "GROUPS\n");
  for(std::set<int>::iterator it = s.begin(); it != s.end(); ++it) fprintf(fp, "N %d 1\n", *it);
  fprintf(fp, "*\n");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "GEOMETRY.txt.cmc");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TOPOLOGY.txt.cmc");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "ATTRIBUTES.txt.cmc");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "MATERIAL.txt.cmc");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "DISPLACEMENTS.txt.cmc");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TEMPERATURES.txt.cmc");
  
  fprintf(fp, "END\n");

  fclose(fp);

  return 1;
}

void generateNozzle(const std::vector<PointData> &points,
                    const std::vector<VertexData> &vertices,
                    const std::vector<SegmentData> &segments,
                    const std::vector<MaterialData> &materials,
                    const std::vector<BoundaryData> &boundaries,
                    double lc, int bf, int tf, int lf)
{
  // bf = 0: constrain inner wall inlet edge
  // bf = 1: constrain baffle outer edges
  // bf = 2: constrain both inner wall inlet and baffle outer edges

  // tf = 0: generate mesh for structural model only
  // tf = 1: generate mesh for both structural and thermal models

  // lf = 0: nonlinear structural model
  // lf = 1: linear structural model

  GmshInitialize();
  GModel *m = new GModel();
  m->setFactory("OpenCASCADE");

  GmshSetOption("General","Terminal", 1.);
  GmshSetOption("Mesh","CharacteristicLengthMin", lc);
  GmshSetOption("Mesh","CharacteristicLengthMax", lc);

  std::vector<double> p1(3), p2(3); // points defining the axis of revolution
  p1[0] = p1[1] = p1[2] = p2[1] = p2[2] = 0; p2[0] = 1;

  std::vector<PointData>::const_iterator pointIt1 = points.begin();
  std::vector<VertexData>::const_iterator vertexIt1 = vertices.begin(), vertexIt2 = vertices.begin()+1;
  GVertex *vertex1 = m->addVertex(pointIt1->xyz[0], pointIt1->xyz[1], pointIt1->xyz[2], lc);

  std::vector<std::pair<GEntity::GeomType,int> > exteriorBoundaryTags;
  std::vector<std::pair<GEntity::GeomType,int> > interiorBoundaryTags;
  std::vector<std::pair<GEntity::GeomType,int> > surfaceTags;
  std::vector<std::pair<GEntity::GeomType,int> > boundaryTags;
  std::vector<std::pair<GEntity::GeomType,int> > cmcBoundaryTags;

  double cth1, cth2;

  // count the number of baffles
  int baffleCount = 0;
  for(std::vector<VertexData>::const_iterator it = vertices.begin(); it != vertices.end(); ++it)
    if(it->wb > 0) baffleCount++;

  // loop over the segments
  int baffleIndex = 0;
  for(std::vector<SegmentData>::const_iterator segmentIt = segments.begin(); segmentIt != segments.end(); ++segmentIt, ++vertexIt2) {

    std::vector<PointData>::const_iterator pointIt2 = points.begin()+vertexIt2->p;
    GVertex *vertex2 = m->addVertex(pointIt2->xyz[0], pointIt2->xyz[1], pointIt2->xyz[2], lc);
    std::vector<std::vector<double> > controlPoints;
    std::vector<PointData>::const_iterator pointIt;
    for(pointIt = pointIt1+1; pointIt != pointIt2; ++pointIt) controlPoints.push_back(pointIt->xyz);
    // inside of inner layer
    GEdge *edge1;
    try {
      edge1 = controlPoints.empty() ? m->addLine(vertex1, vertex2) :  m->addBSpline(vertex1, vertex2, controlPoints); 
    }
    catch(...) {
      std::cerr << "Warning : mesh generation failed in segment " << std::distance(segments.begin(),segmentIt)
                << ", x = [" << vertex1->x() << "," << vertex2->x() << "]\n";
      pointIt1 = pointIt2; vertexIt1 = vertexIt2; vertex1 = vertex2;
      continue;
    } 

    int ns = std::max(1,segmentIt->ns); // number of circumferential segments
    double wb1 = vertexIt1->wb;         // width of baffle
    double wb2 = vertexIt2->wb;         // width of baffle
    int nb1 = vertexIt1->nb;            // number of transfinite points (radial edges of baffle)
    int nb2 = vertexIt2->nb;            // number of transfinite points (radial edges of baffle)
    double ws1 = pointIt1->ws;          // width of stiffeners
    double ws2 = pointIt2->ws;          // width of stiffeners
    int nn = segmentIt->nn;             // number of transfinite points (longitudinal edges)
    int mn = segmentIt->mn;             // number of transfinite points (circumferential edges)
    int sn = segmentIt->sn;             // number of transfinite points (radial edges of stiffener)
    double t1 = std::accumulate(pointIt1->t.begin(), pointIt1->t.end(), 0.); // total thickness of outer structural layer
    double t2 = std::accumulate(pointIt2->t.begin(), pointIt2->t.end(), 0.); // total thickness of outer structural layer
    double tt1 = pointIt1->tt[0];       // thickness of inner thermal insulating layer
    double tt2 = pointIt2->tt[0];       // thickness of inner thermal insulating layer
    int tn_inner = segmentIt->tn[0];    // number of transfinite points through thickness of inner part of insulating layer in thermal mesh
    int tn_outer = segmentIt->tn[1];    // number of transfinite points through thickness of outer part of insulating layer in thermal mesh
    int ln = segmentIt->ln;             // number of transfinite points through each half of the thickness of the load layer in thermal mesh
    double angle = 2*M_PI/ns;

    // middle of inner layer
    cth1 = 1/sqrt(1+pointIt1->dydx*pointIt1->dydx); // cosine of angle between normal and y-axis 
    cth2 = 1/sqrt(1+pointIt2->dydx*pointIt2->dydx); // cosine of angle between normal and y-axis
    std::vector<GEdge*> middleInnerLayer;
    std::vector<std::vector<double> >::iterator it;
    for(unsigned int k=0; k<pointIt1->tt.size()-1; ++k ) {
      GVertex *vertex3 = m->addVertex(pointIt1->xyz[0], pointIt1->xyz[1]+tt1/cth1, pointIt1->xyz[2], lc);
      GVertex *vertex4 = m->addVertex(pointIt2->xyz[0], pointIt2->xyz[1]+tt2/cth2, pointIt2->xyz[2], lc);
      for(it = controlPoints.begin(), pointIt = pointIt1+1; it != controlPoints.end(); ++it, ++pointIt) {
        double cth = 1/sqrt(1+pointIt->dydx*pointIt->dydx); // cosine of angle between normal and y-axis
        (*it)[1] += pointIt->tt[k]/cth;
      }
      middleInnerLayer.push_back(controlPoints.empty() ? m->addLine(vertex3, vertex4) : m->addBSpline(vertex3, vertex4, controlPoints));
      tt1 += pointIt1->tt[k+1];
      tt2 += pointIt2->tt[k+1];
    }
    // outside of inner layer / inside of outer layer
    GVertex *vertex3 = m->addVertex(pointIt1->xyz[0], pointIt1->xyz[1]+tt1/cth1, pointIt1->xyz[2], lc);
    GVertex *vertex4 = m->addVertex(pointIt2->xyz[0], pointIt2->xyz[1]+tt2/cth2, pointIt2->xyz[2], lc);
    for(it = controlPoints.begin(), pointIt = pointIt1+1; it != controlPoints.end(); ++it, ++pointIt) {
      double cth = 1/sqrt(1+pointIt->dydx*pointIt->dydx); // cosine of angle between normal and y-axis
      (*it)[1] += pointIt->tt.back()/cth;
    }
    GEdge *edge2 = controlPoints.empty() ? m->addLine(vertex3, vertex4) : m->addBSpline(vertex3, vertex4, controlPoints);
    // middle of outer layer
    GVertex *vertex5 = m->addVertex(pointIt1->xyz[0], pointIt1->xyz[1]+(tt1+t1/2)/cth1, pointIt1->xyz[2], lc);
    GVertex *vertex6 = m->addVertex(pointIt2->xyz[0], pointIt2->xyz[1]+(tt2+t2/2)/cth2, pointIt2->xyz[2], lc);
    for(it = controlPoints.begin(), pointIt = pointIt1+1; it != controlPoints.end(); ++it, ++pointIt) {
      double cth = 1/sqrt(1+pointIt->dydx*pointIt->dydx); // cosine of angle between normal and y-axis
      (*it)[1] += (std::accumulate(pointIt->t.begin(), pointIt->t.end(), 0.)/2.)/cth;
    }
    GEdge *edge3 = controlPoints.empty() ? m->addLine(vertex5, vertex6) : m->addBSpline(vertex5, vertex6, controlPoints);
    // outside of outer layer
    GVertex *vertex7 = m->addVertex(pointIt1->xyz[0], pointIt1->xyz[1]+(tt1+t1)/cth1, pointIt1->xyz[2], lc);
    GVertex *vertex8 = m->addVertex(pointIt2->xyz[0], pointIt2->xyz[1]+(tt2+t2)/cth2, pointIt2->xyz[2], lc);
    for(it = controlPoints.begin(), pointIt = pointIt1+1; it != controlPoints.end(); ++it, ++pointIt) {
      double cth = 1/sqrt(1+pointIt->dydx*pointIt->dydx); // cosine of angle between normal and y-axis
      (*it)[1] += (std::accumulate(pointIt->t.begin(), pointIt->t.end(), 0.)/2.)/cth;
    }
    GEdge *edge4 = controlPoints.empty() ? m->addLine(vertex7, vertex8) : m->addBSpline(vertex7, vertex8, controlPoints);
    // top of stiffeners
    GVertex *vertex9, *vertex10; GEdge *edge11;
    if(ws1 != 0) {
      vertex9  = m->addVertex(pointIt1->xyz[0], pointIt1->xyz[1]+(tt1+t1/2)/cth1+ws1, pointIt1->xyz[2], lc);
      vertex10 = m->addVertex(pointIt2->xyz[0], pointIt2->xyz[1]+(tt2+t2/2)/cth2+ws2, pointIt2->xyz[2], lc);
      for(it = controlPoints.begin(), pointIt = pointIt1+1; it != controlPoints.end(); ++it, ++pointIt) {
        double cth = 1/sqrt(1+pointIt->dydx*pointIt->dydx); // cosine of angle between normal and y-axis
        (*it)[1] -= (std::accumulate(pointIt->t.begin(), pointIt->t.end(), 0.)/2.)/cth;
        (*it)[1] += pointIt->ws/cth;
      }
      edge11 = controlPoints.empty() ? m->addLine(vertex9, vertex10) : m->addBSpline(vertex9, vertex10, controlPoints);
    }

    GEntity *region1, *region2, *region3, *region4, *face3;
    // loop over the panels
    for(int i=0; i<ns; ++i) {
      if(i > 0) {
        std::list<GEdge*> edges = face3->cast2Face()->edges();
        edge3 = *(++edges.begin()); // the second edge is the one we need to revolve/extrude!
        if(ws1 != 0) {
          double cth = std::cos(2*M_PI/ns), sth = std::sin(2*M_PI/ns);
          // rotate end points
          double x = vertex9->x(), y = vertex9->y(), z = vertex9->z();
          vertex9  = m->addVertex(x, cth*y - sth*z, sth*y + cth*z, lc);
          x = vertex10->x(), y = vertex10->y(), z = vertex10->z();
          vertex10 = m->addVertex(x, cth*y - sth*z, sth*y + cth*z, lc);
          // rotate control points
          for(std::vector<std::vector<double> >::iterator it = controlPoints.begin(); it != controlPoints.end(); ++it) {
            y = (*it)[1], z = (*it)[2];
            (*it)[1] = cth*y - sth*z;
            (*it)[2] = sth*y + cth*z;
          }
          edge11 = controlPoints.empty() ? m->addLine(vertex9, vertex10) : m->addBSpline(vertex9, vertex10, controlPoints);
        }
      }

      // points defining radial unit vector
      std::vector<double> p3(3), p4(3);
      p3[0] = edge3->getBeginVertex()->x(); p3[1] = p3[2] = 0;
      p4[0] = p3[0]; p4[1] = edge3->getBeginVertex()->y(); p4[2] = edge3->getBeginVertex()->z();
      double r = std::sqrt(p4[1]*p4[1]+p4[2]*p4[2]);
      p4[1] /= r; p4[2] /= r;

      // main shell:
      face3 = m->revolve(edge3, p1, p2, angle);
      face3->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
      face3->cast2Face()->meshAttributes.transfiniteArrangement = 1;
      std::list<GEdge*> edges = face3->cast2Face()->edges();
      for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
        (*it)->meshAttributes.method = MESH_TRANSFINITE;
        int edgeIndex = std::distance(edges.begin(),it);
        switch(edgeIndex) {
          case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = mn; break;
          case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = nn; break;
        }
        (*it)->meshAttributes.coeffTransfinite = 1.0;
        surfaceTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
        surfaceTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
        surfaceTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
        if(bf != 1 && segmentIt == segments.begin() && edgeIndex == 0) {
          boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
          boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
          boundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
        }
      }
      face3->addPhysicalEntity(2);
      surfaceTags.push_back(std::make_pair(face3->geomType(),face3->tag()));

      // baffle:
      if(wb1 > 0) {
        if(ws1 == 0) { // no stiffeners
          std::vector<double> p5(p4); p5[1] *= wb1; p5[2] *= wb1; // point defining extrusion vector
          GEntity *baffleEdge = m->extrude(edge3->getBeginVertex(), p3, p5);
          GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
          baffleFace->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
          baffleFace->cast2Face()->meshAttributes.transfiniteArrangement = 1;
          std::list<GEdge*> edges = baffleFace->cast2Face()->edges();
          for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
            (*it)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it);
            switch(edgeIndex) {
              case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = mn; break;
              case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = nb1; break;
            }
            (*it)->meshAttributes.coeffTransfinite = 1.0;
            if(bf != 0 && edgeIndex == 2) {
              boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
              boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
              boundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
            }
          }
          baffleFace->addPhysicalEntity(5+baffleIndex);
        }
        else {
          std::vector<double> p5(p4); p5[1] *= ws1; p5[2] *= ws1; // point defining 1st extrusion vector
          GEntity *baffleEdge = m->extrude(edge3->getBeginVertex(), p3, p5);
          GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
          baffleFace->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
          baffleFace->cast2Face()->meshAttributes.transfiniteArrangement = 1;
          std::list<GEdge*> edges = baffleFace->cast2Face()->edges();
          GEdge* edge5;
          for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
            (*it)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it);
            switch(edgeIndex) {
              case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = mn; break;
              case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = sn; break;
            }
            (*it)->meshAttributes.coeffTransfinite = 1.0;
            if(edgeIndex == 2) {
              if(bf != 0 && ws1 == wb1) { // then, this is the boundary
                boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
              }
              else edge5 = *it;
            }
          }
          baffleFace->addPhysicalEntity(5+baffleIndex);

          if(ws1 < wb1) {
            std::vector<double> p5(p4); p5[1] *= (wb1-ws1); p5[2] *= (wb1-ws1); // point defining 2nd extrusion vector
            GEntity *baffleEdge = m->extrude(edge5->getBeginVertex(), p3, p5);
            GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
            baffleFace->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
            baffleFace->cast2Face()->meshAttributes.transfiniteArrangement = 1;
            std::list<GEdge*> edges = baffleFace->cast2Face()->edges();
            for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
              (*it)->meshAttributes.method = MESH_TRANSFINITE;
              int edgeIndex = std::distance(edges.begin(),it);
              switch(edgeIndex) {
                case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = nb1; break;
              }
              (*it)->meshAttributes.coeffTransfinite = 1.0;
              if(bf != 0 && edgeIndex == 2) {
                boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
              }
            }
            baffleFace->addPhysicalEntity(5+baffleIndex);
          }
          else if(ws1 > wb1) std::cerr << "error: stringer height is greater than baffle height\n"; 
        }
      }

      // stiffener:
      if(ws1 > 0) {
        GEdge *edge12 = m->addLine(edge3->getEndVertex(), edge11->getEndVertex());
        GEdge *edge13 = m->addLine(edge11->getBeginVertex(), edge3->getBeginVertex());
        std::vector<std::vector<GEdge*> > loop(1);
        loop[0].push_back(edge3); loop[0].push_back(edge12); loop[0].push_back(edge11); loop[0].push_back(edge13);
        GEntity *stiffenerFace = m->addPlanarFace(loop);
        stiffenerFace->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
        stiffenerFace->cast2Face()->meshAttributes.transfiniteArrangement = 1;
        std::list<GEdge*> edges = stiffenerFace->cast2Face()->edges();
        for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
          switch(std::distance(edges.begin(),it)) {
            case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = nn; break;
            case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = sn; break;
          }
          (*it)->meshAttributes.coeffTransfinite = 1.0;
        }
        stiffenerFace->addPhysicalEntity(4);
      }

      // final baffle:
      if(wb2 > 0 && segmentIt+1 == segments.end()) {
        int baffleIndex2 = (wb1 > 0) ? baffleIndex+1 : baffleIndex;
        std::vector<double> p3(3), p4(3);
        p3[0] = edge3->getEndVertex()->x(); p3[1] = p3[2] = 0;
        p4[0] = p3[0]; p4[1] = edge3->getEndVertex()->y(); p4[2] = edge3->getEndVertex()->z();
        double r = std::sqrt(p4[1]*p4[1]+p4[2]*p4[2]);
        p4[1] /= r; p4[2] /= r;

        if(ws2 == 0) { // no stiffener
          std::vector<double> p5(p4); p5[1] *= wb2; p5[2] *= wb2; // point defining extrusion vector
          GEntity *baffleEdge = m->extrude(edge3->getEndVertex(), p3, p5);
          GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
          baffleFace->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
          baffleFace->cast2Face()->meshAttributes.transfiniteArrangement = 1;
          std::list<GEdge*> edges = baffleFace->cast2Face()->edges();
          for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
            (*it)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it);
            switch(edgeIndex) {
              case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = mn; break;
              case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = nb2; break;
            }
            (*it)->meshAttributes.coeffTransfinite = 1.0;
            if(bf != 0 && edgeIndex == 2) {
              boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
              boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
              boundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
            }
          }
          baffleFace->addPhysicalEntity(5+baffleIndex2);
        }
        else {
          std::vector<double> p5(p4); p5[1] *= ws2; p5[2] *= ws2; // point defining 1st extrusion vector
          GEntity *baffleEdge = m->extrude(edge3->getEndVertex(), p3, p5);
          GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
          baffleFace->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
          baffleFace->cast2Face()->meshAttributes.transfiniteArrangement = 1;
          std::list<GEdge*> edges = baffleFace->cast2Face()->edges();
          GEdge *edge5;
          for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
            (*it)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it);
            switch(edgeIndex) {
              case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = mn; break;
              case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = sn; break;
            }
            (*it)->meshAttributes.coeffTransfinite = 1.0;
            if(edgeIndex == 2) {
              if(bf !=0 && ws2 == wb2) { // then, this is the boundary
                boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
              }
              else edge5 = *it;
            }
          }
          baffleFace->addPhysicalEntity(5+baffleIndex2);

          if(ws2 < wb2) {
            std::vector<double> p5(p4); p5[1] *= (wb2-ws2); p5[2] *= (wb2-ws2); // point defining 2nd extrusion vector
            GEntity *baffleEdge = m->extrude(edge5->getBeginVertex(), p3, p5);
            GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
            baffleFace->cast2Face()->meshAttributes.method = MESH_TRANSFINITE;
            baffleFace->cast2Face()->meshAttributes.transfiniteArrangement = 1;
            std::list<GEdge*> edges = baffleFace->cast2Face()->edges();
            for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++ it) {
              (*it)->meshAttributes.method = MESH_TRANSFINITE;
              int edgeIndex = std::distance(edges.begin(),it);
              switch(edgeIndex) {
                case 0: case 2: (*it)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it)->meshAttributes.nbPointsTransfinite = nb2; break;
              }
              (*it)->meshAttributes.coeffTransfinite = 1.0;
              if(bf != 0 && edgeIndex == 2) {
                boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(),(*it)->getBeginVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(),(*it)->getEndVertex()->tag()));
                boundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
              }
            }
            baffleFace->addPhysicalEntity(5+baffleIndex2);
          }
          else if(ws2 > wb2) std::cerr << "error: stringer height is greater than baffle height\n";
        }
      }

      // thermal model inside half of inner layer:
      if(tf != 0) {
        GEdge *edge2 = middleInnerLayer.front();
        if(i==0) {
          GEdge *edge5 = m->addLine(edge1->getEndVertex(), edge2->getEndVertex());
          GEdge *edge6 = m->addLine(edge2->getBeginVertex(), edge1->getBeginVertex());
          std::vector<std::vector<GEdge*> > edges(1);
          edges[0].push_back(edge1); edges[0].push_back(edge5); edges[0].push_back(edge2); edges[0].push_back(edge6);
          GFace *solidFace = m->addPlanarFace(edges);
          region1 = m->revolve(solidFace, p1, p2, angle);
        }
        else {
          std::list<GFace*> faces = region1->cast2Region()->faces();
          GFace *solidFace = faces.back();
          region1 = m->revolve(solidFace, p1, p2, angle);
        }

        region1->cast2Region()->meshAttributes.method = MESH_TRANSFINITE;
        region1->cast2Region()->meshAttributes.recombine3D = 1;
  
        std::list<GFace*> faces = region1->cast2Region()->faces();
        for(std::list<GFace*>::iterator it = faces.begin(); it != faces.end(); ++ it) {
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
          (*it)->meshAttributes.recombine = 1;
          int faceIndex = std::distance(faces.begin(),it);
          std::list<GEdge*> edges = (*it)->edges();
          for(std::list<GEdge*>::iterator it2 = edges.begin(); it2 != edges.end(); ++ it2) {
            (*it2)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it2);
            if(faceIndex == 0 || faceIndex == 2) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
              }
            }
            else if(faceIndex == 1 || faceIndex == 3) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = tn_inner; break;
              }
            }
            else if(faceIndex == 4 || faceIndex == 5) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = tn_inner; break;
              }
            }
            (*it2)->meshAttributes.coeffTransfinite = 1.0;
  
            if(faceIndex == 0) { // note: faceIndex 2 is outside, faceIndex 0 is inside
              interiorBoundaryTags.push_back(std::make_pair((*it2)->getBeginVertex()->geomType(),(*it2)->getBeginVertex()->tag()));
              interiorBoundaryTags.push_back(std::make_pair((*it2)->getEndVertex()->geomType(),(*it2)->getEndVertex()->tag()));
              interiorBoundaryTags.push_back(std::make_pair((*it2)->geomType(),(*it2)->tag()));
            }
            if(segmentIt == segments.begin() && faceIndex == 3) {
              cmcBoundaryTags.push_back(std::make_pair((*it2)->getBeginVertex()->geomType(),(*it2)->getBeginVertex()->tag()));
              cmcBoundaryTags.push_back(std::make_pair((*it2)->getEndVertex()->geomType(),(*it2)->getEndVertex()->tag()));
              cmcBoundaryTags.push_back(std::make_pair((*it2)->geomType(),(*it2)->tag()));
            }
          }
          if(faceIndex == 0) {
            interiorBoundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
            (*it)->addPhysicalEntity(0);
          }
          if(segmentIt == segments.begin() && faceIndex == 3) {
            cmcBoundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
          }
        }
  
        region1->addPhysicalEntity(0);
      }

      // thermal model outside half of inner layer:
      if(tf != 0) {
        GEdge *edge1 = middleInnerLayer.front();
        if(i==0) {
          GEdge *edge5 = m->addLine(edge1->getEndVertex(), edge2->getEndVertex());
          GEdge *edge6 = m->addLine(edge2->getBeginVertex(), edge1->getBeginVertex());
          std::vector<std::vector<GEdge*> > edges(1);
          edges[0].push_back(edge1); edges[0].push_back(edge5); edges[0].push_back(edge2); edges[0].push_back(edge6);
          GFace *solidFace = m->addPlanarFace(edges);
          region2 = m->revolve(solidFace, p1, p2, angle);
        }
        else {
          std::list<GFace*> faces = region2->cast2Region()->faces();
          GFace *solidFace = faces.back();
          region2 = m->revolve(solidFace, p1, p2, angle);
        }

        region2->cast2Region()->meshAttributes.method = MESH_TRANSFINITE;
        region2->cast2Region()->meshAttributes.recombine3D = 1;
  
        std::list<GFace*> faces = region2->cast2Region()->faces();
        for(std::list<GFace*>::iterator it = faces.begin(); it != faces.end(); ++ it) {
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
          (*it)->meshAttributes.recombine = 1;
          int faceIndex = std::distance(faces.begin(),it);
          std::list<GEdge*> edges = (*it)->edges();
          for(std::list<GEdge*>::iterator it2 = edges.begin(); it2 != edges.end(); ++ it2) {
            (*it2)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it2);
            if(faceIndex == 0 || faceIndex == 2) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
              }
            }
            else if(faceIndex == 1 || faceIndex == 3) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = tn_outer; break;
              }
            }
            else if(faceIndex == 4 || faceIndex == 5) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = tn_outer; break;
              }
            }
            (*it2)->meshAttributes.coeffTransfinite = 1.0;
          }
          if(faceIndex == 2) { // note: faceIndex 2 is outside, faceIndex 0 is inside
            (*it)->addPhysicalEntity(1);
          }
        }
  
        region2->addPhysicalEntity(999);
      }
  
      // thermal model inside half of outer layer
      if(tf != 0) {
        if(i==0) {
          GEdge *edge7 = m->addLine(edge2->getEndVertex(), edge3->getEndVertex());
          GEdge *edge8 = m->addLine(edge3->getBeginVertex(), edge2->getBeginVertex());
          std::vector<std::vector<GEdge*> > edges(1);
          edges[0].push_back(edge2); edges[0].push_back(edge7); edges[0].push_back(edge3); edges[0].push_back(edge8);
          GFace *solidFace = m->addPlanarFace(edges);
          region3 = m->revolve(solidFace, p1, p2, angle);
        }
        else {
          std::list<GFace*> faces = region3->cast2Region()->faces();
          GFace *solidFace = faces.back();
          region3 = m->revolve(solidFace, p1, p2, angle);
        }

        region3->cast2Region()->meshAttributes.method = MESH_TRANSFINITE;
        region3->cast2Region()->meshAttributes.recombine3D = 1;
  
        std::list<GFace*> faces = region3->cast2Region()->faces();
        for(std::list<GFace*>::iterator it = faces.begin(); it != faces.end(); ++ it) {
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
          (*it)->meshAttributes.recombine = 1;
          int faceIndex = std::distance(faces.begin(),it);
          std::list<GEdge*> edges = (*it)->edges();
          for(std::list<GEdge*>::iterator it2 = edges.begin(); it2 != edges.end(); ++ it2) {
            (*it2)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it2);
            if(faceIndex == 0 || faceIndex == 2) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
              }
            }
            else if(faceIndex == 1 || faceIndex == 3) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = ln; break;
              }
            }
            else if(faceIndex == 4 || faceIndex == 5) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = ln; break;
              }
            }
            (*it2)->meshAttributes.coeffTransfinite = 1.0;
          }
        }
  
        region3->addPhysicalEntity(1);
      }

      // thermal model outside half of outer layer
      if(tf != 0) {
        if(i==0) {
          GEdge *edge9 = m->addLine(edge3->getEndVertex(), edge4->getEndVertex());
          GEdge *edge10 = m->addLine(edge4->getBeginVertex(), edge3->getBeginVertex());
          std::vector<std::vector<GEdge*> > edges(1);
          edges[0].push_back(edge3); edges[0].push_back(edge9); edges[0].push_back(edge4); edges[0].push_back(edge10);
          GFace *solidFace = m->addPlanarFace(edges);
          region4 = m->revolve(solidFace, p1, p2, angle);
        }
        else {
          std::list<GFace*> faces = region4->cast2Region()->faces();
          GFace *solidFace = faces.back();
          region4 = m->revolve(solidFace, p1, p2, angle);
        }

        region4->cast2Region()->meshAttributes.method = MESH_TRANSFINITE;
        region4->cast2Region()->meshAttributes.recombine3D = 1;
  
        std::list<GFace*> faces = region4->cast2Region()->faces();
        for(std::list<GFace*>::iterator it = faces.begin(); it != faces.end(); ++ it) {
          (*it)->meshAttributes.method = MESH_TRANSFINITE;
          (*it)->meshAttributes.recombine = 1;
          int faceIndex = std::distance(faces.begin(),it);
          std::list<GEdge*> edges = (*it)->edges();
          for(std::list<GEdge*>::iterator it2 = edges.begin(); it2 != edges.end(); ++ it2) {
            (*it2)->meshAttributes.method = MESH_TRANSFINITE;
            int edgeIndex = std::distance(edges.begin(),it2);
            if(faceIndex == 0 || faceIndex == 2) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
              }
            }
            else if(faceIndex == 1 || faceIndex == 3) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = mn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = ln; break;
              }
            }
            else if(faceIndex == 4 || faceIndex == 5) {
              switch(edgeIndex) {
                case 0: case 2: (*it2)->meshAttributes.nbPointsTransfinite = nn; break;
                case 1: case 3: (*it2)->meshAttributes.nbPointsTransfinite = ln; break;
              }
            }
            (*it2)->meshAttributes.coeffTransfinite = 1.0;
  
            if(faceIndex == 2) { // note: faceIndex 2 is outside, faceIndex 0 is inside
              exteriorBoundaryTags.push_back(std::make_pair((*it2)->getBeginVertex()->geomType(),(*it2)->getBeginVertex()->tag()));
              exteriorBoundaryTags.push_back(std::make_pair((*it2)->getEndVertex()->geomType(),(*it2)->getEndVertex()->tag()));
              exteriorBoundaryTags.push_back(std::make_pair((*it2)->geomType(),(*it2)->tag()));
            }
          }
          if(faceIndex == 2) {
            exteriorBoundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
            (*it)->addPhysicalEntity(3);
          }
        }
  
        region4->addPhysicalEntity(3);
      }

    }

    pointIt1 = pointIt2;
    vertexIt1 = vertexIt2;
    vertex1 = vertex2;
    if(wb1 > 0) baffleIndex++;
  }

  //m->writeGEO("nozzle.geo");

  m->mesh(3);

  double tolerance = CTX::instance()->geom.tolerance;
  m->removeDuplicateMeshVertices(tolerance);

  m->writeMSH("nozzle.msh");
  writeAEROS(m, points, vertices, segments, materials, boundaries, surfaceTags, boundaryTags, lf, "nozzle.aeros", 2, false, 1.0, (tf==0));
  if(tf != 0) {
    writeAEROH(m, points, vertices, segments, materials, boundaries, interiorBoundaryTags, exteriorBoundaryTags, "nozzle.aeroh", 2);
    writeAEROS2(m, points, vertices, segments, materials, boundaries, cmcBoundaryTags, lf, "nozzle.aeros.cmc", 2);
  }

  delete m;

  GmshFinalize();
}

static PyObject *nozzle_generate(PyObject *self, PyObject *args)
{
  std::ifstream fin("NOZZLE.txt");

  int np, nv, nm; double lc; int bf, tf, nl, nlt, lf;
  fin >> np >> nv >> nm >> lc >> bf >> tf >> nl >> nlt >> lf;
  
  std::vector<PointData> points(np);
  for(int i=0; i<np; ++i) {
    PointData &p = points[i];
    p.xyz.resize(3);
    fin >> p.xyz[0] >> p.xyz[1] >> p.dydx; p.xyz[2] = 0.;
    for(int k=0; k<nl; ++k) { double tk; fin >> tk; p.t.push_back(tk); }
    for(int k=0; k<nlt; ++k) { double tt; fin >> tt; p.tt.push_back(tt); }
    fin >> p.ts >> p.ws;
  }

  std::vector<VertexData> vertices(nv);
  for(int j=0; j<nv; ++j) {
    VertexData &v = vertices[j];
    fin >> v.p >> v.wb >> v.mb >> v.nb >> v.tb;
  }

  std::vector<SegmentData> segments(nv-1);
  for(int j=0; j<nv-1; ++j) {
    SegmentData &s = segments[j];
    for(int k=0; k<nl; ++k) { int mk; fin >> mk; s.m.push_back(mk); }
    fin >> s.ns >> s.ms >> s.nn >> s.mn >> s.sn;
    for(int k=0; k<nlt; ++k) { int mt; fin >> mt; s.mt.push_back(mt); }
    for(int k=0; k<nlt; ++k) { int tn; fin >> tn; s.tn.push_back(tn); }
    fin >> s.ln;
  }

  std::vector<MaterialData> materials(nm);
  for(int k=0; k<nm; ++k) {
    MaterialData &m = materials[k];
    std::string s;
    fin >> s;
    if(s == "ISOTROPIC") {
      m.type = MaterialData::ISOTROPIC;
      fin >> m.E >> m.nu >> m.rho >> m.w >> m.k >> m.h;
    }
    else {
      m.type = MaterialData::ANISOTROPIC;
      fin >> m.E1 >> m.E2 >> m.nu12 >> m.G12 >> m.mu1 >> m.mu2 >> m.rho >> m.w1 >> m.w2 >> m.w12 >> m.k >> m.h;
    }
  }

  fin.close();

  std::ifstream fin2("BOUNDARY.txt");

  int m;
  fin2 >> m;

  std::vector<BoundaryData> boundaries(m); 
  for(int i=0; i<m; ++i) {
    BoundaryData &b = boundaries[i];
    fin2 >> b.x >> b.P >> b.T >> b.Ta;
  }

  fin2.close();

  generateNozzle(points, vertices, segments, materials, boundaries, lc, bf, tf, lf);

  Py_RETURN_NONE;
}

static PyObject *nozzle_convert(PyObject *self, PyObject *args)
{
  {std::ifstream in("TEMP.2");
  std::ofstream out("TEMPERATURES.txt");

  std::string s;
  getline(in,s);
  getline(in,s);
  getline(in,s);

  int node; double x,y,z,t;
  out << "TEMPERATURE\n";
  while (!in.eof()) {
    in >> node >> x >> y >> z >> t;
    out << std::left << std::setw(8) << node << " " << std::setprecision(15) << t << std::endl;
  }

  in.close();
  out.close();}

  {std::ifstream in("TEMP.0");
  std::ofstream out("TEMPERATURES.txt.cmc");

  std::string s;
  getline(in,s);
  getline(in,s);
  getline(in,s);

  int node; double x,y,z,t;
  out << "TEMPERATURE\n";
  while (!in.eof()) {
    in >> node >> x >> y >> z >> t;
    out << std::left << std::setw(8) << node << " " << std::setprecision(15) << t << std::endl;
  }

  in.close();
  out.close();}

  Py_RETURN_NONE;
}

