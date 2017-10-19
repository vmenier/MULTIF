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
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
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
#include "Options.h"
#include "Context.h"
#include "Geo/GModelIO_OCC.h"
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepFeat_SplitShape.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepTools.hxx>
#include <Geom_Plane.hxx>
#include <GeomAPI_IntSS.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2dAPI_PointsToBSpline.hxx>
#include <GeomAPI.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomLProp_SLProps.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopExp.hxx>
#include <TopoDS.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>

struct SectionData {
  double xc, zc, ry, rz;
};

struct ExteriorData {
  double angle, offset, a, b; 
  std::pair<double,double> coord(double dx, double theta) {
    // dx = xexit - x
    double r = a*b/std::sqrt(std::pow(b*std::cos(theta),2) + std::pow(a*std::sin(theta),2));
    double c = offset + dx*std::tan(M_PI*angle/180.);
    double y = r*std::cos(theta);
    double z = r*std::sin(theta) + c;
    return std::make_pair(y,z);
  }
};

struct PointData {
  std::vector<double> xyz, t, tt;
  std::vector<std::pair<double,double>> at1, at2, at3, at4;
  double dydx, ts, ws;
  double T(int l, double theta, bool lflag = false) const { // layer thickness
    // lflag false: l = 0 is the inner (cmc) portion of the thermal insulating layer
    //              l = 1 is the outer (air) portion of the thermal insulating layer
    //              l = 2 is the inner half of all the load bearing layers
    //              l = 3 is the outer half of all the load bearing layers
    // lflag true:  l is the index of the load bearing layer 
    if(theta < 0) theta += 2*M_PI;
    auto interpolate1 = [&](double x1, double y1, double x2, double y2, double x) { return y1 + (y2-y1)*(x-x1)/(x2-x1); };
    auto interpolate2 = [&](const std::vector<std::pair<double,double>> &at, double a) -> double {
      // upper_bound returns an iterator to the first element in the range that is greater than value, or last if no such element is found.
      auto it = std::upper_bound(at.begin(), at.end(), a, [&](const double a, const std::pair<double,double> &b) {return a < b.first;});
      if(it == at.end()) // a is greater than or equal to the last angular coordinate: i.e. back <= a < front+2pi
        return interpolate1((it-1)->first, (it-1)->second, at.front().first+2*M_PI, at.front().second, a);
      else if(it == at.begin()) // a is less than the first angular coordinate: i.e. back-2pi <= a < front]
        return interpolate1(at.back().first-2*M_PI, at.back().second, it->first, it->second, a);
      else
        return interpolate1((it-1)->first, (it-1)->second, it->first, it->second, a);
    };
    if(lflag) {
      switch(l) {
        case 0: // inner load bearing layer
          return at1.empty() ? t[0] : interpolate2(at1, theta);
        case 1: // middle load bearing layer 
          return at2.empty() ? t[1] : interpolate2(at2, theta);
        case 2: // outer load bearing layer
          return at3.empty() ? t[2] : interpolate2(at3, theta); 
        default:
          std::cerr << "Error: thickness of load layer " << l << " is not available\n";
          return 0;
      }
    }
    else {
      switch(l) {
        case 0: // inner (cmc) portion of the thermal insulating layer
          return at4.empty() ? tt[0] : interpolate2(at4, theta);
        case 1: // outer (air) portion of the thermal insulating layer
          return tt[1];
        case 2: case 3: // inner/outer half of the load bearing layer
          return 0.5*( (at1.empty() ? t[0] : interpolate2(at1, theta))
                      +(at2.empty() ? t[1] : interpolate2(at2, theta))
                      +(at3.empty() ? t[2] : interpolate2(at3, theta)) );
        default:
          std::cerr << "Error: thickness of layer " << l << " is not available\n";
          return 0;
      }
    }
  }
  double C(int l, double theta) const { // cumulative layer thickness
    double sum = 0;
    for(int i = 0; i <= l; ++i) sum += T(i, theta);
    return sum;
  }
};

struct VertexData {
  int p, mb;
  double wb;
  int nb;
  double tb;
};

struct SegmentData {
 std::vector<int> m;
 int ns, ms, nn, mn, sn;
 std::vector<int> mt, tn;
 int ln;
};

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

struct BoundaryData {
  double x, P, T, Ta;
};

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
  double Taval = (it == boundaries.begin()) ? it->Ta : ((it == boundaries.end()) ? (it-1)->Ta :
                 ((it-1)->Ta + (it->Ta - (it-1)->Ta)*(x - (it-1)->x)/(it->x - (it-1)->x)));

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
  double x = 0, y = 0, z = 0;
  for(int i = 0; i < 3; i++) {
    x += m->getVertexBDF(i)->x() * scalingFactor;
    y += m->getVertexBDF(i)->y() * scalingFactor;
    z += m->getVertexBDF(i)->z() * scalingFactor;
  } 
  x /= 3;
  y /= 3;
  z /= 3;

  fprintf(fp, "LAYC %d\n", m->getNum());
  // interpolate each layer thickness at centroid from vertex data
  double tval;
  Cmp cmp(points);
  auto it = std::lower_bound(vertices.begin(), vertices.end(), x, cmp);
  int segment_id = std::min<long>(segments.size()-1,std::distance(vertices.begin(), it));
  auto it2 = std::upper_bound(points.begin(), points.end(), x, [&](const double &a, const PointData &b) {return a < b.xyz[0];});
  double zc = (it2 == points.end()) ? points.back().xyz[2] : ((it2 == points.begin()) ? it2->xyz[2]
            : ((it2-1)->xyz[2] + (it2->xyz[2] - (it2-1)->xyz[2])*(x - (it2-1)->xyz[0])/(it2->xyz[0] - (it2-1)->xyz[0])));
  double theta = std::atan2(z-zc, y);
  for(unsigned int k=0; k<it2->t.size(); ++k) {
    tval = (it2 == points.end()) ? points.back().T(k,theta,true) : ((it2 == points.begin()) ? it2->T(k,theta,true)
         : ((it2-1)->T(k,theta,true) + (it2->T(k,theta,true) - (it2-1)->T(k,theta,true))*(x - (it2-1)->xyz[0])/(it2->xyz[0] - (it2-1)->xyz[0])));
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
  double tval = (it == boundaries.begin()) ? it->T : ((it == boundaries.end()) ? (it-1)->T :
                ((it-1)->T + (it->T - (it-1)->T)*(x - (it-1)->x)/(it->x - (it-1)->x)));

  fprintf(fp, "%-8d %-16.9G\n", m->getIndex(), tval);
}

void writeConvec(MElement *m, FILE *fp, double scalingFactor, const std::vector<BoundaryData> &boundaries,
                 const std::vector<PointData> &points, const std::vector<VertexData> &vertices,
                 const std::vector<SegmentData> &segments, const std::vector<MaterialData> &materials)
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
    double x = m->getVertexBDF(i)->x() * scalingFactor;
    std::vector<BoundaryData>::const_iterator it = std::lower_bound(boundaries.begin(), boundaries.end(), x, cmp);
    double Ta = (it == boundaries.begin()) ? it->Ta : ((it == boundaries.end()) ? (it-1)->Ta :
                ((it-1)->Ta + (it->Ta - (it-1)->Ta)*(x - (it-1)->x)/(it->x - (it-1)->x)));

    // find the convection coefficient of the load layer
    Cmp cmp(points);
    std::vector<VertexData>::const_iterator it2 = std::lower_bound(vertices.begin(), vertices.end(), x, cmp);
    double h = (it2 == vertices.end()) ? materials[segments.back().m[0]].h 
                                       : materials[(segments.begin()+std::distance(vertices.begin(),it2-1))->m[0]].h;

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
    fprintf(fp, "failsafe\n");
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
  // stresses and strains in the x and y directions of the material coordinate frame (layers 1 and 3 only)
  fprintf(fp, "stressxx 14 7 \"STRESSXX.1\" 1 NG 2 lower matfrm\n");
  fprintf(fp, "stressxx 14 7 \"STRESSXX.3\" 1 NG 2 upper matfrm\n");
  fprintf(fp, "stressyy 14 7 \"STRESSYY.1\" 1 NG 2 lower matfrm\n");
  fprintf(fp, "stressyy 14 7 \"STRESSYY.3\" 1 NG 2 upper matfrm\n");
  fprintf(fp, "strainxx 14 7 \"STRAINXX.1\" 1 NG 2 lower matfrm\n");
  fprintf(fp, "strainxx 14 7 \"STRAINXX.3\" 1 NG 2 upper matfrm\n");
  fprintf(fp, "strainyy 14 7 \"STRAINYY.1\" 1 NG 2 lower matfrm\n");
  fprintf(fp, "strainyy 14 7 \"STRAINYY.3\" 1 NG 2 upper matfrm\n");
  if(stringerCount > 0) {
    fprintf(fp, "stressxx 14 7 \"STRESSXX.4\" 1 NG 4 median matfrm\n");
    fprintf(fp, "stressyy 14 7 \"STRESSYY.4\" 1 NG 4 median matfrm\n");
    fprintf(fp, "strainxx 14 7 \"STRAINXX.4\" 1 NG 4 median matfrm\n");
    fprintf(fp, "strainyy 14 7 \"STRAINYY.4\" 1 NG 4 median matfrm\n");
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

  // file for computing the mass
  FILE *fp11 = Fopen((name+".mass").c_str(), "w");
  if(!fp11){
    Msg::Error("Unable to open file '%s'", (name+".mass").c_str());
    return 0;
  }

  fprintf(fp11, "MASS\n");
  fprintf(fp11, "*\n");
  fprintf(fp11, "INCLUDE \"%s\"\n*\n", "GEOMETRY.txt");
  fprintf(fp11, "INCLUDE \"%s\"\n*\n", "TOPOLOGY.txt");
  fprintf(fp11, "INCLUDE \"%s\"\n*\n", "ATTRIBUTES.txt");
  fprintf(fp11, "INCLUDE \"%s\"\n*\n", "MATERIAL.txt");
  fprintf(fp11, "INCLUDE \"%s\"\n*\n", "COMPOSITE.txt");
  fprintf(fp11, "INCLUDE \"%s\"\n*\n", "CFRAMES.txt");
  fprintf(fp11, "END\n");

  fclose(fp11);

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
  FILE *fp5 = Fopen("TOPOLOGY.txt.thermal.mass", "w");
  fprintf(fp2, "TOPOLOGY\n");
  fprintf(fp5, "TOPOLOGY\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    bool b = (std::find(exteriorBoundaryTags.begin(), exteriorBoundaryTags.end(),
              std::make_pair(entities[i]->geomType(),entities[i]->tag())) != exteriorBoundaryTags.end());
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
      if(str && std::strcmp(str,"CHEXA") == 0) { 
        writeElem(entities[i]->getMeshElement(j), fp2, 51);
        writeElem(entities[i]->getMeshElement(j), fp5, 17);
      }
      else if(str && b && std::strcmp(str,"CQUAD4") == 0)
        writeElem(entities[i]->getMeshElement(j), fp2, 48); 
    }
  }
  fclose(fp2);
  fclose(fp5);

  // attributes
  FILE *fp3 = Fopen("ATTRIBUTES.txt.thermal", "w");
  FILE *fp8 = Fopen("ATTRIBUTES.txt.thermal.mass", "w");
  fprintf(fp3, "ATTRIBUTES\n");
  fprintf(fp8, "ATTRIBUTES\n");
  for(unsigned int i = 0; i < entities.size(); i++) {
    bool b = (std::find(exteriorBoundaryTags.begin(), exteriorBoundaryTags.end(),
              std::make_pair(entities[i]->geomType(),entities[i]->tag())) != exteriorBoundaryTags.end());
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++) {
      const char *str = entities[i]->getMeshElement(j)->getStringForBDF();
      if(str && (std::strcmp(str,"CHEXA") == 0)) {
        writeAttr(entities[i]->getMeshElement(j), fp3, scalingFactor, entities[i]->physicals[0], points, vertices,
                  segments);
        writeAttr(entities[i]->getMeshElement(j), fp8, scalingFactor, entities[i]->physicals[0], points, vertices,
                  segments);
      }
      else if(str && b && std::strcmp(str,"CQUAD4") == 0) {
        writeAttr(entities[i]->getMeshElement(j), fp3, scalingFactor, entities[i]->physicals[0], points, vertices,
                  segments);
      }
    }
  }
  fclose(fp3);
  fclose(fp8);

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
          writeConvec(entities[i]->getMeshElement(j), fp6, scalingFactor, boundaries, points, vertices, segments, materials);
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

  // file for computing the mass
  FILE *fp9 = Fopen((name+".mass").c_str(), "w");
  if(!fp9){
    Msg::Error("Unable to open file '%s'", (name+".mass").c_str());
    return 0;
  }

  fprintf(fp9, "MASS\n");
  fprintf(fp9, "*\n");
  fprintf(fp9, "INCLUDE \"%s\"\n*\n", "GEOMETRY.txt.thermal");
  fprintf(fp9, "INCLUDE \"%s\"\n*\n", "TOPOLOGY.txt.thermal.mass");
  fprintf(fp9, "INCLUDE \"%s\"\n*\n", "ATTRIBUTES.txt.thermal.mass");
  fprintf(fp9, "INCLUDE \"%s\"\n*\n", "MATERIAL.txt.thermal");
  fprintf(fp9, "END\n");

  fclose(fp9);

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
    fprintf(fp, "failsafe\n");
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
  fprintf(fp, "stressp1 14 7 \"STRESSP1.cmc\" 1\n");
  fprintf(fp, "stressp1 14 7 \"STRESSP1.0\" 1 NG 1\n");
  if(lf == 1) {
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.cmc\" 1 thermal\n");
    fprintf(fp, "stressp1 14 7 \"THERMAL_STRESSP1.0\" 1 NG 1 thermal\n");
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.cmc\" 1 mechanical\n");
    fprintf(fp, "stressp1 14 7 \"MECHANICAL_STRESSP1.0\" 1 NG 1 mechanical\n");
  }
  // 2nd principal stress
  fprintf(fp, "stressp2 14 7 \"STRESSP2.cmc\" 1\n");
  fprintf(fp, "stressp2 14 7 \"STRESSP2.0\" 1 NG 1\n");
  if(lf == 1) {
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.cmc\" 1 thermal\n");
    fprintf(fp, "stressp2 14 7 \"THERMAL_STRESSP2.0\" 1 NG 1 thermal\n");
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.cmc\" 1 mechanical\n");
    fprintf(fp, "stressp2 14 7 \"MECHANICAL_STRESSP2.0\" 1 NG 1 mechanical\n");
  }
  // 3rd principal stress
  fprintf(fp, "stressp3 14 7 \"STRESSP3.cmc\" 1\n");
  fprintf(fp, "stressp3 14 7 \"STRESSP3.0\" 1 NG 1\n");
  if(lf == 1) {
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.cmc\" 1 thermal\n");
    fprintf(fp, "stressp3 14 7 \"THERMAL_STRESSP3.0\" 1 NG 1 thermal\n");
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.cmc\" 1 mechanical\n");
    fprintf(fp, "stressp3 14 7 \"MECHANICAL_STRESSP3.0\" 1 NG 1 mechanical\n");
  }
  // 1st principal strain
  fprintf(fp, "strainp1 14 7 \"STRAINP1.cmc\" 1\n");
  fprintf(fp, "strainp1 14 7 \"STRAINP1.0\" 1 NG 1\n");
  // 2nd principal strain
  fprintf(fp, "strainp2 14 7 \"STRAINP2.cmc\" 1\n");
  fprintf(fp, "strainp2 14 7 \"STRAINP2.0\" 1 NG 1\n");
  // 3rd principal strain
  fprintf(fp, "strainp3 14 7 \"STRAINP3.cmc\" 1\n");
  fprintf(fp, "strainp3 14 7 \"STRAINP3.0\" 1 NG 1\n");

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

  // file for computing the mass
  FILE *fp5 = Fopen((name+".mass").c_str(), "w");
  if(!fp5){
    Msg::Error("Unable to open file '%s'", (name+".mass").c_str());
    return 0;
  }

  fprintf(fp5, "MASS\n");
  fprintf(fp5, "*\n");
  fprintf(fp5, "INCLUDE \"%s\"\n*\n", "GEOMETRY.txt.cmc");
  fprintf(fp5, "INCLUDE \"%s\"\n*\n", "TOPOLOGY.txt.cmc");
  fprintf(fp5, "INCLUDE \"%s\"\n*\n", "ATTRIBUTES.txt.cmc");
  fprintf(fp5, "INCLUDE \"%s\"\n*\n", "MATERIAL.txt.cmc");
  fprintf(fp5, "END\n");

  fclose(fp5);

  return 1;
}

void generateNozzle(std::vector<PointData> &points,
                    const std::vector<VertexData> &vertices,
                    const std::vector<SegmentData> &segments,
                    const std::vector<MaterialData> &materials,
                    const std::vector<BoundaryData> &boundaries,
                    double lc, int tf, int lf, int vf,
                    const std::vector<std::vector<double>> controlpoints[3],
                    const std::vector<double> knots[3],
                    const std::vector<int> mult[3],
                    double theta_in, double theta_out, double baffle_half_width,
                    double top_angle, double top_offset, double top_a, double top_b,
                    double bot_angle, double bot_offset, double bot_a, double bot_b)
{
  // tf = 0: generate mesh for structural model only
  // tf = 1: generate mesh for both structural and thermal models

  // lf = 0: nonlinear structural model
  // lf = 1: linear structural model

  // vf = 0: non-verbose mode
  // vf = 1: verbose mode

  // parameters used by Geom2dAPI_PointsToBSpline
  const int DegMin = 3;                          // default is 3
  const int DegMax = 8;                          // default is 8
  enum Approx_ParametrizationType ParType = Approx_IsoParametric; // default is Approx_ChordLength
  const GeomAbs_Shape Continuity2D = GeomAbs_C2; // default is GeomAbs_C2
  const double Tol2D = 1.0e-15;                  // default is 1.0e-6

  // parameters used by BRepOffsetAPI_ThruSections
  const GeomAbs_Shape Continuity3D = GeomAbs_C2; // default is GeomAbs_C2
  const double Pres3D = 1.0e-6;                  // default is 1.0e-6

  // parameters used by GeomAPI_IntSS
  const double IntSSTol = 1.0e-7;                // recommended value is 1.0e-7

  // parameters used by ConnectEdgesToWires
  const double ConnectEdgesToWiresTol = 1.0e-4;

  // parameters associated with the geometry and mesh
  const bool UseChebyshevNodes = true;
  const int MeshingMethod = MESH_TRANSFINITE;
  const int NbControlPoints = 61;
  const int NbExteriorControlPoints = 41;
  const int NbPanels = 2; // Note: currently this must be set to 2
  const int NbThruSections = 11;
  const int NbInnerThruSections = 11;

  // local variables
  const int NbSegments = segments.size();
  double x_in, x_out;
  std::vector<SectionData> sections(NbThruSections), inner_sections(NbInnerThruSections);
  std::vector<std::pair<GEntity::GeomType,int>> surfaceTags, boundaryTags, interiorBoundaryTags, exteriorBoundaryTags, cmcBoundaryTags;
  ExteriorData top = {top_angle, top_offset, top_a, top_b};
  ExteriorData bot = {bot_angle, bot_offset, bot_a, bot_b};
  GModel *gm;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 1. define some lambda functions
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto make_curve = [&](const std::vector<std::vector<double>> &controlpoints,
                        const std::vector<double> &knots,
                        const std::vector<int> &mult) -> Handle(Geom2d_BSplineCurve) {
    TColgp_Array1OfPnt2d _controlpoints(1, controlpoints.size());
    TColStd_Array1OfReal _weights(1, controlpoints.size());
    for(int i = 0; i < controlpoints.size(); i++) {
      gp_Pnt2d aP(controlpoints[i][0], controlpoints[i][1]);
      _controlpoints.SetValue(i+1, aP);
      _weights.SetValue(i+1, 1.0);
    }
    TColStd_Array1OfReal _knots(1, knots.size());
    for(int i = 0; i < knots.size(); i++) {
      _knots.SetValue(i+1, knots[i]);
    }
    TColStd_Array1OfInteger _mult(1, mult.size());
    int totKnots = 0;
    for(int i = 0; i < mult.size(); i++) {
      _mult.SetValue(i+1, mult[i]);
      totKnots += mult[i];
    }
    const int degree = totKnots - controlpoints.size() - 1;
    return new Geom2d_BSplineCurve(_controlpoints, _weights, _knots, _mult, degree, false);
  };

  auto locate_x_curve_bisection = [&](Handle(Geom2d_Curve) curve, double x, double tol = 2e-14, int maxit = 1000,
                                      bool verbose = false) -> double {
    // solve e->Value(u).X() - x = 0 using bisection method
    double pFirst = curve->FirstParameter(), pLast = curve->LastParameter(), pMiddle;
    for(int i = 0; ; i++) {
      pMiddle = pFirst + (pLast-pFirst)/2;
      gp_Pnt2d p = curve->Value(pMiddle);
      if(verbose) std::cerr << "i = " << i << ", p.X() = " << p.X() << ", x = " << x << std::endl;
      if(std::abs(p.X()-x) < tol) break;
      if(i == maxit) {
        std::cerr << "Warning: bisection solver did not converge in " << maxit << " iterations (residual = " << p.X()-x << ")\n";
        break;
      }
      if(p.X() < x) pFirst = pMiddle;
      else pLast = pMiddle;
    }
    return pMiddle;
  };

  auto locate_y_curve_bisection = [&](Handle(Geom2d_Curve) curve, double y, double tol = 2e-14, int maxit = 1000,
                                      bool verbose = false) -> double {
    // solve e->Value(u).Y() - y = 0 using bisection method
    double pFirst = curve->FirstParameter(), pLast = curve->LastParameter(), pMiddle;
    for(int i = 0; ; i++) {
      pMiddle = pFirst + (pLast-pFirst)/2;
      gp_Pnt2d p = curve->Value(pMiddle);
      if(verbose) std::cerr << "i = " << i << ", p.Y() = " << p.Y() << ", y = " << y << std::endl;
      if(std::abs(p.Y()-y) < tol) break;
      if(i == maxit) {
        std::cerr << "Warning: bisection solver did not converge in " << maxit << " iterations (residual = " << p.Y()-y << ")\n";
        break;
      }
      if(p.Y() < y) pFirst = pMiddle;
      else pLast = pMiddle;
    }
    return pMiddle;
  };

  auto locate_x_surface_bisection = [&](Handle(Geom_Surface) s, double x, double u, double tol = 2e-14, int maxit = 1000,
                                        bool verbose = false) -> double {
    // solve s->Value(u,v).X() - x = 0 for unknown v using bisection method
    double umin, umax, vmin, vmax;
    s->Bounds(umin, umax, vmin, vmax);
    double pFirst = vmin, pLast = vmax, pMiddle;
    for(int i = 0; ; i++) {
      pMiddle = pFirst + (pLast-pFirst)/2;
      gp_Pnt p = s->Value(u, pMiddle);
      if(verbose) std::cerr << "i = " << i << ", p.X() = " << p.X() << ", x = " << x << std::endl;
      if(std::abs(p.X()-x) < tol) break;
      if(i == maxit) {
        std::cerr << "Warning: bisection solver did not converge in " << maxit << " iterations (residual = " << p.X()-x << ")\n";
        break;
      }
      if(p.X() < x) pFirst = pMiddle;
      else pLast = pMiddle;
    }
    return pMiddle;
  };

  auto locate_x_curve_newton = [&](Handle(Geom2d_Curve) e, double x, double u0, int maxit = 20, double tol = 2e-14,
                                   bool verbose = false) -> double {
    // solve e->Value(u).X() - x = 0 using Newton's method with initial guess u0
    double u = u0;
    for(int i = 0; ; ++i) {
      double f = e->Value(u).X() - x;
      if(verbose) std::cerr << "i = " << i << ", f = " << f << std::endl;
      if(std::abs(f) < tol) break;
      if(i == maxit) { std::cerr << "Warning: Newton solver did not converge in " << maxit << " iterations (residual = " << f << ")\n";
                       if(std::isnan(f)) u = u0;
                       break; }
      double dfdu = e->DN(u, 1).X();
      u -= f/dfdu;
    }
    return u;
  };

  auto locate_y_curve_newton = [&](Handle(Geom2d_Curve) e, double y, double u0, int maxit = 20, double tol = 2e-14,
                                   bool verbose = false) -> double {
    // solve e->Value(u).Y() - y = 0 using Newton's method with initial guess u0
    double u = u0;
    for(int i = 0; ; ++i) {
      double f = e->Value(u).Y() - y;
      if(verbose) std::cerr << "i = " << i << ", f = " << f << std::endl;
      if(std::abs(f) < tol) break;
      if(i == maxit) { std::cerr << "Warning: Newton solver did not converge in " << maxit << " iterations (residual = " << f << ")\n";
                       if(std::isnan(f)) u = u0;
                       break; }
      double dfdu = e->DN(u, 1).Y();
      u -= f/dfdu;
    }
    return u;
  };

  auto locate_x_curve = [&](Handle(Geom2d_Curve) e, double x, int maxit = 20, double tol = 2e-14,
                                   bool verbose = false) -> double {
    // solve e->Value(u).X() - x = 0 using Newton's method with initial guess u0 obtained using bisection method
    double u0 = locate_x_curve_bisection(e, x, 0.05, 100, verbose);
    return locate_x_curve_newton(e, x, u0, maxit, tol, verbose);
  };

  auto locate_y_curve = [&](Handle(Geom2d_Curve) e, double y, int maxit = 20, double tol = 2e-14,
                                   bool verbose = false) -> double {
    // solve e->Value(u).Y() - y = 0 using Newton's method with initial guess u0 obtained using bisection method
    double u0 = locate_y_curve_bisection(e, y, 0.05, 100, verbose);
    return locate_y_curve_newton(e, y, u0, maxit, tol, verbose);
  };

  auto make_inner_sections = [&](Handle(Geom2d_Curve) e0, Handle(Geom2d_Curve) e1, Handle(Geom2d_Curve) e2) -> void {
    // compute and store (a) the coordinates of the points on the centerline at the thru-section points: xc, zc
    //                   (b) the radii of the major and minor axes: ry, rz
    // Note: these "inner" sections are used to construct a reference surface from which other surfaces are lofted
    //       to facilitate lofting in the normal direction, the reference surface is constructed by
    //       extrapolating the inner surface of the nozzle by a 10% in the plus and minus x-directions
    double first0  = e0->FirstParameter(), last0 = e0->LastParameter();
    double first1  = e1->FirstParameter(), last1 = e1->LastParameter();
    double first2  = e2->FirstParameter(), last2 = e2->LastParameter();
    double delta_x = (x_out-x_in)/10;
    double x_min   = x_in - delta_x;
    double x_max   = x_out + delta_x;

    for(int i = 0; i < NbInnerThruSections; ++i) {
      // x-coordinate of the i-th thru-section
      double x = UseChebyshevNodes ? ((x_min+x_max)/2 + (x_min-x_max)/2*std::cos(i*M_PI/(NbInnerThruSections-1)))
                 : (x_min + i*(x_max-x_min)/(NbInnerThruSections-1));
      inner_sections[i].xc = x;
      // centerline y-coordinate
      if(x < e0->Value(first0).X()) {
        inner_sections[i].zc = e0->Value(first0).Y() + (x - e0->Value(first0).X())*e0->DN(first0,1).Y();
      }
      else if(x > e0->Value(last0).X()) {
        inner_sections[i].zc = e0->Value(last0).Y() + (x - e0->Value(last0).X())*e0->DN(last0,1).Y();
      }
      else {
        double u0 = locate_x_curve(e0, x);
        inner_sections[i].zc = e0->Value(u0).Y();
      }
      // major axis
      if(x < e1->Value(first1).X()) {
        inner_sections[i].ry = e1->Value(first1).Y() + (x - e1->Value(first1).X())*e1->DN(first1,1).Y();
      }
      else if(x > e1->Value(last1).X()) {
        inner_sections[i].ry = e1->Value(last1).Y() + (x - e1->Value(last1).X())*e1->DN(last1,1).Y();
      }
      else {
        double u1 = locate_x_curve(e1, x);
        inner_sections[i].ry = e1->Value(u1).Y();
      }
      // minor axis
      if(x < e2->Value(first2).X()) {
        inner_sections[i].rz = e2->Value(first2).Y() + (x - e2->Value(first2).X())*e2->DN(first2,1).Y();
      }
      else if(x > e2->Value(last2).X()) {
        inner_sections[i].rz = e2->Value(last2).Y() + (x - e2->Value(last2).X())*e2->DN(last2,1).Y();
      }
      else {
        double u2 = locate_x_curve(e2, x);
        inner_sections[i].rz = e2->Value(u2).Y();
      }
    }
  };

  auto make_sections = [&](Handle(Geom2d_Curve) e0) -> void {
    // compute and store the coordinates of the points on the centerline at the thru-section points: xc, zc
    double first0 = e0->FirstParameter(), last0 = e0->LastParameter();

    for(int i = 0; i < NbThruSections; ++i) {
      // x-coordinate of i-th thru-section
      double x = UseChebyshevNodes ?
                 ((x_in+x_out)/2 + (x_in-x_out)/2*std::cos(i*M_PI/(NbThruSections-1))) : (x_in + i*(x_out-x_in)/(NbThruSections-1));
      sections[i].xc = x;
      // y-coordinate
      double u0 = (i==0) ? first0 : ((i==NbThruSections-1) ? last0 : locate_x_curve(e0, x));
      gp_Pnt2d p0 = e0->Value(u0);
      sections[i].zc = p0.Y();
    }
  };

  auto finalize_points = [&](Handle(Geom2d_Curve) e0) -> void {
    // compute and store the coordinates of the points on the centerline at the points
    for(auto &p : points) {
      double u0 = locate_x_curve(e0, p.xyz[0]);
      gp_Pnt2d p0 = e0->Value(u0);
      p.xyz[2] = p0.Y();
    }
  };

  auto layer_thickness = [&](int l, double x, double theta) -> double {
    // compute the cumulative thickness of the l-th layer at x,theta
    // note: l = 0 is the inner (cmc) portion of the thermal insulating layer
    //       l = 1 is the outer (air) portion of the thermal insulating layer
    //       l = 2 is the inner half of the load bearing layer
    //       l = 3 is the outer half of the load bearing layer
    if(l == -1) return 0.;
    else {
      auto it = std::upper_bound(points.begin(), points.end(), x, [&](const double &a, const PointData &b) {return a < b.xyz[0];});
      if(it == points.end()) return points.back().C(l,theta);
      else if(it == points.begin()) return it->C(l,theta);
      else return (it-1)->C(l,theta) + (it->C(l,theta) - (it-1)->C(l,theta))*(x - (it-1)->xyz[0])/(it->xyz[0] - (it-1)->xyz[0]);
    }
  };

  auto layer_thickness_x_deriv = [&](int l, double x, double theta) -> double {
    // compute the derivative w.r.t x of the of the l-th layer cumulative thickness
    if(l == -1) return 0.;
    else {
      auto it = std::lower_bound(points.begin(), points.end(), x, [&](const PointData &p, const double &b) 
                                 {return p.xyz[0] < (b+std::numeric_limits<double>::epsilon());});
      return (it == points.end() || it == points.begin()) ? 0.0 : ((it->C(l,theta) - (it-1)->C(l,theta))/(it->xyz[0] - (it-1)->xyz[0]));
    }
  };

  auto layer_thickness_a_deriv = [&](int l, double x, double theta) -> double {
    // compute the derivative w.r.t theta of the of the l-th layer cumulative thickness
    const double delta = 1e-6;
    return (layer_thickness(l, x, theta + delta/2) - layer_thickness(l, x, theta - delta/2))/delta;
  };

  auto make_inner_section = [&](int i) -> Handle(Geom2d_Curve) {
    // make a single periodic 2d curve defining the profile of the entire cross-section
    // i is the inner thru-section index
    const double& xc   = inner_sections[i].xc;
    const double& zc   = inner_sections[i].zc;
    const double& ry   = inner_sections[i].ry;
    const double& rz   = inner_sections[i].rz;
    const double alpha = (xc-x_in)/(x_out-x_in);
    double theta_loc = alpha*theta_out + (1-alpha)*theta_in;
    double z_theta = rz*std::cos(theta_loc);
    TColgp_Array1OfPnt2d ctrlPoints(1, NbControlPoints);
    for(int j = 0; j < NbControlPoints; j++) {
      double theta_ell = false/*UseChebyshevNodes*/ ? (-M_PI*std::cos(j*M_PI/(NbControlPoints-1))) : (-M_PI + j*2*M_PI/(NbControlPoints-1));
      double y_ell = -ry*std::sin(theta_ell), z_ell = rz*std::cos(theta_ell);
      double y = y_ell, z = zc + ((z_ell > z_theta) ? z_ell : (1-alpha)*z_ell + alpha*z_theta);
      ctrlPoints.SetValue(j+1, gp_Pnt2d(z, -y)); // local x,y axes of the 2d coordinate system correspond to global z,-y axes
    }
    Geom2dAPI_PointsToBSpline curve_maker;
    curve_maker.Init(ctrlPoints, ParType, DegMin, DegMax, Continuity2D, Tol2D);
    if(!curve_maker.IsDone()) std::cerr << "points to bspline #1 is not done\n";
    Handle(Geom2d_BSplineCurve) curve2d = curve_maker.Curve();
    curve2d->SetPeriodic();
    return curve2d;
  };

  auto make_outer_wire = [&](int i) -> TopoDS_Wire {
    // 1. make 2d curves defining the profiles of the top and bottom surfaces at the i-th thru-section
    const double& xc = sections[i].xc;
    TColgp_Array1OfPnt2d ctrlPoints_top(1, NbExteriorControlPoints), ctrlPoints_bot(1, NbExteriorControlPoints);
    for(int j = 0; j < NbExteriorControlPoints; j++) {
      double theta_top = UseChebyshevNodes ?
                         (0 - M_PI/2*std::cos(j*M_PI/(NbExteriorControlPoints-1))) : (-M_PI/2 + j*M_PI/(NbExteriorControlPoints-1));
      double theta_bot = UseChebyshevNodes ?
                         (M_PI + M_PI/2*std::cos(j*M_PI/(NbExteriorControlPoints-1))) : (3*M_PI/2 - j*M_PI/(NbExteriorControlPoints-1));
      std::pair<double,double> p_top = top.coord(x_out-xc, theta_top+M_PI/2);
      std::pair<double,double> p_bot = bot.coord(x_out-xc, theta_bot+M_PI/2);
      double z_top = p_top.second, z_bot = p_bot.second;
      double y_top = p_top.first,  y_bot = p_bot.first;
      ctrlPoints_top.SetValue(j+1, gp_Pnt2d(z_top, -y_top));
      ctrlPoints_bot.SetValue(j+1, gp_Pnt2d(z_bot, -y_bot));
    }
    Geom2dAPI_PointsToBSpline curve_maker_top;
    curve_maker_top.Init(ctrlPoints_top, ParType, DegMin, DegMax, Continuity2D, Tol2D);
    if(!curve_maker_top.IsDone()) std::cerr << "points to bspline #2 is not done\n";
    Handle(Geom2d_BSplineCurve) curve2d_top = curve_maker_top.Curve();

    Geom2dAPI_PointsToBSpline curve_maker_bot;
    curve_maker_bot.Init(ctrlPoints_bot, ParType, DegMin, DegMax, Continuity2D, Tol2D);
    if(!curve_maker_bot.IsDone()) std::cerr << "points to bspline #3 is not done\n";
    Handle(Geom2d_BSplineCurve) curve2d_bot = curve_maker_bot.Curve();

    // 2. use locate function to find points on curves corresponding to end and mid points
    double p1_top = locate_y_curve(curve2d_top, baffle_half_width);
    double p2_top = locate_y_curve(curve2d_top, 0.);
    double p3_top = locate_y_curve(curve2d_top, -baffle_half_width);
    double p1_bot = locate_y_curve(curve2d_bot, baffle_half_width);
    double p3_bot = locate_y_curve(curve2d_bot, -baffle_half_width);

    // 3. convert 2d curves to 3d
    Handle(Geom_Curve) curve_top = GeomAPI::To3d(curve2d_top, gp_Pln(gp_Pnt(xc, 0., 0.), gp_Dir(1., 0., 0.)));
    Handle(Geom_Curve) curve_bot = GeomAPI::To3d(curve2d_bot, gp_Pln(gp_Pnt(xc, 0., 0.), gp_Dir(1., 0., 0.)));

    // 4. make edges and wire
    BRepBuilderAPI_MakeWire wire_maker;
    TopoDS_Edge edge12_top = BRepBuilderAPI_MakeEdge(curve_top, p1_top, p2_top);
    TopoDS_Edge edge23_top = BRepBuilderAPI_MakeEdge(curve_top, p2_top, p3_top);
    TopoDS_Edge edge13_bot = BRepBuilderAPI_MakeEdge(curve_bot->Reversed(), (1-p1_bot), (1-p3_bot));

    TopTools_IndexedMapOfShape top23_vertexMap, bot13_vertexMap;
    TopExp::MapShapes(edge23_top, TopAbs_VERTEX, top23_vertexMap);
    TopExp::MapShapes(edge13_bot, TopAbs_VERTEX, bot13_vertexMap);
    TopoDS_Edge edge33 = BRepBuilderAPI_MakeEdge(TopoDS::Vertex(bot13_vertexMap(2)), TopoDS::Vertex(top23_vertexMap(1)));
    TopTools_IndexedMapOfShape top12_vertexMap;
    TopExp::MapShapes(edge12_top, TopAbs_VERTEX, top12_vertexMap);
    TopoDS_Edge edge11 = BRepBuilderAPI_MakeEdge(TopoDS::Vertex(top12_vertexMap(2)), TopoDS::Vertex(bot13_vertexMap(1)));

    wire_maker.Add(edge12_top);
    wire_maker.Add(edge11);
    wire_maker.Add(edge13_bot);
    wire_maker.Add(edge33);
    wire_maker.Add(edge23_top);
    if(!wire_maker.IsDone()) std::cerr << "make wire is not done\n";

    return wire_maker.Wire();
  };

  auto make_inner_wire = [&](double xc, Handle(Geom2d_Curve) curve2d) -> TopoDS_Wire {
    // 1. convert 2d curve to a 3d curve in the yz plane
    Handle(Geom_Curve) curve = GeomAPI::To3d(curve2d, gp_Pln(gp_Pnt(xc, 0., 0.), gp_Dir(1., 0., 0.)));

    // 2. make edge and wire
    BRepBuilderAPI_MakeWire wire_maker;
    double pFirst = curve->FirstParameter(), pLast = curve->LastParameter();
    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve, pFirst, pLast);
    wire_maker.Add(edge);
    if(!wire_maker.IsDone()) std::cerr << "make wire is not done\n";
    return wire_maker.Wire();
  };

  auto make_wire = [&](double xc, double zc, Handle(Geom2d_Curve) curve2d) -> TopoDS_Wire {
    // 1. convert 2d curve to a 3d curve in the yz plane
    Handle(Geom_Curve) curve = GeomAPI::To3d(curve2d, gp_Pln(gp_Pnt(xc, 0., 0.), gp_Dir(1., 0., 0.)));

    // 2. make edges and wire
    BRepBuilderAPI_MakeWire wire_maker;
    double pFirst = curve->FirstParameter(), pLast = curve->LastParameter();
    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve, pFirst, pLast);
    wire_maker.Add(edge);
    if(!wire_maker.IsDone()) std::cerr << "make wire is not done\n";

    return wire_maker.Wire();
  };

  auto split_shape_x = [&](const TopoDS_Shape& myShape, double xc) -> TopoDS_Shape {
    // split a surface into pieces to facilitate placement of baffles
    TopoDS_Compound aNewCompound;
    BRep_Builder aBld;
    aBld.MakeCompound(aNewCompound);
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(myShape, TopAbs_FACE, faceMap);
    for(int faceIndex = 1; faceIndex <= faceMap.Extent(); ++faceIndex) {
      const TopoDS_Face& face = TopoDS::Face(faceMap(faceIndex));

      double umin, umax, vmin, vmax;
      BRepTools::UVBounds(face, umin, umax, vmin, vmax);
      Handle_Geom_Surface aSurface1 = BRep_Tool::Surface(face);
      if(aSurface1->Value((umin+umax)/2,vmin).X() < xc && aSurface1->Value((umin+umax)/2,vmax).X() > xc) {
        gp_Pln gp_pln(gp_Pnt(xc, 0., 0.), gp_Dir(1., 0., 0.));
        Handle_Geom_Surface aSurface2 = new Geom_Plane(gp_pln);
        GeomAPI_IntSS aInterSS(aSurface1, aSurface2, IntSSTol);
        if(!aInterSS.IsDone()) std::cerr << "Error: intersection is not done\n";

        BRepFeat_SplitShape aSplitter(face);
        for(int i = 0; i < aInterSS.NbLines(); ++i) {
          Handle(Geom_Curve) aCurve = aInterSS.Line(i+1);
          TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(aCurve);
          aSplitter.Add(edge, face);
        }
        aSplitter.Build();
        if(!aSplitter.IsDone()) std::cerr << "Error: splitting is not done\n";
        TopoDS_Shape splitShape = aSplitter.Shape();
        aBld.Add(aNewCompound, splitShape);
      }
      else {
        aBld.Add(aNewCompound, face);
      }
    }
    return aNewCompound;
  };

  auto split_shape_y = [&](const TopoDS_Shape& myShape) -> TopoDS_Shape {
    // split a surface into pieces to facilitate placement of stringers
    TopoDS_Compound aNewCompound;
    BRep_Builder aBld;
    aBld.MakeCompound(aNewCompound);
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(myShape, TopAbs_FACE, faceMap);
    for(int faceIndex = 1; faceIndex <= faceMap.Extent(); ++faceIndex) {
      const TopoDS_Face& face = TopoDS::Face(faceMap(faceIndex));

      double umin, umax, vmin, vmax;
      BRepTools::UVBounds(face, umin, umax, vmin, vmax);
      Handle_Geom_Surface aSurface1 = BRep_Tool::Surface(face);
      {
        gp_Pln gp_pln(gp_Pnt(0., 0., 0.), gp_Dir(0., 1., 0.));
        Handle_Geom_Surface aSurface2 = new Geom_Plane(gp_pln);
        GeomAPI_IntSS aInterSS(aSurface1, aSurface2, IntSSTol);
        if(!aInterSS.IsDone()) std::cerr << "Error: intersection is not done\n";

        BRepFeat_SplitShape aSplitter(face);
        for(int i = 0; i < aInterSS.NbLines(); ++i) {
          Handle(Geom_Curve) aCurve = aInterSS.Line(i+1);
          TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(aCurve);
          aSplitter.Add(edge, face);
        }
        aSplitter.Build();
        if(!aSplitter.IsDone()) std::cerr << "Error: splitting is not done\n";
        TopoDS_Shape splitShape = aSplitter.Shape();
        aBld.Add(aNewCompound, splitShape);
      }
    }
    return aNewCompound;
  };

  auto make_inner_shell = [&]() -> TopoDS_Shape {
    // construct the inner surface of the nozzle
    // note: this surface is always smoothed (i.e. not ruled)
    std::cout << "Info    : Making the inner shell surface\n";
    BRepOffsetAPI_ThruSections shell_maker(false, false, Pres3D);
    shell_maker.SetContinuity(Continuity3D);
    for(int i = 0; i < NbInnerThruSections; ++i) {
      const double& xc = inner_sections[i].xc;
      TopoDS_Wire wire = make_inner_wire(xc, make_inner_section(i));
      shell_maker.AddWire(wire);
    }
    shell_maker.CheckCompatibility(Standard_True);
    shell_maker.Build();
    if(!shell_maker.IsDone()) std::cerr << "thru sections #1 failed\n";
    return shell_maker.Shape();
  };

  auto make_outer_shell = [&]() -> TopoDS_Shape {
    // construct the exterior geometry
    std::cout << "Info    : Making the outer shell surface\n";
    BRepOffsetAPI_ThruSections shell_maker(false, false, Pres3D);
    shell_maker.SetContinuity(Continuity3D);
    for(int i = 0; i < NbThruSections; ++i) {
      TopoDS_Wire wire = make_outer_wire(i);
      shell_maker.AddWire(wire);
    }
    shell_maker.CheckCompatibility(Standard_True);
    shell_maker.Build();
    if(!shell_maker.IsDone()) std::cerr << "thru sections #2 failed\n";
    TopoDS_Shape myShape = split_shape_y(shell_maker.Shape());
    {
      /*TopTools_IndexedMapOfShape faceMap;
      TopExp::MapShapes(myShape, TopAbs_FACE, faceMap);
      std::cerr << "make_outer_shell #1, faceMap.Extent() = " << faceMap.Extent() << std::endl;*/
      for(int i = 1; i < vertices.size()-1; ++i) {
        myShape = split_shape_x(myShape, points[vertices[i].p].xyz[0]);
        /*TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(myShape, TopAbs_FACE, faceMap);
        std::cerr << "make_outer_shell #2, faceMap.Extent() = " << faceMap.Extent() << std::endl;*/
      }
    }
    return myShape;
  };

  auto locate_x_surface_offset = [&](Handle(Geom_Surface) s, double xc, double zc, double u, double& v0, int l,
                                     int maxit = 20, double tol = 2e-14, bool verbose = false) -> gp_Pnt {
    // solve s->Value(u,v).X()+Normal(u,v).X()*offset - xc = 0 for unknown v using Newton's method with initial guess v0
    double v = v0, offset, umin, umax, vmin, vmax;
    s->Bounds(umin, umax, vmin, vmax);
    GeomLProp_SLProps p(s, 2, 0.01);
    for(int i = 0; ; ++i) {
      p.SetParameters(u, v);
      offset = layer_thickness(l, p.Value().X(), std::atan2(p.Value().Z()-zc, p.Value().Y()));
      double f = p.Value().X()+p.Normal().X()*offset - xc; // assuming outward pointing normal
      if(verbose) std::cerr << "#3: i = " << i << ", |f| = " << std::abs(f) << ", v = " << v << std::endl;
      if(std::abs(f) < tol) { v0 = v; break; }
      if(i == maxit) { std::cerr << "Error: Newton solver #3 did not converge in " << maxit << " iterations (residual = " << std::abs(f) << ")\n";
                       if(std::isnan(f)) v = v0;
                       break; }
      double NormalSquareMagnitude = p.D1U().SquareMagnitude()*p.D1V().SquareMagnitude() - std::pow(p.D1U().Dot(p.D1V()), 2);
      double NormalMagnitude = std::sqrt(NormalSquareMagnitude);
      gp_Vec dNormaldV = (p.DUV().Crossed(p.D1V())) + (p.D1U().Crossed(p.D2V()));
      double dNormalSquareMagnitudedV = 2*p.D1U().Dot(p.DUV())*p.D1V().SquareMagnitude() +
                                        p.D1U().SquareMagnitude()*2*p.D1V().Dot(p.D2V()) -
                                        2*p.D1U().Dot(p.D1V())*(p.D1U().Dot(p.D2V()) + p.DUV().Dot(p.D1V()));
      gp_Vec dUnitNormaldV = 1/NormalMagnitude*dNormaldV - 0.5*dNormalSquareMagnitudedV/NormalSquareMagnitude*p.Normal();
      double dOffsetdX = layer_thickness_x_deriv(l, p.Value().X(), std::atan2(p.Value().Z()-zc, p.Value().Y()));
      double dOffsetdA = layer_thickness_a_deriv(l, p.Value().X(), std::atan2(p.Value().Z()-zc, p.Value().Y()));
      double dOffsetdV = dOffsetdX*p.D1V().X() + dOffsetdA*((zc-p.Value().Z())*p.D1V().Y() + p.Value().Y()*p.D1V().Z())/
                         (p.Value().Y()*p.Value().Y() + (p.Value().Z()-zc)*(p.Value().Z()-zc));
      double dfdv = p.D1V().X() + (dUnitNormaldV.X()*offset + p.Normal().X()*dOffsetdV);
      v -= f/dfdv;
      if(v < vmin || v > vmax) std::cerr << "Warning: v is out-of-bounds: v = " << v << ", vmin = " << vmin << ", vmax = " << vmax << std::endl;
    }
    v0 = v;
    return p.Value().Translated(offset*p.Normal());
  };

  auto locate_xy_surface_offset = [&](Handle(Geom_Surface) s, double xc, double zc, double& u0, double& v0, int l,
                                      int maxit = 20, double tol = 2e-14, bool verbose = false) -> gp_Pnt {
    // solve [s->Value(u,v).X()+Normal(u,v).X()*offset - xc] = [0] for unknowns u,v using Newton's method with initial guess u0,v0
    //       [s->Value(u,v).Y()+Normal(u,v).X()*offset     ]   [0]
    double u = u0, v = v0, offset, umin, umax, vmin, vmax;
    s->Bounds(umin, umax, vmin, vmax);
    GeomLProp_SLProps p(s, 2, 0.01);
    for(int i = 0; ; ++i) {
      p.SetParameters(u, v);
      offset = layer_thickness(l, p.Value().X(), std::atan2(p.Value().Z()-zc, p.Value().Y()));
      double fx = p.Value().X()+p.Normal().X()*offset - xc; // assuming outward pointing normal
      double fy = p.Value().Y()+p.Normal().Y()*offset;
      double fnorm = std::sqrt(fx*fx+fy*fy);
      if(verbose) std::cerr << "#4: i = " << i << ", |f| = " << fnorm << ", u = " << u << ", v = " << v << std::endl;
      if(fnorm < tol) { u0 = u; v0 = v; break; }
      if(i == maxit) { std::cerr << "Error: Newton solver #4 did not converge in " << maxit << " iterations (residual = " << fnorm << ")\n";
                       if(std::isnan(fnorm)) { v = v0; u = u0; }
                       break; }
      double NormalSquareMagnitude = p.D1U().SquareMagnitude()*p.D1V().SquareMagnitude() - std::pow(p.D1U().Dot(p.D1V()), 2);
      double NormalMagnitude = std::sqrt(NormalSquareMagnitude);
      gp_Vec dNormaldU = (p.D2U().Crossed(p.D1V())) + (p.D1U().Crossed(p.DUV()));
      gp_Vec dNormaldV = (p.DUV().Crossed(p.D1V())) + (p.D1U().Crossed(p.D2V()));
      double dNormalSquareMagnitudedU = 2*p.D1V().Dot(p.DUV())*p.D1U().SquareMagnitude() +
                                        p.D1V().SquareMagnitude()*2*p.D1U().Dot(p.D2U()) -
                                        2*p.D1V().Dot(p.D1U())*(p.D1V().Dot(p.D2U()) + p.DUV().Dot(p.D1U()));
      double dNormalSquareMagnitudedV = 2*p.D1U().Dot(p.DUV())*p.D1V().SquareMagnitude() +
                                        p.D1U().SquareMagnitude()*2*p.D1V().Dot(p.D2V()) -
                                        2*p.D1U().Dot(p.D1V())*(p.D1U().Dot(p.D2V()) + p.DUV().Dot(p.D1V()));
      gp_Vec dUnitNormaldU = 1/NormalMagnitude*dNormaldU - 0.5*dNormalSquareMagnitudedU/NormalSquareMagnitude*p.Normal();
      gp_Vec dUnitNormaldV = 1/NormalMagnitude*dNormaldV - 0.5*dNormalSquareMagnitudedV/NormalSquareMagnitude*p.Normal();
      double dOffsetdX = layer_thickness_x_deriv(l, p.Value().X(), std::atan2(p.Value().Z()-zc, p.Value().Y()));
      double dOffsetdA = layer_thickness_a_deriv(l, p.Value().X(), std::atan2(p.Value().Z()-zc, p.Value().Y()));
      double dOffsetdU = dOffsetdX*p.D1U().X() + dOffsetdA*((zc-p.Value().Z())*p.D1U().Y() + p.Value().Y()*p.D1U().Z())/
                         (p.Value().Y()*p.Value().Y() + (p.Value().Z()-zc)*(p.Value().Z()-zc));
      double dOffsetdV = dOffsetdX*p.D1V().X() + dOffsetdA*((zc-p.Value().Z())*p.D1V().Y() + p.Value().Y()*p.D1V().Z())/
                         (p.Value().Y()*p.Value().Y() + (p.Value().Z()-zc)*(p.Value().Z()-zc));
      double a = p.D1U().X() + (dUnitNormaldU.X()*offset + p.Normal().X()*dOffsetdU); // dfxdu
      double b = p.D1V().X() + (dUnitNormaldV.X()*offset + p.Normal().X()*dOffsetdV); // dfxdv
      double c = p.D1U().Y() + (dUnitNormaldU.Y()*offset + p.Normal().Y()*dOffsetdU); // dfydu
      double d = p.D1V().Y() + (dUnitNormaldV.Y()*offset + p.Normal().Y()*dOffsetdV); // dfydv
      double det = a*d - b*c;
      u -= 1/det*( fx*d - fy*b);
      v -= 1/det*(-fx*c + fy*a);
      if(v < vmin || v > vmax) std::cerr << "Warning: v is out-of-bounds: v = " << v << ", vmin = " << vmin << ", vmax = " << vmax << std::endl;
    }
    return p.Value().Translated(offset*p.Normal());
  };

  auto make_lofted_section = [&](double xc, double zc, TopTools_IndexedMapOfShape &faceMap, int l) -> Handle(Geom2d_Curve) {
    // make a single periodic 2d curve defining the profile of the entire cross-section
    // i is the segment index
    double u, umin, umax, v, vmin, vmax;
    TColgp_Array1OfPnt2d ctrlPoints(1, NbControlPoints);
    int ctrlPointIndex = 1;
    assert(faceMap.Extent() == 1);
    const TopoDS_Face &inner_face = TopoDS::Face(faceMap(1));
    Handle_Geom_Surface surface = BRep_Tool::Surface(inner_face);
    BRepTools::UVBounds(inner_face, umin, umax, vmin, vmax);

    double utop = umin + 1*(umax-umin)/2;
    double vtop = locate_x_surface_bisection(surface, xc, utop, 0.01);
    gp_Pnt qtop = locate_xy_surface_offset(surface, xc, zc, utop, vtop, l);

    // set the first control point
    ctrlPoints.SetValue(ctrlPointIndex++, gp_Pnt2d(qtop.Z(), 0.));

    // set the second to the second last control points
    v = vtop;
    for(int k = 1; k < (NbControlPoints-1); ++k) {
      u = false/*UseChebyshevNodes*/ ?
           (utop+(umax-umin)/2 + (umin-umax)/2*std::cos(k*M_PI/(NbControlPoints-1))) : (utop + k*(umax-umin)/(NbControlPoints-1));
      if(u > umax) u -= umax;
      gp_Pnt q = locate_x_surface_offset(surface, xc, zc, u, v, l);
      ctrlPoints.SetValue(ctrlPointIndex++, gp_Pnt2d(q.Z(), -q.Y()));
    }

    // set the last control point (same as first)
    ctrlPoints.SetValue(ctrlPointIndex++, gp_Pnt2d(qtop.Z(), 0.));
   
    Geom2dAPI_PointsToBSpline curve_maker;
    curve_maker.Init(ctrlPoints, ParType, DegMin, DegMax, Continuity2D, Tol2D);
    if(!curve_maker.IsDone()) std::cerr << "points to bspline #4 is not done\n";
    Handle(Geom2d_BSplineCurve) curve2d = curve_maker.Curve();
    curve2d->SetPeriodic();
    return curve2d;
  };

  auto make_lofted_shell = [&](TopoDS_Shape &inner_shell, int l) -> TopoDS_Shape {
    // loft inner shell in the direction of its outward-pointing normal
    std::cout << "Info    : Making lofted shell surface #" << l+1 << std::endl;
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(inner_shell, TopAbs_FACE, faceMap);
    BRepOffsetAPI_ThruSections shell_maker(false, false, Pres3D);
    shell_maker.SetContinuity(Continuity3D);
    for(int i = 0; i < NbThruSections; ++i) {
      const double& xc = sections[i].xc;
      const double& zc = sections[i].zc;
      TopoDS_Wire wire = make_wire(xc, zc, make_lofted_section(xc, zc, faceMap, l));
      shell_maker.AddWire(wire);
    }
    shell_maker.CheckCompatibility(Standard_True);
    shell_maker.Build();
    if(!shell_maker.IsDone()) std::cerr << "thru sections #3 failed\n";
    TopoDS_Shape myShape = split_shape_y(shell_maker.Shape());
    {
      /*TopTools_IndexedMapOfShape faceMap;
      TopExp::MapShapes(myShape, TopAbs_FACE, faceMap);
      std::cerr << "make_lofted_shell #1, faceMap.Extent() = " << faceMap.Extent() << std::endl;*/
      for(int i = 1; i < vertices.size()-1; ++i) {
        myShape = split_shape_x(myShape, points[vertices[i].p].xyz[0]);
        /*TopTools_IndexedMapOfShape faceMap;
        TopExp::MapShapes(myShape, TopAbs_FACE, faceMap);
        std::cerr << "make_lofted_shell #2, i = " << i << ", faceMap.Extent() = " << faceMap.Extent() << std::endl;*/
      }
    }
    return myShape;
  };

  auto add_face_to_model = [&](const TopoDS_Face &face, int physicalTag) -> void {
    // add a face to the gmsh model
    GFace *gface = gm->getOCCInternals()->addFaceToModel(gm, face);
    gface->addPhysicalEntity(physicalTag);
  };

  auto add_lofted_shell_face_to_model = [&](const TopoDS_Face &face, int physicalTag, int nbPointsTransfinite02,
                                            int nbPointsTransfinite13) -> void {
    // add a lofted shell face to the gmsh model and set its mesh parameters
    // nbPointsTransfinite02 is the number of transfinite points on edges 0 and 2
    // nbPointsTransfinite13 is the number of transfinite points on edges 1 and 3
    GFace *gface = gm->getOCCInternals()->addFaceToModel(gm, face);
    gface->addPhysicalEntity(physicalTag);

    if(physicalTag == 2) {
      surfaceTags.push_back(std::make_pair(gface->geomType(), gface->tag()));
    }

    gface->meshAttributes.method = MeshingMethod;
    gface->meshAttributes.transfiniteArrangement = 2;
    std::list<GEdge*> edges = gface->edges();
    for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++it) {
      (*it)->meshAttributes.method = MeshingMethod;
      (*it)->meshAttributes.coeffTransfinite = 1.0;
      int edgeIndex = std::distance(edges.begin(), it);
      if(edgeIndex == 0 || edgeIndex == 2) {
        (*it)->meshAttributes.nbPointsTransfinite = nbPointsTransfinite02;
      }
      if(edgeIndex == 1 || edgeIndex == 3) {
        (*it)->meshAttributes.nbPointsTransfinite = nbPointsTransfinite13;
      }
      if(physicalTag == 2) {
        surfaceTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(), (*it)->getBeginVertex()->tag()));
        surfaceTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(), (*it)->getEndVertex()->tag()));
        surfaceTags.push_back(std::make_pair((*it)->geomType(), (*it)->tag()));
      }
    }
  };

  auto add_stringer_face_to_model = [&](const TopoDS_Face &face, int physicalTag, int nbPointsTransfinite0) -> void {
    // add a stringer face to the gmsh model and set its mesh parameters
    // nbPointsTransfinite0 is the number of transfinite points on edge 0
    GFace *gface = gm->getOCCInternals()->addFaceToModel(gm, face);
    gface->addPhysicalEntity(physicalTag);

    if(physicalTag == 2) {
      surfaceTags.push_back(std::make_pair(gface->geomType(), gface->tag()));
    }

    gface->meshAttributes.method = MESH_UNSTRUCTURED;
    std::list<GEdge*> edges = gface->edges();
    for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++it) {
      (*it)->meshAttributes.method = MeshingMethod;
      (*it)->meshAttributes.coeffTransfinite = 1.0;
      int edgeIndex = std::distance(edges.begin(), it);
      if(edgeIndex == 0) {
        (*it)->meshAttributes.nbPointsTransfinite = nbPointsTransfinite0;
      }
      else {
        double length = (*it)->length((*it)->getLowerBound(), (*it)->getUpperBound());
        (*it)->meshAttributes.nbPointsTransfinite = std::max(2, int(std::ceil(length/lc))+1);
      }
      if(edgeIndex == 2) {
        boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(), (*it)->getBeginVertex()->tag()));
        boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(), (*it)->getEndVertex()->tag()));
        boundaryTags.push_back(std::make_pair((*it)->geomType(), (*it)->tag()));
      }
    }
  };

  auto add_baffle_face_to_model = [&](const TopoDS_Face &face, int physicalTag, int nbPointsTransfinite0) -> void {
    // add a baffle face to the gmsh model and set its mesh parameters
    // nbPointsTransfinite0 is the number of transfinite points on edge 0
    GFace *gface = gm->getOCCInternals()->addFaceToModel(gm, face);
    gface->addPhysicalEntity(physicalTag);

    if(physicalTag == 2) {
      surfaceTags.push_back(std::make_pair(gface->geomType(), gface->tag()));
    }

    gface->meshAttributes.method = MESH_UNSTRUCTURED;
    std::list<GEdge*> edges = gface->edges();
    for(std::list<GEdge*>::iterator it = edges.begin(); it != edges.end(); ++it) {
      (*it)->meshAttributes.method = MeshingMethod;
      (*it)->meshAttributes.coeffTransfinite = 1.0;
      int edgeIndex = std::distance(edges.begin(), it);
      if(edgeIndex == 0) {
        (*it)->meshAttributes.nbPointsTransfinite = nbPointsTransfinite0;
      }
      else {
        double length = (*it)->length((*it)->getLowerBound(), (*it)->getUpperBound());
        (*it)->meshAttributes.nbPointsTransfinite = std::max(2, int(std::ceil(length/lc))+1);
      }
      if(edgeIndex == 2 || edgeIndex == 3 || edgeIndex == 4) {
        boundaryTags.push_back(std::make_pair((*it)->getBeginVertex()->geomType(), (*it)->getBeginVertex()->tag()));
        boundaryTags.push_back(std::make_pair((*it)->getEndVertex()->geomType(), (*it)->getEndVertex()->tag()));
        boundaryTags.push_back(std::make_pair((*it)->geomType(), (*it)->tag()));
      }
    }
  };

  auto add_faces_to_model = [&](const TopoDS_Shape &shape, int physicalTag) -> void {
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(shape, TopAbs_FACE, faceMap);
    for(int i = 0; i < faceMap.Extent(); ++i) {
      const TopoDS_Face &face = TopoDS::Face(faceMap(i+1));
      add_face_to_model(face, physicalTag);
    }
  };

  auto add_lofted_shell_to_model = [&](const TopoDS_Shape &shape, int physicalTag) -> void {
    TopTools_IndexedMapOfShape faceMap;
    TopExp::MapShapes(shape, TopAbs_FACE, faceMap);
    for(int i = 0; i < NbSegments; ++i) {
      for(int j = 0; j < NbPanels; ++j) {
        int faceIndex = j*NbSegments+i+1;
        assert(faceIndex <= faceMap.Extent());
        const TopoDS_Face &face = TopoDS::Face(faceMap(faceIndex));
        int nbPointsTransfinite02 = (i == NbSegments-1 && j == 0) ? segments[i].mn : segments[i].nn;
        int nbPointsTransfinite13 = (i == NbSegments-1 && j == 0) ? segments[i].nn : segments[i].mn;
        add_lofted_shell_face_to_model(face, physicalTag, nbPointsTransfinite02, nbPointsTransfinite13);
      }
    }
  };

  auto add_solid_to_model = [&](const TopoDS_Solid &solid, int physicalTag, int i, int j, int l) -> void {
    // add a solid to the gmsh model and set the mesh parameters
    // i is the segment index, l is the layer index
    gm->getOCCInternals()->buildShapeFromLists(solid);
    gm->getOCCInternals()->buildLists();
    gm->getOCCInternals()->buildGModel(gm);
    GRegion *gregion = gm->getOCCInternals()->getRegionForOCCShape(gm, solid);
    gregion->addPhysicalEntity(physicalTag);

    gregion->cast2Region()->meshAttributes.method = MeshingMethod;
    gregion->cast2Region()->meshAttributes.recombine3D = 1;
    std::list<GFace*> faces = gregion->cast2Region()->faces();
    for(std::list<GFace*>::iterator it = faces.begin(); it != faces.end(); ++ it) {
      (*it)->meshAttributes.method = MeshingMethod;
      (*it)->meshAttributes.recombine = 1;
      int faceIndex = std::distance(faces.begin(),it);
      std::list<GEdge*> edges = (*it)->edges();
      for(std::list<GEdge*>::iterator it2 = edges.begin(); it2 != edges.end(); ++ it2) {
        (*it2)->meshAttributes.method = MeshingMethod;
        (*it2)->meshAttributes.coeffTransfinite = 1.0;
        int edgeIndex = std::distance(edges.begin(),it2);
        if(i == NbSegments-1 && j == 0) {
          if(((faceIndex == 0 || faceIndex == 2) && (edgeIndex == 0 || edgeIndex == 2)) ||
             ((faceIndex == 1 || faceIndex == 5) && (edgeIndex == 0 || edgeIndex == 2))) {
            (*it2)->meshAttributes.nbPointsTransfinite = segments[i].mn; // circumferential edge
          }
          if(((faceIndex == 0 || faceIndex == 2) && (edgeIndex == 1 || edgeIndex == 3)) ||
             ((faceIndex == 3 || faceIndex == 4) && (edgeIndex == 0 || edgeIndex == 2))) {
            (*it2)->meshAttributes.nbPointsTransfinite = segments[i].nn; // longitudinal edge
          }
        }
        else {
          if(((faceIndex == 0 || faceIndex == 2) && (edgeIndex == 0 || edgeIndex == 2)) ||
             ((faceIndex == 1 || faceIndex == 4) && (edgeIndex == 0 || edgeIndex == 2))) {
            (*it2)->meshAttributes.nbPointsTransfinite = segments[i].nn; // longitudinal edge
          }
          if(((faceIndex == 0 || faceIndex == 2) && (edgeIndex == 1 || edgeIndex == 3)) || 
             ((faceIndex == 3 || faceIndex == 5) && (edgeIndex == 0 || edgeIndex == 2))) {
            (*it2)->meshAttributes.nbPointsTransfinite = segments[i].mn; // circumferential edge
          }
        }
        if(((faceIndex == 1 || faceIndex == 4) && (edgeIndex == 1 || edgeIndex == 3)) || 
           ((faceIndex == 3 || faceIndex == 5) && (edgeIndex == 1 || edgeIndex == 3))) {
          (*it2)->meshAttributes.nbPointsTransfinite = (l < 2) ? segments[i].tn[l] : segments[i].ln; // radial edge (i.e. thru thickness)
        }
        if(physicalTag == 0 && faceIndex == 0) {
          interiorBoundaryTags.push_back(std::make_pair((*it2)->getBeginVertex()->geomType(), (*it2)->getBeginVertex()->tag()));
          interiorBoundaryTags.push_back(std::make_pair((*it2)->getEndVertex()->geomType(), (*it2)->getEndVertex()->tag()));
          interiorBoundaryTags.push_back(std::make_pair((*it2)->geomType(), (*it2)->tag()));
        }
        if(physicalTag == 3 && faceIndex == 2) {
          exteriorBoundaryTags.push_back(std::make_pair((*it2)->getBeginVertex()->geomType(), (*it2)->getBeginVertex()->tag()));
          exteriorBoundaryTags.push_back(std::make_pair((*it2)->getEndVertex()->geomType(), (*it2)->getEndVertex()->tag()));
          exteriorBoundaryTags.push_back(std::make_pair((*it2)->geomType(), (*it2)->tag()));
        }
        if(physicalTag == 0 && i == 0 && faceIndex == 3) {
          cmcBoundaryTags.push_back(std::make_pair((*it2)->getBeginVertex()->geomType(), (*it2)->getBeginVertex()->tag()));
          cmcBoundaryTags.push_back(std::make_pair((*it2)->getEndVertex()->geomType(), (*it2)->getEndVertex()->tag()));
          cmcBoundaryTags.push_back(std::make_pair((*it2)->geomType(), (*it2)->tag()));
        }
      }
      if(physicalTag == 0 && faceIndex == 0) { // inside of cmc layer
        interiorBoundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
        (*it)->addPhysicalEntity(physicalTag);
      }
      if(physicalTag == 3 && faceIndex == 2) { // outside of load layer
        exteriorBoundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
        (*it)->addPhysicalEntity(physicalTag);
      }
      if(physicalTag == 1 && faceIndex == 0) { // inside of load layer
        (*it)->addPhysicalEntity(physicalTag);
      }
      if(physicalTag == 0 && i == 0 && faceIndex == 3) { // cmc layer inlet
        cmcBoundaryTags.push_back(std::make_pair((*it)->geomType(),(*it)->tag()));
      }
    }
  };

  auto make_face = [&](const TopoDS_Edge &e1, const TopoDS_Edge &e2) -> TopoDS_Face {
    TopTools_IndexedMapOfShape inner_vertexMap, outer_vertexMap;
    TopExp::MapShapes(e1, TopAbs_VERTEX, inner_vertexMap);
    TopExp::MapShapes(e2, TopAbs_VERTEX, outer_vertexMap);

    Handle_TopTools_HSequenceOfShape edgesWire = new TopTools_HSequenceOfShape;
    Handle_TopTools_HSequenceOfShape myWires = new TopTools_HSequenceOfShape;

    edgesWire->Append(e1);
    edgesWire->Append(e2);
    for(int j = 0; j < 2; ++j) {
      const TopoDS_Vertex &v1 = TopoDS::Vertex(inner_vertexMap(j+1));
      const TopoDS_Vertex &v2 = TopoDS::Vertex(outer_vertexMap(j+1));
      BRepBuilderAPI_MakeEdge edge_maker(v1, v2);
      edgesWire->Append(edge_maker.Edge());
    }

    ShapeAnalysis_FreeBounds::ConnectEdgesToWires(edgesWire, ConnectEdgesToWiresTol, false, myWires);
    assert(myWires->Length() == 1);
    BRepBuilderAPI_MakeFace face_maker(TopoDS::Wire(myWires->Value(1)), true);
    if(!face_maker.IsDone()) std::cerr << "make face failed\n";
    return face_maker.Face();
  };

  auto make_baffle_face = [&](const TopoDS_Edge &e1, const TopoDS_Edge &e2,
                              const TopoDS_Edge &e3, const TopoDS_Edge &e4) -> TopoDS_Face {
    TopTools_IndexedMapOfShape inner_vertexMap, first_outer_vertexMap, last_outer_vertexMap;
    TopExp::MapShapes(e1, TopAbs_VERTEX, inner_vertexMap);
    TopExp::MapShapes(e2, TopAbs_VERTEX, first_outer_vertexMap);
    TopExp::MapShapes(e4, TopAbs_VERTEX, last_outer_vertexMap);

    Handle_TopTools_HSequenceOfShape edgesWire = new TopTools_HSequenceOfShape;
    Handle_TopTools_HSequenceOfShape myWires = new TopTools_HSequenceOfShape;
 
    edgesWire->Append(e1);
    edgesWire->Append(e2);
    edgesWire->Append(e3);
    edgesWire->Append(e4);
    for(int j = 0; j < 2; ++j) {
      const TopoDS_Vertex &v1 = TopoDS::Vertex(inner_vertexMap(j+1));
      const TopoDS_Vertex &v2 = (j == 0) ? TopoDS::Vertex(first_outer_vertexMap(1)) : TopoDS::Vertex(last_outer_vertexMap(2));
      BRepBuilderAPI_MakeEdge edge_maker(v1,v2);
      edgesWire->Append(edge_maker.Edge());
    }
 
    ShapeAnalysis_FreeBounds::ConnectEdgesToWires(edgesWire, ConnectEdgesToWiresTol, false, myWires);
    assert(myWires->Length() == 1);
    BRepBuilderAPI_MakeFace face_maker(TopoDS::Wire(myWires->Value(1)), true);
    if(!face_maker.IsDone()) std::cerr << "make face failed\n";
    return face_maker.Face();
  };

  auto add_stringers_and_baffles_to_model = [&](TopoDS_Shape &inner_shell, TopoDS_Shape &outer_shell) -> void {
    std::cout << "Info    : Making the stringers and baffles\n";
    TopTools_IndexedMapOfShape inner_faceMap, outer_faceMap;
    TopExp::MapShapes(inner_shell, TopAbs_FACE, inner_faceMap);
    TopExp::MapShapes(outer_shell, TopAbs_FACE, outer_faceMap);
    for(int i = 0, baffle_index = 0; i < NbSegments; ++i) {
      for(int j = 0; j < NbPanels; ++j) {
        int innerFaceIndex = j*NbSegments+i+1;
        int outerFaceIndex1 = (3*j)*NbSegments+i+1;
        int outerFaceIndex2 = (3*j+1)*NbSegments+i+1;
        int outerFaceIndex3 = (3*j+2)*NbSegments+i+1;
        assert(innerFaceIndex <= inner_faceMap.Extent()); 
        assert(outerFaceIndex1 <= outer_faceMap.Extent());
        assert(outerFaceIndex2 <= outer_faceMap.Extent());
        assert(outerFaceIndex3 <= outer_faceMap.Extent());
        TopTools_IndexedMapOfShape inner_edgeMap, outer_edgeMap1, outer_edgeMap2, outer_edgeMap3;
        TopExp::MapShapes(inner_faceMap(innerFaceIndex), TopAbs_EDGE, inner_edgeMap);
        TopExp::MapShapes(outer_faceMap(outerFaceIndex1), TopAbs_EDGE, outer_edgeMap1);
        TopExp::MapShapes(outer_faceMap(outerFaceIndex2), TopAbs_EDGE, outer_edgeMap2);
        TopExp::MapShapes(outer_faceMap(outerFaceIndex3), TopAbs_EDGE, outer_edgeMap3);
        if(vertices[i].wb > 0) { // baffle at i-th vertex
          int edgeIndex = (i == NbSegments-1) ? 4 : 2;
          if(i == 0) {
            TopoDS_Face face = make_baffle_face(TopoDS::Edge(inner_edgeMap(edgeIndex)),  TopoDS::Edge(outer_edgeMap1(edgeIndex)),
                                                TopoDS::Edge(outer_edgeMap2(edgeIndex)), TopoDS::Edge(outer_edgeMap3(edgeIndex)));
            add_baffle_face_to_model(face, 5+baffle_index, segments[i].mn);
          }
          else {
            TopoDS_Face face = make_baffle_face(TopoDS::Edge(inner_edgeMap(edgeIndex)),  TopoDS::Edge(outer_edgeMap3(edgeIndex)),
                                                TopoDS::Edge(outer_edgeMap2(edgeIndex)), TopoDS::Edge(outer_edgeMap1(edgeIndex)));
            add_baffle_face_to_model(face, 5+baffle_index, segments[i].mn);
          }
        }
        { // stringer
          int edgeIndex = (i == NbSegments-1) ? 3 : 1;
          if(i == NbSegments-1 && j == 0) { // XXX in first panel of the last segment the edge numberings of the inner and outer faces do not match
            TopoDS_Face face = make_face(TopoDS::Edge(inner_edgeMap(2)), TopoDS::Edge(outer_edgeMap1(edgeIndex)));
            add_stringer_face_to_model(face, 4, segments[i].nn);
          }
          else {
            TopoDS_Face face = make_face(TopoDS::Edge(inner_edgeMap(edgeIndex)), TopoDS::Edge(outer_edgeMap1(edgeIndex)));
            add_stringer_face_to_model(face, 4, segments[i].nn);
          }
        }
      }
      if(vertices[i].wb > 0) baffle_index++;
    }
  };

  auto make_solid = [&](const TopoDS_Face &inner_face, const TopoDS_Face &outer_face, int &err) -> TopoDS_Solid {
    BRepBuilderAPI_Sewing shell_maker;

    // add the inner and outer faces to the set of shapes to be sewed
    shell_maker.Add(inner_face);
    shell_maker.Add(outer_face);

    // construct the faces connecting adjacent edges of inner and outer faces, and add them to the set of shapes to be sewed
    TopTools_IndexedMapOfShape inner_edgeMap, outer_edgeMap;
    TopExp::MapShapes(inner_face, TopAbs_EDGE, inner_edgeMap);
    TopExp::MapShapes(outer_face, TopAbs_EDGE, outer_edgeMap);
    for(int edgeIndex = 1; edgeIndex <= 4; ++edgeIndex) {
      const TopoDS_Edge &inner_edge = TopoDS::Edge(inner_edgeMap(edgeIndex));
      const TopoDS_Edge &outer_edge = TopoDS::Edge(outer_edgeMap(edgeIndex));
      shell_maker.Add(make_face(inner_edge, outer_edge));
    }

    // construct the external boundary of the solid by sewing faces
    shell_maker.Perform();
    if(shell_maker.NbFreeEdges() > 0) {
      std::cerr << "Error: free edges detected in make_solid\n";
      shell_maker.Dump();
    }
    else err = 0;
    TopoDS_Shape shellShape = shell_maker.SewedShape();

    // construct the solid
    BRepBuilderAPI_MakeSolid solid_maker;
    TopTools_IndexedMapOfShape shellMap;
    TopExp::MapShapes(shellShape, TopAbs_SHELL, shellMap);
    for(int shellIndex = 1; shellIndex <= shellMap.Extent(); ++shellIndex) {
      const TopoDS_Shell& shell = TopoDS::Shell(shellMap(shellIndex));
      solid_maker.Add(shell);
    }
    TopoDS_Shape solidShape = solid_maker.Solid();
    return TopoDS::Solid(solidShape);
  };

  auto add_regions_to_model = [&](const TopoDS_Shape &inner_shell, const TopoDS_Shape &outer_shell, int physicalTag, int l) -> void {
    // construct solids using the inner and outer shells and add them to the gmsh model and set the mesh parameters
    std::cout << "Info    : Making solid region #" << physicalTag+1 << std::endl;
    TopTools_IndexedMapOfShape inner_faceMap, outer_faceMap;
    TopExp::MapShapes(inner_shell, TopAbs_FACE, inner_faceMap);
    TopExp::MapShapes(outer_shell, TopAbs_FACE, outer_faceMap);
    for(int i = 0; i < NbSegments; ++i) {
      for(int j = 0; j < NbPanels; ++j) {
        int faceIndex = j*NbSegments+i+1;
        assert(faceIndex <= std::min(inner_faceMap.Extent(), outer_faceMap.Extent()));
        const TopoDS_Face &inner_face = TopoDS::Face(inner_faceMap(faceIndex));
        const TopoDS_Face &outer_face = TopoDS::Face(outer_faceMap(faceIndex));
        int err;
        TopoDS_Solid solid = make_solid(inner_face, outer_face, err);
        if(err == 0) add_solid_to_model(solid, physicalTag, i, j, l);
      }
    }
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 2. do some pre-processing of the input parameters
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Handle(Geom2d_Curve) centerline = make_curve(controlpoints[0], knots[0], mult[0]);
  Handle(Geom2d_Curve) majoraxis  = make_curve(controlpoints[1], knots[1], mult[1]);
  Handle(Geom2d_Curve) minoraxis  = make_curve(controlpoints[2], knots[2], mult[2]);
  x_in = centerline->Value(centerline->FirstParameter()).X();
  x_out = centerline->Value(centerline->LastParameter()).X();
  finalize_points(centerline);
  make_sections(centerline);
  make_inner_sections(centerline, majoraxis, minoraxis);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 3. construct the geometry using OpenCascade
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TopoDS_Shape inner_shell, shell0, shell1, shell2, shell3, shell3_copy, shell4, outer_shell;

  // inside of cmc layer
  inner_shell = make_inner_shell();

  // inside of cmc layer
  if(tf != 0) shell0 = make_lofted_shell(inner_shell, -1);

  // outside of cmc layers / inside of air layer
  if(tf != 0) shell1 = make_lofted_shell(inner_shell, 0);

  // outside of air layer / inside of load layer
  if(tf != 0) shell2 = make_lofted_shell(inner_shell, 1);

  // middle of load layer
  shell3 = make_lofted_shell(inner_shell, 2);
  if(tf != 0) shell3_copy = make_lofted_shell(inner_shell, 2);

  // outside of load layer
  if(tf != 0) shell4 = make_lofted_shell(inner_shell, 3);

  // outside of stringers and baffles / interior aircraft cavity boundary
  outer_shell = make_outer_shell();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 4. construct the Gmsh model and output to file
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  GmshInitialize();
  gm = new GModel();
  gm->setFactory("OpenCASCADE");
  gm->createOCCInternals();

  GmshSetOption("General","Terminal", 1.);
  GmshSetOption("Mesh","CharacteristicLengthMin", lc);
  GmshSetOption("Mesh","CharacteristicLengthMax", lc);

  //DEBUG ONLY: add_faces_to_model(inner_shell, 0);        // add inner surface to gmsh model
  add_lofted_shell_to_model(shell3, 2);                    // add load layer (mid-surface) to gmsh model
  //DEBUG ONLY: add_faces_to_model(outer_shell, 0);        // add outer surface to gmsh model
  add_stringers_and_baffles_to_model(shell3, outer_shell); // add stringers and baffles (mid-surface) to gmsh model

  if(tf != 0) {
    add_regions_to_model(shell0, shell1, 0, 0);            // add cmc layer (volume) to gmsh model
    add_regions_to_model(shell1, shell2, 999, 1);          // add air layer (volume) to gmsh model
    add_regions_to_model(shell2, shell3_copy, 1, 2);       // add inner half of load layer (volume) to gmsh model
    add_regions_to_model(shell3_copy, shell4, 3, 2);       // add outer half of load layer (volume) to gmsh model
  }

  gm->writeGEO("nozzle.geo");

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 5. do the meshing and output mesh to file/s
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

  gm->mesh(3);

  double tolerance = CTX::instance()->geom.tolerance;
  gm->removeDuplicateMeshVertices(tolerance);
  gm->writeMSH("nozzle.msh");
  gm->writeMESH("nozzle.mesh");
  //gm->writeSU2("nozzle.su2", false, 1.0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 6. write the AERO-S input files
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  writeAEROS(gm, points, vertices, segments, materials, boundaries, surfaceTags, boundaryTags, lf, "nozzle.aeros", 2, false, 1.0, (tf==0));
  if(tf != 0) {
    writeAEROH(gm, points, vertices, segments, materials, boundaries, interiorBoundaryTags, exteriorBoundaryTags, "nozzle.aeroh", 2);
    writeAEROS2(gm, points, vertices, segments, materials, boundaries, cmcBoundaryTags, lf, "nozzle.aeros.cmc", 2);
  }

  delete gm;
  GmshFinalize();
}

static PyObject *nozzle_generate(PyObject *self, PyObject *args)
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 1. read NOZZLE.txt
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::ifstream fin("NOZZLE.txt");

  int np, nv, nm; double lc; int bf, tf, nl, nlt, lf, vf, sf;
  fin >> np >> nv >> nm >> lc >> bf >> tf >> nl >> nlt >> lf >> vf >> sf;
  if(vf) {
    std::cout << "Number of points: " << np << std::endl;
    std::cout << "Number of vertices: " << nv << std::endl;
    std::cout << "Number of materials: " << nm << std::endl;
    std::cout << "Characteristic length: " << lc << std::endl;
    std::cout << "Boundary flag: " << bf << std::endl;
    std::cout << "Thermal flag: " << tf << std::endl;
    std::cout << "Number of structural layers: " << nl << std::endl;
    std::cout << "Number of thermal layers: " << nlt << std::endl;
    std::cout << "Linear flag: " << lf << std::endl;
    std::cout << "Verbose flag: " << vf << std::endl;
    std::cout << "Stringer flag: " << sf << std::endl;
  }

  std::vector<PointData> points(np);
  for(int i=0; i<np; ++i) {
    PointData &p = points[i];
    p.xyz.resize(3);
    fin >> p.xyz[0] >> p.xyz[1] >> p.dydx; p.xyz[2] = 0.;
    for(int k=0; k<nl; ++k) { double tk; fin >> tk; p.t.push_back(tk); }
    for(int k=0; k<nlt; ++k) { double tt; fin >> tt; p.tt.push_back(tt); }
    fin >> p.ts >> p.ws;
    if(vf) {
      std::cout << "Coordinates of point " << i << ": " << p.xyz[0] << ", " << p.xyz[1] << std::endl;
      std::cout << "dydx at point " << i << ": " << p.dydx << std::endl;
      std::cout << "Thicknesses of outer (load-bearing) layers at point " << i << ": ";
      for(int k=0; k<nl; ++k) std::cout << p.t[k] << " "; std::cout << std::endl;
      std::cout << "Thicknesses of inner (insulating) layers at point " << i << ": ";
      for(int k=0; k<nlt; ++k) std::cout << p.tt[k] << " "; std::cout << std::endl;
      if(p.ws > 0) std::cout << "Thickness and width of stiffeners at point " << i << ": " << p.ts << ", " << p.ws << std::endl;
    }
  }

  std::vector<VertexData> vertices(nv);
  for(int j=0; j<nv; ++j) {
    VertexData &v = vertices[j];
    fin >> v.p >> v.wb >> v.mb >> v.nb >> v.tb;
    if(vf) {
      std::cout << "Point index of vertex " << j << ": " << v.p << std::endl;
      if(v.wb > 0) std::cout << "baffle width and baffle material id at vertex " << j << ": " << v.wb << ", " << v.mb << std::endl;
      else std::cout << "No baffle at vertex " << j << std::endl;
      std::cout << "Number of transfinite points on radial edges of baffle at vertex " << j << ": " << v.nb << std::endl;
      if(v.wb > 0) std::cout << "Thickness of baffle at vertex " << j << ": " << v.tb << std::endl;
    }
  }

  std::vector<SegmentData> segments(nv-1);
  for(int j=0; j<nv-1; ++j) {
    SegmentData &s = segments[j];
    for(int k=0; k<nl; ++k) { int mk; fin >> mk; s.m.push_back(mk); }
    fin >> s.ns >> s.ms >> s.nn >> s.mn >> s.sn;
    for(int k=0; k<nlt; ++k) { int mt; fin >> mt; s.mt.push_back(mt); }
    for(int k=0; k<nlt; ++k) { int tn; fin >> tn; s.tn.push_back(tn); }
    fin >> s.ln;
    if(vf) {
      std::cout << "Material ids of outer (load-bearing) layers in segment " << j << ": ";
      for(int k=0; k<nl; ++k) std::cout << s.m[k] << " "; std::cout << std::endl;
      if(s.ns > -1) std::cout << "number of panels and stiffener material id in segment " << j
                              << ": " << s.ns  << ", " << s.ms << std::endl;
      else std::cout << "No stiffeners in segment " << j << std::endl;
      std::cout << "Number of transfinite points on longitudinal edges in segment " << j << ": " << s.nn << std::endl;
      std::cout << "Number of transfinite points on circumferential edges in segment " << j << ": " << s.mn << std::endl;
      std::cout << "Number of transfinite points on radial edge of stiffeners in segment " << j << ": " << s.sn << std::endl;
      std::cout << "Material ids of thermal insulating layers in segment " << j << ": ";
      for(int k=0; k<nlt; ++k) std::cout << s.mt[k] << " "; std::cout << std::endl;
      std::cout << "Number of transfinite points through the thickness of the thermal insulating layer in segment " << j << ": ";
      for(int k=0; k<nlt; ++k) std::cout << s.tn[k] << " "; std::cout << std::endl;
      std::cout << "Number of transfinite points through each half of the thickness of the load layer in segment " << j << ": "
                << s.ln << std::endl;
    }
  }

  std::vector<MaterialData> materials(nm);
  for(int k=0; k<nm; ++k) {
    MaterialData &m = materials[k];
    std::string s;
    fin >> s;
    if(s == "ISOTROPIC") {
      m.type = MaterialData::ISOTROPIC;
      fin >> m.E >> m.nu >> m.rho >> m.w >> m.k >> m.h;
      if(vf) std::cout << "Properties of isotropic material " << k << ": " << m.E << ", " << m.nu << ", " << m.rho
                       << ", " << m.w << ", " << m.k << ", " << m.h << std::endl;
    }
    else {
      m.type = MaterialData::ANISOTROPIC;
      fin >> m.E1 >> m.E2 >> m.nu12 >> m.G12 >> m.mu1 >> m.mu2 >> m.rho >> m.w1 >> m.w2 >> m.w12 >> m.k >> m.h;
      if(vf) std::cout << "Properties of anisotropic material " << k << ": " << m.E1 << ", " << m.E2 << ", " << m.nu12 << ", "
                       << m.G12 << ", " << m.mu1 << ", " << m.mu2 << ", " << m.rho << ", " << m.w1 << ", " << m.w2
                       << ", " << m.w12 << ", " << m.k << ", " << m.h << std::endl;
    }
  }

  fin.close();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 2. read BOUNDARY.txt
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::ifstream fin2("BOUNDARY.txt");

  int m;
  fin2 >> m;
  if(vf) std::cout << "Number of boundary locations: " << m << std::endl;

  std::vector<BoundaryData> boundaries(m);
  for(int i=0; i<m; ++i) {
    BoundaryData &b = boundaries[i];
    fin2 >> b.x >> b.P >> b.T >> b.Ta;
    if(vf) std::cout << "Boundary data at location " << i << ": x = " << b.x << ", P = " << b.P << ", T = " << b.T << ", Ta = " << b.Ta << std::endl;
  }

  fin2.close();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 3. read BSPLINE.txt
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::ifstream fin3("BSPLINE.txt");

  int nc[3], nk[3], na1, na2, na3, na4; double sa, ea, bw, ta, to, at, bt, ba, bo, ab, bb;
  fin3 >> nc[0] >> nk[0] >> nc[1] >> nk[1] >> nc[2] >> nk[2] >> sa >> ea >> bw >> ta >> to >> at >> bt >> ba >> bo >> ab >> bb
       >> np >> na1 >> na2 >> na3 >> na4;
  assert(np == points.size());
  if(vf) {
    std::cout << "Number of centerline control points: " << nc[0] << std::endl;
    std::cout << "Number of centerline knots: " << nk[0] << std::endl;
    std::cout << "Number of major axis control points: " << nc[1] << std::endl;
    std::cout << "Number of major axis knots: " << nk[1] << std::endl;
    std::cout << "Number of minor axis control points: " << nc[2] << std::endl;
    std::cout << "Number of minor axis knots: " << nk[2] << std::endl;
    std::cout << "Shovel start angle: " << sa << std::endl;
    std::cout << "Shovel end angle: " << ea << std::endl;
    std::cout << "Baffle half-width: " << bw << std::endl;
    std::cout << "Top angle, offset, semimajor and semiminor axes: " << ta << ", " << to << ", " << at << ", " << bt << std::endl;
    std::cout << "Bottom angle, offset, semimajor and semiminor axes: " << ba << ", " << bo << ", " << ab << ", " << bb << std::endl;
    std::cout << "Number of angular breaks: " << na1 << ", " << na2 << ", " << na3 << ", " << na4 << std::endl;
  }

  std::vector<std::vector<double>> controlpoints[3];
  std::vector<double> knots[3];
  std::vector<int> mult[3];
  for(int j = 0; j < 3; ++j) {
    for(int i = 0; i < nc[j]; ++i) {
      std::vector<double> xy(2);
      fin3 >> xy[0] >> xy[1];
      if(vf) std::cout << "Coordinates of control point " << i << ": " << xy[0] << ", " << xy[1] << std::endl;
      controlpoints[j].push_back(xy);
    }

    for(int i = 0; i < nk[j]; ++i) {
      double u;
      fin3 >> u;
      if(vf) std::cout << "Knot value " << i << ": " << u << std::endl;
      if(knots[j].size() == 0 || u != knots[j].back()) {
        knots[j].push_back(u);
        mult[j].push_back(1);
      }
      else {
        mult[j].back() += 1;
      }
    }
  }
  for(int i = 0; i < np; ++i) {
    for(int j = 0; j < na1; ++j) {
      double xi, ai, ti;
      fin3 >> xi >> ai >> ti;
      assert(xi == points[i].xyz[0]);
      points[i].at1.push_back(std::make_pair(ai*M_PI/180,ti));
      if(vf) std::cout << "Thickness of inner load layer at point " << i+1 << " and angular coordinate " << points[i].at1.back().first
                       << ": " << points[i].at1.back().second << std::endl;
    }
  }
  for(int i = 0; i < np; ++i) {
    for(int j = 0; j < na2; ++j) {
      double xi, ai, ti;
      fin3 >> xi >> ai >> ti;
      assert(xi == points[i].xyz[0]);
      points[i].at2.push_back(std::make_pair(ai*M_PI/180,ti));
      if(vf) std::cout << "Thickness of middle load layer at point " << i+1 << " and angular coordinate " << points[i].at2.back().first
                       << ": " << points[i].at2.back().second << std::endl;
    }
  }
  for(int i = 0; i < np; ++i) {
    for(int j = 0; j < na3; ++j) {
      double xi, ai, ti;
      fin3 >> xi >> ai >> ti;
      assert(xi == points[i].xyz[0]);
      points[i].at3.push_back(std::make_pair(ai*M_PI/180,ti));
      if(vf) std::cout << "Thickness of outer load layer at point " << i+1 << " and angular coordinate " << points[i].at3.back().first
                       << ": " << points[i].at3.back().second << std::endl;
    } 
  }
  for(int i = 0; i < np; ++i) {
    for(int j = 0; j < na4; ++j) {
      double xi, ai, ti;
      fin3 >> xi >> ai >> ti;
      assert(xi == points[i].xyz[0]);
      points[i].at4.push_back(std::make_pair(ai*M_PI/180,ti));
      if(vf) std::cout << "Thickness of thermal layer at point " << i+1 << " and angular coordinate " << points[i].at4.back().first
                       << ": " << points[i].at4.back().second << std::endl;
    } 
  }

  fin3.close();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // 4. build the model
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  generateNozzle(points, vertices, segments, materials, boundaries, lc, tf, lf, vf,
                 controlpoints, knots, mult, sa, ea, bw, ta, to, at, bt, ba, bo, ab, bb);

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

