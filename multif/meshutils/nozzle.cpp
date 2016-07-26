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
    "Generate a mesh for given parameters defined in files NOZZLE.txt and BOUNDARY.txt.";

static PyObject *nozzle_generate(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"generate", nozzle_generate, METH_VARARGS, generate_docstring},
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

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include "MElement.h"
#include "Gmsh.h"
#include "GModel.h"

struct PointData { double x, y; };
struct VertexData { int p, mb; double wb; };
struct SegmentData { int m, ns, ms; double ws; };
struct MaterialData { double E, nu, rho, t, w, Ta; };
struct BoundaryData { double x, P, T; static std::vector<int> tags; };

std::vector<int> BoundaryData::tags;

bool cmp(const BoundaryData &a, const double &b) { return a.x < b; }

void writeNode(MVertex *m, FILE *fp, double scalingFactor)
{
  if(m->getIndex() < 0) return;

  double x1 = m->x() * scalingFactor;
  double y1 = m->y() * scalingFactor;
  double z1 = m->z() * scalingFactor;

  fprintf(fp, "%-8d %-16.9G %-16.9G %-16.9G\n", m->getIndex(), x1, y1, z1);
}

void writeElem(MElement *m, FILE *fp, int elementTagType, int physical)
{
  const char *str = m->getStringForBDF();
  if(!str) return;

  int n = m->getNumVertices();

  int type;
  if(std::strcmp(str,"CHEXA")==0) type = 17;
  else if(std::strcmp(str,"CBAR")==0) type = 6;
  else if(std::strcmp(str,"CPENTA")==0) type = 24;
  else if(std::strcmp(str,"CQUAD4")==0) type = 1515;
  else if(std::strcmp(str,"CTETRA")==0) type = 23;
  else if(std::strcmp(str,"CTRIA3")==0) type = 15;
  else {
    Msg::Error("Element type %s is not supported by AEROS writer", str);
    return;
  }

  if(physical < 0) m->reverse();

  fprintf(fp, "%-8d %-8d", m->getNum(), type);
  for(int i = 0; i < n; i++) {
    fprintf(fp, " %-8d", m->getVertexBDF(i)->getIndex());
  }
  fprintf(fp, "\n");

  if(physical < 0) m->reverse();
}

void writeAttr(MElement *m, FILE *fp, int elementTagType, int elementary, int physical)
{
  int tag = (elementTagType == 3) ? m->getPartition() : (elementTagType == 2) ?
    abs(physical) : elementary;

  fprintf(fp, "%-8d %-8d\n", m->getNum(), tag);
}

void writeDisp(MVertex *m, FILE *fp, double scalingFactor)
{
  if(m->getIndex() < 0) return;

  double x1 = m->x() * scalingFactor;

  if(fabs(x1) < 1e-12) {
    fprintf(fp, "%-8d 1 0.0\n", m->getIndex());
    fprintf(fp, "%-8d 2 0.0\n", m->getIndex());
    fprintf(fp, "%-8d 3 0.0\n", m->getIndex());
    fprintf(fp, "%-8d 4 0.0\n", m->getIndex());
    fprintf(fp, "%-8d 5 0.0\n", m->getIndex());
    fprintf(fp, "%-8d 6 0.0\n", m->getIndex());
  }
}

void writePres(MElement *m, FILE *fp, int elementTagType, int elementary,
               double scalingFactor, const std::vector<BoundaryData> &boundaries)
{
  const char *str = m->getStringForBDF();
  if(!str) return;

  if(std::strcmp(str,"CTRIA3") == 0 &&
     std::find(BoundaryData::tags.begin(), BoundaryData::tags.end(), elementary) != BoundaryData::tags.end()) {

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
}

void writeTemp(MVertex *m, FILE *fp, int elementary, double scalingFactor,
               const std::vector<BoundaryData> &boundaries)
{
  if(m->getIndex() < 0) return;

  if(std::find(BoundaryData::tags.begin(), BoundaryData::tags.end(), elementary) != BoundaryData::tags.end()) {

    double x = m->x() * scalingFactor;

    // interpolate pressure at centroid from boundary data
    std::vector<BoundaryData>::const_iterator it = std::lower_bound(boundaries.begin(), boundaries.end(), x, cmp);
    double tval = (it == boundaries.begin()) ? it->T : ((it-1)->T + (it->T - (it-1)->T)*(x - (it-1)->x)/(it->x - (it-1)->x));

    fprintf(fp, "%-8d %-16.9G\n", m->getIndex(), tval);
  }
}

int writeAEROS(GModel *g,
               const std::vector<MaterialData> &materials,
               const std::vector<BoundaryData> &boundaries,
               const std::string &name,
               int elementTagType=1, bool saveAll=false, double scalingFactor=1.0)
{
  FILE *fp = fopen(name.c_str(), "w");
  if(!fp){
    Msg::Error("Unable to open file '%s'", name.c_str());
    return 0;
  }

  if(g->noPhysicalGroups()) saveAll = true;

  g->indexMeshVertices(saveAll);

  std::vector<GEntity*> entities;
  g->getEntities(entities);

  // nodes
  FILE *fp1 = fopen("GEOMETRY.txt", "w");
  fprintf(fp1, "NODES\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
      writeNode(entities[i]->mesh_vertices[j], fp1, scalingFactor);
  fclose(fp1);

  // elements
  FILE *fp2 = fopen("TOPOLOGY.txt", "w");
  fprintf(fp2, "TOPOLOGY\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++){
      int numPhys = entities[i]->physicals.size();
      if(saveAll || numPhys)
        writeElem(entities[i]->getMeshElement(j),
           fp2, elementTagType,
           numPhys ? entities[i]->physicals[0] : 0);
    }
  fclose(fp2);

  // attributes
  FILE *fp3 = fopen("ATTRIBUTES.txt", "w");
  fprintf(fp3, "ATTRIBUTES\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++){
      int numPhys = entities[i]->physicals.size();
      if(saveAll || numPhys)
        writeAttr(entities[i]->getMeshElement(j),
           fp3, elementTagType, entities[i]->tag(),
           numPhys ? entities[i]->physicals[0] : 0);
    }
  fclose(fp3);

  // material properties
  std::ofstream fout("MATERIAL.txt");
  fout << "MATERIALS\n";
  for(std::vector<MaterialData>::const_iterator it = materials.begin(); it != materials.end(); ++it) {
    fout << std::distance(materials.begin(),it)+1 << " 0.0 " << it->E << " " << it->nu << " " << it->rho
         << " 0.0 0.0 " << it->t << " 0.0 " << it->Ta << " 0.0 " << it->w << " 0.0 0.0 0.0\n";
  }
  fout.close();

  // displacement boundary condition hardcoded to fixed at x=0
  FILE *fp4 = fopen("DISPLACEMENTS.txt", "w");
  fprintf(fp4, "DISPLACEMENTS\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
      writeDisp(entities[i]->mesh_vertices[j], fp4, scalingFactor);
  fclose(fp4);

  // pressure
  FILE *fp5 = fopen("PRESSURES.txt", "w");
  fprintf(fp5, "PRESSURE\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->getNumMeshElements(); j++){
      int numPhys = entities[i]->physicals.size();
      if(saveAll || numPhys)
        writePres(entities[i]->getMeshElement(j),
           fp5, elementTagType, entities[i]->tag(),
           scalingFactor, boundaries);
    }
  fclose(fp5);

  // temperatures
  FILE *fp6 = fopen("TEMPERATURES.txt", "w");
  fprintf(fp6, "TEMPERATURE\n");
  for(unsigned int i = 0; i < entities.size(); i++)
    for(unsigned int j = 0; j < entities[i]->mesh_vertices.size(); j++)
      writeTemp(entities[i]->mesh_vertices[j], fp6, entities[i]->tag(),
                scalingFactor, boundaries);
  fclose(fp6);

  // main file
  fprintf(fp, "STATICS\n");
  fprintf(fp, "sparse\n");
  fprintf(fp, "*\n");
  fprintf(fp, "OUTPUT\n");
  fprintf(fp, "gdisplac \"DISP\" 1\n");
  fprintf(fp, "stressvm \"STRESS\" 1\n");
  fprintf(fp, "*\n");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "GEOMETRY.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TOPOLOGY.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "ATTRIBUTES.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "MATERIAL.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "DISPLACEMENTS.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "PRESSURES.txt");
  fprintf(fp, "INCLUDE \"%s\"\n*\n", "TEMPERATURES.txt");
  fprintf(fp, "END\n");

  fclose(fp);

  return 1;
}

void generateNozzle(const std::vector<std::vector<double> > &points,
                    const std::vector<VertexData> &vertices,
                    const std::vector<SegmentData> &segments,
                    const std::vector<MaterialData> &materials,
                    const std::vector<BoundaryData> &boundaries,
                    double lc)
{
  GmshInitialize();
  GModel *m = new GModel();

  GmshSetOption("General","Terminal", 1.);
  GmshSetOption("Mesh","CharacteristicLengthMin", lc);
  GmshSetOption("Mesh","CharacteristicLengthMax", lc);

  m->setFactory("OpenCASCADE");

  int dimension = 2;

  std::vector<double> p1(3), p2(3); // points defining the axis of revolution
  p1[0] = p1[1] = p1[2] = p2[1] = p2[2] = 0; p2[0] = 1;

  std::vector<std::vector<double> >::const_iterator pointIt1 = points.begin();
  std::vector<VertexData>::const_iterator vertexIt = vertices.begin();
  GVertex *vertex1 = m->addVertex((*pointIt1)[0], (*pointIt1)[1], (*pointIt1)[2], lc);

  // loop over the segments
  for(std::vector<SegmentData>::const_iterator segmentIt = segments.begin(); segmentIt != segments.end(); ++segmentIt, ++vertexIt) {

    std::vector<std::vector<double> >::const_iterator pointIt2 = points.begin()+(vertexIt+1)->p;
    GVertex *vertex2 = m->addVertex((*pointIt2)[0], (*pointIt2)[1], (*pointIt2)[2], lc);
    std::vector<std::vector<double> > controlPoints(pointIt1+1, pointIt2);
    GEdge *edge = m->addBSpline(vertex1, vertex2, controlPoints);

    int mm = segmentIt->m;              // material id of main shell
    int ns = std::max(1,segmentIt->ns); // number of circumferential segments
    double wb = vertexIt->wb;           // width of baffle
    int mb = vertexIt->mb;              // material id of baffle
    double ws = segmentIt->ws;          // width of stiffeners
    int ms = segmentIt->ms;             // material id of stiffeners
    double angle = 2*M_PI/ns;

    GEntity *face;
    for(int i=0; i<ns; ++i) {
      if(i > 0) {
        std::list<GEdge*> edges = face->cast2Face()->edges();
        edge = *(++edges.begin()); // the second edge is the one we need to revolve/extrude!
      }
      //main shell:
      face = m->revolve(edge, p1, p2, angle);
      face->addPhysicalEntity(mm+1);
      BoundaryData::tags.push_back(face->tag());

      std::vector<double> p3(3), p4(3); // points defining radial unit vector
      p3[0] = edge->getBeginVertex()->x(); p3[1] = p3[2] = 0;
      p4[0] = p3[0]; p4[1] = edge->getBeginVertex()->y(); p4[2] = edge->getBeginVertex()->z();
      double r = std::sqrt(p4[1]*p4[1]+p4[2]*p4[2]);
      p4[1] /= r; p4[2] /= r;

      //baffle:
      if(wb > 0) {
        std::vector<double> p5(p4); p5[1] *= wb; p5[2] *= wb; // point defining extrusion vector
        GEntity *baffleEdge = m->extrude(edge->getBeginVertex(), p3, p5);
        GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
        baffleFace->addPhysicalEntity(mb+1);
      }

      //stiffener:
      if(ws > 0) {
        std::vector<double> p6(p4); p6[1] *= ws; p6[2] *= ws; // point defining extrusion vector
        GEntity *stiffenerFace = m->extrude(edge, p3, p6);
        stiffenerFace->addPhysicalEntity(ms+1);
      }

      if(segmentIt+1 == segments.end()) {
        double wb = (vertexIt+1)->wb; // width of baffle
        int mb = (vertexIt+1)->mb;    // material id of baffle
        std::vector<double> p3(3), p4(3);
        p3[0] = edge->getEndVertex()->x(); p3[1] = p3[2] = 0;
        p4[0] = p3[0]; p4[1] = edge->getEndVertex()->y(); p4[2] = edge->getEndVertex()->z();
        double r = std::sqrt(p4[1]*p4[1]+p4[2]*p4[2]);
        p4[1] /= r; p4[2] /= r;

        //final baffle:
        if(wb > 0) {
          std::vector<double> p5(p4); p5[1] *= wb; p5[2] *= wb; // point defining extrusion vector
          GEntity *baffleEdge = m->extrude(edge->getEndVertex(), p3, p5);
          GEntity *baffleFace = m->revolve(baffleEdge, p1, p2, angle);
          baffleFace->addPhysicalEntity(mb+1);
        }
      }
    }

    pointIt1 = pointIt2;
    vertex1 = vertex2;
  }

  m->mesh(dimension);

  double tolerance = lc/50;
  m->removeDuplicateMeshVertices(tolerance);

  m->writeGEO("nozzle.geo");
  m->writeMSH("nozzle.msh");
  writeAEROS(m, materials, boundaries, "nozzle.aeros", 2);

  delete m;
  GmshFinalize();
}

static PyObject *nozzle_generate(PyObject *self, PyObject *args)
{
  std::ifstream fin("NOZZLE.txt");

  int np, nv, nm; double lc;
  fin >> np >> nv >> nm >> lc;
  
  std::vector<std::vector<double> > points;
  for(int i=0; i<np; ++i) {
    std::vector<double> xyz(3);
    fin >> xyz[0] >> xyz[1]; xyz[2] = 0;
    points.push_back(xyz);
  }

  std::vector<VertexData> vertices(nv);
  for(int j=0; j<nv; ++j) {
    VertexData &v = vertices[j];
    fin >> v.p >> v.wb >> v.mb;
  }

  std::vector<SegmentData> segments(nv-1);
  for(int j=0; j<nv-1; ++j) {
    SegmentData &s = segments[j];
    fin >> s.m >> s.ns >> s.ws >> s.ms;
  }

  std::vector<MaterialData> materials(nm);
  for(int k=0; k<nm; ++k) {
    MaterialData &m = materials[k];
    fin >> m.E >> m.nu >> m.rho >> m.t >> m.w >> m.Ta;
  }

  fin.close();

  std::ifstream fin2("BOUNDARY.txt");

  int m;
  fin2 >> m;

  std::vector<BoundaryData> boundaries(m); 
  for(int i=0; i<m; ++i) {
    BoundaryData &b = boundaries[i];
    fin2 >> b.x >> b.P >> b.T;
  }

  fin2.close();

  generateNozzle(points, vertices, segments, materials, boundaries, lc);

  Py_RETURN_NONE;
}
