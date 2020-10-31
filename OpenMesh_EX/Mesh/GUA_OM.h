#pragma once
#ifndef _GUA_OM_H_
#define _GUA_OM_H_

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <windows.h>
#include <gl/gl.h>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "glm/ext.hpp"
using namespace glm;
using namespace std;


#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <utility>
#include <queue>
#include <unordered_map>
#include <deque>
#include <../OpenMesh_EX/Mesh/CustomUtils.h>
using namespace Eigen;


#include <chrono>
using namespace chrono;

struct Face_InnerAngle
{
	double Vertex_Angle[3];
};

class AllPairHarVal_class
{
public:  
	std::vector<double> Har_val;
};

namespace OMT//OpenMesh Triangle mesh
{
	using namespace std;
	/*----------------------------------------------------------------------*/

	/*�w�q�ϥΪ���ǫשM���ݩ�*/
	struct MyTraits : OpenMesh::DefaultTraits
	{
		// let Point and Normal be a vector made from doubles
		typedef OpenMesh::Vec3d Point;
		typedef OpenMesh::Vec3d Normal;

		// add normal property to vertices and faces
		VertexAttributes(OpenMesh::Attributes::Normal);
		FaceAttributes  (OpenMesh::Attributes::Normal);

		// Already defined in OpenMesh::DefaultTraits
		// HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );

		// Uncomment next line to disable attribute PrevHalfedge
		// HalfedgeAttributes( OpenMesh::Attributes::None );
		//
		// or
		//
		// HalfedgeAttributes( 0 );
	};
	/*----------------------------------------------------------------------*/

	/*�w�q�`��type*/
	typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>	    MyMesh	;
	typedef OpenMesh::Vec3d									Vector3d;	//Vec3D type
	typedef MyMesh::Scalar									Scalar	;	//Scalar type
	typedef MyMesh::Point									Point	;	//Point type
	typedef MyMesh::Normal									Normal	;	//Normal type
	typedef MyMesh::VertexHandle							VHandle	;	//VertexHandle type
	typedef MyMesh::HalfedgeHandle							HEHandle;	//HalfedgeHandle type
	typedef MyMesh::EdgeHandle							    EHandle ;	//edgeHandle type
	typedef MyMesh::FaceHandle								FHandle	;	//FaceHandle type
	//-------------Vertex iterators & circulators-------------
	typedef MyMesh::VertexIter								VIter	;	//VertexIter type
	typedef MyMesh::VertexVertexIter						VVIter	;	//VertexVertexIter type
	typedef MyMesh::VertexEdgeIter							VEIter	;	//VertexEdgeIter type
	typedef MyMesh::VertexFaceIter							VFIter	;	//VertexFaceIter type
	typedef MyMesh::EdgeIter								EIter	;	//EdgeIterT	type
	typedef MyMesh::FaceIter								FIter	;	//FaceIter type
	typedef MyMesh::FaceVertexIter							FVIter	;	//FaceVertexIter type
	typedef MyMesh::FaceEdgeIter							FEIter	;	//FaceEdgeIter type
	typedef MyMesh::FaceFaceIter							FFIter	;	//FaceFaceIter type
	typedef MyMesh::VertexOHalfedgeIter						VOHEIter;	//VertexOutHalfEdge type
	typedef MyMesh::ConstVertexVertexIter					CVVIter	;	//ConstVertexVertexIter type
	/*----------------------------------------------------------------------*/

	/*�w�q�B�~��Ƶ��c*/
	using namespace OpenMesh;
	/*----------------------------------------------------------------------*/

	/*�w�q���O*/
	class Model:public MyMesh
	{
	public:
		Model();//constructor
		~Model();//de-constructor
	};
}
/*======================================================================*/

namespace OMP//OpenMesh Polygonal mesh
{
	using namespace std;
	/*----------------------------------------------------------------------*/

	/*�w�q�ϥΪ���ǫשM���ݩ�*/
	struct MyTraits : OpenMesh::DefaultTraits
	{
		// let Point and Normal be a vector made from doubles
		typedef OpenMesh::Vec3d Point;
		typedef OpenMesh::Vec3d Normal;

		// add normal property to vertices and faces
		VertexAttributes(OpenMesh::Attributes::Normal);
		FaceAttributes  (OpenMesh::Attributes::Normal);

		// Already defined in OpenMesh::DefaultTraits
		// HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );

		// Uncomment next line to disable attribute PrevHalfedge
		// HalfedgeAttributes( OpenMesh::Attributes::None );
		//
		// or
		//
		// HalfedgeAttributes( 0 );
	};
	/*----------------------------------------------------------------------*/

	/*�w�q�`��type*/
	typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits>	    MyMesh	;
	typedef OpenMesh::Vec3d									Vector3d;	//Vec3D type
	typedef MyMesh::Scalar									Scalar	;	//Scalar type
	typedef MyMesh::Point									Point	;	//Point type
	typedef MyMesh::Normal									Normal	;	//Normal type
	typedef MyMesh::VertexHandle							VHandle	;	//VertexHandle type
	typedef MyMesh::HalfedgeHandle							HEHandle;	//HalfedgeHandle type
	typedef MyMesh::EdgeHandle							    EHandle ;	//edgeHandle type
	typedef MyMesh::FaceHandle								FHandle	;	//FaceHandle type
	//-------------Vertex iterators & circulators-------------
	typedef MyMesh::VertexIter								VIter	;	//VertexIter type
	typedef MyMesh::VertexVertexIter						VVIter	;	//VertexVertexIter type
	typedef MyMesh::VertexEdgeIter							VEIter	;	//VertexEdgeIter type
	typedef MyMesh::VertexFaceIter							VFIter	;	//VertexFaceIter type
	typedef MyMesh::EdgeIter								EIter	;	//EdgeIterT	type
	typedef MyMesh::FaceIter								FIter	;	//FaceIter type
	typedef MyMesh::FaceVertexIter							FVIter	;	//FaceVertexIter type
	typedef MyMesh::FaceEdgeIter							FEIter	;	//FaceEdgeIter type
	typedef MyMesh::FaceFaceIter							FFIter	;	//FaceFaceIter type
	typedef MyMesh::VertexOHalfedgeIter						VOHEIter;	//VertexOutHalfEdge type
	typedef MyMesh::ConstVertexVertexIter					CVVIter	;	//ConstVertexVertexIter type
	/*----------------------------------------------------------------------*/

	/*�w�q�B�~��Ƶ��c*/
	using namespace OpenMesh;
	//���w�S�O�e�X������Ƶ��c
	struct sp_f 
	{
		FHandle fh;
		float r, g, b;
	};
	//���w�S�O�e�X���I����Ƶ��c
	struct sp_v 
	{
		VHandle vh;
		float r, g, b;
	};
	//���w�t�~�e�X��m����Ƶ��c
	struct sp_p
	{
		Point pt;
		float r, g, b;
	};
	/*----------------------------------------------------------------------*/
	struct srtPath
	{
		std::vector<OMP::VHandle> path;
	};

	/*�w�q���O*/
	class Model:public MyMesh
	{
		//��Ʀ���
	public:
		MyMesh Mesh;													//OpenMesh instance
	private:
		vector< sp_p > sp_p_list;
		vector< sp_v > sp_v_list;
		vector< sp_f > sp_f_list;
		//�禡����
	public:
		Model();//constructor
		~Model();//de-constructor

		bool ReadFile(std::string _fileName);//Ū��mesh���
		bool SaveFile(std::string _fileName);//�x�smesh���

		void Render_solid();			//solid rendering
		void Render_wireframe();		//wireframe rendering

		void RenderSpecifiedPoint();	//�e�X���w��m���I
		void RenderSpecifiedVertex();	//�e�X���w�����I
		void RenderSpecifiedFace();		//�e�X���w����

		void add_sp_p(Point   _p, float _r, float _g, float _b);//���w�B�~�e�X���I
		void add_sp_v(VHandle _v, float _r, float _g, float _b);//���w�B�~�e�X�����I
		void add_sp_f(FHandle _f, float _r, float _g, float _b);//���w�B�~�e�X����
		void clear_sp_p();//�M���B�~�e�X���I
		void clear_sp_v();//�M���B�~�e�X�����I
		void clear_sp_f();//�M���B�~�e�X����

		VHandle addVertex(Point _p);										//�bmodel�W�W�[�s�����I
		FHandle addFace(VHandle _v0, VHandle _v1, VHandle _v2, VHandle _v3);//�bmodel�W�W�[�s����
		void deleteFace(FHandle _f);										//�bmodel�W�R����
		void deleteFace(VHandle _v0, VHandle _v1, VHandle _v2, VHandle _v3);//�bmodel�W�R����

		bool IsVertexVertex( VHandle _vj, VHandle _vi);	//�ˬd�⳻�I�O�_�۾F

		int quad_subdivision1(int _face_id);//1to4 quadrilateral subdivision
		int quad_subdivision2(int _face_id);//1to9 quadrilateral subdivision

	private:
		int findVertex(Point _p);//��X��m_p�O�_�w�g�����I�s�b
	};
}
/*======================================================================*/
class Tri_Mesh:public OMT::Model
{
public:
	Tri_Mesh()
	{
		_deque = LogDeque();
		request_vertex_texcoords2D();
	}

	~Tri_Mesh()
	{
		this->remove_property(cost);
		this->remove_property(QMat);
		this->remove_property(newPoint);
	}

	void loadToBuffer(Tri_Mesh _mesh, std::vector<double> & out_vertices, int & face, std::vector<double> & uv);
	void loadToBufferPatch(std::vector<double> & out_vertices, int & face, std::vector<int> selected, Tri_Mesh & patch);
	void findNearestPoint(Tri_Mesh mesh, std::vector<double> mouse, int face, std::vector<double> &vertex);
	void findNearestVert(Tri_Mesh mesh, std::vector<double> mouse, int face, std::vector<double> &vertex , mat4 MVP , double dis);
	void delVert(VHandle vhandle);
	Matrix4d calculateQ(VertexHandle vhandle);

	void oneRingCollapse(VHandle vhandle);
	bool DetermineConcaveByTwoPoints(std::vector<double> & p1, std::vector<double> & p2, std::vector<double> & vertices);
	void normalizeModel();

	Tri_Mesh simplify(float rate, float threshold = 0);
	Tri_Mesh averageSimplify();

	//-------Edit Flag-------//
    bool                                       Delete_Flag;
	
	int                                        Constrain_num;
	int                                        Boundary_num ;
	OMT::VHandle                               start_vh,end_vh;
	OMT::VHandle                               ExtremeVh[2];
	int                                        PatchType;

	std::vector<OMT::VHandle>                  Pluspt      ;
	std::vector<OMT::VHandle>                  Minuspt     ;
	std::vector<OMT::VHandle>                  Extrme_Pt   ;

	void getUV(std::vector<double> & patchuv, Tri_Mesh patch, float uvRotateAngle);
	void Render_Solid();
	void Render_SolidWireframe();
	void Render_Wireframe();
	void Render_Point();
	void Decimate(int k);
	void Recover(int k);
	void Initialize();
	void Update_Edge(EdgeHandle eh);
	bool equal(Point a, Point b) {
		if ((a - b).length() < 1e-6) {
			return true;
		}
		return false;
	}

private:
	OpenMesh::EPropHandleT<double> cost;
	OpenMesh::VPropHandleT<Matrix4d> QMat;
	OpenMesh::EPropHandleT<Point> newPoint;

	struct CompareCost {
		Tri_Mesh* mesh;
		OpenMesh::EPropHandleT<double>* cost;

		CompareCost(Tri_Mesh* _mesh, OpenMesh::EPropHandleT<double>* _cost) {
			mesh = _mesh;
			cost = _cost;
		}

		bool operator()(const int p1, const int p2) const
		{
			auto edgeHandle1 = mesh->edge_handle(p1);
			auto edgeHandle2 = mesh->edge_handle(p2);
			return mesh->property(*cost, edgeHandle1) < mesh->property(*cost, edgeHandle2);
		}
	};

	struct RollbackInfo {
	public:
		Point p;
		vector<Point> neighborPos;  

		RollbackInfo() {}

		RollbackInfo(int idx, Tri_Mesh* mesh) {
			p = mesh->point(mesh->vertex_handle(idx));
			neighborPos = vector<Point>();
			
			VertexVertexCCWIter vv_it;
			for (vv_it = mesh->vv_ccwbegin(mesh->vertex_handle(idx)); vv_it != mesh->vv_ccwend(mesh->vertex_handle(idx)); ++vv_it) {
				neighborPos.push_back(mesh->point(*vv_it));
			}
		}

	};

	struct DecimationLog {
	public:
		RollbackInfo* collapsePoint;
		RollbackInfo* remainPoint;
		Point remainPointNewPosition;

		DecimationLog(RollbackInfo* info1, RollbackInfo* info2, Point p) {
			collapsePoint = info1;
			remainPoint = info2;
			remainPointNewPosition = p;
		}
	};

	struct LogDeque {
	private:
		deque<DecimationLog*> storage;
		deque<DecimationLog*> decimated;

	public:
		LogDeque() {
			storage = deque<DecimationLog*>();
			decimated = deque<DecimationLog*>();
		}
	
		void pushNewLog(DecimationLog* log) {
			storage.push_front(log);
		}

		void printTop() {
			DecimationLog* log = decimated.front();
		}

		bool canDecimate() {
			return storage.size() > 0;
		}

		bool canRecover() {
			return decimated.size() > 0;
		}

		DecimationLog* toDecimate() {
			if (storage.size() == 0) return nullptr;
			DecimationLog* log = storage.back();
			storage.pop_back();
			decimated.push_front(log);
			// cout << storage.size()<< " " <<decimated.size() << endl;
			return log;
		}

		DecimationLog* toRecover() {
			if (decimated.size() == 0) return nullptr;
			DecimationLog* log = decimated.front();
			decimated.pop_front();
			storage.push_back(log);
			// cout << storage.size() << " " << decimated.size() << endl;
			return log;
		}
	};

	LogDeque _deque;
};

///*======================================================================*/
/*======================================================================*/
bool ReadFile(std::string _fileName, Tri_Mesh *_mesh); //Ū��mesh���
bool SaveFile(std::string _fileName, Tri_Mesh *_mesh); //�x�smesh���
/*��l��view port�]�w�禡*/

#endif

