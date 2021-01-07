#pragma once
#ifndef _GUA_OM_H_
#define _GUA_OM_H_

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

#include <windows.h>
#include <gl/gl.h>
#include <gl/glu.h>
#include <math.h>

// #include "ExpMap/myAxisAlignedBox3d.h"
#ifdef max
	#undef max
#endif
#ifdef min
	#undef min
#endif

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

	/*定義使用的精準度和基本屬性*/
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

	/*定義常用type*/
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

	/*定義額外資料結構*/
	using namespace OpenMesh;
	/*----------------------------------------------------------------------*/

	/*定義類別*/
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

	/*定義使用的精準度和基本屬性*/
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

	/*定義常用type*/
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

	/*定義額外資料結構*/
	using namespace OpenMesh;
	//指定特別畫出面的資料結構
	struct sp_f 
	{
		FHandle fh;
		float r, g, b;
	};
	//指定特別畫出頂點的資料結構
	struct sp_v 
	{
		VHandle vh;
		float r, g, b;
	};
	//指定另外畫出位置的資料結構
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

	/*定義類別*/
	class Model:public MyMesh
	{
		//資料成員
	public:
		MyMesh Mesh;													//OpenMesh instance
	private:
		vector< sp_p > sp_p_list;
		vector< sp_v > sp_v_list;
		vector< sp_f > sp_f_list;
		//函式成員
	public:
		Model();//constructor
		~Model();//de-constructor

		bool ReadFile(std::string _fileName);//讀取mesh資料
		bool SaveFile(std::string _fileName);//儲存mesh資料

		void Render_solid();			//solid rendering
		void Render_wireframe();		//wireframe rendering

		void RenderSpecifiedPoint();	//畫出指定位置的點
		void RenderSpecifiedVertex();	//畫出指定的頂點
		void RenderSpecifiedFace();		//畫出指定的面

		void add_sp_p(Point   _p, float _r, float _g, float _b);//指定額外畫出的點
		void add_sp_v(VHandle _v, float _r, float _g, float _b);//指定額外畫出的頂點
		void add_sp_f(FHandle _f, float _r, float _g, float _b);//指定額外畫出的面
		void clear_sp_p();//清空額外畫出的點
		void clear_sp_v();//清空額外畫出的頂點
		void clear_sp_f();//清空額外畫出的面

		VHandle addVertex(Point _p);										//在model上增加新的頂點
		FHandle addFace(VHandle _v0, VHandle _v1, VHandle _v2, VHandle _v3);//在model上增加新的面
		void deleteFace(FHandle _f);										//在model上刪除面
		void deleteFace(VHandle _v0, VHandle _v1, VHandle _v2, VHandle _v3);//在model上刪除面

		bool IsVertexVertex( VHandle _vj, VHandle _vi);	//檢查兩頂點是否相鄰

		int quad_subdivision1(int _face_id);//1to4 quadrilateral subdivision
		int quad_subdivision2(int _face_id);//1to9 quadrilateral subdivision

	private:
		int findVertex(Point _p);//找出位置_p是否已經有頂點存在
	};
}
/*======================================================================*/
class Tri_Mesh:public OMT::Model
{
public:
	Tri_Mesh()
	{
		add_property(vhandles, "vhandles");
		add_property(UVSets);
		add_property(UVSetIndex);
	}
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


	void Render_Solid();
	void Render_SolidWireframe();
	void Render_Wireframe();
	void Render_Point();

	int Tri_Mesh::findVertex(OMT::Point _p);
	OMT::VPropHandleT< std::vector<OMT::VHandle> >			vhandles;
	OpenMesh::VPropHandleT<std::vector<OMT::Vec2d>>			UVSets;
	OpenMesh::FPropHandleT<int>								UVSetIndex;

	//Exp lib use function
	static const int InvalidID = -1;

	// 特別的向量
	static const OMT::Vec3d ZERO;    // (0,0,0)
	static const OMT::Vec3d UNIT_X;  // (1,0,0)
	static const OMT::Vec3d UNIT_Y;  // (0,1,0)
	static const OMT::Vec3d UNIT_Z;  // (0,0,1)
	static const OMT::Vec3d ONE;     // (1,1,1)

	unsigned int GetVertexCount() {
		return n_vertices();
	};

	void GetVertex(int vID, OMT::Point & vVertex, OMT::Normal * pNormal) const
	{
		OMT::VHandle vh = vertex_handle(vID);
		vVertex = point(vh);
		if (pNormal)
			*pNormal = normal(vh);
	};

	void GetNormal(int vID, OMT::Normal & vNormal) {
		OMT::VHandle vh = vertex_handle(vID);
		vNormal = normal(vh);
	}

	void GetTriangle(int tID, int vTriangle[3])
	{
		FHandle fh = face_handle(tID);
		FVIter fv_it = fv_iter(fh);
		for (int i = 0; fv_it && i < 3; ++fv_it, i++) {
			int index = fv_it->idx();

			vTriangle[i] = index;
		}
	}

	void GetTriangle(int tID, OMT::Point vTriangle[3], OMT::Normal * pNormals = NULL) {

		FHandle fh = face_handle(tID);
		FVIter fv_it = fv_iter(fh);
		for (int i = 0; fv_it && i < 3; ++fv_it, i++) {
			OMT::Point p = point(fv_it);

			vTriangle[i][0] = p[0];
			vTriangle[i][1] = p[1];
			vTriangle[i][2] = p[2];

			if (pNormals) {
				OMT::Normal n = normal(fv_it);

				pNormals[i][0] = n[0];
				pNormals[i][1] = n[1];
				pNormals[i][2] = n[2];
			}
		}
	}

	void GetEdgeLengthStats(double & dMin, double & dMax, double & dAvg)
	{
		dMax = std::numeric_limits<double>::min();
		dMin = std::numeric_limits<double>::max();
		dAvg = 0;

		int nCount = 0;
		FIter curt(faces_begin()), endt(faces_end());
		while (curt != endt) {
			int nID = curt->idx();
			curt++;
			++nCount;

			OMT::Point vVertices[3];
			GetTriangle(nID, vVertices);

			for (int i = 0; i < 3; ++i) {
				double dLen = (vVertices[i] - vVertices[(i + 1) % 3]).norm();
				if (dLen > dMax)
					dMax = dLen;
				if (dLen < dMin)
					dMin = dLen;
				dAvg += dLen / 3.0;
			}
		}
		dAvg /= (double)nCount;
	}

	double GetMaxEdgeLength() {
		double dMaxEdgeLength = std::numeric_limits<double>::min();

		EIter e_it = edges_begin();
		for (e_it; e_it != edges_end(); ++e_it) {
			OMT::HEHandle heh = halfedge_handle(e_it.handle(), 1);

			OMT::Point from_point = point(from_vertex_handle(heh));
			OMT::Point to_point = point(to_vertex_handle(heh));
			OMT::Vec3d temp_vector = to_point - from_point;

			double length = temp_vector.norm();
			if (length > dMaxEdgeLength)
				dMaxEdgeLength = length;
		}

		//std::cout << "dMaxEdgeLength = " << dMaxEdgeLength << std::endl;
		return dMaxEdgeLength;
	}

	/*void Union(myAxisAlignedBox3d & dest, const OMT::Point & point)
	{
		for (int k = 0; k < 3; ++k) {
			if (point[k] < dest.Min[k])
				dest.Min[k] = point[k];
			if (point[k] > dest.Max[k])
				dest.Max[k] = point[k];
		}
	}*/

	bool HasUVSet(int nSetID)
	{
		return nSetID < property(UVSets, vertices_begin()).size();
	}

	int AppendUVSet()
	{
		//std::cout << "start AppendUVSet" << std::endl;

		for (VIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
		{
			property(UVSets, v_it).push_back(OMT::Vec2d(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()));
		}

		//std::cout << "end Append UVSet" << std::endl;
		return property(UVSets, vertices_begin()).size() - 1;
	}

	void InitializeUVSet(int nSetID)
	{
		for (VIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
		{
			property(UVSets, v_it)[nSetID] = OMT::Vec2d(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
		}

	}

	bool GetUV(int vID, int nSetID, OMT::Vec2d & vUV)
	{
		OMT::VHandle vh = vertex_handle(vID);
		if (!vh.is_valid())
			return false;

		vUV[0] = property(UVSets, vh)[nSetID][0];
		vUV[1] = property(UVSets, vh)[nSetID][1];

		return true;
	}

	bool GetTriangleUV(int tID, int nSetID, OMT::Vec2d vUV[3])
	{
		int nTri[3];
		GetTriangle(tID, nTri);

		for (int i = 0; i < 3; i++) {
			OMT::VHandle vh = vertex_handle(nTri[i]);

			if (!vh.is_valid())
				return false;

			vUV[i][0] = property(UVSets, vh)[nSetID][0];
			vUV[i][1] = property(UVSets, vh)[nSetID][1];
		}

		return true;
	}

	void SetUV(int vID, int nSetID, const OMT::Vec2d & vUV)
	{
		OMT::VHandle vh = vertex_handle(vID);

		property(UVSets, vh)[nSetID][0] = vUV[0];
		property(UVSets, vh)[nSetID][1] = vUV[1];
	}

	/*void findNearestVertex(OMT::Point p)
	{
		float length = std::numeric_limits<double>::max();
		for(OMT::VIter v_it = vertices_begin(); v_it != vertices_end(); v_it++)
	}

	void Lenght*/

	void ClearUVSet(int nSetID)
	{
		for (VIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
		{
			property(UVSets, v_it)[nSetID] = OMT::Vec2d(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
		}
	}

	OMT::Point Tri_Mesh::ClosestPtPointTriangle(FHandle fh, OMT::Point p, double &bu, double &bv, double &bw, double scalar)
	{
		int index = 0;
		OMT::Point fPoints[3];

		for (OMT::FVIter fv_it = fv_iter(fh); fv_it; ++fv_it)
		{
			OMT::Point fp = point(fv_it.handle());
			fPoints[index] = fp * scalar;
			++index;
		}

		OMT::Point center = (fPoints[0] + fPoints[1] + fPoints[2]) / 3.0;
		OMT::Vec3d centerP = p - center;
		//projection to triangle plane
		OMT::Normal fn = normal(fh);
		double len = dot(centerP, fn);
		OMT::Point projPoint = p - len * fn;

		//calculate distance
		OMT::Vec3d ab = OMT::Vec3d(fPoints[1] - fPoints[0]);
		OMT::Vec3d ac = OMT::Vec3d(fPoints[2] - fPoints[0]);
		OMT::Vec3d bc = OMT::Vec3d(fPoints[2] - fPoints[1]);
		OMT::Vec3d ap = OMT::Vec3d(projPoint - fPoints[0]);
		OMT::Vec3d bp = OMT::Vec3d(projPoint - fPoints[1]);
		OMT::Vec3d cp = OMT::Vec3d(projPoint - fPoints[2]);

		double d1 = dot(ab, ap);
		double d2 = dot(ac, ap);
		double d3 = dot(ab, bp);
		double d4 = dot(ac, bp);
		double d5 = dot(ab, cp);
		double d6 = dot(ac, cp);

		double va = d3 * d6 - d5 * d4;
		double vb = d5 * d2 - d1 * d6;
		double vc = d1 * d4 - d3 * d2;

		double snom = d1;
		double sdenom = -d3;
		double tnom = d2;
		double tdenom = -d6;
		double unom = d4 - d3;
		double udenom = d5 - d6;

		if (d1 <= 0 && d2 <= 0)
		{
			bu = 1;
			bv = 0;
			bw = 0;
			return fPoints[0];
		}

		if (d3 >= 0 && d4 <= d3)
		{
			bu = 0;
			bv = 1;
			bw = 0;
			return fPoints[1];
		}

		if (vc <= 0 && d1 >= 0 && d3 <= 0)
		{
			double v = d1 / (d1 - d3);
			bu = 1 - v;
			bv = v;
			bw = 0;

			return fPoints[0] + v * ab;
		}

		if (d6 >= 0 && d5 <= d6)
		{
			bu = 0;
			bv = 0;
			bw = 1;

			return fPoints[2];
		}

		if (vb <= 0 && d2 >= 0 && d6 <= 0)
		{
			double w = d2 / (d2 - d6);
			bu = 1 - w;
			bv = 0;
			bw = w;

			return fPoints[0] + w * ac;
		}

		if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
		{
			float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
			bu = 0;
			bv = 1 - w;
			bw = w;

			return fPoints[1] + w * bc;
		}

		double denom = 1 / (va + vb + vc);
		double v = vb * denom;
		double w = vc * denom;

		bu = 1 - v - w;
		bv = v;
		bw = w;

		return fPoints[0] + ab * v + ac * w;
	}
private:
};

///*======================================================================*/
/*======================================================================*/
bool ReadFile(std::string _fileName, Tri_Mesh *_mesh); //讀取mesh資料
bool SaveFile(std::string _fileName, Tri_Mesh *_mesh); //儲存mesh資料
/*初始化view port設定函式*/

#endif

