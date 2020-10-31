#include "GUA_OM.h"
#include <algorithm>
#include <iostream>


using namespace std;

class TriVertex
{
public:
	double x;
	double y;
	double z;
	TriVertex()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	TriVertex(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	~TriVertex()
	{
	}

	bool operator ==(const TriVertex& b)
	{
		if (this->x == b.x && this->y == b.y && this->z == b.z)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	friend ostream &operator<<(ostream& os, const TriVertex& tv)
	{
		os << "(" << tv.x << ", " << tv.y << ", " << tv.z << ")";
		return os;
	}

};
class TriLine
{
public:
	TriVertex p1;
	TriVertex p2;
	double vector[3] ;

	TriLine(TriVertex &p1, TriVertex &p2)
	{
		this->p1 = p1;
		this->p2 = p2;
		vector[0] = p2.x - p1.x;
		vector[1] = p2.y - p1.y;
		vector[2] = p2.z - p1.z;
	}
	~TriLine()
	{
	}
	friend ostream &operator<<(ostream& os, const TriLine& line)
	{
		os << "[" << line.p1 << " -> " << line.p2 <<  "]";
		return os;
	}
	double static cross(TriLine & line1, TriLine& line2)
	{
		return line1.vector[0] * line2.vector[0] + line1.vector[1] * line2.vector[1] + line1.vector[2] * line2.vector[2];
	}

};


namespace OMT {
	/*======================================================================*/
	Model::Model() {
		request_vertex_status();
		request_edge_status();
		request_face_status();
	}
	Model::~Model() {
		release_vertex_status();
		release_edge_status();
		release_face_status();
	}
}
/*======================================================================*/
namespace OMP {
	Model::Model() {
		Mesh.request_vertex_status();
		Mesh.request_edge_status();
		Mesh.request_face_status();
	}
	Model::~Model() {
		Mesh.release_vertex_status();
		Mesh.release_edge_status();
		Mesh.release_face_status();
	}
	/*======================================================================*/
	bool Model::ReadFile(std::string _fileName) {
		bool isRead = false;
		OpenMesh::IO::Options opt;
		if (OpenMesh::IO::read_mesh(Mesh, _fileName, opt)) {
			//read mesh from filename OK!
			isRead = true;
		}
		if (isRead) {
			// If the file did not provide vertex normals and mesh has vertex normal ,then calculate them
			if (!opt.check(OpenMesh::IO::Options::VertexNormal) && Mesh.has_vertex_normals()) {
				Mesh.update_normals();
			}
		}
		return isRead;
	}
	bool Model::SaveFile(std::string _fileName) {
		bool isSave = false;
		OpenMesh::IO::Options opt;
		if (OpenMesh::IO::write_mesh(Mesh, _fileName, opt)) {
			//read mesh from filename OK!
			isSave = true;
		}
		return isSave;
	}
	/*======================================================================*/
	void Model::Render_solid() {
		FIter f_it;
		FVIter	fv_it;
		glPushAttrib(GL_LIGHTING_BIT);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glEnable(GL_DEPTH_TEST);
		glPolygonOffset(2.0, 2.0);
		glBegin(GL_POLYGON);
		//glColor4f(1.0, 0.5, 1.0, 0.5);
		for (f_it = Mesh.faces_begin(); f_it != Mesh.faces_end(); ++f_it) {
			for (fv_it = Mesh.fv_iter(f_it); fv_it; ++fv_it) {
				glNormal3dv(Mesh.normal(fv_it.handle()).data());
				glVertex3dv(Mesh.point(fv_it.handle()).data());
			}
		}
		glEnd();
		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	void Model::Render_wireframe() {
		MyMesh::HalfedgeHandle _hedge;
		EIter e_it = Mesh.edges_begin();

		glDisable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);
		glColor3f(0.0, 0.0, 0.0);
		glLineWidth(1);
		glBegin(GL_LINES);
		for (e_it = Mesh.edges_begin(); e_it != Mesh.edges_end(); ++e_it) {
			_hedge = Mesh.halfedge_handle(e_it.handle(), 1);

			glVertex3dv(Mesh.point(Mesh.from_vertex_handle(_hedge)).data());
			glVertex3dv(Mesh.point(Mesh.to_vertex_handle(_hedge)).data());
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
	/*======================================================================*/
	void Model::RenderSpecifiedPoint() {
		glPushAttrib(GL_LIGHTING_BIT);
		glDisable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);
		glPointSize(5.0f);
		glBegin(GL_POINTS);
		vector<sp_p>::iterator p_itr = sp_p_list.begin();
		for (p_itr; p_itr != sp_p_list.end(); ++p_itr) {
			glColor3f(p_itr->r, p_itr->g, p_itr->b);
			glVertex3dv(p_itr->pt.data());
		}
		glEnd();
		glEnable(GL_LIGHTING);
		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	void Model::RenderSpecifiedVertex() {
		glPushAttrib(GL_LIGHTING_BIT);
		glDisable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);
		glPointSize(5.0f);
		glBegin(GL_POINTS);
		vector< sp_v >::iterator v_itr = sp_v_list.begin();
		for (v_itr; v_itr != sp_v_list.end(); ++v_itr) {
			glColor3f(v_itr->r, v_itr->g, v_itr->b);
			glVertex3dv(Mesh.point(v_itr->vh).data());
		}
		glEnd();
		glEnable(GL_LIGHTING);
		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	void Model::RenderSpecifiedFace() {
		glDisable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glPushAttrib(GL_LIGHTING_BIT);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(0.5, 1.0);
		glDisable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);
		glBegin(GL_QUADS);
		FVIter fv_itr;
		vector< sp_f >::iterator f_itr;
		for (f_itr = sp_f_list.begin(); f_itr != sp_f_list.end(); ++f_itr) {
			glColor3f(f_itr->r, f_itr->g, f_itr->b);
			for (fv_itr = Mesh.fv_iter(f_itr->fh); fv_itr; ++fv_itr) {
				glNormal3dv(Mesh.normal(fv_itr.handle()).data());
				glVertex3dv(Mesh.point(fv_itr.handle()).data());
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
		glDisable(GL_POLYGON_OFFSET_FILL);
		glPolygonMode(GL_FRONT, GL_FILL);
		glEnable(GL_CULL_FACE);
	}

	/*======================================================================*/
	void Model::add_sp_p(Point   _p, float _r, float _g, float _b) {
		sp_p input_data;
		input_data.pt = _p;
		input_data.r = _r;
		input_data.g = _g;
		input_data.b = _b;
		sp_p_list.push_back(input_data);
	}
	void Model::add_sp_v(VHandle _v, float _r, float _g, float _b) {
		sp_v input_data;
		input_data.vh = _v;
		input_data.r = _r;
		input_data.g = _g;
		input_data.b = _b;
		sp_v_list.push_back(input_data);
	}
	void Model::add_sp_f(FHandle _f, float _r, float _g, float _b) {
		sp_f input_data;
		input_data.fh = _f;
		input_data.r = _r;
		input_data.g = _g;
		input_data.b = _b;
		sp_f_list.push_back(input_data);
	}
	void Model::clear_sp_p() {
		sp_p_list.clear();
	}
	void Model::clear_sp_v() {
		sp_v_list.clear();
	}
	void Model::clear_sp_f() {
		sp_f_list.clear();
	}
	/*======================================================================*/
	VHandle Model::addVertex(Point _p) {
		int find_result = findVertex(_p);
		if (find_result != -1) {
			return Mesh.vertex_handle(find_result);
		}
		else {
			return Mesh.add_vertex(_p);
		}
	}
	FHandle Model::addFace(VHandle _v0, VHandle _v1, VHandle _v2, VHandle _v3) {
		vector<VHandle> face_vhandles;

		face_vhandles.clear();
		face_vhandles.push_back(_v0);
		face_vhandles.push_back(_v1);
		face_vhandles.push_back(_v2);
		face_vhandles.push_back(_v3);
		return Mesh.add_face(face_vhandles);
	}
	void Model::deleteFace(FHandle _f) {
		Mesh.delete_face(_f);
		Mesh.garbage_collection();
	}
	void Model::deleteFace(VHandle _v0, VHandle _v1, VHandle _v2, VHandle _v3) {
		/*
		v1				v0
		*<--------------*
		|				|
		|				|
		|				|
		|		f		|
		|				|
		|				|
		|				|
		* --------------*
		v2				v3
		*/

		HEHandle v0v1 = Mesh.find_halfedge(_v0, _v1);
		if (v0v1.is_valid()) {
			FHandle fh = Mesh.face_handle(v0v1);
			if (fh.is_valid()) {
				Mesh.delete_face(fh);
				Mesh.garbage_collection();
			}
		}
	}
	/*======================================================================*/
	bool Model::IsVertexVertex(VHandle _vj, VHandle _vi) {
		for (VVIter vvit = Mesh.vv_iter(_vi); vvit; ++vvit)
			if (vvit.handle() == _vj)
				return true;
		return false;
	}
	/*======================================================================*/
	int Model::quad_subdivision1(int _face_id) {
		/*----------------------------------------------------------------------*/
		//這段是為了解決index問題
		VHandle vq, vw, vt, vr;
		vq = addVertex(Point(0, 0, 100));
		vw = addVertex(Point(1, 0, 100));
		vt = addVertex(Point(1, 1, 100));
		vr = addVertex(Point(0, 1, 100));
		addFace(vq, vw, vt, vr);
		/*----------------------------------------------------------------------*/
		/*收集需要subdivision的face*/
		//初始table
		bool *table = new bool[Mesh.n_faces()];
		for (unsigned i = 0; i < Mesh.n_faces(); i++) {
			table[i] = false;
		}

		vector< FHandle > candidate_faces, last_found_face;
		last_found_face.push_back(Mesh.face_handle(_face_id));
		table[_face_id] = true;

		while (last_found_face.size() != 0) {
			vector< FHandle > new_found_faces;
			for (vector< FHandle >::iterator crnt_f = last_found_face.begin(); crnt_f != last_found_face.end(); ++crnt_f) {
				for (FFIter ff_itr = Mesh.ff_iter(*crnt_f); ff_itr; ++ff_itr) {
					if (table[ff_itr.handle().idx()] != true) {
						new_found_faces.push_back(ff_itr.handle());
						table[ff_itr.handle().idx()] = true;
					}
				}
			}
			for (vector< FHandle >::iterator f_itr = last_found_face.begin(); f_itr != last_found_face.end(); ++f_itr) {
				candidate_faces.push_back(*f_itr);
			}
			last_found_face.assign(new_found_faces.begin(), new_found_faces.end());
		}
		delete table;
		/*----------------------------------------------------------------------*/
		/*對candidate faces做subdivision*/
		/*
			v0		vd		v3
			*-------*-------*
			|		|		|
			|		|		|
			|	  ve|		|
		  va*-------*-------*vc
			|		|		|
			|		|		|
			|		|		|
			*-------*-------*
			v1		vb		v2
		*/
		for (vector< FHandle >::iterator face_itr = candidate_faces.begin(); face_itr != candidate_faces.end(); ++face_itr) {
			VHandle v[4], va, vb, vc, vd, ve;
			FVIter fv_itr = Mesh.fv_iter(*face_itr);
			for (int i = 0; fv_itr; ++fv_itr) {
				v[i++] = fv_itr.handle();
			}

			deleteFace(v[0], v[1], v[2], v[3]);
			va = addVertex((Mesh.point(v[0]) + Mesh.point(v[1])) / 2);
			vb = addVertex((Mesh.point(v[1]) + Mesh.point(v[2])) / 2);
			vc = addVertex((Mesh.point(v[2]) + Mesh.point(v[3])) / 2);
			vd = addVertex((Mesh.point(v[3]) + Mesh.point(v[0])) / 2);
			ve = addVertex((Mesh.point(v[0]) + Mesh.point(v[1]) + Mesh.point(v[2]) + Mesh.point(v[3])) / 4);
			addFace(vd, v[0], va, ve);
			addFace(va, v[1], vb, ve);
			addFace(vb, v[2], vc, ve);
			addFace(vc, v[3], vd, ve);
		}
		/*----------------------------------------------------------------------*/
		deleteFace(vq, vw, vt, vr);//這行只是為了解決index問題
		/*----------------------------------------------------------------------*/
		return 0;
	}
	int Model::quad_subdivision2(int _face_id) {
		/*----------------------------------------------------------------------*/
		//這段是為了解決index問題
		VHandle vq, vw, vt, vr;
		vq = addVertex(Point(0, 0, 100));
		vw = addVertex(Point(1, 0, 100));
		vt = addVertex(Point(1, 1, 100));
		vr = addVertex(Point(0, 1, 100));
		addFace(vq, vw, vt, vr);
		/*----------------------------------------------------------------------*/
		/*收集需要subdivision的face*/
		//初始table
		bool *table = new bool[Mesh.n_faces()];
		for (unsigned i = 0; i < Mesh.n_faces(); i++) {
			table[i] = false;
		}

		vector< FHandle > candidate_faces, last_found_face;
		last_found_face.push_back(Mesh.face_handle(_face_id));
		table[_face_id] = true;

		while (last_found_face.size() != 0) {
			vector< FHandle > new_found_faces;
			for (vector< FHandle >::iterator crnt_f = last_found_face.begin(); crnt_f != last_found_face.end(); ++crnt_f) {
				for (FFIter ff_itr = Mesh.ff_iter(*crnt_f); ff_itr; ++ff_itr) {
					if (table[ff_itr.handle().idx()] != true) {
						new_found_faces.push_back(ff_itr.handle());
						table[ff_itr.handle().idx()] = true;
					}
				}
			}
			for (vector< FHandle >::iterator f_itr = last_found_face.begin(); f_itr != last_found_face.end(); ++f_itr) {
				candidate_faces.push_back(*f_itr);
			}
			last_found_face.assign(new_found_faces.begin(), new_found_faces.end());
		}
		delete table;
		/*----------------------------------------------------------------------*/
		/*對candidate faces做subdivision*/
		/*
			v0		vh		vg		v3
			*-------*-------*-------*
			|		|		|		|
			|		|		|		|
			|	  vi|	  vl|		|
		 va *-------*-------*-------*vf
			|		|		|		|
			|		|		|		|
			|	  vj|	  vk|		|
		 vb *-------*-------*-------*ve
			|		|		|		|
			|		|		|		|
			|		|		|		|
			*-------*-------*-------*
			v1		vc		vd		v2
		*/
		for (vector< FHandle >::iterator face_itr = candidate_faces.begin(); face_itr != candidate_faces.end(); ++face_itr) {
			VHandle v[4], va, vb, vc, vd, ve, vf, vg, vh, vi, vj, vk, vl;
			FVIter fv_itr = Mesh.fv_iter(*face_itr);
			for (int i = 0; fv_itr; ++fv_itr) {
				v[i++] = fv_itr.handle();
			}

			deleteFace(v[0], v[1], v[2], v[3]);
			va = addVertex((Mesh.point(v[0]) * 2 + Mesh.point(v[1])) / 3);
			vb = addVertex((Mesh.point(v[0]) + Mesh.point(v[1]) * 2) / 3);
			vc = addVertex((Mesh.point(v[1]) * 2 + Mesh.point(v[2])) / 3);
			vd = addVertex((Mesh.point(v[1]) + Mesh.point(v[2]) * 2) / 3);
			ve = addVertex((Mesh.point(v[2]) * 2 + Mesh.point(v[3])) / 3);
			vf = addVertex((Mesh.point(v[2]) + Mesh.point(v[3]) * 2) / 3);
			vg = addVertex((Mesh.point(v[3]) * 2 + Mesh.point(v[0])) / 3);
			vh = addVertex((Mesh.point(v[3]) + Mesh.point(v[0]) * 2) / 3);

			vi = addVertex((Mesh.point(vh) * 2 + Mesh.point(vc)) / 3);
			vj = addVertex((Mesh.point(vh) + Mesh.point(vc) * 2) / 3);
			vk = addVertex((Mesh.point(vd) * 2 + Mesh.point(vg)) / 3);
			vl = addVertex((Mesh.point(vd) + Mesh.point(vg) * 2) / 3);
			addFace(v[0], va, vi, vh);
			addFace(va, vb, vj, vi);
			addFace(vb, v[1], vc, vj);
			addFace(vc, vd, vk, vj);
			addFace(vd, v[2], ve, vk);
			addFace(ve, vf, vl, vk);
			addFace(vf, v[3], vg, vl);
			addFace(vg, vh, vi, vl);
			addFace(vi, vj, vk, vl);
		}
		/*----------------------------------------------------------------------*/
		deleteFace(vq, vw, vt, vr);//這行只是為了解決index問題
		/*----------------------------------------------------------------------*/
		return 0;
	}
	/*======================================================================*/
	int Model::findVertex(Point _p) {
		for (VIter v_itr = Mesh.vertices_begin(); v_itr != Mesh.vertices_end(); ++v_itr)
			if (Mesh.point(v_itr) == _p)
				return v_itr.handle().idx();
		return -1;
	}
	/*======================================================================*/
};
/*======================================================================*/

const OpenMesh::VertexHandle Tri_Mesh::InvalidVertexHandle;

void Tri_Mesh::Render_Solid() {
	FIter f_it;
	FVIter	fv_it;
	//glPushAttrib(GL_LIGHTING_BIT);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_LIGHTING);
	glPolygonOffset(2.0, 2.0);
	glBegin(GL_TRIANGLES);
	glColor4f(0.81, 0.74, 0.33, 0.3);
	for (f_it = faces_begin(); f_it != faces_end(); ++f_it) {
		for (fv_it = fv_iter(f_it); fv_it; ++fv_it) {
			glNormal3dv(normal(fv_it.handle()).data());
			glVertex3dv(point(fv_it.handle()).data());
		}
	}
	glEnd();

	glDisable(GL_POLYGON_OFFSET_FILL);
}

void Tri_Mesh::Render_SolidWireframe() {
	FIter f_it;
	FVIter	fv_it;

	glDisable(GL_LIGHTING);
	glPushAttrib(GL_LIGHTING_BIT);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_DEPTH_TEST);
	glPolygonOffset(2.0, 2.0);
	glBegin(GL_TRIANGLES);
	glColor4f(1.0, 0.96, 0.49, 1.0);
	for (f_it = faces_begin(); f_it != faces_end(); ++f_it) {
		for (fv_it = fv_iter(f_it); fv_it; ++fv_it) {
			//glNormal3dv(normal(fv_it.handle()));
			glVertex3dv(point(fv_it.handle()).data());
		}
	}
	glEnd();

	//glDisable(GL_POLYGON_OFFSET_FILL);

	glPushAttrib(GL_LIGHTING_BIT);
	glDisable(GL_LIGHTING);
	glLineWidth(1.0);
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
	for (OMT::EIter e_it = edges_begin(); e_it != edges_end(); ++e_it) {
		OMT::HEHandle _hedge = halfedge_handle(e_it.handle(), 1);

		OMT::Point curVertex = point(from_vertex_handle(_hedge));
		glVertex3dv(curVertex.data());

		curVertex = point(to_vertex_handle(_hedge));
		glVertex3dv(curVertex.data());
	}
	glEnd();
	glPopAttrib();
}
//cache mesh to buffer
void Tri_Mesh::loadToBuffer(Tri_Mesh _mesh, std::vector<double> & out_vertices, int & face, std::vector<double> & uv) {
	//to buffer , then shader
	FIter f_it;
	FVIter	fv_it;
	out_vertices.clear();
	face = 0;
	/*std::cout << "mesh support : " << std::endl;
	std::cout << "  " << "texcoords" << ": " << ((_mesh.has_vertex_texcoords1D()) ? "yes\n" : "no\n") << std::endl;
	std::cout << "  " << "texcoords" << ": " << ((_mesh.has_vertex_texcoords2D()) ? "yes\n" : "no\n") << std::endl;
	std::cout << "  " << "texcoords" << ": " << ((_mesh.has_vertex_texcoords3D()) ? "yes\n" : "no\n") << std::endl;
	std::cout << "  " << "normal" << ": " << ((_mesh.has_vertex_normals()) ? "yes\n" : "no\n") << std::endl;*/

	for (f_it = faces_begin(); f_it != faces_end(); ++f_it) {
		face++;
		fv_it = fv_iter(f_it);
		//force 3vert or infinite loop??
		for (int i = 0; i < 3; i++, ++fv_it) {
			// 每個點有三個vertexes
			out_vertices.push_back(*(point(fv_it.handle()).data()));
			out_vertices.push_back(*(point(fv_it.handle()).data() + 1));
			out_vertices.push_back(*(point(fv_it.handle()).data() + 2));
			//std::cout << "x = " << *(point(fv_it.handle()).data()) << " y = " << *(point(fv_it.handle()).data() + 1) << " z = " << *(point(fv_it.handle()).data() + 2) << std::endl;
			uv.push_back(_mesh.texcoord2D(fv_it.handle())[0]);
			uv.push_back(_mesh.texcoord2D(fv_it.handle())[1]);
			//std::cout << "s = " << _mesh.texcoord2D(fv_it.handle())[0] << " t = " << _mesh.texcoord2D(fv_it.handle())[0] << std::endl;
		}
		//cout << "face: " << out_vertices[out_vertices.size() - 3] << " " << out_vertices[out_vertices.size() - 2] << " " << out_vertices[out_vertices.size() - 1] << endl;
	}

}

void Tri_Mesh::delVert(VHandle vhandle) {
	std::vector<VHandle> neighborVert = std::vector<VHandle>();
	VVCWIter vv_it;
	//find one ring
	for (vv_it = vv_cwiter(vhandle); vv_it.is_valid(); ++vv_it) {
		neighborVert.push_back(*vv_it);
		cout << point(*vv_it) << endl;
	}
	//delete
	delete_vertex(vhandle);
	//repair
	cout << "should add " << neighborVert.size() - 2 << " face" << endl;

	for (int index = 1; index < neighborVert.size() - 1; index++) {
		cout << "add." << endl;
		//important.  clockwise or get complex edge error
		add_face(neighborVert[0], neighborVert[index + 1], neighborVert[index]);
	}

	garbage_collection();
	return;
}

//find choosed vert via choosed face
void Tri_Mesh::findNearestVert(Tri_Mesh mesh, std::vector<double> mouse, int face, std::vector<double> &vertex, mat4 MVP, double dis) {
	return;
	FIter f_it;
	FVIter fv_it;
	VHandle min;
	int isFaceMatch = 0;
	for (f_it = faces_begin(); f_it != faces_end() && !isFaceMatch; ++f_it) {
		//find chosen face
		if (f_it.handle().idx() != face) continue;
		else min = fv_iter(f_it).handle();
		// compare three vert
		for (fv_it = ++fv_iter(f_it); fv_it; ++fv_it) {
			//put point into vec to mult MVP
			vec3 tempPoint, curPoint;
			for (int i = 0; i < 3; i++) {
				tempPoint[i] = point(min)[i];
				curPoint[i] = point(fv_it.handle())[i];
			}
			//fix range issue 
			vec4 tempPosOnScreen = MVP * vec4(tempPoint, 1) / pow(dis, 1.35) / 0.7;
			vec4 curPosOnScreen = MVP * vec4(curPoint, 1) / pow(dis, 1.35) / 0.7;
			//compare mouse dis via iterator
			double temp = pow((tempPosOnScreen[0] - mouse[0]), 2) + pow((tempPosOnScreen[1] - mouse[1]), 2);
			double current = pow((curPosOnScreen[0] - mouse[0]), 2) + pow((curPosOnScreen[1] - mouse[1]), 2);
			if (temp > current) min = fv_it.handle();
			//cout << "cur iter point scr position: " << curPosOnScreen[0] << "," << curPosOnScreen[1] << std::endl;
			isFaceMatch = 1;
		}
	}
	for (int i = 0; i < 3 && isFaceMatch; i++) vertex.push_back(point(min)[i]);
	
	//delVert(min);
	oneRingCollapse(min);

	//if (isFaceMatch) printf("selected point is : %f %f %f\n", vertex[0], vertex[1], vertex[2]);
}

///////render func /////////
void Tri_Mesh::Render_Wireframe() {
	//glPushAttrib(GL_LIGHTING_BIT);	
	glDisable(GL_LIGHTING);
	glLineWidth(1.0);

	glColor3f(0.0, 0.0, 0.0);

	glBegin(GL_LINES);
	for (OMT::EIter e_it = edges_begin(); e_it != edges_end(); ++e_it) {
		OMT::HEHandle _hedge = halfedge_handle(e_it.handle(), 1);

		OMT::Point curVertex = point(from_vertex_handle(_hedge));
		glVertex3dv(curVertex.data());

		curVertex = point(to_vertex_handle(_hedge));
		glVertex3dv(curVertex.data());
	}
	glEnd();

}

void Tri_Mesh::Render_Point() {
	glPointSize(8.0);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_POINTS);
	for (OMT::VIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {
		glVertex3dv(point(v_it).data());
	}
	glEnd();
}
///////// file  func ////////////
bool ReadFile(std::string _fileName, Tri_Mesh *_mesh) {
	bool isRead = false;
	OpenMesh::IO::Options opt;
	opt += OpenMesh::IO::Options::VertexTexCoord;
	if (OpenMesh::IO::read_mesh(*_mesh, _fileName, opt)) {
		//read mesh from filename OK!
		isRead = true;
	}
	if (isRead) {
		// If the file did not provide vertex normals and mesh has vertex normal ,then calculate them
		if (!opt.check(OpenMesh::IO::Options::VertexNormal) && _mesh->has_vertex_normals()) {
			_mesh->update_normals();
		}
	}

	printf("opt.check( OpenMesh::IO::Options::VertexNormal ) = %d\n", opt.check(OpenMesh::IO::Options::VertexNormal));
	printf("_mesh->has_vertex_normals() = %d\n", _mesh->has_vertex_normals());

	printf("opt.check( OpenMesh::IO::Options::VertexTexCoord ) = %d\n", opt.check(OpenMesh::IO::Options::VertexTexCoord));
	printf("_mesh->has_vertex_texcoord2D() = %d\n", _mesh->has_vertex_texcoords2D());

	return isRead;
}

bool SaveFile(std::string _fileName, Tri_Mesh *_mesh) {
	bool isSave = false;
	OpenMesh::IO::Options opt = OpenMesh::IO::Options::VertexTexCoord;
	opt += OpenMesh::IO::Options::VertexNormal;
	if (OpenMesh::IO::write_mesh(*_mesh, _fileName, opt))isSave = true;
	return isSave;
}

//////////// not using funcs////////////
//@Deprecated func
void Tri_Mesh::findNearestPoint(Tri_Mesh mesh, std::vector<double> mouse, int face, std::vector<double> &vertex) {
	FIter f_it;
	FVIter fv_it;
	VHandle min;
	int isFaceMatch = 0;
	for (f_it = faces_begin(); f_it != faces_end() && !isFaceMatch; ++f_it) {

		if (f_it.handle().idx() != face) continue;
		else min = fv_iter(f_it).handle();

		for (fv_it = ++fv_iter(f_it); fv_it; ++fv_it) { // 每個面有三個點
			//printf("in for\n");
			double temp = std::sqrtl(std::pow((point(min)[0] - mouse[0]), 2) + std::pow((point(min)[1] - mouse[1]), 2) + std::pow((point(min)[2] - mouse[2]), 2));
			double current = std::sqrtl(std::pow((point(fv_it.handle())[0] - mouse[0]), 2) + std::pow((point(fv_it.handle())[1] - mouse[1]), 2) + std::pow((point(fv_it.handle())[2] - mouse[2]), 2));
			printf(" temp %f  v.s. current %f\n", temp, current);
			std::cout << "cur iter point : " << point(fv_it.handle())[0] << "," << point(fv_it.handle())[1] << "," << point(fv_it.handle())[2] << std::endl;
			if (temp > current) min = fv_it.handle();
			isFaceMatch = 1;
		}
	}
	//printf("here?\n");

	for (int i = 0; i < 3 && isFaceMatch; i++) {
		vertex.push_back(point(min)[i]);
	}
	if (isFaceMatch) printf("selected point is : %f %f %f\n", vertex[0], vertex[1], vertex[2]);
}
//get uv coordinate
void Tri_Mesh::getUV(std::vector<double> & patchuv, Tri_Mesh patch, float uvRotateAngle) {
	VVIter vv_it;
	VIter v_it;
	EIter e_it;
	HHandle heh;
	FVIter fv_it;
	FIter f_it;
	//Step1 : 先找到第一個邊界
	for (e_it = edges_begin(); e_it != edges_end(); e_it++) {
		printf("enter edge iterator...\n");
		bool isBoundary = patch.is_boundary(*e_it);
		if (isBoundary) {
			heh = patch.halfedge_handle(e_it.handle(), 1);
			break;
		}
	}
	// Step2 : 沿著該邊界找尋下一個邊
	double perimeter = 0;
	std::vector<double> segLength;
	std::vector<Tri_Mesh::VertexHandle> vhs; // 儲存排序好的邊界點
	HHandle hehNext = heh;
	do {
		Point from = patch.point(patch.from_vertex_handle(hehNext));
		Point to = patch.point(patch.to_vertex_handle(hehNext));
		perimeter += (from - to).length(); // v0 - v? 的長度
		printf("perimeter = %f\n", perimeter);
		segLength.push_back(perimeter); // 存入vector中，以便之後做texcoord
		vhs.push_back(patch.from_vertex_handle(hehNext)); // 把邊界上的點一一存入
		hehNext = patch.next_halfedge_handle(hehNext); // 可以直接依序往下一個heh走
	} while (heh != hehNext);

	//Step3 : 找完所有的邊界點和距離後，可以將邊界已知點換成texcoord
	float rd = (225 + uvRotateAngle) * M_PI / 180.0;
	float initDist = 0;
	Tri_Mesh::TexCoord2D st(0, 0);
	float R = std::sqrt(2) / 2.0; // 根號2/2
	st[0] = R * cos(rd) + 0.5;
	st[1] = R * sin(rd) + 0.5;

	if (st[0] > 1) {
		st[0] = 1;
		st[1] = tan(rd) / 2 + 0.5;
	}
	if (st[0] < 0) {
		st[0] = 0;
		st[1] = 0.5 - tan(rd) / 2;
	}
	if (st[1] > 1) {
		st[0] = tan(M_PI_2 - rd) / 2 + 0.5;
		st[1] = 1;
	}

	if (st[1] < 0) {
		st[0] = 0.5 - tan(M_PI_2 - rd) / 2;
		st[1] = 0;
	}

	if (uvRotateAngle <= 90) {
		initDist = st.length();
	}
	else if (uvRotateAngle <= 180) {
		initDist = 1 + st[1];
	}
	else if (uvRotateAngle <= 270) {
		initDist = 3 - st[0];
	}
	else {
		initDist = 4 - st[1];
	}

	patch.request_vertex_texcoords2D();
	patch.set_texcoord2D(vhs[0], st);
	perimeter /= 4.0;
	for (int i = 1; i < vhs.size(); ++i) {
		double curLen = segLength[i - 1] / perimeter + initDist; // 意即L0n/Ltotal*4
		if (curLen > 4) {
			curLen -= 4;
		}
		if (curLen <= 1) {
			st[0] = curLen;
			st[1] = 0;
		}
		else if (curLen <= 2) {
			st[0] = 1;
			st[1] = curLen - 1;
		}
		else if (curLen <= 3) {
			st[0] = 3 - curLen;
			st[1] = 1;
		}
		else {
			st[0] = 0;
			st[1] = 4 - curLen;
		}
		patch.set_texcoord2D(vhs[i], st); // 先求出邊界點的uv座標，做為已知點
	}

	//Step4 : 接下來，算出非邊界的權重(使用均質公式)
	OpenMesh::HPropHandleT<double> heWeight;
	OpenMesh::VPropHandleT<int> row;
	patch.add_property(heWeight, "heWeight"); //加入property：權重
	patch.add_property(row, "row");

	for (e_it = patch.edges_begin(); e_it != patch.edges_end(); ++e_it) {
		bool isBound = patch.is_boundary(*e_it);
		if (!isBound) { // 如果不是邊界，就要算他的權重，兩個half_edge都要
			GLdouble angle1, angle2, w;
			Tri_Mesh::HalfedgeHandle _heh = patch.halfedge_handle(e_it.handle(), 0);
			Point pFrom = patch.point(patch.from_vertex_handle(_heh)); // 起點
			Point pTo = patch.point(patch.to_vertex_handle(_heh)); // 終點
			//問題：opposite_he_opposite_vh, opposite_vh抱錯
			Point p1 = patch.point(patch.opposite_vh(_heh)); // 線上對點
			Point p2 = patch.point(patch.opposite_he_opposite_vh(_heh)); // 線下對點
			//auto _vh = patch.opposite_vh(_heh);
			double edgeLen = (pFrom - pTo).length(); //線長度
			OpenMesh::Vec3d v1 = (OpenMesh::Vec3d)(pTo - pFrom);
			v1.normalize();
			OpenMesh::Vec3d v2 = (OpenMesh::Vec3d)(p1 - pFrom);
			v2.normalize();
			angle1 = std::acos(OpenMesh::dot(v1, v2)); //線上點與線的arccos，取得角度

			v2 = (OpenMesh::Vec3d)(p2 - pFrom);
			v2.normalize();

			angle2 = std::acos(OpenMesh::dot(v1, v2)); //線下點與線的arccos，取得角度
			//均質公式
			w = (std::tan(angle1 / 2.0f) + std::tan(angle2 / 2.0f)) / edgeLen;
			patch.property(heWeight, _heh) = w;  //將其中一個edge算好丟入邊的性質中

			// 再來計算反方向的halfedge weight
			v1 = -v1;
			v2 = (OpenMesh::Vec3d)(p1 - pTo);
			v2.normalize();
			angle1 = std::acos(OpenMesh::dot(v1, v2));//線上點與線的arccos，取得角度

			v2 = (OpenMesh::Vec3d)(p2 - pTo);
			v2.normalize();
			angle2 = std::acos(OpenMesh::dot(v1, v2));//線下點與線的arccos，取得角度

			w = (std::tan(angle1 / 2.0f) + std::tan(angle2 / 2.0f)) / edgeLen;
			patch.property(heWeight, patch.opposite_halfedge_handle(_heh)) = w;//將反方向edge算好丟入邊的性質中
		}
	}
	//Step5 : 找出矩陣的大小(NxN)，總點數-邊界點
	int count = 0;
	for (v_it = patch.vertices_begin(); v_it != patch.vertices_end(); v_it++) { // debug：patch.vertices_end()要加patch
		if (patch.is_boundary(*v_it)) patch.property(row, *v_it) = -1;
		else patch.property(row, *v_it) = count++;
	}
	//Step6 : 填寫矩陣
	typedef Eigen::SparseMatrix<double> SpMat;
	SpMat A(count, count);
	Eigen::VectorXd BX(count);
	Eigen::VectorXd BY(count);//x 和 y 各解一次
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > linearSolver;
	BX.setZero();
	BY.setZero();
	// fiil matrix
	for (v_it = patch.vertices_begin(); v_it != patch.vertices_end(); ++v_it) {
		if (!patch.is_boundary(*v_it)) { // 未知點是我們要解的東西
			int i = patch.property(row, *v_it);
			double totalWeight = 0;
			for (vv_it = patch.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
				HHandle _heh = patch.find_halfedge(*v_it, *vv_it);
				double w = patch.property(heWeight, _heh);

				if (patch.is_boundary(*vv_it)) { //是邊界的為已知，填入B矩陣
					TexCoord2D texCoord = patch.texcoord2D(*vv_it);
					BX[i] += w * texCoord[0]; // 單行全相加
					BY[i] += w * texCoord[1];
				}
				else { // 不是邊界的為未知點，填入A矩陣
					int j = patch.property(row, *vv_it);
					A.insert(i, j) = -w;
				}
				totalWeight += w;
			}
			A.insert(i, i) = totalWeight; // 這裡的解法不用讓所有的w除上大W
		}
	}
	A.makeCompressed();
	//Step7 : 解開未知的uv座標們
	SpMat At = A.transpose();
	linearSolver.compute(At * A);

	Eigen::VectorXd TX = linearSolver.solve(At * BX);
	Eigen::VectorXd TY = linearSolver.solve(At * BY);
	for (v_it = patch.vertices_begin(); v_it != patch.vertices_end(); ++v_it) {
		if (!patch.is_boundary(*v_it)) {
			int i = patch.property(row, *v_it);
			patch.set_texcoord2D(*v_it, TexCoord2D(TX[i], TY[i]));
		}
	}
	for (f_it = faces_begin(); f_it != faces_end(); ++f_it) {
		for (fv_it = fv_iter(f_it); fv_it; ++fv_it) {
			patchuv.push_back(patch.texcoord2D(*fv_it)[0]);
			patchuv.push_back(patch.texcoord2D(*fv_it)[1]);
			//patchuv.push_back(*(point(fv_it.handle()).data()));
			//patchuv.push_back(*(point(fv_it.handle()).data() + 1));
		}
	}
	//printf("patchuv : \n");
	//printf("patchuv.size() : %d\n", patchuv.size());
	//for (int i = 0; i < patchuv.size(); i += 2) printf("s = %f , t = %f\n", patchuv.at(i), patchuv.at(i + 1));
}
//save selected patch for right up window
void Tri_Mesh::loadToBufferPatch(std::vector<double> & out_vertices, int & face, std::vector<int> selectedFace, Tri_Mesh & patch) {
	FIter f_it;
	FVIter fv_it;
	VIter v_it;
	face = 0;
	std::vector<Tri_Mesh::VertexHandle> vhandle;
	std::vector<Tri_Mesh::VertexHandle> face_vhandles;
	int verticesSeq[3] = { 0,0,0 };
	int verticesSeqIndex = 0;
	int isPatchHasPoint = 0;
	std::map<int, int> mapVerticesToVHandle;
	for (f_it = faces_begin(); f_it != faces_end(); ++f_it) {
		// check selected face
		int isFaceSelected = 0;
		for (int i = 0; i < selectedFace.size(); i++) {
			if (f_it.handle().idx() == selectedFace[i]) { // found face
				isFaceSelected = 1;
				break;
			}
			else if (f_it.handle().idx() > selectedFace[i] && i == selectedFace.size() - 1) isFaceSelected = 0; // face not selected
		}
		if (isFaceSelected == 0) continue;// goto next face if current not selected
		face++;
		for (fv_it = fv_iter(f_it); fv_it; ++fv_it) { // get verts of selected face
			out_vertices.push_back(*(point(fv_it.handle()).data()));
			out_vertices.push_back(*(point(fv_it.handle()).data() + 1));
			out_vertices.push_back(*(point(fv_it.handle()).data() + 2));
			//check repeat vert, if not => add to vhandle
			isPatchHasPoint = 0;
			int isSame = 0;
			for (v_it = patch.vertices_begin(); v_it != patch.vertices_end(); v_it++) {
				isPatchHasPoint = 1;
				//check repeat vert
				if ((out_vertices.at(out_vertices.size() - 3) == patch.point(v_it.handle())[0] && out_vertices.at(out_vertices.size() - 2) == patch.point(v_it.handle())[1] && out_vertices.at(out_vertices.size() - 1) == patch.point(v_it.handle())[2])) {
					verticesSeq[verticesSeqIndex++] = out_vertices.size() - 3;
					mapVerticesToVHandle[out_vertices.size() - 3] = v_it.handle().idx();
					isSame = 1;
					break;
				}
			}
			if (isSame) continue;
			if (isPatchHasPoint == 0) { // first vert
				//printf("first point...\n");
				vhandle.push_back(patch.add_vertex(Tri_Mesh::Point(out_vertices.at(0), out_vertices.at(1), out_vertices.at(2))));
				verticesSeq[verticesSeqIndex++] = out_vertices.size() - 3; // record out_vert position
				mapVerticesToVHandle[out_vertices.size() - 3] = patch.vertices_begin().handle().idx(); // mapping id between patch & out_vert
				//printf("mapVToH = %d : %d\n", out_vertices.size() - 3, patch.vertices_begin().handle().idx());
			}
			else { // new point
				vhandle.push_back(patch.add_vertex(Tri_Mesh::Point(out_vertices.at(out_vertices.size() - 3), out_vertices.at(out_vertices.size() - 2), out_vertices.at(out_vertices.size() - 1))));
				verticesSeq[verticesSeqIndex++] = out_vertices.size() - 3;
				mapVerticesToVHandle[out_vertices.size() - 3] = (v_it).handle().idx();
				//printf("mapVToH = %d : %d\n", out_vertices.size() - 3, (v_it).handle().idx());
			}
		}
		//add face to patch from vhandle
		int num;
		face_vhandles.clear();
		for (int i = 0; i < 3; i++) {
			num = mapVerticesToVHandle.find(verticesSeq[i])->second;
			face_vhandles.push_back(vhandle[num]);
		}
		//std::cout << "add face to face_vhandles" << std::endl;
		patch.add_face(face_vhandles);
		verticesSeqIndex = 0;
	}
}
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

Tri_Mesh Tri_Mesh::simplify(float rate, float threshold) {
	Tri_Mesh* simplified = this;
	OpenMesh::EPropHandleT<double> cost;
	OpenMesh::VPropHandleT<Matrix4d> QMat;
	OpenMesh::EPropHandleT<Point> newPoint;
	simplified->add_property(cost);
	simplified->add_property(QMat);
	simplified->add_property(newPoint);


	auto t0 = high_resolution_clock::now();
	VIter v_it;
	EIter e_it;
	int vertexCount;
	// calculate the Q matrix for all vertices
	for (vertexCount = 0, v_it = simplified->vertices_begin(); v_it != simplified->vertices_end(); ++v_it, vertexCount++) {
		Matrix4d q = calculateQ(v_it.handle());
		simplified->property(QMat, *v_it) = q;
	}
	cout << "done preprocessing in " << duration<double>(high_resolution_clock::now() - t0).count() << " s" << endl;
	t0 = high_resolution_clock::now();

	auto update_edge = [&](EdgeHandle& eh) {
		VertexHandle to = to_vertex_handle(halfedge_handle(eh, 0));
		VertexHandle from = from_vertex_handle(halfedge_handle(eh, 0));

		Matrix4d newQ = simplified->property(QMat, to) + simplified->property(QMat, from);
		//cout << "newQ: " << newQ << endl;
		//cout << "oldQ: " << calculateQ(to) + calculateQ(from) << endl;;
		//Matrix4d newQ = calculateQ(to) + calculateQ(from);
		Matrix4d m = Matrix4d(newQ);

		m(3, 0) = 0.0f;
		m(3, 1) = 0.0f;
		m(3, 2) = 0.0f;
		m(3, 3) = 1.0f;

		//VectorXd x = m.fullPivLu().solve(Vector4d(0, 0, 0, 1));
		Vector4d x = m.inverse() * (Vector4d(0, 0, 0, 1));
		Vector4d newV;
		if (m.determinant() == 0.0f) {
			Point temp = (point(to) + point(from)) / 2;
			newV = Vector4d(temp[0], temp[1], temp[2], 1);
		}
		else {
			newV = x;
		}
		
		Point newP = Point(newV(0), newV(1), newV(2));	

		auto cur_cost = (newV.transpose() * newQ * newV)(0, 0);
		simplified->property(cost, eh) = cur_cost;
		simplified->property(newPoint, eh) = newP;
	};

	CompareCost compare = CompareCost(simplified, &cost);

	std::set<int, CompareCost> pq(compare);


	if (threshold == 0) {
		for (e_it = simplified->edges_begin(); e_it != simplified->edges_end(); ++e_it) {
			update_edge(e_it.handle());
			pq.insert(e_it.handle().idx());
		}
	}

	cout << "done pushing all edges in " << duration<double>(high_resolution_clock::now() - t0).count() << " s" << endl;
	t0 = high_resolution_clock::now();

	// collapse the edge with smallest cost
	// if the connected vertices form a concave polygon, ignore this edge
	// repeat until the vertex number is lower than the target number
	int targetVertexCount = vertexCount * rate;
	while (vertexCount > targetVertexCount) {
		if (pq.size() == 0) break;
		auto top = pq.begin();
		EdgeHandle eh = simplified->edge_handle(*top);

		VertexHandle from = from_vertex_handle(halfedge_handle(eh, 0));
		VertexHandle remain = to_vertex_handle(halfedge_handle(eh, 0));
		
		pq.erase(*top);

		if (!from.is_valid() || !remain.is_valid() || !eh.is_valid()) {
			continue;
		}

		if (is_collapse_ok(halfedge_handle(eh, 0))) {
			RollbackInfo* info1 = new RollbackInfo(from.idx(), this);
			RollbackInfo* info2 = new RollbackInfo(remain.idx(), this);
			
			// remove connected edges in pq
			VertexEdgeIter ve_it = simplified->ve_iter(from);
			for (; ve_it.is_valid(); ++ve_it) {
				pq.erase(ve_it.handle().idx());
			}

			ve_it = simplified->ve_iter(remain);
			for (; ve_it.is_valid(); ++ve_it) {
				pq.erase(ve_it.handle().idx());
			}

			// collapse
			simplified->collapse(halfedge_handle(eh, 0));

			// new point position
			point(remain) = simplified->property(newPoint, eh);
			//cout << "FINAL POINT: " << point(remain)<< " score: "<< simplified->property(cost, eh) << endl;
			from.invalidate();

			DecimationLog* log = new DecimationLog(info1, info2, point(remain));
			_deque.pushNewLog(log);

			// update the cost of connected edges and push them back to the pq
			ve_it = simplified->ve_iter(remain);
			for (; ve_it.is_valid(); ++ve_it) {
				update_edge(ve_it.handle());
				pq.insert(ve_it.handle().idx());
			}
			cout << "done an iteration in " << duration<double>(high_resolution_clock::now() - t0).count() << " s" << endl;
			t0 = high_resolution_clock::now();
			
			vertexCount--;
			_deque.toDecimate();
		}

	}

	cout << "FINISH" << endl;
	simplified->remove_property(cost);
	simplified->remove_property(QMat);
	simplified->remove_property(newPoint);
	simplified->garbage_collection();
	return *simplified;
}


Matrix4d Tri_Mesh::calculateQ(VertexHandle vhandle) {
	VVIter vv_it;
	Matrix4d q = Matrix4d();
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			q(i, j) = 0.0;
		}
	}
	Point v = point(vhandle);
	for (VFIter vf_it = vf_iter(vhandle); vf_it.is_valid(); ++vf_it) {
		Point p = normal(vf_it);
		
		double a = p[0];
		double b = p[1];
		double c = p[2];
		double d = -a * v[0] - b * v[1] - c * v[2];
		//cout << a << ", " << b << ", " << c << ", " << d << endl;
		q(0, 0) += a * a;
		q(0, 1) += a * b;
		q(0, 2) += a * c;
		q(0, 3) += a * d;

		q(1, 0) += b * a;
		q(1, 1) += b * b;
		q(1, 2) += b * c;
		q(1, 3) += b * d;

		q(2, 0) += c * a;
		q(2, 1) += c * b;
		q(2, 2) += c * c;
		q(2, 3) += c * d;

		q(3, 0) += d * a;
		q(3, 1) += d * b;
		q(3, 2) += d * c;
		q(3, 3) += d * d;
	}
	
	//p.normalize();
	//cout << "face normal: " << p << endl;
	//cout << q << endl << endl;
	return q;
}


Tri_Mesh Tri_Mesh::averageSimplify() {
	Tri_Mesh simplified = Tri_Mesh(*this);
	//only update once
	OpenMesh::VPropHandleT<bool> updated;
	simplified.add_property(updated);

	VIter v_it;
	EIter e_it;
	int vertexCount;
	// label all vert not updated
	for (vertexCount = 0, v_it = simplified.vertices_begin(); v_it != simplified.vertices_end(); ++v_it, vertexCount++) {
		simplified.property(updated, *v_it) = false;
	}
	for (vertexCount = 0, v_it = simplified.vertices_begin(); v_it != simplified.vertices_end(); ++v_it, vertexCount++) {
		if (simplified.property(updated, *v_it)) continue;
		cout << "simplify." << endl;
		VHandle pickedVert = v_it.handle();
		simplified.property(updated, *v_it) = true;
		/*
		float newVertX = point(pickedVert)[0];
		float newVertY = point(pickedVert)[1];
		float newVertZ = point(pickedVert)[2];
		int clusterCount = 1;
		std::vector<VHandle> clusterVert = std::vector<VHandle>();
		*/
		oneRingCollapse(pickedVert);
	}
	cout << "simplify complete" << endl;
	simplified.garbage_collection();
	return simplified;
}

void Tri_Mesh::oneRingCollapse(VHandle vhandle) {
	vector<HHandle> edges = vector<HHandle>();
	VIHIter vv_it;
	//find one ring
	for (vv_it = vih_iter(vhandle); vv_it; ++vv_it) {
		//cout << "collect." << endl;
		edges.push_back(vv_it.handle());
	}
	//collapse one ring
	if (edges.size() > 0) {
		for (size_t i = 0; i <= edges.size() - 1; i++) {
			if (i >= edges.size() || i < 0) { cout << "out of bound" << endl; break; }
			//if (!edges[i].is_valid) { cout << "unvalid edge" << endl; break; }
			//cout << i << "collapse." << endl;
			collapse(edges[i]);
		}
	}
	garbage_collection();
	return;
}


bool Tri_Mesh::DetermineConcaveByTwoPoints(std::vector<double> & p1, std::vector<double> & p2, std::vector<double> & vertices)
{
	std::vector<int> faceBuffer;
	//std::vector<TriVertex> vertexBuffer;
	std::vector<TriLine> edgeBuffer;
	TriVertex center1 = TriVertex(p1[0], p1[1], p1[2]);
	TriVertex center2 = TriVertex(p2[0], p2[1], p2[2]);
	for (int faceID = 0; faceID < vertices.size() / 9; faceID++)
	{
		for (int i = 0; i < 3; i++)
		{
			// 每個面有三個vertexes
			double vertexX = vertices[faceID * 9 + 3 * i];
			double vertexY = vertices[faceID * 9 + 3 * i + 1];
			double vertexZ = vertices[faceID * 9 + 3 * i + 2];
			if (p1[0] == vertexX && p1[1] == vertexY && p1[2] == vertexZ)
			{
				//比對是否屬於該平面
				if (std::find(faceBuffer.begin(), faceBuffer.end(), faceID) == faceBuffer.end())
				{
					faceBuffer.push_back(faceID);
					//加入該平面其餘2點
					TriVertex temp[2];
					for (int k = 1; k < 3; k++)
					{
						int offset = (i + k) % 3;
						double X = vertices[faceID * 9 + 3 * offset];
						double Y = vertices[faceID * 9 + 3 * offset + 1];
						double Z = vertices[faceID * 9 + 3 * offset + 2];
						temp[k - 1] = TriVertex(X, Y, Z);
					}
					if (!(temp[0] == center2 || temp[1] == center2))
					{
						edgeBuffer.push_back(TriLine(temp[0], temp[1]));
					}
				}
			}
			if (p2[0] == vertexX && p2[1] == vertexY && p2[2] == vertexZ)
			{
				if (std::find(faceBuffer.begin(), faceBuffer.end(), faceID) == faceBuffer.end())
				{
					faceBuffer.push_back(faceID);
					//加入該平面其餘2點
					TriVertex temp[2];
					for (int k = 1; k < 3; k++)
					{
						int offset = (i + k) % 3;
						double X = vertices[faceID * 9 + 3 * offset];
						double Y = vertices[faceID * 9 + 3 * offset + 1];
						double Z = vertices[faceID * 9 + 3 * offset + 2];
						temp[k - 1] = TriVertex(X, Y, Z);
					}
					if (!(temp[0] == center1 || temp[1] == center1))
					{
						edgeBuffer.push_back(TriLine(temp[0], temp[1]));
					}
				}
			}
		}
	}
	//face buffer
	/*cout << "faceBuffer: ";
	for (int i = 0; i < faceBuffer.size(); i++)
	{
	cout << faceBuffer[i] << ", ";
	}
	cout << endl;*/

	//edge buffer
	//cout << edgeBuffer.size() << endl;
	for (int i = 0; i < edgeBuffer.size(); i++)
	{
		//cout << edgeBuffer[i] << endl;
		glPointSize(8.0);
		glColor3f(1.0, 0.0, 0.0);
		glBegin(GL_POINTS);
		for (OMT::VIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {
			GLdouble po[3] = { edgeBuffer[i].p1.x, edgeBuffer[i].p1.y, edgeBuffer[i].p1.z };
			glVertex3dv(po);
		}
		glEnd();
	}
	int rounds = 0;
	int index = 0;
	while (rounds < edgeBuffer.size())
	{
		rounds++;
		TriLine line1 = edgeBuffer[index];
		//find next line 
		for (int i = 0; i < edgeBuffer.size(); i++)
		{
			if (edgeBuffer[i].p1 == line1.p2)
			{
				index = i;
				break;
			}
		}
		TriLine line2 = edgeBuffer[index];
		double crossValue = TriLine::cross(line1, line2);
		//cout << "crossValue: " << crossValue << endl;
		if (crossValue < 0)
			return false;
	}
	return true;
}


void Tri_Mesh::normalizeModel()
{
	double maxX = *(point(vertices_begin().handle()).data());
	double minX = *(point(vertices_begin().handle()).data());
	double maxY = *(point(vertices_begin().handle()).data() + 1);
	double minY = *(point(vertices_begin().handle()).data() + 1);
	double maxZ = *(point(vertices_begin().handle()).data() + 2);
	double minZ = *(point(vertices_begin().handle()).data() + 2);

	for (VertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
	{
		double posX = *(point(v_it.handle()).data());
		double posY = *(point(v_it.handle()).data() + 1);
		double posZ = *(point(v_it.handle()).data() + 2);
		maxX = std::max(posX, maxX);
		minX = std::min(posX, minX);
		maxY = std::max(posY, maxY);
		minY = std::min(posY, minY);
		maxZ = std::max(posZ, maxZ);
		minZ = std::min(posZ, minZ);
	}

	double scalar = std::max(maxX - minX, maxY - minY);
	scalar = std::max(scalar, maxZ - minZ);
	double center[3] = { (maxX + minX) / 2, (maxY + minY) / 2, (maxZ + minZ) / 2 };

	
	for (VertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
	{
		double posX = *(point(v_it.handle()).data());
		double posY = *(point(v_it.handle()).data() + 1);
		double posZ = *(point(v_it.handle()).data() + 2);
		Point newPoint((posX - center[0]) / scalar, (posY - center[1]) / scalar, (posZ - center[2]) / scalar);
		set_point(*v_it, newPoint);
		//cout << newPoint << endl;
		//cout << "Vertex: " << newPoint << endl;
	}
}


// Get Decimation Logs from the deque and decimate for k times.
// Return if k-times decimation is finished or the log is empty.
void Tri_Mesh::Decimate(int k) {
	while (k--) {
		DecimationLog* log = _deque.toDecimate();
		if (log == nullptr) {
			cout << "NULL" << endl;
			return;
		}

		VertexIter v_it;
		VertexEdgeIter ve_it;
		bool found = false;
		Point p1 = log->collapsePoint->p;
		Point p2 = log->remainPoint->p;

		// Find the corresponding edge to collapse
		for (v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {
			if (equal(point(v_it.handle()), p1)) {
				for (ve_it = ve_iter(v_it); ve_it.is_valid(); ++ve_it) {
					VertexHandle v1 = from_vertex_handle(halfedge_handle(ve_it, 0));
					VertexHandle v2 = to_vertex_handle(halfedge_handle(ve_it, 0));

					// Determine the halfedge direction
					if (equal(point(v1), p1) && equal(point(v2), p2)) {
						found = true;
						//cout << "collapse" << endl;
						collapse(halfedge_handle(ve_it, 0));
						point(v2) = log->remainPointNewPosition;
						break;
					}
					else if ((equal(point(v1), p2) && equal(point(v2), p1))) {
						found = true;
						//cout << "collapse" << endl;
						collapse(halfedge_handle(ve_it, 1));
						point(v1) = log->remainPointNewPosition;
						break;
					}
				}
			}
			if (found) break;
		}
	}
	cout << "FINISH" << endl;
}


// Get Recovering Logs from the deque and recover for k times.
// Return if k-times recovering is finished or the log is empty.
void Tri_Mesh::Recover(int k) {
	while (k--) {
		DecimationLog* log = _deque.toRecover();
		if (log == nullptr) {
			cout << "NULL" << endl;
			return;
		}

		vector<vector<VertexHandle>> allFaces;
		VertexIter v_it;
		Point target = log->remainPointNewPosition;
		for (v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {

			// Find the remain vertex after halfedge collapse
			if (equal(point(v_it.handle()), target)) {

				// Create two new vertices
				VertexHandle recoverV = add_vertex(log->collapsePoint->p);
				VertexHandle remain = add_vertex(log->remainPoint->p);

				// Recover all faces formed with the collapsed point
				int n = log->collapsePoint->neighborPos.size();
				for (int i = 0; i < n; i++) {
					vector<VertexHandle> face;
					face.push_back(recoverV);

					// There is no vertex at original position of the remain point, 
					// so we initialize the vertex handle with this value,
					// or we will miss the remain vertex
					double minL = (point(remain)- log->collapsePoint->neighborPos[i]).length();
					VertexHandle hd = remain, hd1 = remain;

					// Find which points we want to form a face
					for (VVIter vv_it = vv_iter(v_it); vv_it.is_valid(); ++vv_it) {
						double l = (point(vv_it) - log->collapsePoint->neighborPos[i]).length();
						if (l < minL) {
							minL = l;
							hd = vv_it.handle();
						}
					}
					face.push_back(hd);

					minL = (point(remain) - log->collapsePoint->neighborPos[(i + 1) % n]).length();
					for (VVIter vv_it = vv_iter(v_it); vv_it.is_valid(); ++vv_it) {
						double l = (point(vv_it) - log->collapsePoint->neighborPos[(i + 1) % n]).length();
						if (l < minL) {
							minL = l;
							hd1 = vv_it.handle();
						}
					}
					face.push_back(hd1);

					if (face.size() != 3) {
						cout << "ERROR SIZE " << face.size() << endl;
						continue;
					}
					else {
						// save the result to a temporary vector
						allFaces.push_back(face);
					}
				}


				// Recover all faces formed with the remain point
				n = log->remainPoint->neighborPos.size();
				for (int i = 0; i < n; i++) {
					vector<VertexHandle> face;
					face.push_back(remain);

					// There is no edge connected with the recovered vertex, 
					// so we initialize the vertex handle with this value,
					// or we will miss the recovered vertex
					double minL = (point(recoverV) - log->remainPoint->neighborPos[i]).length();
					VertexHandle hd = recoverV, hd1 = recoverV;

					// Find which points we want to form a face
					for (VVIter vv_it = vv_iter(v_it); vv_it.is_valid(); ++vv_it) {
						double l = (point(vv_it) - log->remainPoint->neighborPos[i]).length();
						if (l < minL) {
							minL = l;
							hd = vv_it.handle();
						}
					}
					face.push_back(hd);

					minL = (point(recoverV) - log->remainPoint->neighborPos[(i + 1) % n]).length();
					for (VVIter vv_it = vv_iter(v_it); vv_it.is_valid(); ++vv_it) {
						double l = (point(vv_it) - log->remainPoint->neighborPos[(i + 1) % n]).length();
						if (l < minL) {
							minL = l;
							hd1 = vv_it.handle();
						}
					}
					face.push_back(hd1);

					if (face.size() != 3) {
						cout << "ERROR SIZE" << endl;
						continue;
					}
					else {
						// Save the result to a temporary vector
						allFaces.push_back(face);
					}
				}

				break;

			}
		}

		// Delete the original remain vertex before we add the face back,
		// or we will get some errors
		delete_vertex(v_it);
		for (auto face: allFaces) {
			add_face(face);
		}
	}

	
	garbage_collection();
	cout << "FINISH" << endl;
}