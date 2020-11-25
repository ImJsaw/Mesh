#include "GUA_OM.h"
#include <algorithm>
#include <iostream>


using namespace std;

class TriVertex {
public:
	double x;
	double y;
	double z;
	TriVertex() {
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	TriVertex(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	~TriVertex() {}

	bool operator ==(const TriVertex& b) {
		if (this->x == b.x && this->y == b.y && this->z == b.z) {
			return true;
		}
		else {
			return false;
		}
	}
	friend ostream &operator<<(ostream& os, const TriVertex& tv) {
		os << "(" << tv.x << ", " << tv.y << ", " << tv.z << ")";
		return os;
	}

};
class TriLine {
public:
	TriVertex p1;
	TriVertex p2;
	double vector[3];

	TriLine(TriVertex &p1, TriVertex &p2) {
		this->p1 = p1;
		this->p2 = p2;
		vector[0] = p2.x - p1.x;
		vector[1] = p2.y - p1.y;
		vector[2] = p2.z - p1.z;
	}
	~TriLine() {}
	friend ostream &operator<<(ostream& os, const TriLine& line) {
		os << "[" << line.p1 << " -> " << line.p2 << "]";
		return os;
	}
	double static cross(TriLine & line1, TriLine& line2) {
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
		//�o�q�O���F�ѨMindex���D
		VHandle vq, vw, vt, vr;
		vq = addVertex(Point(0, 0, 100));
		vw = addVertex(Point(1, 0, 100));
		vt = addVertex(Point(1, 1, 100));
		vr = addVertex(Point(0, 1, 100));
		addFace(vq, vw, vt, vr);
		/*----------------------------------------------------------------------*/
		/*�����ݭnsubdivision��face*/
		//��ltable
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
		/*��candidate faces��subdivision*/
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
		deleteFace(vq, vw, vt, vr);//�o��u�O���F�ѨMindex���D
								   /*----------------------------------------------------------------------*/
		return 0;
	}
	int Model::quad_subdivision2(int _face_id) {
		/*----------------------------------------------------------------------*/
		//�o�q�O���F�ѨMindex���D
		VHandle vq, vw, vt, vr;
		vq = addVertex(Point(0, 0, 100));
		vw = addVertex(Point(1, 0, 100));
		vt = addVertex(Point(1, 1, 100));
		vr = addVertex(Point(0, 1, 100));
		addFace(vq, vw, vt, vr);
		/*----------------------------------------------------------------------*/
		/*�����ݭnsubdivision��face*/
		//��ltable
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
		/*��candidate faces��subdivision*/
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
		deleteFace(vq, vw, vt, vr);//�o��u�O���F�ѨMindex���D
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
			// �C���I���T��vertexes
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

		for (fv_it = ++fv_iter(f_it); fv_it; ++fv_it) { // �C�ӭ����T���I
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
	//Step1 : �����Ĥ@�����
	for (e_it = edges_begin(); e_it != edges_end(); e_it++) {
		printf("enter edge iterator...\n");
		bool isBoundary = patch.is_boundary(*e_it);
		if (isBoundary) {
			heh = patch.halfedge_handle(e_it.handle(), 1);
			break;
		}
	}
	// Step2 : �u�۸���ɧ�M�U�@����
	double perimeter = 0;
	std::vector<double> segLength;
	std::vector<Tri_Mesh::VertexHandle> vhs; // �x�s�ƧǦn������I
	HHandle hehNext = heh;
	do {
		Point from = patch.point(patch.from_vertex_handle(hehNext));
		Point to = patch.point(patch.to_vertex_handle(hehNext));
		perimeter += (from - to).length(); // v0 - v? ������
		printf("perimeter = %f\n", perimeter);
		segLength.push_back(perimeter); // �s�Jvector���A�H�K���ᰵtexcoord
		vhs.push_back(patch.from_vertex_handle(hehNext)); // ����ɤW���I�@�@�s�J
		hehNext = patch.next_halfedge_handle(hehNext); // �i�H�����̧ǩ��U�@��heh��
	} while (heh != hehNext);

	//Step3 : �䧹�Ҧ�������I�M�Z����A�i�H�N��ɤw���I����texcoord
	float rd = (225 + uvRotateAngle) * M_PI / 180.0;
	float initDist = 0;
	Tri_Mesh::TexCoord2D st(0, 0);
	float R = std::sqrt(2) / 2.0; // �ڸ�2/2
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
		double curLen = segLength[i - 1] / perimeter + initDist; // �N�YL0n/Ltotal*4
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
		patch.set_texcoord2D(vhs[i], st); // ���D�X����I��uv�y�СA�����w���I
	}

	//Step4 : ���U�ӡA��X�D��ɪ��v��(�ϥΧ��褽��)
	OpenMesh::HPropHandleT<double> heWeight;
	OpenMesh::VPropHandleT<int> row;
	patch.add_property(heWeight, "heWeight"); //�[�Jproperty�G�v��
	patch.add_property(row, "row");

	for (e_it = patch.edges_begin(); e_it != patch.edges_end(); ++e_it) {
		bool isBound = patch.is_boundary(*e_it);
		if (!isBound) { // �p�G���O��ɡA�N�n��L���v���A���half_edge���n
			GLdouble angle1, angle2, w;
			Tri_Mesh::HalfedgeHandle _heh = patch.halfedge_handle(e_it.handle(), 0);
			Point pFrom = patch.point(patch.from_vertex_handle(_heh)); // �_�I
			Point pTo = patch.point(patch.to_vertex_handle(_heh)); // ���I
																   //���D�Gopposite_he_opposite_vh, opposite_vh���
			Point p1 = patch.point(patch.opposite_vh(_heh)); // �u�W���I
			Point p2 = patch.point(patch.opposite_he_opposite_vh(_heh)); // �u�U���I
																		 //auto _vh = patch.opposite_vh(_heh);
			double edgeLen = (pFrom - pTo).length(); //�u����
			OpenMesh::Vec3d v1 = (OpenMesh::Vec3d)(pTo - pFrom);
			v1.normalize();
			OpenMesh::Vec3d v2 = (OpenMesh::Vec3d)(p1 - pFrom);
			v2.normalize();
			angle1 = std::acos(OpenMesh::dot(v1, v2)); //�u�W�I�P�u��arccos�A���o����

			v2 = (OpenMesh::Vec3d)(p2 - pFrom);
			v2.normalize();

			angle2 = std::acos(OpenMesh::dot(v1, v2)); //�u�U�I�P�u��arccos�A���o����
													   //���褽��
			w = (std::tan(angle1 / 2.0f) + std::tan(angle2 / 2.0f)) / edgeLen;
			patch.property(heWeight, _heh) = w;  //�N�䤤�@��edge��n��J�䪺�ʽ褤

												 // �A�ӭp��Ϥ�V��halfedge weight
			v1 = -v1;
			v2 = (OpenMesh::Vec3d)(p1 - pTo);
			v2.normalize();
			angle1 = std::acos(OpenMesh::dot(v1, v2));//�u�W�I�P�u��arccos�A���o����

			v2 = (OpenMesh::Vec3d)(p2 - pTo);
			v2.normalize();
			angle2 = std::acos(OpenMesh::dot(v1, v2));//�u�U�I�P�u��arccos�A���o����

			w = (std::tan(angle1 / 2.0f) + std::tan(angle2 / 2.0f)) / edgeLen;
			patch.property(heWeight, patch.opposite_halfedge_handle(_heh)) = w;//�N�Ϥ�Vedge��n��J�䪺�ʽ褤
		}
	}
	//Step5 : ��X�x�}���j�p(NxN)�A�`�I��-����I
	int count = 0;
	for (v_it = patch.vertices_begin(); v_it != patch.vertices_end(); v_it++) { // debug�Gpatch.vertices_end()�n�[patch
		if (patch.is_boundary(*v_it)) patch.property(row, *v_it) = -1;
		else patch.property(row, *v_it) = count++;
	}
	//Step6 : ��g�x�}
	typedef Eigen::SparseMatrix<double> SpMat;
	SpMat A(count, count);
	Eigen::VectorXd BX(count);
	Eigen::VectorXd BY(count);//x �M y �U�Ѥ@��
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > linearSolver;
	BX.setZero();
	BY.setZero();
	// fiil matrix
	for (v_it = patch.vertices_begin(); v_it != patch.vertices_end(); ++v_it) {
		if (!patch.is_boundary(*v_it)) { // �����I�O�ڭ̭n�Ѫ��F��
			int i = patch.property(row, *v_it);
			double totalWeight = 0;
			for (vv_it = patch.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
				HHandle _heh = patch.find_halfedge(*v_it, *vv_it);
				double w = patch.property(heWeight, _heh);

				if (patch.is_boundary(*vv_it)) { //�O��ɪ����w���A��JB�x�}
					TexCoord2D texCoord = patch.texcoord2D(*vv_it);
					BX[i] += w * texCoord[0]; // �����ۥ[
					BY[i] += w * texCoord[1];
				}
				else { // ���O��ɪ��������I�A��JA�x�}
					int j = patch.property(row, *vv_it);
					A.insert(i, j) = -w;
				}
				totalWeight += w;
			}
			A.insert(i, i) = totalWeight; // �o�̪��Ѫk�������Ҧ���w���W�jW
		}
	}
	A.makeCompressed();
	//Step7 : �Ѷ}������uv�y�Э�
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

void Tri_Mesh::Update_Edge(EdgeHandle eh) {
	VertexHandle to = to_vertex_handle(halfedge_handle(eh, 0));
	VertexHandle from = from_vertex_handle(halfedge_handle(eh, 0));

	Matrix4d newQ = this->property(QMat, to) + this->property(QMat, from);
	Matrix4d m = Matrix4d(newQ);

	m(3, 0) = 0.0f;
	m(3, 1) = 0.0f;
	m(3, 2) = 0.0f;
	m(3, 3) = 1.0f;
	Vector4d newV;
	if (m.determinant() == 0.0f) {
		Point temp = (point(to) + point(from)) / 2;
		newV = Vector4d(temp[0], temp[1], temp[2], 1);
	}
	else {
		Vector4d x = m.inverse() * (Vector4d(0, 0, 0, 1));
		newV = x;
	}
	
	Point newP = Point(newV(0), newV(1), newV(2));

	//double cur_cost = (newV.transpose() * newQ * newV)(0, 0);
	double newcost = 0;
	double VtQ[4];
	
	for (int i = 0; i < 4; i++)
	{
		VtQ[i] = newV(0) * newQ(0, i) + newV(1) * newQ(1, i) + newV(2) * newQ(2, i) + newV(3) * newQ(3, i);
	}
	
	for (int i = 0; i < 4; i++)
	{
		newcost += VtQ[i] * newV(i);
	}
	//cout << cur_cost << ", " << newcost << endl;
	//this->property(cost, eh) = cur_cost;
	this->property(cost, eh) = newcost;
	this->property(newPoint, eh) = newP;

	
}

void Tri_Mesh::simplify(float rate) {
	clock_t t1;
	t1 = clock();
	cout << (double)(t1) / CLOCKS_PER_SEC << endl;
	int vertexCount = this->n_vertices();
	int targetVertexCount = vertexCount * rate;

	while (vertexCount > targetVertexCount && this->_deque.canDecimate()) {
		this->Decimate(1);
		vertexCount--;
	}
	if (vertexCount <= targetVertexCount) return;

	VIter v_it;
	EIter e_it;

	CompareCost compare = CompareCost(this, &cost);
	std::set<int, CompareCost> pq(compare);

	for (e_it = this->edges_begin(); e_it != this->edges_end(); ++e_it) {
		Update_Edge(e_it.handle());
		pq.insert(e_it.handle().idx());
	}

	// collapse the edge with smallest cost
	// if the connected vertices form a concave polygon, ignore this edge
	// repeat until the vertex number is lower than the target number
	while (vertexCount > targetVertexCount) {
		if (pq.size() == 0) break;

		auto top = pq.begin();
		EdgeHandle eh = this->edge_handle(*top);

		VertexHandle from = from_vertex_handle(halfedge_handle(eh, 0));
		VertexHandle remain = to_vertex_handle(halfedge_handle(eh, 0));
		Point np = this->property(newPoint, eh);

		// pop the edge with the smallest cost
		pq.erase(*top);

		if (!from.is_valid() || !remain.is_valid() || !eh.is_valid()) {
			continue;
		}

		if (is_collapse_ok(halfedge_handle(eh, 0)) && DetermineConcaveByTwoPoints(&from, &remain, &np)) {
			Point np = this->property(newPoint, eh);
			RollbackInfo* info1 = new RollbackInfo(from.idx(), this);
			RollbackInfo* info2 = new RollbackInfo(remain.idx(), this);

			// remove connected edges in pq
			VertexEdgeIter ve_it = this->ve_iter(from);
			for (; ve_it.is_valid(); ++ve_it) {
				pq.erase(ve_it.handle().idx());
			}

			ve_it = this->ve_iter(remain);
			for (; ve_it.is_valid(); ++ve_it) {
				pq.erase(ve_it.handle().idx());
			}

			// collapse
			this->collapse(halfedge_handle(eh, 0));

			// new point position
			point(remain) = this->property(newPoint, eh);
			from.invalidate();

			// add new log to the deque
			DecimationLog* log = new DecimationLog(info1, info2, point(remain));
			_deque.pushNewLog(log);

			// update the cost of connected edges and push them back to the pq
			ve_it = this->ve_iter(remain);
			for (; ve_it.is_valid(); ++ve_it) {
				Update_Edge(ve_it.handle());
				pq.insert(ve_it.handle().idx());
			}

			vertexCount--;
			_deque.toDecimate();
		}
	}

	cout << "FINISH" << endl;
	this->garbage_collection();
	return;
}


Matrix4d Tri_Mesh::calculateQ(VertexHandle vhandle) {
	VVIter vv_it;
	Matrix4d q = Matrix4d();

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			q(i, j) = 0.0;
		}
	}
	Point v = point(vhandle);
	//#pragma omp parallel for
	for (VFIter vf_it = vf_iter(vhandle); vf_it.is_valid(); ++vf_it) {
		Point p = normal(vf_it);

		double a = p[0];
		double b = p[1];
		double c = p[2];
		double d = -a * v[0] - b * v[1] - c * v[2];

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

	return q;
}

Matrix4d Tri_Mesh::calculateFQ(VertexHandle vhandle) {
	VVIter vv_it;
	Matrix4d q = Matrix4d();
	//init matrix
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			q(i, j) = 0.0;
		}
	}
	
	Point vi = point(vhandle);

	for (VVIter vv_it = vv_iter(vhandle); vv_it.is_valid(); ++vv_it) {
		Point vj = point(vv_it);
		// a = normalized edge vector of edge ij, b = a * vi
		Point a = (vj - vi).normalized();
		Point b = a % vi;
		MatrixXd kij = MatrixXd(3, 4);
		//init matrix
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				kij(i, j) = 0.0;
			}
		}

		kij(0, 1) -= a[2];
		kij(0, 2) += a[1];
		kij(0, 3) -= b[0];

		kij(1, 0) += a[2];
		kij(1, 2) -= a[0];
		kij(1, 3) -= b[1];

		kij(2, 0) -= a[1];
		kij(2, 1) += a[0];
		kij(2, 3) -= b[2];
		Matrix4d tmpQ= kij.transpose() * kij;
		q += tmpQ;
	}

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


bool Tri_Mesh::DetermineConcaveByTwoPoints(VertexHandle *p1, VertexHandle *p2, Point* np) {
	Point v1 = point((*p1));
	Point v2 = point((*p2));
	for (VertexFaceIter vf_it = vf_begin(*p1); vf_it != vf_end(*p1); ++vf_it) {
		vector<Point> oldVector;
		vector<Point> newVector;
		for (FaceVertexIter fv_it = fv_iter(vf_it); fv_it; fv_it++) {
			Point p = point(fv_it);
			if (p == v1) {
				continue;
			}
			else {
				oldVector.push_back((v1 - p).normalize());
				newVector.push_back(((*np) - p).normalize());
			}
		}
		double oldCross[3] = { oldVector[0][1] * oldVector[1][2] - oldVector[0][2] * oldVector[1][1],
			oldVector[0][0] * oldVector[1][2] - oldVector[0][2] * oldVector[1][0],
			oldVector[0][0] * oldVector[1][1] - oldVector[0][1] * oldVector[1][0] };
		double newCross[3] = { newVector[0][1] * newVector[1][2] - newVector[0][2] * newVector[1][1],
			newVector[0][0] * newVector[1][2] - newVector[0][2] * newVector[1][0],
			newVector[0][0] * newVector[1][1] - newVector[0][1] * newVector[1][0] };
		double dot = oldCross[0] * newCross[0] + oldCross[1] * newCross[1] + oldCross[2] * newCross[2];
		if (dot < 0) {
			return false;
		}
	}


	for (VertexFaceIter vf_it = vf_begin(*p2); vf_it != vf_end(*p2); ++vf_it) {
		vector<Point> oldVector;
		vector<Point> newVector;

		for (FaceVertexIter fv_it = fv_iter(vf_it); fv_it; fv_it++) {
			Point p = point(fv_it);
			if (p == v2) {
				continue;
			}
			else {
				oldVector.push_back((p - v2).normalize());
				newVector.push_back((p - (*np)).normalize());
			}
		}
		double oldCross[3] = { oldVector[0][1] * oldVector[1][2] - oldVector[0][2] * oldVector[1][1],
			oldVector[0][0] * oldVector[1][2] - oldVector[0][2] * oldVector[1][0],
			oldVector[0][0] * oldVector[1][1] - oldVector[0][1] * oldVector[1][0] };
		double newCross[3] = { newVector[0][1] * newVector[1][2] - newVector[0][2] * newVector[1][1],
			newVector[0][0] * newVector[1][2] - newVector[0][2] * newVector[1][0],
			newVector[0][0] * newVector[1][1] - newVector[0][1] * newVector[1][0] };
		double dot = oldCross[0] * newCross[0] + oldCross[1] * newCross[1] + oldCross[2] * newCross[2];
		if (dot < 0) {
			return false;
		}
	}

	return true;
}


void Tri_Mesh::normalizeModel() {
	double maxX = *(point(vertices_begin().handle()).data());
	double minX = *(point(vertices_begin().handle()).data());
	double maxY = *(point(vertices_begin().handle()).data() + 1);
	double minY = *(point(vertices_begin().handle()).data() + 1);
	double maxZ = *(point(vertices_begin().handle()).data() + 2);
	double minZ = *(point(vertices_begin().handle()).data() + 2);

	for (VertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {
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


	for (VertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {
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
			garbage_collection();
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

	garbage_collection();
	cout << "FINISH" << endl;
}


// Get Recovering Logs from the deque and recover for k times.
// Return if k-times recovering is finished or the log is empty.
void Tri_Mesh::Recover(int k) {
	while (k--) {
		DecimationLog* log = _deque.toRecover();
		if (log == nullptr) {
			cout << "NULL" << endl;
			garbage_collection();
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
					double minL = (point(remain) - log->collapsePoint->neighborPos[i]).length();
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
		omerr().disable();
		for (auto face : allFaces) {
			add_face(face);
		}
	}


	garbage_collection();
	cout << "FINISH" << endl;
}

void Tri_Mesh::Initialize() {
	this->normalizeModel();
	this->add_property(QEMcost);
	this->add_property(cost);
	this->add_property(QMat);
	this->add_property(FMat);
	this->add_property(newPoint);

	// calculate the Q matrix for all vertices
	for (VIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it) {
		Matrix4d q = calculateQ(v_it.handle());
		this->property(QMat, *v_it) = q;
		this->property(FMat, *v_it) = calculateFQ(v_it.handle());
	}
	//calc all QEMCost  at init

	cout << "calc QEM" << endl;
	for (HIter h_it = halfedges_begin(); h_it != halfedges_end(); ++h_it) {
		this->property(QEMcost, *h_it) = qemCost(h_it.handle());
	}
}

bool cmp(pair<int, double> &a, pair<int, double> &b) { return a.second < b.second; }

double Tri_Mesh::qemCost(const HHandle hhandle) {
	VertexHandle from = from_vertex_handle(hhandle);
	VertexHandle to = to_vertex_handle(hhandle);
	double wa = 1;
	double wb = 20;
	//fij + fjj
	double fa = getFa(from,to) + getFa(to,to);
	double fb;
	double oneRingLenSigma=0;

	for (VEIter ve_it = ve_iter(from); ve_it.is_valid(); ++ve_it) {
		oneRingLenSigma += calc_edge_length(ve_it.handle());
	}
	fb = calc_edge_length(hhandle) * oneRingLenSigma;
	return (wa * fa + wb * fb);
}

double Tri_Mesh::getFa(const VHandle f,const VHandle vHandle) {
	//pt * Q * p
	Point v = point(vHandle);
	Vector4d P = Vector4d(v[0], v[1], v[2], 1);
	Matrix4d newQ = this->property(FMat, f);
	double VtQ[4];
	double Fa = 0;
	for (int i = 0; i < 4; i++) {
		VtQ[i] = P(0) * newQ(0, i) + P(1) * newQ(1, i) + P(2) * newQ(2, i) + P(3) * newQ(3, i);
	}

	for (int i = 0; i < 4; i++) {
		Fa += VtQ[i] * P(i);
	}
	return Fa;
}

void Tri_Mesh::face2Edge(int faces) {
	
	//single time collapse amount
	double collapseForce = 1;
	HIter h_it;
	std::vector<pair<int, double>> costMap;
	
	std::vector<int> droped = vector<int>();

	CompareSkeCost compare = CompareSkeCost(this, &QEMcost);
	std::set<int, CompareSkeCost> pq(compare);
	//log all QEM cost precalc
	for (HIter h_it = halfedges_begin(); h_it != halfedges_end(); ++h_it) {
		pq.insert(h_it.handle().idx());
	}
	int curFace = this->n_faces();
	while (curFace > faces) {

		//cout << "remove." << endl;
		if (pq.size() == 0) break;
		auto top = pq.begin();
		HHandle hh = this->halfedge_handle(*top);

		VertexHandle from = from_vertex_handle(hh);
		VertexHandle to = to_vertex_handle(hh);
		// pop the edge with the smallest cost
		pq.erase(pq.begin());
		if (!from.is_valid() || !to.is_valid() || !hh.is_valid()) {
			cout << "invalid!" << endl;
			continue;
		}

		if (is_collapse_ok(hh)) {
		//if (is_collapse_ok(hh)) {

			VVIter vv_it;

			//cout << "collapse" << endl;
			vector<int> toUpdate = vector<int>();
			// log one ring vert, need to update cost
			for (vv_it = vv_iter(from); vv_it.is_valid(); ++vv_it) {
				toUpdate.push_back(vv_it.handle().idx());
			}
			toUpdate.push_back(to.idx());
			//update Q
			this->property(FMat, to) = this->property(FMat, from) + this->property(FMat, to);

			//cout << "collapsed" << endl;
			// collapse
			this->collapse(hh);
			from.invalidate();


			for (VOHIter voh_it = voh_iter(to); voh_it; ++voh_it) {
				auto index = pq.find(voh_it.handle().idx());
				if(index != pq.end())
					pq.erase(index);
				//pq.erase(pq.find(voh_it.handle().idx()));
			}
			// update the cost of connected edges and push them back to the pq
			//cout << "update one ring & collapse point cost" << endl;
			for (int i = 0; i < toUpdate.size(); i++) {
				int idx = toUpdate[i];
				for (VOHIter voh_it = voh_iter(vertex_handle(idx)); voh_it; ++voh_it) {
					auto index = pq.find(voh_it.handle().idx());
					if (index != pq.end())
						pq.erase(index);
					//pq.erase(pq.find(voh_it.handle().idx()));
				}
			}


			for (VOHIter voh_it = voh_iter(to); voh_it; ++voh_it) {
				//cout << "cur half edge" << voh_it.handle().idx() << endl;
				this->property(QEMcost, *voh_it) = qemCost(voh_it.handle());
				//cout << "re insert" << endl;
				//re insert to update sort
				pq.insert(voh_it.handle().idx());
			}
			// update the cost of connected edges and push them back to the pq
			//cout << "update one ring & collapse point cost" << endl;
			for (int i = 0; i < toUpdate.size(); i++) {
				int idx = toUpdate[i];
				//recalc one ring matrix
				this->property(FMat, vertex_handle(idx)) = calculateFQ(vertex_handle(idx));
				//cout << "cur vert" << idx << endl;
				for (VOHIter voh_it = voh_iter(vertex_handle(idx)); voh_it; ++voh_it) {
					//cout << "cur half edge" << voh_it.handle().idx() << endl;
					this->property(QEMcost, *voh_it) = qemCost(voh_it.handle());
					//cout << "re insert" << endl;
					//re insert to update sort
					pq.insert(voh_it.handle().idx());
				}
			}
			curFace--;
		}

	}

	cout << "FINISH" << endl;
	this->garbage_collection();
	return;
}

void Tri_Mesh::getSkeleton() {
	cout << "get skeleton" << endl;
	const int smoothTime = 5;
	//model avg area
	double avgArea;
	//get avg face area
	int faceCount = 0;
	double areaSum = 0;
	const int vertices = this->n_vertices();
	OpenMesh::FPropHandleT<double> oriArea;
	OpenMesh::FPropHandleT<double> curArea;
	this->add_property(oriArea);
	this->add_property(curArea);
	w0 = 1.0;
	//WH init
	_WH.resize(vertices);

	for (FIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it) {
		faceCount++;
		double area = getArea(f_it.handle());
		property(oriArea, f_it.handle()) = area;
		property(curArea, f_it.handle()) = area;
		areaSum += area;
	}
	avgArea = areaSum / (double)faceCount;

	for (int i = 0; i < _WH.size(); i++)
		_WH[i] = w0;
	//WL init
	_WL = 5;

	cout << "avg area:" << avgArea << "from " << faceCount << "face, area total" << areaSum << endl;

	for (int i = 0; i < smoothTime; i++) {
		cout << "--" << endl << "smooth" << i << endl;
		//update area
		areaSum = 0;
		for (FIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it) {
			double area = getArea(f_it.handle());
			property(curArea, f_it.handle()) = area;
			areaSum += area;
		}
		cout << "area ratio : " << areaSum / avgArea / faceCount << endl;
		if (areaSum / avgArea / faceCount < 0.4) {
			cout << "ratio :" << areaSum / avgArea / faceCount << ", stop" << endl;
			break;
		}

		vector<VectorXd> newV = getNewVert(_WL, _WH);
		cout << "updating vert.." << endl;

		//extract new vert from matrix.
		for (VIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it) {
			int index = v_it->idx();
			Point p = Point();
			p[0] = newV[0][index];
			p[1] = newV[1][index];
			p[2] = newV[2][index];
			point(v_it.handle()) = p;
			//this->set_point(*v_it, p);
			//cout << "update point " << index << "to " << p << endl;
		}


		//cout << "update." << endl;
		//update WH, WL
		_WH.resize(vertices);
		for (VIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it) {
			int index = v_it->idx();
			double originArea = 0;
			double currentArea = 0;
			//collect one ring
			for (VFIter vf_it = vf_iter(v_it.handle()); vf_it; ++vf_it) {
				originArea += property(oriArea, vf_it.handle());
				currentArea += property(curArea, vf_it.handle());
			}
			if (currentArea == 0 || originArea == 0) {
				cout << "!!!!!!!!! area zero error!!!!!!" << endl;
				cout << "cur" << currentArea << "ori" << originArea << endl;
			}
			else {
				_WH[index] = sqrt(originArea / currentArea);
			}
		}
		_WL *= _SL;
	}
	garbage_collection();
	return;
}

vector<VectorXd> Tri_Mesh::getNewVert(double WL, vector<double> WH) {
	// (WL) L   V' = 0
	// WH		V' = WH V

	SparseMatrix<double> L = calculateL()*WL;
	SparseMatrix<double> WHMatrix(L.rows(), L.cols());

	const int col = L.cols();
	for (int colIndex = 0; colIndex < col; colIndex++) {
		//cout << "assign WH " << Left.rows()+ colIndex << "," << colIndex << ":"<< WH[colIndex] << endl;
		//assign WH
		WHMatrix.insert(colIndex, colIndex) = WH[colIndex];
	}
	SparseMatrix<double> Left(L.rows() + WHMatrix.rows(), L.cols());
	Left.reserve(L.nonZeros() + WHMatrix.nonZeros());
	for (Index c = 0; c < L.cols(); ++c) {
		Left.startVec(c); // Important: Must be called once for each column before inserting!
		for (SparseMatrix<double>::InnerIterator itL(L, c); itL; ++itL)
			Left.insertBack(itL.row(), c) = itL.value();
		for (SparseMatrix<double>::InnerIterator itC(WHMatrix, c); itC; ++itC)
			Left.insertBack(itC.row() + L.rows(), c) = itC.value();
	}
	Left.finalize();
	Left.makeCompressed();
	cout << "compressed" << endl;
	//right matrix
	const int vertices = this->n_vertices();
	vector<VectorXd> right;
	right.resize(3);
	right[0].resize(vertices * 2);
	right[1].resize(vertices * 2);
	right[2].resize(vertices * 2);

	for (int i = 0; i < vertices; i++) {
		right[0][i] = 0;
		right[1][i] = 0;
		right[2][i] = 0;
	}

	//iterate all vert, assign B
	for (VIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it) {
		int index = v_it->idx();
		//first  n row 0 , offset n
		right[0][vertices + index] = point(v_it.handle())[0] * WH[index];
		right[1][vertices + index] = point(v_it.handle())[1] * WH[index];
		right[2][vertices + index] = point(v_it.handle())[2] * WH[index];
	}

	cout << "solving linear system..." << endl;
	//solve linear system.
	SimplicialLDLT<SparseMatrix<double>> linearSolver;
	linearSolver.compute(Left.transpose() * Left);
	vector<VectorXd> newVert;
	newVert.resize(3);
	newVert[0] = linearSolver.solve(Left.transpose() * right[0]);
	newVert[1] = linearSolver.solve(Left.transpose() * right[1]);
	newVert[2] = linearSolver.solve(Left.transpose() * right[2]);
	return newVert;
}

SparseMatrix<double> Tri_Mesh::calculateL() {
	const int edge_size = this->n_edges();
	SparseMatrix<double> L = prepareLaplacian();
	VIter v_it;
	VOHIter ve_it;
	//i
	//iterate every vert
	for (v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it) {
		int index = v_it->idx();
		float sum = 0.0;
		int count = 0;
		//iterate all edge out from vert
		for (ve_it = voh_iter(v_it.handle()); ve_it; ++ve_it) {
			count++;
			//cot => i,j
			float cot = getCot(ve_it.handle());
			const int toID = to_vertex_handle(ve_it.handle()).idx();
			//cout << " L[" << index << "][" << toID<< "]="<<cot << endl;
			L.coeffRef(index, toID) = cot;
			//L[index][toID] = cot;
			sum -= cot;
		}
		// i,i
		//cout << " L[" << index << "][" << index << "]=" << sum / count << endl;
		L.coeffRef(index, index) = sum;
		//L[index][index] = sum / count;
	}
	return L;
}

SparseMatrix<double> Tri_Mesh::prepareLaplacian() {
	const int vertices = this->n_vertices();
	SparseMatrix<double> A(vertices, vertices);
	VectorXi sizes(vertices);
	vector<vector<VHandle>> neighbors(vertices);
	for (int i = 0; i < vertices; i++) {
		//assume max 10 connected edge per vert
		neighbors[i].reserve(10);
		for (VHandle vv : this->vv_range(VHandle(i))) {
			neighbors[i].push_back(vv);
		}
		//plus one for L(i,i)
		sizes[i] = neighbors[i].size() + 1;
	}

	//assign memory for element
	A.reserve(sizes);
	for (int i = 0; i < vertices; i++) {
		for (int j = 0; j < neighbors[i].size(); j++) {
			A.insert(i, neighbors[i][j].idx());
		}
		A.insert(i, i);
	}
	return A;
}

double Tri_Mesh::getArea(const FHandle f) {
	FVIter fv_it = fv_iter(f);
	const Point& P = point(fv_it);  ++fv_it;
	const Point& Q = point(fv_it);  ++fv_it;
	const Point& R = point(fv_it);
	return ((Q - P) % (R - P)).norm() * 0.5f;
}

float Tri_Mesh::getCot(const HHandle eh) {
	const HHandle he01 = eh;
	const HHandle he10 = opposite_halfedge_handle(he01);;
	const HHandle he12 = next_halfedge_handle(he01);
	const HHandle he03 = next_halfedge_handle(he10);
	const VHandle v0 = from_vertex_handle(he01);
	const VHandle v1 = to_vertex_handle(he01);
	const VHandle v2 = to_vertex_handle(he12);
	const VHandle v3 = to_vertex_handle(he03);
	const Point p0 = point(v0);
	const Point p1 = point(v1);
	const Point p2 = point(v2);
	const Point p3 = point(v3);

	float cot = 0.0;
	if (!is_boundary(he01)) {
		const Point v21 = p2 - p1;
		const Point v20 = p2 - p0;
		cot += OpenMesh::dot(v21, v20) / OpenMesh::cross(v21, v20).norm();
	}
	if (!is_boundary(he10)) {
		const Point v31 = p3 - p1;
		const Point v30 = p3 - p0;
		cot += OpenMesh::dot(v31, v30) / OpenMesh::cross(v31, v30).norm();
	}

	return cot;
}
