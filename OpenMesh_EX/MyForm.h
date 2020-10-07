#pragma once

#include <vgl.h>
#include "DotNetUtilities.h"
#include "Mesh/GUA_OM.h"
#include "Mesh/DP.h"
#include <LoadShaders.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "glm/ext.hpp"
using namespace glm;
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
using namespace std;

#define FACE_SIZE 5000
int OBJ_NUM = 13;
int objptr = 0;
Tri_Mesh *mesh;
Tri_Mesh *patch; // selected patch

GLint SCR_WIDTH = 817;
GLint SCR_HEIGHT = 541;

GLuint quadVAO;
GLuint quadVBO;
unsigned int framebuffer;
unsigned int textureColorbuffer;
unsigned int programFrame;
unsigned int programUV;
unsigned int programImg;
glm::vec4 pixel;
int facesid[FACE_SIZE];
vector<int> facesid2;
std::vector<double> modelCenter;

float quadVertices1[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
	// positions   // texCoords
	-1.0f,  1.0f,  0.0f, 1.0f,
	-1.0f, -1.0f,  0.0f, 0.0f,
	1.0f, -1.0f,  1.0f, 0.0f,

	-1.0f,  1.0f,  0.0f, 1.0f,
	1.0f, -1.0f,  1.0f, 0.0f,
	1.0f,  1.0f,  1.0f, 1.0f
};

bool isLoad = false;
bool wholeModel = false;
bool line = true;
std::vector<double> vertices[13];

std::vector<double> meshUV[13];
std::vector<double> verticesPatch; // patch的點，給vbo用
std::vector<double> patchUV; // patch的uv座標，給vbo用
std::vector<double> selectedVertices;

unsigned int checkerBoardImg[3]; // image on model
//unsigned int imgptr[13] = { 1, 0, 2, 1, 0, 2, 1, 0, 2, 0, 1, 2, 0 };
unsigned int imgptr[4][13] = { 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 2, 0, 1,
							   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
							   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };

unsigned int img = 0; // one patch
unsigned int design = 0;
unsigned int checkerBoardImg2[3]; // show the picture
double rotateAngle = 0.0f;

float eyeAngleX = 0.0;
float eyeAngleY = 0.0;
float translateX = 0.0;
float translateY = 0.0;
float eyedistanceuv = 2.0;

float eyedistance = 2.0;
#define DOR(angle) (angle*3.1415/180);
int prevMouseX, prevMouseY;

mat4 MVPuv;
mat4 ProjectionUV;
mat4 ViewMatrixUV;

GLuint VBO;
GLuint meshVBO;
GLuint VAO;
GLuint VBOuv;
GLuint VAOuv;
GLuint VBOi;
GLuint VAOi;
int face[13];
int facePatch;

GLuint program;

GLuint mvpID;
GLuint ColorID;

mat4 Projection;
mat4 ViewMatrix;
mat4 Model;
mat4 MVP;

ShaderInfo meshShaders[] = {
	{ GL_VERTEX_SHADER, "meshShader.vp" },//vertex shader
	{ GL_FRAGMENT_SHADER, "meshShader.fp" },//fragment shader
	{ GL_NONE, NULL } };

typedef struct _TextureData {
	_TextureData() : width(0), height(0), data(0) {}
	int width;
	int height;
	unsigned char* data;
} TextureData;

TextureData Load_png(const char* path, bool mirroredY = true) {
	TextureData texture;
	int n;
	stbi_uc *data = stbi_load(path, &texture.width, &texture.height, &n, 4);
	if (data != NULL) {
		texture.data = new unsigned char[texture.width * texture.height * 4 * sizeof(unsigned char)];
		memcpy(texture.data, data, texture.width * texture.height * 4 * sizeof(unsigned char));
		// vertical-mirror image data
		if (mirroredY) {
			for (size_t i = 0; i < texture.width; i++) {
				for (size_t j = 0; j < texture.height / 2; j++) {
					for (size_t k = 0; k < 4; k++) {
						std::swap(texture.data[(j * texture.width + i) * 4 + k], texture.data[((texture.height - j - 1) * texture.width + i) * 4 + k]);
					}
				}
			}
		}
		stbi_image_free(data);
		printf("texture load complete at path : %s\n", path);
	}
	else {
		printf("texture load un-complete at path : %s\n", path);
	}
	return texture;
}

//not used
// GLCamera  read_depth
bool read_mouse(int x, int y, point &p, mat4 mv, mat4 proj) {

	GLdouble M[16], P[16]; GLint V[4];

	//glGetDoublev(GL_MODELVIEW_MATRIX, M);
	//glGetDoublev(GL_PROJECTION_MATRIX, P);
	glGetIntegerv(GL_VIEWPORT, V);

	/*
	ViewMatrix = lookAt(
		glm::vec3(eyedistance*cos(horizonAngle)*cos(verticleAngle), eyedistance*sin(verticleAngle), eyedistance*sin(horizonAngle)*cos(verticleAngle)),
		glm::vec3(0, 0, 0), // and looks at the origin
		glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
	);
	*/


	mat4 Model = translate(translateX, translateY, 0.0f);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			M[i] = Model[i][j];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			P[i] = proj[i][j];
		}
	}

	static const float dx[] =
	{ 0, 1,-1,-1, 1, 3,-3, 0, 0, 6,-6,-6, 6, 25,-25,  0,  0 };
	static const float dy[] =
	{ 0, 1, 1,-1,-1, 0, 0, 3,-3, 6, 6,-6,-6,  0,  0, 25,-25 };
	const float scale = 0.01f;
	const int displacements = sizeof(dx) / sizeof(float);

	int xmin = V[0], xmax = V[0] + V[2] - 1, ymin = V[1], ymax = V[1] + V[3] - 1;

	for (int i = 0; i < displacements; i++) {
		int xx = std::min(std::max(x + int(dx[i] * scale*V[2]), xmin), xmax);
		int yy = std::min(std::max(y + int(dy[i] * scale*V[3]), ymin), ymax);

		float d;
		glReadPixels(xx, yy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &d);

		static float maxd = 0.0f;
		if (!maxd) {
			glScissor(xx, yy, 1, 1);
			glEnable(GL_SCISSOR_TEST);
			glClearDepth(1);
			glClear(GL_DEPTH_BUFFER_BIT);
			glReadPixels(xx, yy, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &maxd);
			if (maxd) {
				glClearDepth(d / maxd);
				glClear(GL_DEPTH_BUFFER_BIT);
			}
			glDisable(GL_SCISSOR_TEST);
			glClearDepth(1);
			if (!maxd)
				return false;
		}

		d /= maxd;
		if (d > 0.0001f && d < 0.9999f) {
			GLdouble X, Y, Z;
			gluUnProject(xx, yy, d, M, P, V, &X, &Y, &Z);
			p = point((float)X, (float)Y, (float)Z);
			return true;
		}
	}
	return false;
}

namespace OpenMesh_EX {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// MyForm 的摘要
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form {
	public:
		MyForm(void) {
			//constructer
			InitializeComponent();
			std::cout << "construct" << std::endl;
			pixel.r = 0.0f;
			pixel.g = 0.0f;
			pixel.b = 0.0f;
			pixel.a = 0.0f;
			for (int i = 0; i < FACE_SIZE; i++) facesid[i] = -1;
			Projection = glm::perspective(100.0f, 4.0f / 3.0f, 0.001f, 1000.0f);
		}

	protected:
		/// <summary>
		/// 清除任何使用中的資源。
		/// </summary>
		~MyForm() {
			if (components) {
				delete components;
			}
		}

	private: System::Windows::Forms::MenuStrip^  menuStrip1;
	private: System::Windows::Forms::ToolStripMenuItem^  fileToolStripMenuItem;
	private: System::Windows::Forms::ToolStripMenuItem^  loadModelToolStripMenuItem;
	private: System::Windows::Forms::OpenFileDialog^  openModelDialog;
	private: System::Windows::Forms::SaveFileDialog^  saveModelDialog;
	private: System::Windows::Forms::ToolStripMenuItem^  saveModelToolStripMenuItem;
	private: System::Windows::Forms::TableLayoutPanel^  tableLayoutPanel1;
	private: HKOGLPanel::HKOGLPanelControl^  hkoglPanelControl1;
	private: System::Windows::Forms::ToolStripMenuItem^  errorQuadricToolStripMenuItem;
	protected:

	private:
		/// <summary>
		/// 設計工具所需的變數。
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// 此為設計工具支援所需的方法 - 請勿使用程式碼編輯器修改
		/// 這個方法的內容。
		/// </summary>
		void InitializeComponent(void) {
			HKOGLPanel::HKCOGLPanelCameraSetting^  hkcoglPanelCameraSetting1 = (gcnew HKOGLPanel::HKCOGLPanelCameraSetting());
			HKOGLPanel::HKCOGLPanelPixelFormat^  hkcoglPanelPixelFormat1 = (gcnew HKOGLPanel::HKCOGLPanelPixelFormat());
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->fileToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->loadModelToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->saveModelToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->errorQuadricToolStripMenuItem = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->openModelDialog = (gcnew System::Windows::Forms::OpenFileDialog());
			this->saveModelDialog = (gcnew System::Windows::Forms::SaveFileDialog());
			this->hkoglPanelControl1 = (gcnew HKOGLPanel::HKOGLPanelControl());
			this->tableLayoutPanel1 = (gcnew System::Windows::Forms::TableLayoutPanel());
			this->menuStrip1->SuspendLayout();
			this->tableLayoutPanel1->SuspendLayout();
			this->SuspendLayout();
			// 
			// menuStrip1
			// 
			this->menuStrip1->ImageScalingSize = System::Drawing::Size(20, 20);
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->fileToolStripMenuItem,
					this->errorQuadricToolStripMenuItem
			});
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Size = System::Drawing::Size(1051, 24);
			this->menuStrip1->TabIndex = 1;
			this->menuStrip1->Text = L"menuStrip1";
			this->menuStrip1->ItemClicked += gcnew System::Windows::Forms::ToolStripItemClickedEventHandler(this, &MyForm::menuStrip1_ItemClicked);
			// 
			// fileToolStripMenuItem
			// 
			this->fileToolStripMenuItem->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(2) {
				this->loadModelToolStripMenuItem,
					this->saveModelToolStripMenuItem
			});
			this->fileToolStripMenuItem->Name = L"fileToolStripMenuItem";
			this->fileToolStripMenuItem->Size = System::Drawing::Size(38, 20);
			this->fileToolStripMenuItem->Text = L"File";
			// 
			// loadModelToolStripMenuItem
			// 
			this->loadModelToolStripMenuItem->Name = L"loadModelToolStripMenuItem";
			this->loadModelToolStripMenuItem->Size = System::Drawing::Size(144, 22);
			this->loadModelToolStripMenuItem->Text = L"Load Model";
			this->loadModelToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::loadModelToolStripMenuItem_Click);
			// 
			// saveModelToolStripMenuItem
			// 
			this->saveModelToolStripMenuItem->Name = L"saveModelToolStripMenuItem";
			this->saveModelToolStripMenuItem->Size = System::Drawing::Size(144, 22);
			this->saveModelToolStripMenuItem->Text = L"Save Model";
			this->saveModelToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::saveModelToolStripMenuItem_Click);
			// 
			// errorQuadricToolStripMenuItem
			// 
			this->errorQuadricToolStripMenuItem->Name = L"errorQuadricToolStripMenuItem";
			this->errorQuadricToolStripMenuItem->Size = System::Drawing::Size(94, 20);
			this->errorQuadricToolStripMenuItem->Text = L"Error Quadric";
			this->errorQuadricToolStripMenuItem->Click += gcnew System::EventHandler(this, &MyForm::errorQuadricToolStripMenuItem_Click);
			// 
			// openModelDialog
			// 
			this->openModelDialog->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::openModelDialog_FileOk);
			// 
			// saveModelDialog
			// 
			this->saveModelDialog->DefaultExt = L"obj";
			this->saveModelDialog->FileOk += gcnew System::ComponentModel::CancelEventHandler(this, &MyForm::saveModelDialog_FileOk);
			// 
			// hkoglPanelControl1
			// 
			hkcoglPanelCameraSetting1->Far = 1000;
			hkcoglPanelCameraSetting1->Fov = 45;
			hkcoglPanelCameraSetting1->Near = -1000;
			hkcoglPanelCameraSetting1->Type = HKOGLPanel::HKCOGLPanelCameraSetting::CAMERATYPE::ORTHOGRAPHIC;
			this->hkoglPanelControl1->Camera_Setting = hkcoglPanelCameraSetting1;
			this->hkoglPanelControl1->Dock = System::Windows::Forms::DockStyle::Fill;
			this->hkoglPanelControl1->Location = System::Drawing::Point(3, 3);
			this->hkoglPanelControl1->Name = L"hkoglPanelControl1";
			hkcoglPanelPixelFormat1->Accumu_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Alpha_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			hkcoglPanelPixelFormat1->Stencil_Buffer_Bits = HKOGLPanel::HKCOGLPanelPixelFormat::PIXELBITS::BITS_0;
			this->hkoglPanelControl1->Pixel_Format = hkcoglPanelPixelFormat1;
			this->hkoglPanelControl1->Size = System::Drawing::Size(1044, 494);
			this->hkoglPanelControl1->TabIndex = 2;
			this->hkoglPanelControl1->Load += gcnew System::EventHandler(this, &MyForm::hkoglPanelControl1_Load);
			this->hkoglPanelControl1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MyForm::hkoglPanelControl1_Paint);
			this->hkoglPanelControl1->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseDown);
			this->hkoglPanelControl1->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseMove);
			this->hkoglPanelControl1->MouseWheel += gcnew System::Windows::Forms::MouseEventHandler(this, &MyForm::hkoglPanelControl1_MouseWheel);
			// 
			// tableLayoutPanel1
			// 
			this->tableLayoutPanel1->ColumnCount = 1;
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent,
				61.05442F)));
			this->tableLayoutPanel1->ColumnStyles->Add((gcnew System::Windows::Forms::ColumnStyle(System::Windows::Forms::SizeType::Percent,
				38.94558F)));
			this->tableLayoutPanel1->Controls->Add(this->hkoglPanelControl1, 0, 0);
			this->tableLayoutPanel1->Location = System::Drawing::Point(0, 22);
			this->tableLayoutPanel1->Margin = System::Windows::Forms::Padding(2);
			this->tableLayoutPanel1->Name = L"tableLayoutPanel1";
			this->tableLayoutPanel1->RowCount = 1;
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 50)));
			this->tableLayoutPanel1->RowStyles->Add((gcnew System::Windows::Forms::RowStyle(System::Windows::Forms::SizeType::Percent, 50)));
			this->tableLayoutPanel1->Size = System::Drawing::Size(1050, 500);
			this->tableLayoutPanel1->TabIndex = 3;
			this->tableLayoutPanel1->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MyForm::tableLayoutPanel1_Paint);
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1051, 521);
			this->Controls->Add(this->tableLayoutPanel1);
			this->Controls->Add(this->menuStrip1);
			this->MainMenuStrip = this->menuStrip1;
			this->Name = L"MyForm";
			this->Text = L"OpenMesh_EX";
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			this->tableLayoutPanel1->ResumeLayout(false);
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion

		//init
	private: System::Void hkoglPanelControl1_Load(System::Object^  sender, System::EventArgs^  e) {

		SCR_WIDTH = this->hkoglPanelControl1->Size.Width;
		SCR_HEIGHT = this->hkoglPanelControl1->Size.Height;

		//////openGL init//////////
		glewExperimental = GL_TRUE; //置於glewInit()之前
		if (glewInit()) {
			std::cerr << "Unable to initialize GLEW ... exiting" << std::endl;//c error
			exit(EXIT_FAILURE);
		}
		else std::cout << "initialize GLEW success" << std::endl;//c error
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

		//VAO
		glGenVertexArrays(1, &VAO);
		glBindVertexArray(VAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &meshVBO);

		program = LoadShaders(meshShaders);//讀取shader
		glUseProgram(program);//uniform參數數值前必須先use shader
		mvpID = glGetUniformLocation(program, "MVP");
		ColorID = glGetUniformLocation(program, "Color");

		ViewMatrix = glm::lookAt(
			glm::vec3(0, 5, 5), // Camera is at (0,10,25), in World Space
			glm::vec3(0, 0, 0), // and looks at the origin
			glm::vec3(0, 1, 0)  // Head is up (set to 0,1,0 to look upside-down)
		);


		for (int i = 0; i < 3; i++) {
			std::string ProjectName;
			if (i == 0) ProjectName = "circle2.jpg";
			else if (i == 1) ProjectName = "circle3.jpg";
			else if (i == 2) ProjectName = "circle4.jpg";

			TextureData tdata = Load_png((ProjectName).c_str(), true);

			//Generate empty texture
			glGenTextures(1, &checkerBoardImg[i]);

			glBindTexture(GL_TEXTURE_2D, checkerBoardImg[i]);

			//Do texture setting
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, tdata.width, tdata.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, tdata.data);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glBindTexture(GL_TEXTURE_2D, 0);
			//-----------------------------------------------------------
			glUniform1i(glGetUniformLocation(program, "sprite"), 0);
		}

		//clear canvas
		glClearColor(0.0, 0.0, 0.0, 1);

		//use frameBuffer to store face id
		ShaderInfo shaderframe[] = {
		{ GL_VERTEX_SHADER, "framebuffer.vp" },//vertex shader
		{ GL_FRAGMENT_SHADER, "framebuffer.fp" },//fragment shader
		{ GL_NONE, NULL } };
		programFrame = LoadShaders(shaderframe);
		glUseProgram(programFrame);//uniform參數數值前必須先use shader

		// screen quad VAO
		glGenVertexArrays(1, &quadVAO);
		glGenBuffers(1, &quadVBO);

		glUseProgram(programFrame);
		glUniform1i(glGetUniformLocation(programFrame, "screenTexture"), 0);

		// framebuffer configuration
		// -------------------------
		glGenFramebuffers(1, &framebuffer);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
		// create a color attachment texture

		glGenTextures(1, &textureColorbuffer);
		glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
		//GL_RGBA32F for store value > 1
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, hkoglPanelControl1->Width, hkoglPanelControl1->Height, 0, GL_RGBA, GL_FLOAT, NULL);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);

		GLenum dr[2] = { GL_COLOR_ATTACHMENT0 ,GL_DEPTH_ATTACHMENT };
		glDrawBuffers(2, dr);

		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
			cout << "Framebuffer is not complete!" << endl;
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT)
			cout << "Framebuffer is not complete attach!" << endl;

		//bind to normal
		glBindTexture(GL_TEXTURE_2D, 0);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}

	//display
	private: System::Void hkoglPanelControl1_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) {
		//std::cout << "refresh" << std::endl;
		SCR_WIDTH = this->hkoglPanelControl1->Size.Width;
		SCR_HEIGHT = this->hkoglPanelControl1->Size.Height;

		///////////////draw origin obj///////////

		for (int i = 0; i < objptr; i++) {

			//assign mesh data & UV data to VBO
			if (mesh != NULL) {
				//std::cout << "refresh mesh not null" << std::endl;
				glBindBuffer(GL_ARRAY_BUFFER, VBO);
				glBufferData(GL_ARRAY_BUFFER, (vertices[i]).size() * sizeof(double), &(vertices[i][0]), GL_STATIC_DRAW);
				glBindBuffer(GL_ARRAY_BUFFER, meshVBO);
				glBufferData(GL_ARRAY_BUFFER, (meshUV[i]).size() * sizeof(double), &(meshUV[i][0]), GL_STATIC_DRAW);
			}
			glEnable(GL_DEPTH_TEST);

			glBindVertexArray(VAO);
			glUseProgram(program);//uniform參數數值前必須先use shader

			//set MVP matrix
			float horizonAngle = DOR(eyeAngleX);
			float verticleAngle = DOR(eyeAngleY);
			//ViewMatrix = lookAt(
			//	glm::vec3(eyedistance*cos(horizonAngle)*cos(verticleAngle), eyedistance*sin(verticleAngle), eyedistance*sin(horizonAngle)*cos(verticleAngle)),
			//	glm::vec3(0, 0, 0), // and looks at the origin
			//	glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
			//);
			ViewMatrix = lookAt(
				glm::vec3(modelCenter[0], modelCenter[1], modelCenter[2] + eyedistance ),
				glm::vec3(modelCenter[0], modelCenter[1], modelCenter[2]), // and looks at the origin
				glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
			);
			mat4 Model = translate(translateX, translateY, 0.0f) * rotate(horizonAngle, float(modelCenter[0]), 1.0f, float(modelCenter[2]));
			//MVP = Model * Projection * ViewMatrix;//translate via screen viewport, so model last

			MVP = Projection * ViewMatrix * Model;

			glUniformMatrix4fv(mvpID, 1, GL_FALSE, &MVP[0][0]);

			//VBO setting
			glBindBuffer(GL_ARRAY_BUFFER, VBO);
			// 1rst attribute buffer : vertices
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0,				//location
				3,				//vec3
				GL_DOUBLE,			//type
				GL_FALSE,			//not normalized

				0,		//strip
				0);//buffer offset
			glBindBuffer(GL_ARRAY_BUFFER, meshVBO); // this attribute comes from a different vertex buffer
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 2, GL_DOUBLE, GL_FALSE, 0, (void*)0); // 從instanceVBO傳入的
			glBindBuffer(GL_ARRAY_BUFFER, 0);


			////////////////////ready to draw////////////////

			//clear canvas
			glClearDepth(1.0);
			glClearColor(0.0, 0.0, 0.0, 1.0);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			//draw faceID to frameBuffer
			if (isLoad) {
				glBindFramebuffer(GL_DRAW_FRAMEBUFFER, framebuffer);
				glEnable(GL_DEPTH_TEST); // enable depth testing (is disabled for rendering screen-space quad)
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				//set rgb to -1 to check is storing face id
				glm::vec3 color = glm::vec3(-1.0, 0.0, 0.0);
				glUniform3fv(ColorID, 1, &color[0]);
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, checkerBoardImg[imgptr[design][i]]);
				glDrawArrays(GL_TRIANGLES, 0, face[i] * 3);
			}
			glClearColor(1.0, 1.0, 1.0, 1.0);
			glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

			//draw to screen
			if (isLoad) {
				//line
				if (line) {
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					glm::vec3 color = glm::vec3(0.0, 0.0, 0.0);
					glUniform3fv(ColorID, 1, &color[0]);
					glDrawArrays(GL_TRIANGLES, 0, face[i] * 3);
				}
				//face
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				glm::vec3 color = glm::vec3(1.0, 0.85, 0.5);
				glUniform3fv(ColorID, 1, &color[0]);
				glDrawArrays(GL_TRIANGLES, 0, face[i] * 3);

			}
		}

		//----------------------------
		//draw current picked vert
		//---------------------------

		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		//glBindVertexArray(VAO);
		glUseProgram(program);//uniform參數數值前必須先use shader
		glUniformMatrix4fv(mvpID, 1, GL_FALSE, &MVP[0][0]);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);

		if (selectedVertices.size() != 0) {
			glBindBuffer(GL_ARRAY_BUFFER, VBO);
			glBufferData(GL_ARRAY_BUFFER, selectedVertices.size() * sizeof(double), &selectedVertices[0], GL_STATIC_DRAW);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0,				//location
				3,				//vec3
				GL_DOUBLE,			//type
				GL_FALSE,			//not normalized
				0,				//strip
				0);//buffer offset
			glDisable(GL_DEPTH_TEST);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glm::vec3 color = glm::vec3(1.0, 0.0, 1.0);
			glUniform3fv(ColorID, 1, &color[0]);
			glPointSize(10.0);
			glDrawArrays(GL_POINTS, 0, 1);
		}
	}

	//mouseClick
	private: System::Void hkoglPanelControl1_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
		if (e->Button == System::Windows::Forms::MouseButtons::Left) {
			//leftClick	
			//record mouse position for drag event
			prevMouseX = e->X;
			prevMouseY = e->Y;
			//cout << "click : " << e->X << "," << e->Y << endl;
			std::vector<double> mousePosition;
			double mouseOnScreenX = (double)e->X * 2 / (double)SCR_WIDTH - 1;
			double mouseOnScreenY = (double)(SCR_HEIGHT - e->Y) * 2 / (double)SCR_HEIGHT - 1;
			mousePosition.push_back(mouseOnScreenX);
			mousePosition.push_back(mouseOnScreenY);
			//cout << "click SCR : " << mouseOnScreenX << "," << mouseOnScreenY << endl;
			//read face
			glBindFramebuffer(GL_READ_FRAMEBUFFER, framebuffer);
			glReadBuffer(GL_COLOR_ATTACHMENT0);
			glReadPixels(e->X, hkoglPanelControl1->Height - e->Y, 1, 1, GL_RGBA, GL_FLOAT, &pixel);
			cout << "face id : " << pixel.r << endl;
			if (pixel.r > 0) {
				selectedVertices.clear();
				mesh->findNearestVert(*mesh, mousePosition, pixel.r - 1, selectedVertices, MVP, eyedistance);
				//update mesh
				
				mesh->loadToBuffer(*mesh, vertices[0], face[0], meshUV[0], modelCenter);
				std::cout << "meshUV.size() : " << meshUV[objptr].size() << "vertices.size()" << vertices[objptr].size() << endl;
				std::cout << "face" << face[objptr] << std::endl;
			}


			/*
			//printf("mouse x = %d mouse y = %d\n", e->X, hkoglPanelControl1->Height - e->Y);
			if (isLoad) {
				//detect same face already stored
				for (int i = 0; i < facesid2.size(); i++) {
					if (facesid2[i] == int(pixel.r) - 1) break;
					if (pixel.r != 0 && (facesid2[i] > int(pixel.r) - 1 || i == facesid2.size() - 1)) {
						//when id > curID mean no repeat 'cause vector sorted
						facesid2.push_back(int(pixel.r) - 1);
						break;
					}
				}
				//first value
				if (pixel.r != 0 && facesid2.size() == 0) facesid2.push_back(int(pixel.r) - 1);
				std::sort(facesid2.begin(), facesid2.end());

				cout << "selected faceID: ";
				for (int i = 0; i < facesid2.size(); i++) cout << facesid2[i] << " ";
				cout << endl;
				//cout << endl << "selected face count : " << facesid2.size() << endl;
			}
			*/
			glReadBuffer(GL_NONE);
			glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

		}
		if (e->Button == System::Windows::Forms::MouseButtons::Middle) {
			//record mouse position for drag event
			prevMouseX = e->X;
			prevMouseY = e->Y;
		}
		hkoglPanelControl1->Invalidate();
	}

	//mouseDrag
	private: System::Void hkoglPanelControl1_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
		if (e->Button == System::Windows::Forms::MouseButtons::Left) {
			//std::cout << "left" << std::endl;
			//eyeAngleX += (e->X - prevMouseX)*0.05;
			//eyeAngleY += (e->Y - prevMouseY)*0.05;
			eyeAngleX += (e->X - prevMouseX)*10;
			eyeAngleY += (e->Y - prevMouseY)*10;
			//record mouse position for drag event
			prevMouseX = e->X;
			prevMouseY = e->Y;
		}
		if (e->Button == System::Windows::Forms::MouseButtons::Middle) {
			//std::cout << "middle" << std::endl;
			translateX += (e->X - prevMouseX)*0.002;
			translateY -= (e->Y - prevMouseY)*0.002;
			//record mouse position for drag event
			prevMouseX = e->X;
			prevMouseY = e->Y;
		}
		hkoglPanelControl1->Invalidate();
	}

	//mouseWheel
	private: System::Void hkoglPanelControl1_MouseWheel(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) {
		if (e->Delta < 0) eyedistance += 0.2;
		else {
			eyedistance -= 0.2;
			//if (eyedistance < 0.4) eyedistance = 0.4;
			//std::cout << "wheel up, distance : "  << eyedistance << std::endl;
		}
		hkoglPanelControl1->Invalidate();
	}

	//click "openModel"
	private: System::Void loadModelToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) {
		openModelDialog->Filter = "Model(*.obj)|*obj";
		openModelDialog->Multiselect = false;
		openModelDialog->ShowDialog();
	}

	//check openModel dialog
	private: System::Void openModelDialog_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e) {
		std::string filename;
		MarshalString(openModelDialog->FileName, filename);
		//del old mesh on screen
		if (mesh == NULL) {
			mesh = new Tri_Mesh;
			for (int i = 0; i < OBJ_NUM; i++) {
				vertices[i].clear();
				face[i] = 0;
				meshUV[i].clear();
			}
		}
		else if (mesh != NULL) delete mesh;
		mesh = new Tri_Mesh;

		if (ReadFile(filename, mesh)) std::cout << filename << std::endl;
		isLoad = true;

		mesh->loadToBuffer(*mesh, vertices[0], face[0], meshUV[0], modelCenter);
		std::cout << "meshUV.size() : " << meshUV[0].size() << "vertices.size()" << vertices[0].size() << endl;
		std::cout << "face" << face[0] << std::endl;
		objptr++;
		hkoglPanelControl1->Invalidate();
	}

	//saveObj menu open
	private: System::Void saveModelToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) {
		saveModelDialog->Filter = "Model(*.obj)|*obj";
		saveModelDialog->ShowDialog();
	}

	//check saveObj in dialog 
	private: System::Void saveModelDialog_FileOk(System::Object^  sender, System::ComponentModel::CancelEventArgs^  e) {
		std::string filename;
		MarshalString(saveModelDialog->FileName, filename);
		if (SaveFile(filename, mesh)) std::cout << filename << std::endl;
	}

	private: System::Void tableLayoutPanel1_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) {}
	private: System::Void menuStrip1_ItemClicked(System::Object^  sender, System::Windows::Forms::ToolStripItemClickedEventArgs^  e) {}

	// button2 is for 模型上的線on/off
	private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
		if (line) line = false;
		else line = true;
	}

	private: System::Void errorQuadricToolStripMenuItem_Click(System::Object^  sender, System::EventArgs^  e) {
		if (mesh != NULL) {
			Tri_Mesh t = mesh->simplify(0.7);
			//Tri_Mesh t = mesh->averageSimplify();
			mesh = &t;
			hkoglPanelControl1->Invalidate();
		}
	}
};
}
