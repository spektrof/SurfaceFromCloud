#ifndef GLRENDER_H
#define GLRENDER_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>
#include <QGLBuffer>

#include "camera.h"
#include "surfacegenerator.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram);
QT_FORWARD_DECLARE_CLASS(QOpenGLTexture)

class GLRender : public QOpenGLWidget, protected QOpenGLFunctions
{
	Q_OBJECT

#define PROGRAM_VERTEX_ATTRIBUTE 0
#define PROGRAM_TEXCOORD_ATTRIBUTE 1

public:
	explicit GLRender(QWidget *parent = 0);
	~GLRender();

	static void qNormalizeAngle(int &angle);
	void rotateBy(int xAngle, int yAngle, int zAngle);

	void setClearColor(const QColor &color);

	void getCloudCalculationTime(float& load_time, float& scal_time, float& all, float& join, float& sort);
	void getCloudPoints(int& data);
	void attemptUpdate(const int& interval_ms);
	void attemptChangeDrawType(const unsigned int& ind);
	void attemptChangePowerCrustDrawType(const unsigned int& ind);
	void attemptChangeThreadCapacity(const unsigned int& cloud_t, const unsigned int& surface_t);

	void resetCloud() { sfg.resetCloud(); }	//TODO new name mybe
	void goNode(const int&);
	void goCell(const int&);
	void nextNeighb();
	void nextCell() { sfg.nextCInd(); }

	void swapStartFlag();
	//void showNearest();
	//void showKNearest();

	void calculate_surface()
	{
		ProcessReturns res;
		switch (res = sfg.surface_process())
		{
		case PROCESS_DONE:
		{
			makeObjectFromCloud();
			update();
			break;
		}
		case EMPTY_POINTS_SET:
		{
			emit empty_points();
			break;
		}
		case NOT_ENABLED_COMPONENT:
		{
			emit not_enable_component();
			break;
		}
		default:
		{
			emit unknown_error();
		}
		}
	}

	void setMeshAccuracy(const int& val)
	{
		sfg.setMeshAccuracy(val);
	}

	bool newFileInput(char* filename, const int& ma, bool repeat)
	{
		bool ret = sfg.newFileInput(filename, ma, repeat);
		makeObjectFromCloud();
		update();
		return ret;
	}

	bool newRandomInput(const unsigned int& n, const int& ma, bool repeat)
	{
		bool ret = sfg.newRandomInput(n, ma, repeat);
		makeObjectFromCloud();
		update();
		return ret;
	}

	bool addNewOutlierFilter(const unsigned int& type, const float& kn_px, const float& py = 0.0f, const float& pz = 0.0f, const float& nx = 0.0f, const float& ny = 0.0f, const float& nz = 0.0f)
	{
		return sfg.addNewOutlierFilter(type, kn_px, py, pz, nx, ny, nz);
	}

	bool addNewSmootherFilter(const unsigned int& type, const float& first, const float& second = 0.0f)
	{
		return sfg.addNewSmootherFilter(type, first, second);
	}

	bool addNewSimplifierFilter(const unsigned int& type, const float& first, const float& second = 0.0f)
	{
		return sfg.addNewSimplifierFilter(type, first, second);
	}

	bool addNewSurfaceReconstructor(const unsigned int& type, const float& first = 0.0f, const float& second = 0.0f, const float& third = 0.0f)
	{
		return sfg.addNewSurfaceReconstructor(type, first, second, third);
	}

	void SetFilterProperty(const unsigned int& id, const double& first, const double& second = 0.0f, const double& third = 0.0f)
	{
		qDebug() << " SETTING PROP : " << first << " " << second << "\n";
		sfg.SetFilterProperty(id, first, second);
	}

	void SetFilterProperty(const unsigned int& id, const double& x, const double& y, const double& z, const double& nx, const double& ny, const double& nz)
	{
		sfg.SetFilterProperty(id, x, y, z, nx, ny, nz);
	}

	void SetSurfaceProperty(const unsigned int& id, const double& first, const double& second = 0.0f, const double& third = 0.0f)
	{
		sfg.SetSurfaceProperty(id, first, second, third);
	}

	void deleteFilter(const unsigned int& id)
	{
		sfg.deleteFilter(id);
	}

signals:
	void clicked();		//dont need click signal
	void step_finished();
	void data_ended();
	void empty_points();
	void not_enable_component();
	void unknown_error();

protected:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int width, int height) override;
	void mousePressEvent(QMouseEvent *event) override;
	void mouseMoveEvent(QMouseEvent *event) override;
	void mouseReleaseEvent(QMouseEvent *event) override;
	void wheelEvent(QWheelEvent *event) override;
	void keyPressEvent(QKeyEvent *event) override;
	void keyReleaseEvent(QKeyEvent *event) override;

	void initBaseShaders();
	void initSphereShaders();

private:
	void makeObject();
	void makeObjectFromCloud();
	void makeSphereObject();
	void makeSphereRaycastObject();
	void makeObjectOfPowerCrustPart(GLPaintFormat&);
	QVector3D GLRender::GetUV(float u, float v);
	void renderObject(const QMatrix4x4& mvp, const QMatrix4x4& m);
	void renderPowerCrustObject(const QMatrix4x4& mvp, const QMatrix4x4& m);

	void DrawCells(QOpenGLBuffer&, QOpenGLBuffer&, const QMatrix4x4&, QVector3D&, bool transparent = true);
	void DrawSphere(QOpenGLBuffer&, QOpenGLBuffer&, const QMatrix4x4&, std::vector<std::pair<QVector3D, float>>&,const int&, QVector3D&, bool transparent = true);
	void DrawSegments(QOpenGLBuffer& point_buffer, QOpenGLBuffer& index_buffer, const int& starter_index, const int& last_index, const QMatrix4x4& mvp, QVector3D& col, bool transparent = true);

	enum AdditionalBind
	{
		TEXTURE,
		ALL_MATRICIES,
		UNKNOWN
	};

	void shaderBinding(QOpenGLShaderProgram* shader, std::vector<AdditionalBind> texture_bind, QMatrix4x4 mvp);

	QColor clearColor;
	QPoint lastPos;
	int xRot;
	int yRot;
	int zRot;
	QOpenGLTexture *textures[6];
	QOpenGLShaderProgram *base;
	QOpenGLShaderProgram *sphere;
	QOpenGLBuffer vbo;		//TODO: delete this later, just an example for drawing
	QOpenGLBuffer cloud_vbo;
	QOpenGLBuffer elementbuffer;
	QOpenGLBuffer m_quad_vb;

	QOpenGLBuffer sphere_vbo;
	QOpenGLBuffer sphere_indicies;

	std::vector<std::pair<QVector3D, float>> center_radius_pairs;
	std::vector<QVector3D> colors;
	std::vector<QVector3D> points_colors;	//TODO: jobban összevonni
	std::vector<unsigned int> point_draw_parts;
	std::vector<unsigned int> center_draw_parts;

	SurfaceGenerator sfg;

	Camera m_camera;

	int counter;	//TODO: change this with fps counter

};

#endif