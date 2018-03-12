#ifndef GLRENDER_H
#define GLRENDER_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>
#include <QGLBuffer>

#include "surfacegenerator.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram);
QT_FORWARD_DECLARE_CLASS(QOpenGLTexture)

class GLRender : public QOpenGLWidget, protected QOpenGLFunctions
{
	Q_OBJECT

public:
	explicit GLRender(QWidget *parent = 0);
	~GLRender();

	static void qNormalizeAngle(int &angle);
	void rotateBy(int xAngle, int yAngle, int zAngle);

	void setClearColor(const QColor &color);

	void getCloudCalculationTime(float& load_time, float& scal_time, float& all, float& join, float& sort);
	void getCloudPoints(int& data);
	void attemptUpdate();
	void attemptChangeDrawType(const unsigned int& ind);
	void attemptChangeThreadCapacity(const unsigned int& cloud_t, const unsigned int& surface_t);

	void resetCloud() { sfg.resetCloud(); }	//TODO new name mybe
	void goNode(const int&);

	void swapStartFlag();
	void showNearest();
	void showKNearest();

	void calculate_surface()
	{
		ProcessReturns res;
		switch (res = sfg.surface_process())
		{
		case PROCESS_DONE:
		{
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

	bool addNewSurfaceReconstructor(const unsigned int& type)
	{
		return false;
	//	sfc.
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

private:
	void makeObject();
	void makeObjectFromCloud();
	void renderObject();

	QColor clearColor;
	QPoint lastPos;
	int xRot;
	int yRot;
	int zRot;
	QOpenGLTexture *textures[6];
	QOpenGLShaderProgram *program;
	QOpenGLBuffer vbo;		//TODO: delete this later, just an example for drawing
	QOpenGLBuffer cloud_vbo;
	QOpenGLBuffer elementbuffer;

	SurfaceGenerator sfg;

	int counter;	//TODO: change this with fps counter
	float zoom;

	size_t drawing_size;
};

#endif