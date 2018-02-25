#ifndef GLRENDER_H
#define GLRENDER_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>

#include "cloudplugin.h"
#include "types.h"

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

	void getCloudCalculationTime(float& fps, float& join, float& sort, float& all);
	void getCloudPoints(int& data, int& box);
	void attemptUpdate();
	void attemptChangeDrawType(const unsigned int& ind);
	void attemptChangeThreadCapacity(const unsigned int& cloud_t, const unsigned int& surface_t);

	void resetCloud() { cloud->resetCloud(); }	//TODO new name mybe
	void goNode(const int&);

	void changeComputation();

signals:
	void clicked();		//dont need click signal
	void step_finished();

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

	CloudContainer* cloud;

	int counter;	//TODO: change this with fps counter
	float zoom;
};

#endif