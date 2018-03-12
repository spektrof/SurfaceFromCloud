#include "glrender.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLTexture>
#include <QMouseEvent>

#ifdef ENABLE_DEBUG
	#include <QDebug>
#endif

GLRender::GLRender(QWidget *parent)
	: QOpenGLWidget(parent),
	clearColor(Qt::black),
	xRot(0), yRot(0), zRot(0),
	counter(0),
	zoom(1.0f),
	drawing_size(0),
	sfg(SurfaceGenerator()),
	program(0)
{
	memset(textures, 0, sizeof(textures));
	elementbuffer = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);

}

GLRender::~GLRender()
{
	makeCurrent();
	vbo.destroy();
	cloud_vbo.destroy();

	for (int i = 0; i < 6; ++i)
		delete textures[i];
	delete program;
	doneCurrent();
}

void GLRender::qNormalizeAngle(int &angle)
{
	while (angle < 0)
		angle += 360 * 16;
	while (angle > 360 * 16)
		angle -= 360 * 16;
}

void GLRender::rotateBy(int xAngle, int yAngle, int zAngle)
{
	xRot += xAngle;
	yRot += yAngle;
	zRot += zAngle;

	qNormalizeAngle(xRot);
	qNormalizeAngle(yRot);
	qNormalizeAngle(zRot);

	update();
}

void GLRender::setClearColor(const QColor &color)
{
	clearColor = color;
	update();
}

void GLRender::initializeGL()
{
	initializeOpenGLFunctions();

	//makeObject();		//TODO: REMOVE THESE later
	makeObjectFromCloud();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

#define PROGRAM_VERTEX_ATTRIBUTE 0
#define PROGRAM_TEXCOORD_ATTRIBUTE 1

	QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
	const char *vsrc =
		"attribute highp vec4 vertex;\n"
		"attribute mediump vec4 texCoord;\n"
		"varying mediump vec4 texc;\n"
		"uniform mediump mat4 matrix;\n"
		"void main(void)\n"
		"{\n"
		"    gl_Position = matrix * vertex;\n"
		"    texc = texCoord;\n"
		"}\n";
	vshader->compileSourceCode(vsrc);

	QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
	const char *fsrc =
		"uniform sampler2D texture;\n"
		"varying mediump vec4 texc;\n"
		"void main(void)\n"
		"{\n"
		"    gl_FragColor = vec4(1,0,0,0.5) /*texture2D(texture, texc.st)*/;\n"
		"}\n";
	fshader->compileSourceCode(fsrc);

	program = new QOpenGLShaderProgram;
	program->addShader(vshader);
	program->addShader(fshader);
	program->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
	program->bindAttributeLocation("texCoord", PROGRAM_TEXCOORD_ATTRIBUTE);
	program->link();

	program->bind();
	program->setUniformValue("texture", 0);
}

void GLRender::attemptUpdate()
{
	switch(sfg.process())
	{
		case PROCESS_DONE:
		{
			makeObjectFromCloud();
			emit step_finished();					//calculation bar refresh
			
			update();

			break;
		}
		case NO_INPUT_DATA:
		{
			emit data_ended();		//TODO: can be other error
			break;
		}
		default:
			break;
	}
}

void GLRender::paintGL()
{
	switch (sfg.getDrawObject())
	{
		case POINTSS:
			cloud_vbo.bind();
			break;
		case BOX:
		{
			elementbuffer.bind();
			break;
		}
		case SURFACE:
		{
			break;
		}
		default:
			break;
	}

	counter++;		//TODO: this will be the fps claculator
	glClearColor(clearColor.redF(), clearColor.greenF(), clearColor.blueF(), clearColor.alphaF());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	QMatrix4x4 m;
	m.ortho(-2.5f / zoom, +2.5f / zoom, +2.5f / zoom, -2.5f / zoom, -10.0f, 85.0f);

	m.translate(0.0f, 0.0f, -50.0f );
	m.rotate(xRot / 16.0f, 1.0f, 0.0f, 0.0f);
	m.rotate(yRot / 16.0f, 0.0f, 1.0f, 0.0f);
	m.rotate(zRot / 16.0f, 0.0f, 0.0f, 1.0f);

	program->setUniformValue("matrix", m);
	program->enableAttributeArray(PROGRAM_VERTEX_ATTRIBUTE);
	program->enableAttributeArray(PROGRAM_TEXCOORD_ATTRIBUTE);
	program->setAttributeBuffer(PROGRAM_VERTEX_ATTRIBUTE, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
	program->setAttributeBuffer(PROGRAM_TEXCOORD_ATTRIBUTE, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));

	renderObject();
}

void GLRender::resizeGL(int width, int height)
{
	int side = qMin(width, height);
	glViewport((width - side) / 2, (height - side) / 2, side, side);
}

void GLRender::mousePressEvent(QMouseEvent *event)
{
	lastPos = event->pos();
}

void GLRender::mouseMoveEvent(QMouseEvent *event)
{
	int dx = event->x() - lastPos.x();
	int dy = event->y() - lastPos.y();

	if (event->buttons() & Qt::LeftButton) {
		rotateBy(/*8 * dy*/0, 8 * dx, 0);
	}
	else if (event->buttons() & Qt::RightButton) {
		rotateBy(8 * dy, 0, /*8 * dx*/0);
	}
	lastPos = event->pos();
}

void GLRender::mouseReleaseEvent(QMouseEvent * /* event */)
{
	emit clicked();
}

void GLRender::wheelEvent(QWheelEvent *event)
{
	float magnitude = event->angleDelta().y() / 240.0f;

	zoom = zoom + magnitude > 0 ? zoom + magnitude : zoom;
	update();
}

//we wont need this - only test basics
void GLRender::makeObject()
{
	static const int coords[6][4][3] = {
		{ { +1, -1, -1 },{ -1, -1, -1 },{ -1, +1, -1 },{ +1, +1, -1 } },
	{ { +1, +1, -1 },{ -1, +1, -1 },{ -1, +1, +1 },{ +1, +1, +1 } },
	{ { +1, -1, +1 },{ +1, -1, -1 },{ +1, +1, -1 },{ +1, +1, +1 } },
	{ { -1, -1, -1 },{ -1, -1, +1 },{ -1, +1, +1 },{ -1, +1, -1 } },
	{ { +1, -1, +1 },{ -1, -1, +1 },{ -1, -1, -1 },{ +1, -1, -1 } },
	{ { -1, -1, +1 },{ +1, -1, +1 },{ +1, +1, +1 },{ -1, +1, +1 } }
	};

	for (int j = 0; j < 6; ++j)
		textures[j] = new QOpenGLTexture(QImage(QString("images/side%1.png").arg(j + 1)).mirrored());

	QVector<GLfloat> vertData;
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 4; ++j) {
			// vertex position
			vertData.append(0.2 * coords[i][j][0]);
			vertData.append(0.2 * coords[i][j][1]);
			vertData.append(0.2 * coords[i][j][2]);
			// texture coordinate
			vertData.append(j == 0 || j == 3);
			vertData.append(j == 0 || j == 1);
		}
	}

	vbo.create();
	vbo.bind();
	vbo.allocate(vertData.constData(), vertData.count() * sizeof(GLfloat));
}

void GLRender::makeObjectFromCloud()
{
	GLPaintFormat paintData = sfg.getPaintData();
	drawing_size = 0;

	if (paintData.p.empty() && paintData.box_points.empty()) { return; }

	switch (sfg.getDrawObject())
	{
		case POINTSS:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();

			QVector<GLfloat> vertData;
			for (auto it : paintData.p)
			{
				// vertex position
				vertData.append(it.x());
				vertData.append(it.y());
				vertData.append(it.z());
				// texture coordinate
				vertData.append(0);
				vertData.append(1);
				//TODO
			}

			drawing_size = paintData.p.size();

			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(vertData.constData(), vertData.count() * sizeof(GLfloat));
			return;
		}
		case BOX:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.box_points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
				// texture coordinate
				boxData.append(0);
				boxData.append(1);
				//TODO
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			//-------------
			drawing_size = paintData.ix.size();

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			return;
		}
		case SURFACE:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.box_points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
				// texture coordinate
				boxData.append(0);
				boxData.append(1);
				//TODO
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			//-------------
			drawing_size = paintData.ix.size();

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));

#ifdef ENABLE_DEBUG
			qDebug() << *std::max_element(paintData.ix.begin(), paintData.ix.end()) << "\n" << boxData.size();
#endif
			return;
		}
		default:
			return;
	}

	
}

void GLRender::renderObject()
{
	//TODO FPS COUNTER WITH INDEXING
	//http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-9-vbo-indexing/

	switch (sfg.getDrawObject())
	{
	case POINTSS:
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
		cloud_vbo.bind();
		glDrawArrays(GL_POINTS, 0, drawing_size);
		return;
	case BOX:
	{
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
		cloud_vbo.bind();
		elementbuffer.bind();

		glDrawElements(
			GL_TRIANGLES,      // mode
			drawing_size,    // count
			GL_UNSIGNED_INT,   // type
			(void*)0           // element array buffer offset
		);
		return;
	}
	case SURFACE:		//crash!!
	{
		if (drawing_size == 0) return;

#ifdef ENABLE_DEBUG
		qDebug() << "Surface draw\n";
#endif

		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
		//cloud_vbo.bind();
		elementbuffer.bind();

#ifdef ENABLE_DEBUG
		qDebug() << "\tBinded\n";
#endif

		glDrawElements(
			GL_TRIANGLES,      // mode
			drawing_size,    // count
			GL_UNSIGNED_INT,   // type
			(void*)0           // element array buffer offset
		);

#ifdef ENABLE_DEBUG
		qDebug() << "\tDrawed\n";
#endif

		return;
	}
	default:
		return;
	}
}

void GLRender::getCloudCalculationTime(float& c_load_time, float& s_cal_time, float& all_time, float& join_time, float& sort_time)
{
	sfg.getTimes(c_load_time, s_cal_time, all_time, join_time, sort_time);
}

void GLRender::getCloudPoints(int& data)
{
	data = drawing_size;
}

void GLRender::attemptChangeDrawType(const unsigned int& ind)
{ 
	if (sfg.getDrawObject() - ind)
	{
		sfg.setDrawObject(ind);
		makeObjectFromCloud(); 
		update(); 
	} 
}

void GLRender::attemptChangeThreadCapacity(const unsigned int& cloud_t, const unsigned int& surface_t)
{
	if (sfg.getThreadCapacity() - cloud_t)
	{
		//TODO: wait until it finish a session
		sfg.changeThreadCapacity(cloud_t);
	}

	/*if (surf->getThreadCapacity() - surface_t)
	{
		//TODO: wait until it finish a session
		cloud->changeThreadCapacity(surface_t);
	}*/
}

void GLRender::goNode(const int& direction)
{
	switch (direction)
	{
	case 0:
	{
		sfg.goParentNode();
		break;
	}
	case 1:
	{
		sfg.goLeftNode();
		break;
	}
	case 2:
	{
		sfg.goRightNode();
		break;
	}
	default:
		break;
	}

	makeObjectFromCloud();
	update();
}

void GLRender::swapStartFlag()
{
	sfg.swapStartFlag();
}

void GLRender::showNearest()
{
	/*cloud->nearest();
	makeObjectFromCloud();
	update();*/
}

void GLRender::showKNearest()
{
	/*cloud->knearest();
	makeObjectFromCloud();
	update();*/
}