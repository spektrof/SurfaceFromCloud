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
	sfg(SurfaceGenerator()),
	p_base(0), pc_base(0), pt_base(0), sphere(0)
{
	memset(textures, 0, sizeof(textures));
	elementbuffer = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
	m_quad_vb = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	cloud_vbo = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	sphere_vbo = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
	sphere_indicies = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);

	setFocusPolicy(Qt::StrongFocus);

	memset(kinect_textures, 0, sizeof(kinect_textures));
	k_ind = 0;
}

GLRender::~GLRender()
{
	makeCurrent();
	sphere_vbo.destroy();
	cloud_vbo.destroy();

	for (int i = 0; i < 6; ++i)
		delete textures[i];

	delete p_base;
	delete pt_base;
	delete pc_base;
	delete sphere;
	doneCurrent();
}

void GLRender::load_texture(std::vector<std::string>& file_list)
{
	const char* colorf[4] = { "images/red.png","images/green.png","images/purple.png","images/blue.png" };
	for (int j = 0; j < file_list.size(); ++j)
	{
	//	qDebug() << file_list[j].c_str();
		kinect_textures[j] = new QOpenGLTexture(QImage(QString(file_list[j].c_str())));
		//kinect_textures[j] = new QOpenGLTexture(QImage(QString(colorf[j])));
	}
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

	makeObjectFromCloud();	//CHECK THIS TODO
	makeSphereRaycastObject();
	makeSphereObject();
	load_texture(sfg.get_file_name_list());

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);

	initBaseShaders();
	initSphereShaders();

	//m_camera.SetProj(-2.5f / zoom, 2.5f / zoom, 2.5f / zoom, -2.5f / zoom, -10.0f, 85.0f);
}

void GLRender::initBaseShaders()
{
	/**************************
	*	Point with uniform col *
	**************************/
	QOpenGLShader *p_vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
	const char *p_vsrc =
		"attribute highp vec4 vertex;\n"
		"uniform mediump mat4 matrix;\n"
		"void main(void)\n"
		"{\n"
		"    gl_Position = matrix * vertex;\n"
		"}\n";
	p_vshader->compileSourceCode(p_vsrc);

	QOpenGLShader *p_fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
	const char *p_fsrc =
		"uniform vec3 color = vec3(1,0,0);"
		"void main(void)\n"
		"{\n"
		"    gl_FragColor = vec4(color,0.5);\n"
		"}\n";
	p_fshader->compileSourceCode(p_fsrc);

	p_base = new QOpenGLShaderProgram;
	p_base->addShader(p_vshader);
	p_base->addShader(p_fshader);
	p_base->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
	p_base->link();

	/**************************
	*	Point, color pairs*
	**************************/
	QOpenGLShader *pc_vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
	const char *pc_vsrc =
		"attribute highp vec4 vertex;\n"
		"attribute mediump vec3 vertex_color;\n"
		"varying mediump vec3 color;\n"
		"uniform mediump mat4 matrix;\n"
		"void main(void)\n"
		"{\n"
		"    gl_Position = matrix * vertex;\n"
		"    color = vertex_color;\n"
		"}\n";
	pc_vshader->compileSourceCode(pc_vsrc);

	QOpenGLShader *pc_fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
	const char *pc_fsrc =
		"varying mediump vec3 color;\n"
		"void main(void)\n"
		"{\n"
		"    gl_FragColor = vec4(color,0.5);\n"
		"}\n";
	pc_fshader->compileSourceCode(pc_fsrc);

	pc_base = new QOpenGLShaderProgram;
	pc_base->addShader(pc_vshader);
	pc_base->addShader(pc_fshader);
	pc_base->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
	pc_base->bindAttributeLocation("vertex_color", PROGRAM_COLOR_ATTRIBUTE);
	pc_base->link();

	/**************************
	*	Point texture coord pair*
	**************************/

	QOpenGLShader *pt_vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
	const char *pt_vsrc =
		"attribute highp vec4 vertex;\n"
		"attribute mediump vec4 texCoord;\n"
		"varying mediump vec4 texc;\n"
		"uniform mediump mat4 matrix;\n"
		"void main(void)\n"
		"{\n"
		"    gl_Position = matrix * vertex;\n"
		"    texc = texCoord;\n"
		"}\n";
	pt_vshader->compileSourceCode(pt_vsrc);

	QOpenGLShader *pt_fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
	const char *pt_fsrc =
		"uniform sampler2D texture;\n"
		"varying mediump vec4 texc;\n"
		"uniform vec3 color = vec3(1,0,0);"
		"void main(void)\n"
		"{\n"
		"    gl_FragColor = texture2D(texture, texc.st);\n"
		"}\n";
	pt_fshader->compileSourceCode(pt_fsrc);

	pt_base = new QOpenGLShaderProgram;
	pt_base->addShader(pt_vshader);
	pt_base->addShader(pt_fshader);
	pt_base->bindAttributeLocation("vertex", PROGRAM_VERTEX_ATTRIBUTE);
	pt_base->bindAttributeLocation("texCoord", PROGRAM_TEXCOORD_ATTRIBUTE);
	pt_base->link();
}

void GLRender::initSphereShaders()
{
	QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, this);
	vshader->compileSourceFile("src/sphere.vert");

	QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, this);
	fshader->compileSourceFile("src/sphere.frag");

	sphere = new QOpenGLShaderProgram;
	sphere->addShader(vshader);
	sphere->addShader(fshader);
	sphere->bindAttributeLocation("vs_in_pos", PROGRAM_VERTEX_ATTRIBUTE);
	sphere->link();
}

void GLRender::attemptUpdate(const int& interval_ms)
{
	switch(sfg.process())
	{
		case PROCESS_DONE:
		{
			makeObjectFromCloud();
			emit step_finished();					//calculation bar refresh
		//	qDebug() << "ATTEMPTUPDATE";
			update();

			break;
		}
		case NO_INPUT_DATA:
		{
			emit data_ended();		//TODO: can be other error
			break;
		}
		default:
			if (m_camera.Update((float)interval_ms / 1000.0f))
				update();
			break;
	}
}

void GLRender::paintGL()
{
	glClearColor(clearColor.redF(), clearColor.greenF(), clearColor.blueF(), clearColor.alphaF());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	QMatrix4x4 m;

	m.setToIdentity();
	m.translate(0.0f, 0.0f, 0.0f);
	m.rotate(xRot / 16.0f, 1.0f, 0.0f, 0.0f);
	m.rotate(yRot / 16.0f, 0.0f, 1.0f, 0.0f);
	m.rotate(zRot / 16.0f, 0.0f, 0.0f, 1.0f);

	QMatrix4x4 mvp = m_camera.GetViewProj() * m;
	renderObject(mvp,m);
}

void GLRender::resizeGL(int width, int height)
{
	//TODO: check halado graf
	int side = qMin(width, height);
	glViewport((width - side) / 2, (height - side) / 2, side, side);

	//glViewport(0, 0, width, height);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//glScalef(height *1. / width, 1.0, 1.0);
	//m_camera.Resize(side, side);
	//glMatrixMode(GL_MODELVIEW);
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
}

QVector3D GLRender::GetUV(float u, float v)
{
	// origó középpontú, egységsugarú gömb parametrikus alakja: http://hu.wikipedia.org/wiki/G%C3%B6mb#Egyenletek 
	// figyeljünk:	matematikában sokszor a Z tengely mutat felfelé, de nálunk az Y, tehát a legtöbb képlethez képest nálunk
	//				az Y és Z koordináták felcserélve szerepelnek
	u *= 2 * 3.1415f;
	v *= 3.1415f;
	float cu = cosf(u), su = sinf(u), cv = cosf(v), sv = sinf(v);
	float r = 1;

	return QVector3D(r*cu*sv, r*cv, r*su*sv);
}

void GLRender::makeSphereObject()
{
	QVector<GLfloat> vertData;
	vertData.clear();

	const int N = 100;
	for (int i = 0; i <= N; ++i)
		for (int j = 0; j <= N; ++j)
		{
			float u = i / (float)N;
			float v = j / (float)N;

			QVector3D tmp = GetUV(u, v);
			vertData.append(tmp.x());
			vertData.append(tmp.y());
			vertData.append(tmp.z());
		}

	// indexpuffer adatai: NxM négyszög = 2xNxM háromszög = háromszöglista esetén 3x2xNxM index
	std::vector<GLuint> indicies;
	for (int i = 0; i<N; ++i)
		for (int j = 0; j<N; ++j)
		{
			indicies.push_back(  i +       j *      (N + 1) );
			indicies.push_back( (i + 1) +  j *      (N + 1) );
			indicies.push_back(  i +      (j + 1) * (N + 1) );
			indicies.push_back( (i + 1) +  j *      (N + 1) );
			indicies.push_back( (i + 1) + (j + 1) * (N + 1) );
			indicies.push_back(  i +      (j + 1) * (N + 1)	);

		}

	sphere_vbo.create();
	sphere_vbo.bind();
	sphere_vbo.allocate(vertData.constData(), vertData.count() * sizeof(GLfloat));
	sphere_vbo.release();
	//-------------

	sphere_indicies.create();
	sphere_indicies.bind();
	sphere_indicies.allocate(&indicies[0], indicies.size() * sizeof(GLuint));
	sphere_indicies.release();
}

void GLRender::makeSphereRaycastObject()
{
	QVector<GLfloat> vertData;
	vertData.append(-1); vertData.append(-1); vertData.append(0);
	vertData.append(1); vertData.append(-1); vertData.append(0);
	vertData.append(-1); vertData.append(1); vertData.append(0);
	vertData.append(1); vertData.append(1); vertData.append(0);
	
	
	m_quad_vb.create();
	m_quad_vb.bind();
	m_quad_vb.allocate(vertData.constData(), vertData.count() * sizeof(GLfloat));
	/*
	if (voronoi_cell_centers.empty() && voronoi_polar_ball_centers.empty()) return;
		qDebug() << voronoi_cell_centers.size() << " surf sphere drawing" << "\n";

		std::vector<AdditionalBind> attr;
		attr.push_back(ALL_MATRICIES);

		m_quad_vb.bind();
		shaderBinding(sphere, attr, m);

		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);

		std::vector<QVector3D> polar_ball_centers;
		std::vector<GLfloat> polar_ball_radius;

		int counter = 0;
		for (auto it : voronoi_cell_centers)
		{
			counter++;
			polar_ball_centers.push_back(it);
			polar_ball_radius.push_back(0.05f);
			//qDebug() << it.first.x() << " " << it.first.y() << " " << it.first.z() << "\n";
			//qDebug() << it.second << "\n";

			sphere->setUniformValue("size", (GLint)polar_ball_centers.size());
			sphere->setUniformValueArray("centers", &polar_ball_centers[0], polar_ball_centers.size());
			sphere->setUniformValueArray("radius", &polar_ball_radius[0], polar_ball_radius.size(), 1);
			sphere->setUniformValue("color", QVector3D(148.0f / 255.0f, 0, 211.0f / 255.0f));

			polar_ball_radius.clear();
			polar_ball_centers.clear();

			glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);		//mxa 100-as részekben rajzoljunk gömböket TODO: késõbb (sok gömbnél)
			if (counter == 8) break;
		}

		qDebug() << voronoi_polar_ball_centers.size() << " polar ball drawing" << "\n";
		
		
		for (auto it : voronoi_polar_ball_centers)
		{
			counter++;
			polar_ball_centers.push_back(it.first);
			polar_ball_radius.push_back(it.second);
		//	qDebug() << it.first.x() << " " << it.first.y() << " " << it.first.z() << "\n";
		//	qDebug() << it.second << "\n";

			sphere->setUniformValue("size", (GLint)polar_ball_centers.size());
			sphere->setUniformValueArray("centers", &polar_ball_centers[0], polar_ball_centers.size());
			sphere->setUniformValueArray("radius", &polar_ball_radius[0], polar_ball_radius.size(), 1);
			sphere->setUniformValue("color", QVector3D(0, 0, 1));

			polar_ball_radius.clear();
			polar_ball_centers.clear();

			glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);		//mxa 100-as részekben rajzoljunk gömböket TODO: késõbb (sok gömbnél)
			if (counter == 12) break;
		}

		m_quad_vb.release();
		sphere->release();
		*/
}

void GLRender::makeObjectFromCloud()
{
	GLPaintFormat paintData = sfg.getPaintData();

	cloud_vbo.destroy();
	elementbuffer.destroy();
	m_quad_vb.destroy();

	if (paintData.points.empty() && paintData.centers_with_radius.empty()) {
		center_radius_pairs.clear();
		colors			   .clear();
		points_colors	   .clear();
		point_draw_parts   .clear();
		center_draw_parts  .clear();
		return; 
	}
	
	switch (sfg.getDrawObject())
	{
		case POINTSS:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			QVector<GLfloat> vertData;
			for (auto it : paintData.points)
			{
				// vertex position
				vertData.append(it.x());
				vertData.append(it.y());
				vertData.append(it.z());
			}

		//	qDebug() << paintData.points.size() << " point incoming" << "\n";
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(vertData.constData(), vertData.count() * sizeof(GLfloat));
			cloud_vbo.release();

			colors = paintData.points_col;
			return;
		}
		case BOX:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
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
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();
			return;
		}
		case SURFACE:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (size_t it = 0; it < paintData.points.size(); ++it)
			{
				// vertex position
				boxData.append(paintData.points[it].x());
				boxData.append(paintData.points[it].y());
				boxData.append(paintData.points[it].z());
				// texture coordinate
				boxData.append(paintData.points_col[it].x());
				boxData.append(paintData.points_col[it].y());
				boxData.append(paintData.points_col[it].z());
			}

			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			return;
		}
		case ALGORITHM:
			makeObjectOfPowerCrustPart(paintData);
		default:
			return;
	}

}

void GLRender::makeObjectOfPowerCrustPart(GLPaintFormat& paintData)
{
	switch (sfg.getPowerCrustDrawObject())
	{
		case	    RESULT								   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;

			bool contain_uv = paintData.points_col.empty();
			for (size_t it = 0; it < paintData.points.size(); ++it)
			{
				// vertex position
				boxData.append(paintData.points[it].x());
				boxData.append(paintData.points[it].y());
				boxData.append(paintData.points[it].z());
				// texture coordinate
				if (contain_uv)
				{
					boxData.append(paintData.uv_coords[it].first);
					boxData.append(paintData.uv_coords[it].second);
				}
				else
				{
					boxData.append(paintData.points_col[it].x());
					boxData.append(paintData.points_col[it].y());
					boxData.append(paintData.points_col[it].z());
				}
			}

			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			points_colors = paintData.points_col;
			point_draw_parts = paintData.point_part_lengths;

			return;
		}
		case		DELAUNEY							   :
		case		VORONOI_DIAGRAM						   :
		case		POWER_DIAGRAM						   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			colors = paintData.col;

			return;
		}
		case		VORONOI_WITH_SURF_POINT				   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;

			return;
		}
		case		POWER_DIAGRAM_BY_CELL:
		case		VORONOI_BY_CELL_FULL				   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			/*voronoi_cell_centers.clear();
			for (auto it : paintData.surf_centers)
				voronoi_cell_centers.push_back(QVector3D(it.x(), it.y(), it.z()));

			voronoi_polar_ball_centers.clear();
			for (auto it : paintData.polar_ball_centers)
				voronoi_polar_ball_centers.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));
*/
			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;

			return;
		}
		case		POWER_DIAGRAM_CELL_WITH_NEIGHBOURS	   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			/*voronoi_cell_centers.clear();

			voronoi_polar_ball_centers.clear();
			for (auto it : paintData.polar_ball_centers)
				voronoi_polar_ball_centers.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			common_point_centers.clear();
			qDebug() << "COMMON POINTS\n";

			for (auto it : paintData.common_points)
			{
				common_point_centers.push_back(QVector3D(it.x(), it.y(), it.z()));
				qDebug() << it.x() << "," << it.y() << "," << it.z();
			}
			qDebug() << "\\/\\/\\/\\/\\/\\/\\/\\/\n";*/
			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;
			point_draw_parts = paintData.point_part_lengths;

			return;
		}
		case		POLES								   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			/*voronoi_cell_centers.clear();
			for (auto it : paintData.surf_centers)
				voronoi_cell_centers.push_back(QVector3D(it.x(), it.y(), it.z()));

			voronoi_polar_ball_centers.clear();
			for (auto it : paintData.polar_ball_centers)
				voronoi_polar_ball_centers.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));
				*/
			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;
			return;
		}
		
		case		INNER_POLES							   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			/*voronoi_polar_ball_centers.clear();
			for (auto it : paintData.polar_ball_centers)
				voronoi_polar_ball_centers.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			*/
			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;

			qDebug() << " WE HAVE " << center_radius_pairs.size() << " INNER pole\n";
			return;
		}
		case		OUTER_POLES							   :
		{
			/*if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			voronoi_polar_ball_centers.clear();
			for (auto it : paintData.polar_ball_centers)
				voronoi_polar_ball_centers.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));
*/
			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;

			qDebug() << " WE HAVE " << center_radius_pairs.size() << " OUTER pole\n";

			return;
		}
		case		UNKNOWN_POLES						   :
		{
			/*if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			voronoi_polar_ball_centers.clear();
			for (auto it : paintData.polar_ball_centers)
				voronoi_polar_ball_centers.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));
*/

			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;

		//	qDebug() << " WE HAVE " << center_radius_pairs.size() << " UNKNOWN pole\n";

			return;
		}
		case		INNER_OUTER_POLES					   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			colors = paintData.col;
			points_colors = paintData.points_col;
			point_draw_parts = paintData.point_part_lengths;

			qDebug() << " WE HAVE " << center_radius_pairs.size() << " INNER and Outer pole\n";

			return;
		}
		case		MEDIAL_AXIS							   :
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------
			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			points_colors = paintData.points_col;
			point_draw_parts = paintData.point_part_lengths;

			return;
		}
		case VORONOI_WITH_CELLDUAL:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			center_radius_pairs.clear();
			for (auto it : paintData.centers_with_radius)
				center_radius_pairs.push_back(std::pair<QVector3D, float>(QVector3D(it.first.x(), it.first.y(), it.first.z()), it.second));

			colors = paintData.col;
			points_colors = paintData.points_col;
			center_draw_parts = paintData.center_part_lengths;
			point_draw_parts = paintData.point_part_lengths;
			return;
		}
		case POWER_SHAPE:
		{
			if (cloud_vbo.isCreated()) cloud_vbo.release();
			if (elementbuffer.isCreated()) elementbuffer.release();

			QVector<GLfloat> boxData;
			for (auto it : paintData.points)
			{
				// vertex position
				boxData.append(it.x());
				boxData.append(it.y());
				boxData.append(it.z());
			}
			cloud_vbo.create();
			cloud_vbo.bind();
			cloud_vbo.allocate(boxData.constData(), boxData.count() * sizeof(GLfloat));
			cloud_vbo.release();
			//-------------

			elementbuffer.create();
			elementbuffer.bind();
			elementbuffer.allocate(&paintData.ix[0], paintData.ix.size() * sizeof(GLuint));
			elementbuffer.release();

			points_colors = paintData.points_col;
			point_draw_parts = paintData.point_part_lengths;

			return;
		}
		default:
			return;
	
	}
}

void GLRender::shaderBinding(QOpenGLShaderProgram* shader, std::vector<AdditionalBind> additional_bind, QMatrix4x4 mvp)
{
	shader->bind();
	shader->setUniformValue("matrix", mvp);		//TODO replace name
	shader->enableAttributeArray(PROGRAM_VERTEX_ATTRIBUTE);
	shader->setAttributeBuffer(PROGRAM_VERTEX_ATTRIBUTE, GL_FLOAT, 0, 3, 3 * sizeof(GLfloat));

	for (auto it : additional_bind)
	{
		switch(it)
		{
		case COLOR:
		{
			shader->setAttributeBuffer(PROGRAM_VERTEX_ATTRIBUTE, GL_FLOAT, 0, 3, 6 * sizeof(GLfloat));
			shader->enableAttributeArray(PROGRAM_COLOR_ATTRIBUTE);
			shader->setAttributeBuffer(PROGRAM_COLOR_ATTRIBUTE, GL_FLOAT, 3 * sizeof(GLfloat), 3, 6 * sizeof(GLfloat));
			return;
		}
		case TEXTURE:
		{
			kinect_textures[k_ind]->bind();
		//	qDebug() << k_ind;
			shader->setAttributeBuffer(PROGRAM_VERTEX_ATTRIBUTE, GL_FLOAT, 0, 3, 5 * sizeof(GLfloat));
			shader->enableAttributeArray(PROGRAM_TEXCOORD_ATTRIBUTE);
			shader->setAttributeBuffer(PROGRAM_TEXCOORD_ATTRIBUTE, GL_FLOAT, 3 * sizeof(GLfloat), 2, 5 * sizeof(GLfloat));
			return;
		}
		case ALL_MATRICIES:
		{
			shader->setUniformValue("model", mvp);
			shader->setUniformValue("modelI", mvp.inverted());

			shader->setUniformValue("viewProj", m_camera.GetViewProj());
			shader->setUniformValue("viewIprojI", (m_camera.GetProj() * m_camera.GetViewMatrix()).inverted());
			shader->setUniformValue("view", m_camera.GetViewMatrix());
		}
		default:
			return;
		}
	}
}

void GLRender::renderObject(const QMatrix4x4& mvp, const QMatrix4x4& m)
{
	//TODO FPS COUNTER WITH INDEXING
	//http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-9-vbo-indexing/

	switch (sfg.getDrawObject())
	{
		case POINTSS:
		{
			if (!cloud_vbo.isCreated()) return;

			cloud_vbo.bind();
		//	qDebug() << cloud_vbo.size() / sizeof(GLfloat) / 3 << " point drawing" << "\n";
			shaderBinding(p_base, std::vector<AdditionalBind>(), mvp);
			p_base->setUniformValue("color", colors[0]);

			glDisable(GL_BLEND);
			glDepthMask(GL_TRUE);

			glDrawArrays(GL_POINTS, 0, cloud_vbo.size() / sizeof(GLfloat));
			cloud_vbo.release();
			p_base->release();
			return;
		}
		case BOX:
		{
			/*glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glDepthMask(GL_FALSE);
			cloud_vbo.bind();
			elementbuffer.bind();

			glDrawElements(
				GL_TRIANGLES,      // mode
				drawing_size,    // count
				GL_UNSIGNED_INT,   // type
				(void*)0           // element array buffer offset
			);*/
			return;
		}
		case SURFACE:
		{
			if (!elementbuffer.isCreated()) return;

			cloud_vbo.bind();
			elementbuffer.bind();
		//	qDebug() << elementbuffer.size() / sizeof(GLfloat) << " surface indicies drawing" << "\n";

			std::vector<AdditionalBind> attr;
			attr.push_back(COLOR);
			shaderBinding(pc_base, attr, mvp);

			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glDepthMask(GL_FALSE);

			glDrawElements(
				GL_TRIANGLES,								// mode
				elementbuffer.size() / sizeof(GLfloat),     // count
				GL_UNSIGNED_INT,						    // type
				(void*)0									// element array buffer offset
			);

			elementbuffer.release();
			cloud_vbo.release();
			p_base->release();
			return;
		}
		case ALGORITHM:
		{
			renderPowerCrustObject(mvp, m);
			return;
		}
		case TEST:
		{
			qDebug() << "Sphere mesh binded\n";

			sphere_vbo.bind();
			sphere_indicies.bind();
			qDebug() << sphere_indicies.size() / sizeof(GLfloat) << " sphere drawing" << "\n";
			shaderBinding(p_base, std::vector<AdditionalBind>(), mvp);

			glDisable(GL_BLEND);
			glDepthMask(GL_FALSE);

			glDrawElements(
				GL_TRIANGLES,      // mode
				sphere_indicies.size() / sizeof(GLfloat),    // count
				GL_UNSIGNED_INT,   // type
				(void*)0           // element array buffer offset
			);

			sphere_indicies.release();
			sphere_vbo.release();

			p_base->release();

			return;
		}
		default:
			return;
			//LINE DRAW
			/*glDisable(GL_BLEND);
			glDepthMask(GL_TRUE);
			cloud_vbo.bind();

			glDrawArrays(
				GL_LINES,      // mode
				0,
				drawing_size    // count
			);*/
	}
}

void GLRender::DrawCells(QOpenGLShaderProgram* shader, QOpenGLBuffer& point_buffer, QOpenGLBuffer& index_buffer, const int& starter_index, const int& last_index, const QMatrix4x4& mvp, QVector3D& col, bool transparent)
{
	if (!elementbuffer.isCreated()) return;

	cloud_vbo.bind();
	elementbuffer.bind();
//	qDebug() << elementbuffer.size() / sizeof(GLfloat) << " cell drawing" << "\n";
//	qDebug() << "start : " << starter_index << " " << last_index;
	std::vector<AdditionalBind> binds;
	if (shader == pc_base) binds.push_back(COLOR);
	if (shader == pt_base) binds.push_back(TEXTURE);

	shaderBinding(shader, binds, mvp);
	if (shader == p_base) shader->setUniformValue("color", col);

	if (transparent)
	{
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
	}
	else
	{
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
	}

	glDrawElements(
		GL_TRIANGLES,      // mode
		last_index,    // count
		GL_UNSIGNED_INT,   // type
		(void*)(starter_index*sizeof(unsigned int))           // element array buffer offset
	);

	elementbuffer.release();
	cloud_vbo.release();
	shader->release();
}

void GLRender::DrawSphere(QOpenGLShaderProgram* shader, QOpenGLBuffer& sphere_vbo, QOpenGLBuffer& sphere_indicies, const QMatrix4x4& mvp, std::vector<std::pair<QVector3D, float>>& centers_with_radius, const int& size, QVector3D& col, bool transparent)
{
	if (centers_with_radius.empty()) return;

	sphere_vbo.bind();
	sphere_indicies.bind();

	std::vector<AdditionalBind> binds;
	if (shader == pc_base) binds.push_back(COLOR);
	if (shader == pt_base) binds.push_back(TEXTURE);

	shaderBinding(shader, binds, mvp);

	if (transparent)
	{
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
	}
	else
	{
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
	}

	QMatrix4x4 transforms;

//	qDebug() << size << "sphere drawing" << "\n";

	for (std::vector<std::pair<QVector3D, float>>::iterator it = centers_with_radius.begin(); it != centers_with_radius.begin() + size; ++it)
	{
		transforms.setToIdentity();
		transforms.translate(it->first);
		transforms.scale(it->second);
		shader->setUniformValue("matrix", mvp * transforms);		//TODO replace name
		shader->setUniformValue("color", col);

		glDrawElements(
			GL_TRIANGLES,      // mode
			sphere_indicies.size() / sizeof(GLfloat),    // count
			GL_UNSIGNED_INT,   // type
			(void*)0           // element array buffer offset
		);

	}

	sphere_indicies.release();
	sphere_vbo.release();

	shader->release();
}

void GLRender::DrawSegments(QOpenGLShaderProgram* shader, QOpenGLBuffer& point_buffer, QOpenGLBuffer& index_buffer, const int& starter_index, const int& last_index, const QMatrix4x4& mvp, QVector3D& col, bool transparent)
{
	if (!elementbuffer.isCreated()) return;

	cloud_vbo.bind();
//	qDebug() << cloud_vbo.size() / sizeof(GLfloat) << " poi" << "\n";
	elementbuffer.bind();
//	qDebug() << starter_index << " start ind" << "\n";
//	qDebug() << last_index << " last ind" << "\n";

	std::vector<AdditionalBind> binds;
	if (shader == pc_base) binds.push_back(COLOR);
	if (shader == pt_base) binds.push_back(TEXTURE);

	shaderBinding(shader, binds, mvp);

	if (shader == p_base) shader->setUniformValue("color", col);


	if (transparent)
	{
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
	}
	else
	{
		glDisable(GL_BLEND);
		glDepthMask(GL_TRUE);
	}

	glDrawElements(
		GL_LINES,      // mode
		last_index,    // count
		GL_UNSIGNED_INT,   // type
		(void*) (starter_index * sizeof(unsigned int))           // element array buffer offset
	);

	elementbuffer.release();
	cloud_vbo.release();
	shader->release();
}

void GLRender::renderPowerCrustObject(const QMatrix4x4& mvp, const QMatrix4x4& m)
{
	switch (sfg.getPowerCrustDrawObject())
	{
		case	RESULT:
		{
			elementbuffer.bind();
			if (elementbuffer.size() == 0) return;
			//qDebug() << point_draw_parts.size();
			if (point_draw_parts.empty()) return;

			int start = 0;
	
			if (points_colors.empty())		//we use uv
			{
				k_ind = 0;
				DrawCells(pt_base, cloud_vbo, elementbuffer, start, point_draw_parts[0], mvp, QVector3D(0,0,0), false);
				k_ind = 1;
				start += point_draw_parts[0];
				DrawCells(pt_base, cloud_vbo, elementbuffer, start, point_draw_parts[1], mvp, QVector3D(0, 0, 0), false);
				k_ind = 2;
				start += point_draw_parts[1];
				DrawCells(pt_base, cloud_vbo, elementbuffer, start, point_draw_parts[2], mvp, QVector3D(0, 0, 0), false);
				k_ind = 3;
				start += point_draw_parts[2];
				DrawCells(pt_base, cloud_vbo, elementbuffer, start, point_draw_parts[3], mvp, QVector3D(0, 0, 0), false);
			}
			else
			{
				DrawCells(pc_base, cloud_vbo, elementbuffer, start, elementbuffer.size() / sizeof(GLfloat), mvp, points_colors[0], false);
			}
			return;
		}
		case		DELAUNEY:
		{
			elementbuffer.bind();
			if (elementbuffer.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, elementbuffer.size() / sizeof(GLfloat), mvp, colors[0], false);
			return;
		}
		case		POWER_DIAGRAM:
		case		VORONOI_DIAGRAM:
		{
			elementbuffer.bind();
			if (elementbuffer.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, elementbuffer.size() / sizeof(GLfloat), mvp, colors[0],false);
			return;
		}
		case		VORONOI_WITH_SURF_POINT:
		{
			//DrawCells(cloud_vbo, elementbuffer, mvp, points_colors[0]);
			elementbuffer.bind();
			if (elementbuffer.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, elementbuffer.size() / sizeof(GLfloat), mvp, points_colors[0]);
			//-------------------------------
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, center_radius_pairs, center_draw_parts[0], colors[0]);
			return;
		}
		case		VORONOI_BY_CELL_FULL:
		{
			/*if (!elementbuffer.isCreated()) return;

			cloud_vbo.bind();
			elementbuffer.bind();
			qDebug() << elementbuffer.size() / sizeof(GLfloat) << " cell drawing" << "\n";

			shaderBinding(base, std::vector<AdditionalBind>(), mvp);
			base->setUniformValue("color", QVector3D(1, 0, 0));

			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glDepthMask(GL_FALSE);

			glDrawElements(
				GL_TRIANGLES,      // mode
				elementbuffer.size() / sizeof(GLfloat),    // count
				GL_UNSIGNED_INT,   // type
				(void*)0           // element array buffer offset
			);

			elementbuffer.release();
			cloud_vbo.release();
			base->release();*/
			//DrawCells(cloud_vbo, elementbuffer, mvp, points_colors[0]);

			elementbuffer.bind();
			if (elementbuffer.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, elementbuffer.size() / sizeof(GLfloat), mvp, points_colors[0]);
			//-------------------------------

			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, center_radius_pairs, center_draw_parts[0], colors[0]);

			/*qDebug() << voronoi_polar_ball_centers.size() << " polar sphere drawing" << "\n";

			sphere_vbo.bind();
			sphere_indicies.bind();

			shaderBinding(base, std::vector<AdditionalBind>(), mvp);

			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable(GL_BLEND);
			glDepthMask(GL_FALSE);

			QMatrix4x4 transforms;

			for (auto it : voronoi_polar_ball_centers)
			{
				transforms.setToIdentity();
				transforms.translate(it.first);
				transforms.scale(it.second);
				base->setUniformValue("matrix", m_camera.GetViewProj() * m * transforms);
				base->setUniformValue("color", QVector3D(0, 0, 1));

				glDrawElements(
					GL_TRIANGLES,      // mode
					sphere_indicies.size() / sizeof(GLfloat),    // count
					GL_UNSIGNED_INT,   // type
					(void*)0           // element array buffer offset
				);
			}*/

			std::vector<std::pair<QVector3D, float>> poles(center_radius_pairs.begin() + center_draw_parts[0], center_radius_pairs.end());
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, poles, center_draw_parts[1], colors[1]);

			/*qDebug() << voronoi_cell_centers.size() << " voronoi sphere drawing" << "\n";

			for (auto it : voronoi_cell_centers)
			{
				transforms.setToIdentity();
				transforms.translate(it);
				transforms.scale(0.05f);
				base->setUniformValue("matrix", m_camera.GetViewProj() * m * transforms);		//TODO replace name
				base->setUniformValue("color", QVector3D(148.0f / 255.0f, 0, 211.0f / 255.0f));

				glDrawElements(
					GL_TRIANGLES,      // mode
					sphere_indicies.size() / sizeof(GLfloat),    // count
					GL_UNSIGNED_INT,   // type
					(void*)0           // element array buffer offset
				);

			}

			sphere_indicies.release();
			sphere_vbo.release();

			base->release();
			*/
			return;
		}
		case		POWER_DIAGRAM_BY_CELL:
		{
		//	DrawCells(cloud_vbo, elementbuffer, mvp, points_colors[0]);
			elementbuffer.bind();
			if (elementbuffer.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, elementbuffer.size() / sizeof(GLfloat), mvp, points_colors[0]);
			//-------------------------------
			std::vector<std::pair<QVector3D, float>> polar_balls = center_radius_pairs;
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, polar_balls, center_draw_parts[0], colors[0]);
			return;
	}
		case		POWER_DIAGRAM_CELL_WITH_NEIGHBOURS:
		{
		//	DrawCells(cloud_vbo, elementbuffer, mvp, points_colors[0]);
			if (point_draw_parts.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, point_draw_parts[0], mvp, points_colors[0]);
			DrawSegments(p_base, cloud_vbo, elementbuffer, point_draw_parts[0], point_draw_parts[1], mvp, points_colors[1]);
				//-------------------------------

			std::vector<std::pair<QVector3D, float>> polar_balls = center_radius_pairs;
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, polar_balls, center_draw_parts[0], colors[0]);

			std::vector<std::pair<QVector3D, float>> common_points(center_radius_pairs.begin() + center_draw_parts[0], center_radius_pairs.end());
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, common_points, center_draw_parts[1], colors[1]);
			return;
		}
		case		POLES:
		{
			if (center_radius_pairs.size() == 0) return;

			std::vector<std::pair<QVector3D, float>> surface_centers = center_radius_pairs;
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, surface_centers, center_draw_parts[0], colors[0]);

			std::vector<std::pair<QVector3D, float>> polar_balls(center_radius_pairs.begin() + center_draw_parts[0], center_radius_pairs.end());
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, polar_balls, center_draw_parts[1], colors[1]);

			return;
		}
		case		INNER_POLES:
		case		OUTER_POLES:
		case		UNKNOWN_POLES:
		{
			if (center_radius_pairs.size() == 0) return;

			std::vector<std::pair<QVector3D, float>> poles = center_radius_pairs;
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, poles, center_draw_parts[0], colors[0]);
			return;
		}
		case  INNER_OUTER_POLES:
		{
			if (point_draw_parts.size() == 0) return;
			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, point_draw_parts[0], mvp, points_colors[0]);
			DrawSegments(p_base, cloud_vbo, elementbuffer, point_draw_parts[0], point_draw_parts[1], mvp, points_colors[1]);

			return;
		}
		case VORONOI_WITH_CELLDUAL:
		{
			if (point_draw_parts.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, point_draw_parts[0], mvp, points_colors[0]);
			DrawSegments(p_base, cloud_vbo, elementbuffer, point_draw_parts[0], point_draw_parts[1], mvp, points_colors[1]);

			std::vector<std::pair<QVector3D, float>> surface_centers = center_radius_pairs;
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, surface_centers, center_draw_parts[0], colors[0]);

			std::vector<std::pair<QVector3D, float>> polar_balls(center_radius_pairs.begin() + center_draw_parts[0], center_radius_pairs.end());
			DrawSphere(p_base, sphere_vbo, sphere_indicies, mvp, polar_balls, center_draw_parts[1], colors[1]);

			return;
		}
		case		MEDIAL_AXIS:
		{
		//	qDebug() << point_draw_parts.size();
			if (point_draw_parts.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, point_draw_parts[0], mvp, points_colors[0], false);

			return;
		}
		case POWER_SHAPE:
		{
		//	qDebug() << point_draw_parts.size();
			if (point_draw_parts.size() == 0) return;

			DrawSegments(p_base, cloud_vbo, elementbuffer, 0, point_draw_parts[0], mvp, points_colors[0], false);

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
	data = 0;
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

void GLRender::attemptChangePowerCrustDrawType(const unsigned int& ind)
{
	if (sfg.getPowerCrustDrawObject() - ind)
	{
		sfg.setPowerCrustDrawObject(ind);
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

void GLRender::goCell(const int& direction)
{
	switch (direction)
	{
	case 0:
	{
		sfg.goPreviousCell();
		break;
	}
	case 1:
	{
		sfg.goNextCell();
		break;
	}
	default:
		break;
	}

	makeObjectFromCloud();	//todo: akkor is make ha nem ez az aktuális
	update();
}

void GLRender::nextNeighb()
{
	sfg.nextNeighb();
	makeObjectFromCloud();
	update();
}

void GLRender::swapStartFlag()
{
	sfg.swapStartFlag();
	//TODO: stop??

	makeObjectFromCloud();
	update();
}

//---------------------------------------

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
	m_camera.wheelEvent(event);
	update();
}

void GLRender::keyPressEvent(QKeyEvent *event)
{
	switch (event->key())
	{
	default:
		m_camera.KeyboardDown(event);
	}
}

void GLRender::keyReleaseEvent(QKeyEvent *event)
{
	switch (event->key())
	{
	default:
		m_camera.KeyboardUp(event);
	}
}
/*void GLRender::showNearest()
{
cloud->nearest();
makeObjectFromCloud();
update();
}

void GLRender::showKNearest()
{
cloud->knearest();
makeObjectFromCloud();
update();
}
*/