#include <QtWidgets>

#include "glrender.h"
#include "graphicinterface.h"

GraphicInterface::GraphicInterface()
{
	ui.glrender = new GLRender();

	timer = new QTimer(this);
	timer->start(20);

	fps = new QLabel("FPS: ",this);
	_fps = new QLCDNumber(4,this);
	join = new QLabel("Join(ms): ", this);
	_join = new QLCDNumber(6, this);
	sort = new QLabel("Sort(ms): ", this);
	_sort = new QLCDNumber(7, this);
	all = new QLabel("All(ms): ", this);
	_all = new QLCDNumber(5, this);
	poi = new QLabel("Poi(dat): ", this);
	_poi = new QLCDNumber(5, this);
	poi2 = new QLabel("Poi(box): ", this);
	_poi2 = new QLCDNumber(6, this);
	
	ui.setupUi(this);
	setCalculationBar();

	//---------------------------------
	ui.glrender->setClearColor(QColor(255, 255, 0));
	ui.glrender->attemptChangeDrawType(ui.display->currentIndex());
	setCore(ui.core->value());

	//---------------------------------
	setupSignals();
}

void GraphicInterface::setupSignals()
{
	connect(timer, &QTimer::timeout, this, &GraphicInterface::attemptUpdate);

	connect(ui.actionExit, &QAction::triggered, this, &GraphicInterface::stop);
	connect(ui.actionExport, &QAction::triggered, this, &GraphicInterface::exportSurface);
	connect(ui.actionInfo, &QAction::triggered, this, &GraphicInterface::info);
	connect(ui.actionNew, &QAction::triggered, this, &GraphicInterface::newConnection);
	connect(ui.actionOpen, &QAction::triggered, this, &GraphicInterface::open);
	connect(ui.actionSave, &QAction::triggered, this, &GraphicInterface::save);
	connect(ui.actionSave_as, &QAction::triggered, this, &GraphicInterface::saveas);

	connect(ui.compButton, SIGNAL(clicked()), this, SLOT(compChanged()));

	connect(ui.left_child, SIGNAL(clicked()), this, SLOT(goLeftInTree()));
	connect(ui.right_child, SIGNAL(clicked()), this, SLOT(goRightInTree()));
	connect(ui.parent, SIGNAL(clicked()), this, SLOT(goParentInTree()));

	connect(ui.core, SIGNAL(valueChanged(int)), this, SLOT(setCore(const int&)));
	connect(ui.thread, SIGNAL(valueChanged(int)), this, SLOT(setThread(const int&)));

	connect(ui.display, SIGNAL(currentIndexChanged(int )), this, SLOT(setDisplayType(const int &)));

	connect(ui.glrender, &GLRender::step_finished, this, &GraphicInterface::refreshCalculationBar);
	connect(ui.glrender, &GLRender::clicked, this, &GraphicInterface::rotateOneStep);
}

void GraphicInterface::setCalculationBar()
{
	ui.statusBar->addWidget(fps);
	ui.statusBar->addWidget(_fps);
	_fps->display(0.00);

	ui.statusBar->addWidget(all);
	ui.statusBar->addWidget(_all);
	_all->display(0.00);

	ui.statusBar->addWidget(sort);
	ui.statusBar->addWidget(_sort);
	_sort->display(0.00);

	ui.statusBar->addWidget(join);
	ui.statusBar->addWidget(_join);
	_join->display(0.00);

	ui.statusBar->addWidget(poi);
	ui.statusBar->addWidget(_poi);
	_poi->display(0);

	ui.statusBar->addWidget(poi2);
	ui.statusBar->addWidget(_poi2);
	_poi2->display(0);
}

void GraphicInterface::refreshCalculationBar()
{
	float all_t, join_t, sort_t;
	float fps_t;
	int poi, poi2;
	ui.glrender->getCloudCalculationTime(fps_t, join_t, sort_t, all_t);
	ui.glrender->getCloudPoints(poi,poi2);

	_fps->display(fps_t);
	_all->display(all_t);
	_sort->display(sort_t);
	_join->display(join_t);
	_poi->display(poi);
	_poi2->display(poi2);
}

void GraphicInterface::stop()
{
	this->close();
}

void GraphicInterface::exportSurface()
{
	int a = 2;
}

void GraphicInterface::info()
{
	int a = 2;
}

void GraphicInterface::open()
{
	int a = 2;
}

void GraphicInterface::save()
{
	int a = 2;
}

void GraphicInterface::saveas()
{
	int a = 2;
}

void GraphicInterface::newConnection()
{
	ui.compButton->setText("Stop");
	ui.glrender->resetCloud();
}

//TODO: replace with lambda?
void GraphicInterface::goLeftInTree()
{
	ui.glrender->goNode(1);
}	
void GraphicInterface::goRightInTree()
{
	ui.glrender->goNode(2);
}
void GraphicInterface::goParentInTree()
{
	ui.glrender->goNode(0);
}
void GraphicInterface::rotateOneStep()
{
	//ui.glrender->rotateBy(+2 * 16, +2 * 16, -1 * 16);
}

void GraphicInterface::attemptUpdate()
{
	ui.glrender->attemptUpdate();
}

void GraphicInterface::setCore(const int& core)
{
	//TODO: rethink later
	//main thread: render, 1-2 thread for surface making and the remains for cloud operations
	unsigned int threads = (unsigned int) std::ceil((float)core * 1.5);
	unsigned int surf_t = 1 + threads / 8;
	unsigned int cloud_t = threads - 1 - surf_t;

	ui.glrender->attemptChangeThreadCapacity(cloud_t, surf_t);
}

void GraphicInterface::setThread(const int& i)
{
	int a = 2;
}

void GraphicInterface::setDisplayType(const int &ind)
{
	ui.glrender->attemptChangeDrawType(ind);
}

void GraphicInterface::compChanged()
{
	ui.compButton->setText(ui.compButton->text() == "Stop" ? "Start" : "Stop");
	ui.glrender->changeComputation();
}