/*
MSc thesis work: Surface generating from pointcloud -possibly in real time -
Peter Lukacs
*/

#include <QApplication>

#include "graphicinterface.h"

int main(int argc, char *argv[])
{
	Q_INIT_RESOURCE(sfc);

	QApplication app(argc, argv);

	GraphicInterface window;
	window.show();
	return app.exec();
}