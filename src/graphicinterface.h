#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets>
#include "qnspinbox.h"

#include "ui_sfc.h"

class GLWidget;
class GLNSpinBox;

class GraphicInterface : public QMainWindow
{
	Q_OBJECT

public:
	GraphicInterface(/*QWidget *parent = Q_NULLPTR*/);

private slots:
void rotateOneStep();
void attemptUpdate();
void goLeftInTree();
void goRightInTree();
void goParentInTree();
void setCore(const int&);
void setThread(const int&);
void setDisplayType(const int &text);
void stop();
void exportSurface();
void info();
void open();
void save();
void saveas();
void newConnection();
void compChanged();

protected:
	void setupSignals();

	void setCalculationBar();
	void refreshCalculationBar();
private:
	Ui::sfcGUI ui;

	QTimer *timer;

	QLabel *fps;
	QLabel *join;
	QLabel *sort;
	QLabel *all;
	QLabel *poi;
	QLabel *poi2;
	QLCDNumber *_fps;
	QLCDNumber *_join;
	QLCDNumber *_sort;
	QLCDNumber *_all;
	QLCDNumber *_poi;
	QLCDNumber *_poi2;
};

#endif