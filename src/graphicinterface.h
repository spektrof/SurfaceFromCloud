#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets>
#include "qnspinbox.h"

#include "ui_sfc.h"
#include "ui_new.h"
#include "ui_new_outl_filter.h"
#include "ui_new_simpl_filter.h"
#include "ui_new_smoth_filter.h"
#include "ui_new_surf_rec.h"
#include "ui_filter_widget.h"

class GLWidget;
class GLNSpinBox;

class GraphicInterface : public QMainWindow
{
	Q_OBJECT

public:
	GraphicInterface(/*QWidget *parent = Q_NULLPTR*/);

private slots:
void newConnection();
void exportSurface();
void info();
void open();
void save();
void saveas();

void acceptNew();
void acceptNewOutlierFilter();
void acceptNewSimplFilter();
void acceptNewSmothFilter();
void acceptNewSurf();
void deleteFilter();

void cancelDialog();
void openFileDial();

void sourceChanged(const int &);
void setMaxPoi(const int &);
void setCore(const int&);
void swapStartFlag();
void noMoreData();

//void showNearest();
//void showKNearest();

protected:
	void setupSignals();
	void setupFilterActions();

	void openOutlierFilterPropertyWindow();
	void openSimplifierPropertyWindow();
	void openSmootherFilterPropertyWindow();
	void openSurfacePropertyWindow();

	void setCalculationBar();
	void refreshCalculationBar();
private:
	Ui::sfcGUI ui;
	Ui::SfcNew ui_new;
	Ui::sfcOutF ui_out_filter;
	Ui::sfcSimpF ui_simpl_filter;
	Ui::sfcSmoF ui_smo_filter;
	Ui::sfcSurf ui_surf;

	QDialog* dialog;

	std::vector<Ui::sfcWidgetF*> ui_filters;
	QWidget* widget;
	unsigned int active_filter_item;

	unsigned int hovered_filter_item;
	unsigned int hovered_surface_item;

	QTimer *timer;

	QLabel *c_load;
	QLabel *join;
	QLabel *sort;
	QLabel *all;
	QLabel *poi;
	QLabel *s_cal;
	QLCDNumber *_c_load;
	QLCDNumber *_join;
	QLCDNumber *_sort;
	QLCDNumber *_all;
	QLCDNumber *_poi;
	QLCDNumber * _s_cal;
};

#endif