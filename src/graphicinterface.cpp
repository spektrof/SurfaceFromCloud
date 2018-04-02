#include <QtWidgets>
#include <QFileInfo>
#include <QMessageBox>
#include <qdialog.h>

#include "graphicinterface.h"

GraphicInterface::GraphicInterface()
{
	dialog = NULL;
	timer = new QTimer(this);
	timer->start(20);

	c_load = new QLabel("Load: ",this);
	_c_load = new QLCDNumber(4,this);
	join = new QLabel("Join(ms): ", this);
	_join = new QLCDNumber(4, this);
	sort = new QLabel("Sort(ms): ", this);
	_sort = new QLCDNumber(4, this);
	all = new QLabel("All(ms): ", this);
	_all = new QLCDNumber(5, this);
	poi = new QLabel("Poi(dat): ", this);
	_poi = new QLCDNumber(8, this);
	s_cal = new QLabel("Surf cal: ", this);
	_s_cal = new QLCDNumber(6, this);
	
	ui.setupUi(this);

	setupFilterActions();

	setCalculationBar();
	ui.max_poi->setText("1");

	//---------------------------------
	ui.glrender->setClearColor(QColor(255, 255, 0));
	ui.glrender->attemptChangeDrawType(ui.display->currentIndex());
//	ui.glrender->attemptChangePowerCrustDrawType(ui.power_c_draw->currentIndex());
	ui.forkdtree->setVisible(false);
	ui.powercrust->setVisible(false);
	ui.poisson->setVisible(false);
	ui.forkdtree->setVisible(false);
	setCore(ui.core->value());

	hovered_surface_item = 2;
	acceptNewSurf();

	//---------------------------------
	setupSignals();

	active_filter_item = 0;
}

void GraphicInterface::setupSignals()
{
	connect(timer, &QTimer::timeout, this, [this]() { this->ui.glrender->attemptUpdate(this->timer->interval()); } );

	connect(ui.actionExit, &QAction::triggered, this, [this]() { this->close(); });
	connect(ui.actionExport, &QAction::triggered, this, &GraphicInterface::exportSurface);
	connect(ui.actionInfo, &QAction::triggered, this, &GraphicInterface::info);
	connect(ui.actionNew, &QAction::triggered, this, &GraphicInterface::newConnection);
	connect(ui.actionOpen, &QAction::triggered, this, &GraphicInterface::open);
	connect(ui.actionSave, &QAction::triggered, this, &GraphicInterface::save);
	connect(ui.actionSave_as, &QAction::triggered, this, &GraphicInterface::saveas);

	connect(ui.compButton, SIGNAL(clicked()), this, SLOT(swapStartFlag()));
	connect(ui.surf, &QPushButton::clicked, this, [this]() { this->ui.glrender->calculate_surface(); this->refreshCalculationBar(); });
	//connect(ui.nearest, SIGNAL(clicked()), this, SLOT(showNearest()));
	//connect(ui.knearest, SIGNAL(clicked()), this, SLOT(showKNearest()));

	connect(ui.left_child, &QPushButton::clicked, this, [this]() { this->ui.glrender->goNode(1); });
	connect(ui.right_child, &QPushButton::clicked, this, [this]() { this->ui.glrender->goNode(2); });
	connect(ui.parent, &QPushButton::clicked, this, [this]() { this->ui.glrender->goNode(0); });
	connect(ui.parent, &QPushButton::clicked, this, [this]() { this->ui.glrender->nextCell(); });

	connect(ui.core, SIGNAL(valueChanged(int)), this, SLOT(setCore(const int&)));
	
	connect(ui.display, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [this](const int & ind) { this->ui.glrender->attemptChangeDrawType(ind); });

	connect(ui.max_slider, SIGNAL(valueChanged(int)), this, SLOT(setMaxPoi(const int &)));

	connect(ui.glrender, &GLRender::step_finished, this, &GraphicInterface::refreshCalculationBar);
	connect(ui.glrender, &GLRender::data_ended, this, &GraphicInterface::noMoreData);
	
	connect(ui.glrender, &GLRender::empty_points, this, [this]() {	QMessageBox::warning(this, "Oops", "Press Start before Surface!"); });
	connect(ui.glrender, &GLRender::not_enable_component, this, [this]() {	QMessageBox::warning(this, "Oops", "No Surface component found, recompile with x FLAG!"); });
	connect(ui.glrender, &GLRender::unknown_error, this, [this]() {	QMessageBox::warning(this, "Oops", "UNKNOWN error during surface calculation!"); });

	//------------------
	connect(ui.filters, &QToolBox::currentChanged, this, [this](const int & id) { this->active_filter_item = id; } );


	//------------------- POWERCRUST

	connect(ui.next_neigh, &QPushButton::clicked, this, [this]() { this->ui.glrender->nextNeighb(); });
	connect(ui.next_cell, &QPushButton::clicked, this, [this]() { this->ui.glrender->goCell(1); });
	connect(ui.prev_cell, &QPushButton::clicked, this, [this]() { this->ui.glrender->goCell(0); });
	connect(ui.power_c_draw, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged), [this](const int & ind) { this->ui.glrender->attemptChangePowerCrustDrawType(ind); });

	//------------------- POISSON

	connect(ui.pois_angle, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val) {this->ui.glrender->SetSurfaceProperty(0, ui.pois_angle->value(), ui.pois_maxtri->value(), ui.pois_surfapprerr->value()); });
	connect(ui.pois_maxtri, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val) {this->ui.glrender->SetSurfaceProperty(0, ui.pois_angle->value(), ui.pois_maxtri->value(), ui.pois_surfapprerr->value()); });
	connect(ui.pois_surfapprerr, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val) {this->ui.glrender->SetSurfaceProperty(0, ui.pois_angle->value(), ui.pois_maxtri->value(), ui.pois_surfapprerr->value()); });

}

void GraphicInterface::setupFilterActions()
{
	ui.filters->removeItem(0);

	ui.actionSimplifiers->setDisabled(true);
	ui.actionSmoothers->setDisabled(true);
	ui.actionOutlier_removals->setDisabled(true);

	//Outlier Removals
	connect(ui.actionCGAL_based_outl, &QAction::triggered, this, &GraphicInterface::openOutlierFilterPropertyWindow);
	connect(ui.actionCGAL_based_outl, &QAction::hover, this, [this]() { { this->hovered_filter_item = 1; }});
	connect(ui.actionPlane_base_out, &QAction::triggered, this, &GraphicInterface::openOutlierFilterPropertyWindow);
	connect(ui.actionPlane_base_out, &QAction::hover, this, [this]() { { this->hovered_filter_item = 2; }});

	//Simplifier
	connect(ui.actionCGAL_based_Grid, &QAction::triggered, this, &GraphicInterface::openSimplifierPropertyWindow);
	connect(ui.actionCGAL_based_Grid, &QAction::hover, this, [this]() { { this->hovered_filter_item = 4; }});
	connect(ui.actionCGAL_based_Hieararchy, &QAction::triggered, this, &GraphicInterface::openSimplifierPropertyWindow);
	connect(ui.actionCGAL_based_Hieararchy, &QAction::hover, this, [this]() { { this->hovered_filter_item = 5; }});
	connect(ui.actionGrid_mine, &QAction::triggered, this, &GraphicInterface::openSimplifierPropertyWindow);
	connect(ui.actionGrid_mine, &QAction::hover, this, [this]() { { this->hovered_filter_item = 3; }});
	connect(ui.actionCGAL_WLOP_algorithm, &QAction::triggered, this, &GraphicInterface::openSimplifierPropertyWindow);
	connect(ui.actionCGAL_WLOP_algorithm, &QAction::hover, this, [this]() { { this->hovered_filter_item = 6; }});

	//Smoother
	connect(ui.actionCGAL_based_Jet_smoother, &QAction::triggered, this, &GraphicInterface::openSmootherFilterPropertyWindow);
	connect(ui.actionCGAL_based_Jet_smoother, &QAction::hover, this, [this]() { { this->hovered_filter_item = 7; }});

	//Surface Rec
	connect(ui.actionCGAL_based_Poisson, &QAction::triggered, this, &GraphicInterface::openSurfacePropertyWindow);
	connect(ui.actionCGAL_based_Poisson, &QAction::hover, this, [this]() { { this->hovered_surface_item = 1; }});
	connect(ui.actionVoronoi_based, &QAction::triggered, this, &GraphicInterface::openSurfacePropertyWindow);
	connect(ui.actionVoronoi_based, &QAction::hover, this, [this]() { { this->hovered_surface_item = 2; }});
	connect(ui.actionPartial_diff, &QAction::triggered, this, &GraphicInterface::openSurfacePropertyWindow);
	connect(ui.actionPartial_diff, &QAction::hover, this, [this]() { { this->hovered_surface_item = 3; }});
}

void GraphicInterface::openOutlierFilterPropertyWindow()
{
	dialog = new QDialog(this);
	ui_out_filter.setupUi(dialog);

	connect(ui_out_filter.ok, SIGNAL(clicked()), this, SLOT(acceptNewOutlierFilter()));
	connect(ui_out_filter.cancel, SIGNAL(clicked()), this, SLOT(cancelDialog()));

	switch (hovered_filter_item)
	{
		case 1:
		{
			ui_out_filter.out->raise();
			ui_out_filter.name->setText("CGAL Outlier removal");
			break;
		}
		case 2:
		{
			ui_out_filter.out_plane->raise();
			ui_out_filter.name_pl->setText("Plane based outlier removal");
			break;
		}
	}

	dialog->show();
}

void GraphicInterface::acceptNewOutlierFilter()
{
	QString name;

	switch (hovered_filter_item)
	{
		case 1:
		{
			int kn = ui_out_filter.kn->value();
			double outl_lim = ui_out_filter.outl_lim->value();
			if (!ui.glrender->addNewOutlierFilter(1, kn, outl_lim))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			name = ui_out_filter.name->text();

			widget = new QWidget();
			ui_filters.push_back(new Ui::sfcWidgetF());
			ui_filters.back()->setupUi(widget);
			ui_filters.back()->outlier->raise();
			ui_filters.back()->kn->setValue(kn);

			connect(ui_filters.back()->delete_out, SIGNAL(clicked()), this, SLOT(deleteFilter()));
			connect(ui_filters.back()->kn, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), this, [this](const int& val) {this->ui.glrender->SetFilterProperty(active_filter_item, val, this->ui_filters[active_filter_item]->outl_limit->value()); });
			connect(ui_filters.back()->outl_limit, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val) {this->ui.glrender->SetFilterProperty(active_filter_item, this->ui_filters[active_filter_item]->kn->value(), val); });
		
			break;
		}
		case 2:
		{
			float px = ui_out_filter.px->value();
			float py = ui_out_filter.py->value();
			float pz = ui_out_filter.pz->value();
			float nx = ui_out_filter.nx->value();
			float ny = ui_out_filter.ny->value();
			float nz = ui_out_filter.nz->value();

			if (!ui.glrender->addNewOutlierFilter(2, px, py, pz, nx, ny, nz))
			{	
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			name = ui_out_filter.name_pl->text();

			widget = new QWidget();
			ui_filters.push_back(new Ui::sfcWidgetF());
			ui_filters.back()->setupUi(widget);
			ui_filters.back()->outlierPlane->raise();

			ui_filters.back()->px->setValue(px);
			ui_filters.back()->py->setValue(py);
			ui_filters.back()->pz->setValue(pz);
			ui_filters.back()->nx->setValue(nx);
			ui_filters.back()->ny->setValue(ny);
			ui_filters.back()->nz->setValue(nz);

			connect(ui_filters.back()->delete_pl, SIGNAL(clicked()), this, SLOT(deleteFilter()));
			connect(ui_filters.back()->px, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val) 
			{this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->px->value(), ui_filters[active_filter_item]->py->value(), ui_filters[active_filter_item]->pz->value(),
																	  ui_filters[active_filter_item]->nx->value(), ui_filters[active_filter_item]->ny->value(), ui_filters[active_filter_item]->nz->value()); });
			connect(ui_filters.back()->py, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->px->value(), ui_filters[active_filter_item]->py->value(), ui_filters[active_filter_item]->pz->value(),
																	  ui_filters[active_filter_item]->nx->value(), ui_filters[active_filter_item]->ny->value(), ui_filters[active_filter_item]->nz->value()); });
			connect(ui_filters.back()->pz, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->px->value(), ui_filters[active_filter_item]->py->value(), ui_filters[active_filter_item]->pz->value(),
																	  ui_filters[active_filter_item]->nx->value(), ui_filters[active_filter_item]->ny->value(), ui_filters[active_filter_item]->nz->value()); });
			connect(ui_filters.back()->nx, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->px->value(), ui_filters[active_filter_item]->py->value(), ui_filters[active_filter_item]->pz->value(),
																	  ui_filters[active_filter_item]->nx->value(), ui_filters[active_filter_item]->ny->value(), ui_filters[active_filter_item]->nz->value()); });
			connect(ui_filters.back()->ny, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->px->value(), ui_filters[active_filter_item]->py->value(), ui_filters[active_filter_item]->pz->value(),
																	  ui_filters[active_filter_item]->nx->value(), ui_filters[active_filter_item]->ny->value(), ui_filters[active_filter_item]->nz->value()); });
			connect(ui_filters.back()->nz, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->px->value(), ui_filters[active_filter_item]->py->value(), ui_filters[active_filter_item]->pz->value(),
																	  ui_filters[active_filter_item]->nx->value(), ui_filters[active_filter_item]->ny->value(), ui_filters[active_filter_item]->nz->value()); });

			break;
		}
	}

	ui.filters->addItem(widget,name);
	active_filter_item = ui_filters.size() - 1;
	ui.filters->setCurrentIndex(active_filter_item);
	//ui_filters.back()->scrollArea->setWidgetResizable(true);

	dialog->close();
	delete dialog;
}

void GraphicInterface::openSimplifierPropertyWindow()
{
	dialog = new QDialog(this);
	ui_simpl_filter.setupUi(dialog);

	connect(ui_simpl_filter.ok, SIGNAL(clicked()), this, SLOT(acceptNewSimplFilter()));
	connect(ui_simpl_filter.cancel, SIGNAL(clicked()), this, SLOT(cancelDialog()));

	switch (hovered_filter_item)
	{
		case 4:
		{
			ui_simpl_filter.grid->raise();
			ui_simpl_filter.name->setText("CGAL Grid simplifier");
			break;
		}
		case 3:
		{
			ui_simpl_filter.grid->raise();
			ui_simpl_filter.name->setText("Grid simplifier by me");
			break;
		}
		case 5:
		{
			ui_simpl_filter.hor->raise();
			ui_simpl_filter.name_h->setText("CGAL Hieararchical simplifier");
			break;
		}
		case 6:
		{
			ui_simpl_filter.wlop_name->setText("WLOP simplifier");
			ui_simpl_filter.wlop->raise();
			break;
		}
	}

	dialog->show();
}

void GraphicInterface::acceptNewSimplFilter()
{
	QString name;

	switch (hovered_filter_item)
	{
		case 4:
		{
			float cell_s = ui_simpl_filter.c_size->value();

			if (!ui.glrender->addNewSimplifierFilter(1, cell_s))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			name = ui_simpl_filter.name->text();

			widget = new QWidget();
			ui_filters.push_back(new Ui::sfcWidgetF());
			ui_filters.back()->setupUi(widget);
			ui_filters.back()->gridCgal->raise();

			ui_filters.back()->c_size_cg->setValue(cell_s);

			connect(ui_filters.back()->delete_gcg, SIGNAL(clicked()), this, SLOT(deleteFilter()));
			connect(ui_filters.back()->c_size_cg, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{ this->ui.glrender->SetFilterProperty(active_filter_item, val);});

			break;
		}
		case 3:
		{
			float cell_s = ui_simpl_filter.c_size->value();

			if (!ui.glrender->addNewSimplifierFilter(2, cell_s))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			name = ui_simpl_filter.name->text();

			widget = new QWidget();
			ui_filters.push_back(new Ui::sfcWidgetF());
			ui_filters.back()->setupUi(widget);
			ui_filters.back()->gridBme->raise();
			ui_filters.back()->c_size->setValue(cell_s);

			connect(ui_filters.back()->delete_gbm, SIGNAL(clicked()), this, SLOT(deleteFilter()));
			connect(ui_filters.back()->c_size, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{ this->ui.glrender->SetFilterProperty(active_filter_item, val); });

			break;
		}
		case 5:
		{
			float mcs = ui_simpl_filter.mcs->value();
			float msv = ui_simpl_filter.msv->value();

			if (!ui.glrender->addNewSimplifierFilter(3, mcs, msv))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			name = ui_simpl_filter.name_h->text();

			widget = new QWidget();
			ui_filters.push_back(new Ui::sfcWidgetF());
			ui_filters.back()->setupUi(widget);
			ui_filters.back()->horCgal->raise();
			ui_filters.back()->mcs->setValue(mcs);
			ui_filters.back()->msv->setValue(msv);

			connect(ui_filters.back()->delete_hor, SIGNAL(clicked()), this, SLOT(deleteFilter()));
			connect(ui_filters.back()->mcs, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{ this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->mcs->value(), ui_filters[active_filter_item]->msv->value()); });
			connect(ui_filters.back()->msv, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{ this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->mcs->value(), ui_filters[active_filter_item]->msv->value()); });


			break;
		}
		case 6:
		{
			float rp = ui_simpl_filter.rp->value();
			float nr = ui_simpl_filter.nr->value();

			if (!ui.glrender->addNewSimplifierFilter(4, rp, nr))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			name = ui_simpl_filter.wlop_name->text();

			widget = new QWidget();
			ui_filters.push_back(new Ui::sfcWidgetF());
			ui_filters.back()->setupUi(widget);
			ui_filters.back()->wlopCgal->raise();
			ui_filters.back()->rp->setValue(rp);
			ui_filters.back()->nr->setValue(nr);

			connect(ui_filters.back()->delete_wlop, SIGNAL(clicked()), this, SLOT(deleteFilter()));
			connect(ui_filters.back()->rp, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{ this->ui.glrender->SetFilterProperty(active_filter_item,  val, ui_filters[active_filter_item]->nr->value()); });
			connect(ui_filters.back()->nr, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), this, [this](const double& val)
			{ this->ui.glrender->SetFilterProperty(active_filter_item, ui_filters[active_filter_item]->rp->value(), val); });

			break;
		}

	}
	ui.filters->addItem(widget, name);
	active_filter_item = ui_filters.size() - 1;
	ui.filters->setCurrentIndex(active_filter_item);

	dialog->close();
	delete dialog;
}

void GraphicInterface::openSmootherFilterPropertyWindow()
{
	dialog = new QDialog(this);
	ui_smo_filter.setupUi(dialog);

	connect(ui_smo_filter.ok, SIGNAL(clicked()), this, SLOT(acceptNewSmothFilter()));
	connect(ui_smo_filter.cancel, SIGNAL(clicked()), this, SLOT(cancelDialog()));

	switch (hovered_filter_item)
	{
		case 7:
		{
			ui_smo_filter.name->setText("Jet smoother");
			ui_smo_filter.jet->raise();
			break;
		}
	}

	dialog->show();
}

void GraphicInterface::acceptNewSmothFilter()
{
	QString name;

	switch (hovered_filter_item)
	{
		case 7:
		{
			int kn = ui_smo_filter.kn_jet->value();

			if (!ui.glrender->addNewSmootherFilter(1, kn))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}
			name = ui_smo_filter.name->text();

			widget = new QWidget();
			ui_filters.push_back(new Ui::sfcWidgetF());
			ui_filters.back()->setupUi(widget);
			ui_filters.back()->jetCgal->raise();
			ui_filters.back()->kn_jet->setValue(kn);

			connect(ui_filters.back()->delete_jet, SIGNAL(clicked()), this, SLOT(deleteFilter()));
			connect(ui_filters.back()->kn_jet, static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged), this, [this](const int& val) {this->ui.glrender->SetFilterProperty(active_filter_item, val); });

			break;
		}
	}
	ui.filters->addItem(widget, name);
	active_filter_item = ui_filters.size() - 1;
	ui.filters->setCurrentIndex(active_filter_item);

	dialog->close();
	delete dialog;
}

void GraphicInterface::openSurfacePropertyWindow()
{
	dialog = new QDialog(this);
	ui_surf.setupUi(dialog);

	connect(ui_surf.ok, SIGNAL(clicked()), this, SLOT(acceptNewSurf()));
	connect(ui_surf.cancel, SIGNAL(clicked()), this, SLOT(cancelDialog()));

	switch (hovered_surface_item)
	{
		case 1:
		{
			ui_surf.poisson->raise();
			break;
		}
		case 2:
		{
			ui_surf.powercrust->raise();
			break;
		}
		default:
			delete dialog;
			return;
	}

	dialog->show();
}

void GraphicInterface::acceptNewSurf()
{
	switch (hovered_surface_item)
	{
		case 1:
		{
			float angle = ui_surf.angle->value();
			float max_tri = ui_surf.trisize->value();
			float surf_appr_err = ui_surf.apprerr->value();

			if (!ui.glrender->addNewSurfaceReconstructor(0, angle, max_tri, surf_appr_err))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			ui.pois_angle->setValue(angle);
			ui.pois_maxtri->setValue(max_tri);
			ui.pois_surfapprerr->setValue(surf_appr_err);

			ui.powercrust->setVisible(false);
			ui.poisson->setVisible(true);
			ui.poisson->raise();

			break;
		}

		case 2:
		{
			if (!ui.glrender->addNewSurfaceReconstructor(1))
			{
				QMessageBox::warning(this, "Oops", "Something happened");
				return;
			}

			ui.power_c_draw->setCurrentIndex(0);
			ui.glrender->attemptChangePowerCrustDrawType(ui.power_c_draw->currentIndex());

			ui.powercrust->setVisible(true);
			ui.poisson->setVisible(false);
			ui.powercrust->raise();

			break;
		}
	}
	
	if (dialog == NULL) return;	//TODO REMOVE THIS LATER

	dialog->close();
	delete dialog;
}

void GraphicInterface::deleteFilter()
{
	ui.glrender->deleteFilter(active_filter_item);
	Ui::sfcWidgetF* tmp = ui_filters[active_filter_item];
	unsigned int act_tmp = active_filter_item;

	ui_filters.erase(ui_filters.begin() + active_filter_item, ui_filters.begin() + active_filter_item + 1);
	ui.filters->removeItem(active_filter_item);

	active_filter_item = act_tmp == 0 ? act_tmp : act_tmp -1;
	if (ui_filters.size()!=0) ui.filters->setCurrentIndex(active_filter_item);

	delete tmp;
}

//-----------------------------------
void GraphicInterface::setCalculationBar()
{
	ui.statusBar->addWidget(c_load);
	ui.statusBar->addWidget(_c_load);
	_c_load->display(0.00);

	ui.statusBar->addWidget(poi);
	ui.statusBar->addWidget(_poi);
	_poi->display(0);

	ui.statusBar->addWidget(all);
	ui.statusBar->addWidget(_all);
	_all->display(0.00);

	ui.statusBar->addWidget(s_cal);
	ui.statusBar->addWidget(_s_cal);
	_s_cal->display(0);

	ui.statusBar->addWidget(sort);
	ui.statusBar->addWidget(_sort);
	_sort->display(0.00);

	ui.statusBar->addWidget(join);
	ui.statusBar->addWidget(_join);
	_join->display(0.00);

	
}

void GraphicInterface::refreshCalculationBar()
{
	float all_t, join_t, sort_t, c_load, s_cal;
	int poi;
	ui.glrender->getCloudCalculationTime(c_load, s_cal, all_t, join_t, sort_t);
	ui.glrender->getCloudPoints(poi);

	_c_load->display(c_load);
	_all->display(all_t);
	_sort->display(sort_t);
	_join->display(join_t);
	_poi->display(poi);
	_s_cal->display(s_cal);
}

void GraphicInterface::exportSurface()
{
	if (ui.compButton->text() == "Stop")
	{
		ui.glrender->swapStartFlag();
	}

	QString file_name = QFileDialog::getOpenFileName(this, "Export a file");
}

void GraphicInterface::info()
{

}

void GraphicInterface::open()
{
	if (ui.compButton->text() == "Stop")
	{
		ui.glrender->swapStartFlag();
	}

	QString file_name = QFileDialog::getOpenFileName(this, "Open a file");

}

void GraphicInterface::save()
{
	if (ui.compButton->text() == "Stop")
	{
		ui.glrender->swapStartFlag();
	}

	QString file_name = QFileDialog::getOpenFileName(this, "Open a file");
}

void GraphicInterface::saveas()
{
	if (ui.compButton->text() == "Stop")
	{
		ui.glrender->swapStartFlag();
	}

	QString file_name = QFileDialog::getOpenFileName(this, "Open a file");
}

void GraphicInterface::setMaxPoi(const int& val)
{
	const int ma = pow(2, val);
	ui.max_poi->setText(QString::number( ma ));
	ui.glrender->setMeshAccuracy(ma);
}

void GraphicInterface::newConnection()
{
	if (ui.compButton->text() == "Stop")
		ui.glrender->swapStartFlag();
	ui.compButton->setText("Start");

	dialog = new QDialog(this);

	ui_new.setupUi(dialog);

	connect(ui_new.accept_but, SIGNAL(clicked()), this, SLOT(acceptNew()));
	connect(ui_new.cancel_but, SIGNAL(clicked()), this, SLOT(cancelDialog()));
	connect(ui_new.br_file, SIGNAL(clicked()), this, SLOT(openFileDial()));
	connect(ui_new.source, SIGNAL(currentIndexChanged(int)), this, SLOT(sourceChanged(const int &)));

	ui_new.fromrandom->setVisible(false);
	ui_new.random->setValidator(new QIntValidator(0, 99999999, this));

	dialog->show();
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

//mybe put to lambda as well
void GraphicInterface::swapStartFlag()
{
	ui.compButton->setText(ui.compButton->text() == "Stop" ? "Start" : "Stop");
	ui.glrender->swapStartFlag();
}

void GraphicInterface::acceptNew()
{
	if (ui_new.source->currentIndex() == 0 && !QFileInfo(ui_new.file->text()).exists())
	{
		QMessageBox::warning(this,"Wrong filename", "The file is not exist.\nPlease choose a valid file!");
		return;
	}
	else if (ui_new.source->currentIndex() == 1 && ui_new.random->text() == "")
	{
		QMessageBox::warning(this, "Wrong details", "Please give the number of points!");
		return;
	}

	dialog->close();

	ui.compButton->setText("Start");

	if (ui_new.source->currentIndex() == 0 && !ui.glrender->newFileInput( ui_new.file->text().toLatin1().data(), ui.max_poi->text().toInt(), ui_new.repeat->isChecked()))
	{
		ui.compButton->setVisible(false);
		QMessageBox::warning(this, "Wrong file", "The file is not suitable for our format!");
		delete dialog;
		//ui.compButton->setEnabled(false);
		return;
	}
	else if (ui_new.source->currentIndex() == 1 && !ui.glrender->newRandomInput( ui_new.random->text().toInt(), ui.max_poi->text().toInt(), ui_new.repeat->isChecked()))
	{
		ui.compButton->setVisible(false);
		QMessageBox::warning(this, "Wrong", "Something bad happened!");
		delete dialog;
		return;
	}

	delete dialog;
	//ui.compButton->setEnabled(true);
	ui.compButton->setVisible(true);

	while (ui.filters->count() != 0) ui.filters->removeItem(0);
	active_filter_item = 0;
}
//mybe put to lambda.. but too much existance
void GraphicInterface::cancelDialog()
{
	dialog->close();
	delete dialog;
}

void GraphicInterface::openFileDial()	//it seems that this have memory leak but actually that caching thing for later using.
{
	QString file_name = QFileDialog::getOpenFileName(this,"Open a file");

	ui_new.file->setText(file_name);
}

void GraphicInterface::sourceChanged(const int & id)
{
	if (id == 0)
	{
		ui_new.fromrandom->setVisible(false);
		ui_new.fromfile->setVisible(true);
	}
	else
	{
		ui_new.fromrandom->setVisible(true);
		ui_new.fromfile->setVisible(false);
	}
}

void GraphicInterface::noMoreData()
{
	ui.compButton->setText("Start");
	ui.glrender->swapStartFlag();
	ui.glrender->resetCloud();
}

/*void GraphicInterface::showNearest()
{
ui.nearest->setText(ui.nearest->text() == "Previous" ? "Nearest" : "Previous");
ui.glrender->showNearest();
}
void GraphicInterface::showKNearest()
{
ui.knearest->setText(ui.knearest->text() == "Previous" ? "KNearest" : "Previous");
//ui.glrender->changeComputation();
}*/
