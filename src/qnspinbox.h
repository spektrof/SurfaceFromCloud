#pragma once
/*
QSpinBox valuchanged(int) event called twice because of timerEvent function
If we override that then it solves the problem and call the right slot only once.
Ref: http://www.qtcentre.org/archive/index.php/t-43078.html
*/

#include <QtWidgets>

class QNSpinBox : public QSpinBox
{
public:
	QNSpinBox(QWidget* parent = NULL) : QSpinBox(parent) {}
protected:
	void timerEvent(QTimerEvent *event) override { event->accept(); }
};