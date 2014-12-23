#ifndef MANDELBROTWINDOW_H
#define MANDELBROTWINDOW_H

#include <QMainWindow>

namespace Ui {
class MandelbrotWindow;
}

class MandelbrotWindow : public QMainWindow
{
	Q_OBJECT

public:
	explicit MandelbrotWindow(QWidget *parent = 0);
	~MandelbrotWindow();

private:
	void paintEvent(QPaintEvent *);
	Ui::MandelbrotWindow *ui;
};

#endif // MANDELBROTWINDOW_H
