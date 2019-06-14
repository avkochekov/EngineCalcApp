#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <QVector>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:

    void on__btn_CalcKinematic_clicked();

    void on__btn_CalcTemp_clicked();

    void on__btn_CalcEng_clicked();

    void on__btn_CalcGroup_1_clicked();

    void on__btn_CalcGroup0_Apply_clicked();

    void on__btn_CalcGroup_2_clicked();

    void on__btn_CalcGroup_3_clicked();

    void on__btn_CalcDinam_clicked();

//    void SetTableVision(int ccount, int rcount, QStringList hlabels);

//    void PlotGraphic(QCustomPlot<int> plot, QVector<double> xData, QVector<double> yData, QString xLabel, QString yLabel, int yMin, int yMax);

//    void on_radioButton_7_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
