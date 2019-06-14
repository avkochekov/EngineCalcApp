#ifndef KINEMATIC_H
#define KINEMATIC_H

#include <QWidget>

namespace Ui {
class Kinematic;
}

class Kinematic : public QWidget
{
    Q_OBJECT
    
public:
    explicit Kinematic(QWidget *parent = 0);
    ~Kinematic();
    
private slots:
    void on_pushButton_clicked();

private:
    Ui::Kinematic *ui;
};

#endif // KINEMATIC_H
