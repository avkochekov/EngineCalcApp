#include "kinematic.h"
#include "ui_kinematic.h"

#include "math.h"

Kinematic::Kinematic(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Kinematic)
{
    ui->setupUi(this);
}

Kinematic::~Kinematic()
{
    delete ui;
}

void Kinematic::on_pushButton_clicked()
{

}
