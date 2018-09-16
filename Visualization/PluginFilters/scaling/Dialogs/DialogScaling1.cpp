#include "DialogScaling1.h"
#include "ui_DialogScaling1.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogScaling1::DialogScaling1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogScaling1)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogScaling1::~DialogScaling1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogScaling1::setScale(double x, double y, double z)
{
  ui->lineEdit_X->setText(QString::number(x));
  ui->lineEdit_Y->setText(QString::number(y));
  ui->lineEdit_Z->setText(QString::number(z));
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogScaling1::getScale(double &x, double &y, double &z)
{
  x = ui->lineEdit_X->text().toDouble();
  y = ui->lineEdit_Y->text().toDouble();
  z = ui->lineEdit_Z->text().toDouble();
}
////////////////////////////////////////////////////////////////////////////////
