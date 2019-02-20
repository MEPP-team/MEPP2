#include "DialogHelloworld1.h"
#include "ui_DialogHelloworld1.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogHelloworld1::DialogHelloworld1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogHelloworld1)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogHelloworld1::~DialogHelloworld1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogHelloworld1::setHelloworld(double x, double y, double z)
{
  ui->lineEdit_X->setText(QString::number(x));
  ui->lineEdit_Y->setText(QString::number(y));
  ui->lineEdit_Z->setText(QString::number(z));
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogHelloworld1::getHelloworld(double &x, double &y, double &z)
{
  x = ui->lineEdit_X->text().toDouble();
  y = ui->lineEdit_Y->text().toDouble();
  z = ui->lineEdit_Z->text().toDouble();
}
////////////////////////////////////////////////////////////////////////////////
