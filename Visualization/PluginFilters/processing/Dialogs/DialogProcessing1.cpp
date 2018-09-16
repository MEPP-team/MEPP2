#include "DialogProcessing1.h"
#include "ui_DialogProcessing1.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogProcessing1::DialogProcessing1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogProcessing1)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogProcessing1::~DialogProcessing1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogProcessing1::setProcess(double x, double y, double z)
{
  ui->lineEdit_X->setText(QString::number(x));
  ui->lineEdit_Y->setText(QString::number(y));
  ui->lineEdit_Z->setText(QString::number(z));
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogProcessing1::getProcess(double &x, double &y, double &z)
{
  x = ui->lineEdit_X->text().toDouble();
  y = ui->lineEdit_Y->text().toDouble();
  z = ui->lineEdit_Z->text().toDouble();
}
////////////////////////////////////////////////////////////////////////////////
