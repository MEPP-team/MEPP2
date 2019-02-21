#include "helloworld_dialog.h"
#include "ui_helloworld_dialog.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::HelloworldDialog::HelloworldDialog(QWidget *parent)
    : QDialog(parent), ui(new Ui::HelloworldDialog)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::HelloworldDialog::~HelloworldDialog() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::HelloworldDialog::setParameters(double x, double y, double z)
{
  ui->lineEdit_X->setText(QString::number(x));
  ui->lineEdit_Y->setText(QString::number(y));
  ui->lineEdit_Z->setText(QString::number(z));
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::HelloworldDialog::getParameters(double &x, double &y, double &z)
{
  x = ui->lineEdit_X->text().toDouble();
  y = ui->lineEdit_Y->text().toDouble();
  z = ui->lineEdit_Z->text().toDouble();
}
////////////////////////////////////////////////////////////////////////////////
