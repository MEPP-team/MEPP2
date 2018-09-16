#include "DialogCurvature1.h"
#include "ui_DialogCurvature1.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogCurvature1::DialogCurvature1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogCurvature1)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogCurvature1::~DialogCurvature1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogCurvature1::setCurvature(bool geod, double radius)
{
  Qt::CheckState state =
      geod ? Qt::CheckState::Checked : Qt::CheckState::Unchecked;
  ui->checkBox_geod->setCheckState(state);
  ui->lineEdit_radius->setText(QString::number(radius));
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogCurvature1::getCurvature(bool &geod, double &radius)
{
  geod = ui->checkBox_geod->isChecked();
  radius = ui->lineEdit_radius->text().toDouble();
}
////////////////////////////////////////////////////////////////////////////////
