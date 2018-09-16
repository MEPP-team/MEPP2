#include "DialogJnd1.h"
#include "ui_DialogJnd1.h"
#include <cmath> // M_PI
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogJnd1::DialogJnd1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogJnd1)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogJnd1::~DialogJnd1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogJnd1::getProcess(int &screenH,
                             int &screenW,
                             double &screenS,
                             int &sceneH,
                             double &sceneFov,
                             double &userD,
                             int &nbLights,
                             float &log_disp,
                             bool &use_log,
                             bool &force_jnd)
{
  screenH = ui->spinBox_screenH->value();
  screenW = ui->spinBox_screenW->value();
  screenS = ui->doubleSpinBox_screenS->value();
  sceneH = ui->spinBox_sceneH->value();
  sceneFov = M_PI * ui->doubleSpinBox_sceneFOV->value();
  userD = ui->doubleSpinBox_userDist->value();
  nbLights = ui->spinBox_nbLights->value();
  use_log = ui->checkBox_log->isChecked();
  force_jnd = ui->checkBox_force->isChecked();
  log_disp = ui->doubleSpinBox_logD->value();
}
////////////////////////////////////////////////////////////////////////////////
