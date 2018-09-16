#include "DialogDecompressionValence1.h"
#include "ui_DialogDecompressionValence1.h"

#include <QFileDialog>

////////////////////////////////////////////////////////////////////////////////
FEVV::DialogDecompressionValence1::DialogDecompressionValence1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogDecompressionValence1)
{
  ui->setupUi(this);
  connect(
      ui->pushButtonBrowse, SIGNAL(clicked()), this, SLOT(selectFilename()));
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogDecompressionValence1::~DialogDecompressionValence1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogDecompressionValence1::setDecompressionValenceParams(
    const std::string &p3dFilePath,
    bool write_info,
    bool write_intermediate_meshes,
    bool display_intermediate_meshes)
{
  Qt::CheckState write_info_state =
      write_info ? Qt::CheckState::Checked : Qt::CheckState::Unchecked;
  Qt::CheckState write_intermediate_meshes_state =
      write_intermediate_meshes ? Qt::CheckState::Checked
                                : Qt::CheckState::Unchecked;
  Qt::CheckState display_intermediate_meshes_state =
      display_intermediate_meshes ? Qt::CheckState::Checked
                                  : Qt::CheckState::Unchecked;

  ui->p3dFilePath->setText(QString::fromStdString(p3dFilePath));
  ui->write_information->setCheckState(write_info_state);
  ui->write_intermediate_meshes->setCheckState(write_intermediate_meshes_state);
  ui->display_intermediate_meshes->setCheckState(
      display_intermediate_meshes_state);
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogDecompressionValence1::getDecompressionValenceParams(
    std::string &p3dFilePath,
    bool &write_info,
    bool &write_intermediate_meshes,
    bool &display_intermediate_meshes)
{
  p3dFilePath = ui->p3dFilePath->text().toStdString();
  // remove starting and finishing  quotes if present
  if(p3dFilePath.front() == '"' && p3dFilePath.back() == '"')
  {
    p3dFilePath.pop_back();
    p3dFilePath.erase(0, 1);
  }
#ifdef WIN32
  // replace '\' in windows path with '/'
  std::replace(p3dFilePath.begin(), p3dFilePath.end(), '\\', '/');
#endif

  write_info = ui->write_information->isChecked();
  write_intermediate_meshes = ui->write_intermediate_meshes->isChecked();
  display_intermediate_meshes = ui->display_intermediate_meshes->isChecked();
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogDecompressionValence1::selectFilename()
{
  QString filename =
      QFileDialog::getOpenFileName(this,
                                   tr("Open File"),
                                   ui->p3dFilePath->text(),
                                   tr("P3D files (*.p3d);;All files (*)"),
                                   0,
                                   QFileDialog::DontUseNativeDialog);

  if(!filename.isEmpty())
    ui->p3dFilePath->setText(filename);
}
////////////////////////////////////////////////////////////////////////////////
