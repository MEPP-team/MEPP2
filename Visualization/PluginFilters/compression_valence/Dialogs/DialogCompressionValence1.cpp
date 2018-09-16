#include "DialogCompressionValence1.h"
#include "ui_DialogCompressionValence1.h"

#include <QFileDialog>

////////////////////////////////////////////////////////////////////////////////
FEVV::DialogCompressionValence1::DialogCompressionValence1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogCompressionValence1)
{
  ui->setupUi(this);
  connect(
      ui->pushButtonBrowse, SIGNAL(clicked()), this, SLOT(selectFilename()));
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogCompressionValence1::~DialogCompressionValence1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogCompressionValence1::setCompressionValenceParams(
    const std::string &p3dFilePath,
    bool with_adaptative_quantization,
    int quantization_bits,
    int max_vertices)
{
  Qt::CheckState state = with_adaptative_quantization
                             ? Qt::CheckState::Checked
                             : Qt::CheckState::Unchecked;

  ui->p3dFilePath->setText(QString::fromStdString(p3dFilePath));
  ui->with_adaptative_quantization->setCheckState(state);
  ui->quantization_bits->setText(QString::number(quantization_bits));
  ui->max_vertices->setText(QString::number(max_vertices));
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogCompressionValence1::getCompressionValenceParams(
    std::string &p3dFilePath,
    bool &with_adaptative_quantization,
    int &quantization_bits,
    int &max_vertices)
{
  p3dFilePath = ui->p3dFilePath->text().toStdString();
  with_adaptative_quantization = ui->with_adaptative_quantization->isChecked();
  quantization_bits = ui->quantization_bits->text().toUInt();
  max_vertices = ui->max_vertices->text().toUInt();
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogCompressionValence1::selectFilename()
{
  QString filename =
      QFileDialog::getSaveFileName(this,
                                   tr("Open File"),
                                   ui->p3dFilePath->text(),
                                   tr("P3D files (*.p3d);;All files (*)"),
                                   0,
                                   QFileDialog::DontUseNativeDialog);

  if(!filename.isEmpty())
    ui->p3dFilePath->setText(filename);
}
////////////////////////////////////////////////////////////////////////////////
