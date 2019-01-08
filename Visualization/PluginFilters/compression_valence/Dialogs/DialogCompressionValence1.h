#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <QDialog>
////////////////////////////////////////////////////////////////////////////////
namespace Ui {
class DialogCompressionValence1;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class DialogCompressionValence1 : public QDialog
{
  Q_OBJECT

public:
  explicit DialogCompressionValence1(QWidget *parent = 0);
  ~DialogCompressionValence1();

  void setCompressionValenceParams(const std::string &p3dFilePath,
                                   bool with_adaptative_quantization,
                                   int quantization_bits,
                                   int max_verticesdouble);
  void getCompressionValenceParams(std::string &p3dFilePath,
                                   bool &with_adaptative_quantization,
                                   int &quantization_bits,
                                   int &max_vertices);

private slots:
  void selectFilename();

private:
  Ui::DialogCompressionValence1 *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
