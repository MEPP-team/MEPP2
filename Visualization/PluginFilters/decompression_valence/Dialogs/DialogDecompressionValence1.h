#ifndef DialogDecompressionValence1_H
#define DialogDecompressionValence1_H
////////////////////////////////////////////////////////////////////////////////
#include <QDialog>
////////////////////////////////////////////////////////////////////////////////
namespace Ui {
class DialogDecompressionValence1;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class DialogDecompressionValence1 : public QDialog
{
  Q_OBJECT

public:
  explicit DialogDecompressionValence1(QWidget *parent = 0);
  ~DialogDecompressionValence1();

  void setDecompressionValenceParams(const std::string &p3dFilePath,
                                     bool write_info,
                                     bool write_intermediate_meshes,
                                     bool display_intermediate_meshes);
  void getDecompressionValenceParams(std::string &p3dFilePath,
                                     bool &write_info,
                                     bool &write_intermediate_meshes,
                                     bool &display_intermediate_meshes);

private slots:
  void selectFilename();

private:
  Ui::DialogDecompressionValence1 *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
#endif // DialogDecompressionValence1_H
