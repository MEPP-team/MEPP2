#ifndef DialogJnd1_H
#define DialogJnd1_H
////////////////////////////////////////////////////////////////////////////////
#include <QDialog>
////////////////////////////////////////////////////////////////////////////////
namespace Ui {
class DialogJnd1;
}
////////////////////////////////////////////////////////////////////////////////
namespace FEVV {

class DialogJnd1 : public QDialog
{
  Q_OBJECT

public:
  explicit DialogJnd1(QWidget *parent = 0);
  ~DialogJnd1();

  void getProcess(int &screenH,
                  int &screenW,
                  double &screenS,
                  int &sceneH,
                  double &sceneFov,
                  double &userD,
                  int &nbLights,
                  float &log_disp,
                  bool &use_log,
                  bool &force_jnd);

private:
  Ui::DialogJnd1 *ui;
};

} // namespace FEVV

////////////////////////////////////////////////////////////////////////////////
#endif // DialogJnd1_H
