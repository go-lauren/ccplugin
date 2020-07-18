#ifndef FPARAMSDLG_H
#define FPARAMSDLG_H

#include <QDialog>

QT_BEGIN_NAMESPACE
namespace Ui { class fParamsDlg; }
QT_END_NAMESPACE

class fParamsDlg : public QDialog
{
    Q_OBJECT

public:
    fParamsDlg(QWidget *parent = nullptr);
    ~fParamsDlg();

    QString getMeshPath();
    QString getFramesPath();
    QString getSupportingAreasPath();

private slots:
    void on_browseMesh_clicked();

    void on_browseFrames_clicked();

    void on_browseSupportingAreas_clicked();

private:
    Ui::fParamsDlg *ui;

    QString browseFile();
};
#endif // FPARAMSDLG_H
