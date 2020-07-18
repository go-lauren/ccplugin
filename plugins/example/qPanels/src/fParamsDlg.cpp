#include "fParamsDlg.h"
#include "ui_fParamsDlg.h"

#include <QFileDialog>

fParamsDlg::fParamsDlg(QWidget *parent)
    : QDialog(parent)
    , ui(new Ui::fParamsDlg)
{
    ui->setupUi(this);
}

fParamsDlg::~fParamsDlg()
{
    delete ui;
}

QString fParamsDlg::browseFile()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open the file");
    return filename;
}

void fParamsDlg::on_browseMesh_clicked()
{
    QString filename = browseFile();
    ui->meshLabel->setText(filename);
}

void fParamsDlg::on_browseFrames_clicked()
{
    QString filename = browseFile();
    ui->framesLabel->setText(filename);
}

void fParamsDlg::on_browseSupportingAreas_clicked()
{
    QString filename = browseFile();
    ui->supportingAreasLabel->setText(filename);
}

QString fParamsDlg::getMeshPath()
{
    return ui->meshLabel->text();
}

QString fParamsDlg::getFramesPath()
{
    return ui->framesLabel->text();
}

QString fParamsDlg::getSupportingAreasPath()
{
    return ui->supportingAreasLabel->text();
}
