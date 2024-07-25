#include "MeshParamDialog.h"
#include <QApplication>
#include <QDesktopWidget>

#include "..\src\Common\CommonDefinitions.h"

MeshParamDialog::MeshParamDialog(QWidget* parent /* = 0 */)
	:QDialog(parent)
{
	initDialog();
}

MeshParamDialog::~MeshParamDialog()
{
}

void MeshParamDialog::initDialog()
{
	createWidget();
	createLayout();
}

void MeshParamDialog::createWidget()
{
	create_Basic_Operation_Information_Widget();
}

void MeshParamDialog::createLayout()
{
	tabWidget = new QTabWidget();
	tabWidget->addTab(view_BOI, "QP");

	QGridLayout *layout = new QGridLayout();
	layout->addWidget(tabWidget, 0, 0, 1, 1);
	setLayout(layout);
}

void MeshParamDialog::create_Basic_Operation_Information_Widget()
{
	load_ref_mesh_AM = new QPushButton("Load Ref Mesh");
	do_remehsing_AM = new QPushButton("Do Re meshing");
	target_edge_length_AM = new QLabel("Target Edge Length:");
	target_edge_length_line_AM = new QLineEdit("1.0");
	target_edge_length_line_AM->setValidator(new QDoubleValidator(0.0, 1000000.0, 10));
	sample_ratio_AM = new QLabel("Edge Length Ratio:");
	sample_ratio_line_AM = new QLineEdit("1.5");
	sample_ratio_line_AM->setValidator(new QDoubleValidator(1.0, 1000.0, 10));

	QGridLayout* LCOT_layout = new QGridLayout(); int LCOT_layout_index = 0;
	LCOT_layout->addWidget(load_ref_mesh_AM, LCOT_layout_index, 0, 1, 2); LCOT_layout_index += 1;
	LCOT_layout->addWidget(target_edge_length_AM, LCOT_layout_index, 0, 1, 1);
	LCOT_layout->addWidget(target_edge_length_line_AM, LCOT_layout_index, 1, 1, 1); LCOT_layout_index += 1;
	LCOT_layout->addWidget(sample_ratio_AM, LCOT_layout_index, 0, 1, 1);
	LCOT_layout->addWidget(sample_ratio_line_AM, LCOT_layout_index, 1, 1, 1); LCOT_layout_index += 1;
	LCOT_layout->addWidget(do_remehsing_AM, LCOT_layout_index, 0, 1, 2); LCOT_layout_index += 1;

	Basic_Operation_And_Information = new QWidget();
	Basic_Operation_And_Information->setLayout(LCOT_layout);

	view_BOI = new QScrollArea;
	view_BOI->setFocusPolicy(Qt::NoFocus);
	view_BOI->setFrameStyle(QFrame::NoFrame);
	view_BOI->setWidget(Basic_Operation_And_Information);
	view_BOI->setWidgetResizable(true);

	connect(load_ref_mesh_AM, SIGNAL(clicked()), SIGNAL(load_ref_mesh_AM_signal()));
	connect(do_remehsing_AM, SIGNAL(clicked()), SIGNAL(do_remehsing_AM_signal()));
}
