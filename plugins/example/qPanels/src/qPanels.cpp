//##########################################################################
//#                                                                        #
//#                CLOUDCOMPARE PLUGIN: ExamplePlugin                      #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 of the License.               #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#                             COPYRIGHT: XXX                             #
//#                                                                        #
//##########################################################################

// First:
//	Replace all occurrences of 'ExamplePlugin' by your own plugin class name in this file.
//	This includes the resource path to info.json in the constructor.

// Second:
//	Open ExamplePlugin.qrc, change the "prefix" and the icon filename for your plugin.
//	Change the name of the file to <yourPluginName>.qrc

// Third:
//	Open the info.json file and fill in the information about the plugin.
//	 "type" should be one of: "Standard", "GL", or "I/O" (required)
//	 "name" is the name of the plugin (required)
//	 "icon" is the Qt resource path to the plugin's icon (from the .qrc file)
//	 "description" is used as a tootip if the plugin has actions and is displayed in the plugin dialog
//	 "authors", "maintainers", and "references" show up in the plugin dialog as well

#include <QtGui>

#include "qPanels.h"

#include "ActionA.h"
#include "fParamsDlg.h"
#include "facade.h"
#include "tools.h"
#include <vector>
// Default constructor:
//	- pass the Qt resource path to the info.json file (from <yourPluginName>.qrc file) 
//  - constructor should mainly be used to initialize actions and other members
qPanels::qPanels ( QObject *parent )
	: QObject( parent )
	, ccStdPluginInterface( ":/CC/plugin/qPanels/info.json" )
	, m_action( nullptr )
{
}

// This method should enable or disable your plugin actions
// depending on the currently selected entities ('selectedEntities').
void qPanels::onNewSelection( const ccHObject::Container &selectedEntities )
{
	if ( m_action == nullptr )
	{
		return;
	}
	
	int cloud = 0;
	int polys = 0;
	for ( ccHObject *object : selectedEntities ) {
		switch( object->getClassID()) {
			case CC_TYPES::POINT_CLOUD:
				cloud++;
				break;
			case CC_TYPES::POLY_LINE:
				polys++;
				break;
			default:
				m_action->setEnabled( false );
				return;
		}
	}

	m_action->setEnabled( true );
}

// This method returns all the 'actions' your plugin can perform.
// getActions() will be called only once, when plugin is loaded.
QList<QAction *> qPanels::getActions()
{
	// default action (if it has not been already created, this is the moment to do it)
	if ( !m_action )
	{
		// Here we use the default plugin name, description, and icon,
		// but each action should have its own.
		m_action = new QAction( getName(), this );
		m_action->setToolTip( getDescription() );
		m_action->setIcon( getIcon() );
		
		// Connect appropriate signal
		connect( m_action, &QAction::triggered, this, [this]()
		{
			doAction( m_app );
		});
	}

	return { m_action };
}


void qPanels::doAction( ccMainAppInterface *appInterface )
{
	if ( appInterface == nullptr )
	{
		// The application interface should have already been initialized when the plugin is loaded
		Q_ASSERT( false );
		
		return;
	}
	
	fParamsDlg dlg;
	if (!dlg.exec()) {
		appInterface->dispToConsole( "[Panels] Returned false from exec!", ccMainAppInterface::STD_CONSOLE_MESSAGE );
		return;
	}
	appInterface->dispToConsole( "[Panels] Returned true from exec!", ccMainAppInterface::STD_CONSOLE_MESSAGE );

	QString mesh_file = dlg.getMeshPath();
	QString frames_file = dlg.getFramesPath();
	QString sa_file = dlg.getSupportingAreasPath();
	double width = 21.741828;
	double height = 22.034304;
	double min_width = 0.5;
	double min_height = 0.5;
	double max_width = 13.5;
	double max_height = 5.5;
	double delta = 0.01;

	FacadeSolver solver = FacadeSolver();
	vector<Frame> frame;
	vector<SupportingArea> v;
	vector<SupportingArea> h;

	Mesh mesh;

	mesh.loadFromFile(mesh_file.toStdString());
	loadFramesFromFile(frame, frames_file.toStdString());
	loadSupportingAreasFromFile(v, h, sa_file.toStdString());

	solver.init(width, height, min_height, min_width, max_height, max_width, delta, frame, v, h, mesh);
	solver.solve();
	solver.printSolution();
}