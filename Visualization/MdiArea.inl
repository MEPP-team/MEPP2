// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include <QMimeData>

#include <QUrl>

#include <iostream>

#include "Visualization/SimpleWindow.h"

inline FEVV::MdiArea::MdiArea(QWidget *parent) : QMdiArea(parent)
{
	bType = bNone;
}

inline void
FEVV::MdiArea::dragEnterEvent(QDragEnterEvent *event)
{
	//std::cout << "dragEnterEvent called" << std::endl;
	if (event->mimeData()->hasUrls()) // hasFormat("text/uri-list")
	{
// ---------------------------------------------------
#ifdef __linux__
		// new
		m_mw->activateWindow();
		m_mw->raise();
		// new

		bType=bLeft;
		if (event->mouseButtons() & Qt::RightButton)
			bType=bRight;
#endif
// ---------------------------------------------------

		//std::cout << "acceptProposedAction" << std::endl;
		event->acceptProposedAction();
	}
}

inline void
FEVV::MdiArea::dropEvent(QDropEvent *event)
{
	//std::cout << "dropEvent called" << std::endl;

// ---------------------------------------------------
#ifdef __APPLE__
	// new
	m_mw->activateWindow();
	m_mw->raise();
	// new

	bType=bLeft;
	if (event->keyboardModifiers() & Qt::MetaModifier)
		bType=bRight;
#endif
#ifdef _MSC_VER
	// new
	m_mw->activateWindow();
	m_mw->raise();
	// new

	bType=bLeft;
	if (event->mouseButtons() & Qt::RightButton)
		bType=bRight;
#endif
// ---------------------------------------------------

	QList<QUrl> urls;
	urls = event->mimeData()->urls();

	m_mw->drag_files.clear();

	for (int i=0; i<urls.size(); i++)
	{
		QFileInfo fi(urls[i].toLocalFile());
		QString ext = fi.suffix();

		if ( (ext.toLower()=="obj") || (ext.toLower()=="off") || (ext.toLower()=="coff") || (ext.toLower()=="ply") || (ext.toLower()=="msh")
#ifdef FEVV_USE_VTK
				|| (ext.toLower()=="vtk") || (ext.toLower()=="vtp") || (ext.toLower()=="vtu")
#endif
#ifdef FEVV_USE_FBX
				|| (ext.toLower()=="fbx")
#endif
// Point clouds
#ifdef FEVV_USE_CGAL
				|| (ext.toLower()=="xyz")// CGALPOINTSET
#endif 
#ifdef FEVV_USE_PCL
				|| (ext.toLower()=="pcd")// PCLPOINTCLOUD
#endif 
				)
		{
			QString sFile = urls[i].toLocalFile();
			m_mw->drag_files << sFile;

			//std::cout << "File: " << sFile.toStdString() << std::endl;
		}
	}

	if (! m_mw->drag_files.isEmpty())
	{
		m_mw->drag=true;
		m_mw->recent=false; // here only for 'protection'

#ifdef _MSC_VER
		if (bType==bRight)
			m_mw->shift_drag = true;
		else
			m_mw->shift_drag = false;
#else
		if (event->keyboardModifiers() & Qt::ShiftModifier)
			m_mw->shift_drag = true;
#endif
		if (event->keyboardModifiers() & Qt::AltModifier)
			m_mw->alt_drag = true;
		if (event->keyboardModifiers() & Qt::ControlModifier)
			m_mw->ctrl_drag = true;

		//std::cout << "acceptProposedAction" << std::endl;
		event->acceptProposedAction();

		m_mw->on_actionOpen_triggered();
	}
}