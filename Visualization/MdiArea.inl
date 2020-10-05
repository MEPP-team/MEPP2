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
//#include "mdiarea.hxx"
//#include "mainwindow.hxx"

#include <QMimeData>

#include <iostream>

#include "Visualization/SimpleWindow.h"

inline FEVV::MdiArea::MdiArea(QWidget *parent) : QMdiArea(parent)
{
	bType = bNone;
}

#if(0)
inline void
FEVV::MdiArea::paintEvent(QPaintEvent *paintEvent)
{
	QMdiArea::paintEvent(paintEvent);

	// and now only paint your image here
	QPainter painter(viewport());
 
	painter.fillRect(paintEvent->rect(), QColor(23, 74, 124));
	//painter.drawImage(paintEvent->rect()/*QPoint(0, 0)*/, QImage("./mepp_background.bmp"));
 
	painter.end();
}
#endif

inline void
FEVV::MdiArea::dragEnterEvent(QDragEnterEvent *event)
{
    std::cout << "dragEnterEvent called" << std::endl;
	if (event->mimeData()->hasUrls())	//hasFormat("text/uri-list")
	{
		event->acceptProposedAction();
        std::cout << "acceptProposedAction" << std::endl;

// ---------------------------------------------------
#ifdef __linux__
	bType=bLeft;
	if (event->mouseButtons() & Qt::RightButton)
		bType=bRight;
#endif
// ---------------------------------------------------
	}
}

inline void
FEVV::MdiArea::dropEvent(QDropEvent *event)
{
    std::cout << "dropEvent called" << std::endl;
	int res = 0;
	//Viewer *viewer = NULL;
	//QList<QUrl> urls = event->mimeData()->urls();

	// test
	if (m_mw->activeMdiChild() == 0)
		std::cout << "No MdiChild" << std::endl;
	else
		std::cout << "MdiChild" << std::endl;
	// test

// ---------------------------------------------------
#ifdef __APPLE__
	bType=bLeft;
	if (event->keyboardModifiers() & Qt::MetaModifier)
		bType=bRight;
#endif
#ifdef _MSC_VER
	bType=bLeft;
	if (event->mouseButtons() & Qt::RightButton)
		bType=bRight;
#endif
// ---------------------------------------------------

#if(0)
	for (int i=0; i<urls.size(); i++)
	{
		QFileInfo fi(urls[i].toLocalFile());
		QString ext = fi.suffix();

		if ( (ext.toLower()=="off") || (ext.toLower()=="obj") || (ext.toLower()=="smf") || (ext.toLower()=="ply") || (ext.toLower()=="x3d")
#ifdef WITH_ASSIMP
				|| (ext.toLower()=="dae") || (ext.toLower()=="3ds") || (ext.toLower()=="lwo")
#endif
				)
		{
			if (m_mw->activeMdiChild() == 0)
				res = m_mw->loadFile(urls[i].toLocalFile(), Normal, NULL);
			else
			{
				viewer = qobject_cast<Viewer *>(m_mw->activeMdiChild());

				if (bType == bLeft)
					res = m_mw->loadFile(urls[i].toLocalFile(), Normal, NULL);
				else if (bType == bRight)
				{
	#ifdef __linux__
					if (event->keyboardModifiers() & Qt::MetaModifier)
	#else
					if (event->keyboardModifiers() & Qt::AltModifier)
	#endif
						res = m_mw->addFile(viewer, urls[i].toLocalFile(), Time, NULL);
					else
						res = m_mw->addFile(viewer, urls[i].toLocalFile(), Space, NULL);
				}
			}

			if (res)
				break;

			m_mw->writeSettings();
			m_mw->readSettings();
		}
	}
	if (viewer)
		viewer->recreateListsAndUpdateGL();
#endif

	event->acceptProposedAction();
}