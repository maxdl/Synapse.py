#
#    Module      : wxAboutDialog.py
#    Date        : June 3, 2010
#    Description : an About dialog
#
#    Copyright 2010 Max Larsson <m.d.larsson@medisin.uio.no>
#
#    This file is part of Synapse.
#
#    Synapse is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Synapse is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Synapse.  If not, see <http://www.gnu.org/licenses/>.

import wx
import gui
import version

class wxAboutDialog(gui.AboutDialog):
    def __init__(self, parent):
        gui.AboutDialog.__init__(self, parent)
        self.TitleLabel.SetLabel(version.title)
        self.IconBitmap.SetBitmap(wx.Bitmap(version.icon, wx.BITMAP_TYPE_ANY))
        self.VersionLabel.SetLabel("Version %s" % version.version)
        self.LastModLabel.SetLabel("Last modified %s %s, %s." % version.date)
        self.CopyrightLabel.SetLabel("Copyright %s %s." % (version.date[2], 
                                                           version.author))
        self.LicenseLabel.SetLabel("Released under the terms of the GPLv3"
                                        " license.")
        self.EmailHyperlink.SetLabel("%s" % version.email)
        self.EmailHyperlink.SetURL("mailto://%s" % version.email)
        self.WebHyperlink.SetLabel("http://%s" % version.homepage)
        self.WebHyperlink.SetURL("http://%s" % version.homepage)
        self.SetIcon(wx.Icon(version.icon, wx.BITMAP_TYPE_ICO))              
        self.Fit()

    def OnClose( self, event ):
        self.Destroy()
    
    
