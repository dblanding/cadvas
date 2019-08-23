#!/usr/bin/env python
#
# cadvas.py     1/11/04
# A 2D CAD application using Python and the Gnome Canvas
# The latest  version of this file can be found at:
# http://members.localnet.com/~blanding/cadvas
#
# Author: Doug Blanding   <doug dot blanding at kodak dot com>
#
# CADvas is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# CADvas is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with CADvas; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

import os
import pickle
import math
import gtk
import gnome.canvas

def intersection(line1, line2):
    """Return intersection (x,y) of 2 lines expressed in (a,b,c) form."""
    a,b,c = line1
    d,e,f = line2
    i = b*f-c*e
    j = c*d-a*f
    k = a*e-b*d
    if k:
        return (i/k, j/k)
    else:
        return None

def cnvrt_2pts_to_coef(pt1, pt2):
    """Return (a,b,c) coefficients of line defined by 2 (x,y) pts."""
    x1, y1 = pt1
    x2, y2 = pt2
    a = y2 - y1
    b = x1 - x2
    c = x2*y1-x1*y2
    return (a, b, c)

def closer(p0, p1, p2):
    """Return closer of p1 or p2 to point p0."""
    d1 = (p1[0] - p0[0])**2 + (p1[1] - p0[1])**2
    d2 = (p2[0] - p0[0])**2 + (p2[1] - p0[1])**2
    if d1 < d2: return p1
    else: return p2

def midpoint(p1, p2, f=.5):
    """Return point part way (f=.5 by def) between points p1 and p2."""
    return (((p2[0]-p1[0])*f)+p1[0], ((p2[1]-p1[1])*f)+p1[1])

def arc_end_pts(param):
    """Return end pts (p1, p2) of arc specified as (x, y, r, a0, a1),
    where x,y = cntr pt, r=radius, and a0, a1 are starting and ending
    angles, expressed in radians."""
    x, y, r, a0, a1 = param
    x0 = x + r*math.cos(a0)
    y0 = y + r*math.sin(a0)
    x1 = x + r*math.cos(a1)
    y1 = y + r*math.sin(a1)
    return ((x0, y0), (x1, y1))

def proj_pt_on_line(line, pt):
    """Return point which is the projection of pt on line."""
    a, b, c = line
    x, y = pt
    denom = a**2 + b**2
    xp = (b**2*x - a*b*y -a*c)/denom
    yp = (a**2*y - a*b*x -b*c)/denom
    return (xp, yp)

def perp_line(line, pt):
    """Return coeff of newline thru pt and perpend to line."""
    a, b, c = line
    x, y = pt
    cnew = a*y - b*x
    return (b, -a, cnew)

def para_line(line, pt):
    """Return coeff of newline thru pt and parallel to line."""
    a, b, c = line
    x, y = pt
    cnew = -(a*x + b*y)
    return (a, b, cnew)

def para_lines(line, d):
    """Return 2 parallel lines straddling line, offset d."""
    a, b, c = line
    c1 = math.sqrt(a**2 + b**2)*d
    cline1 = (a, b, c + c1)
    cline2 = (a, b, c - c1)
    return (cline1, cline2)

def pt2pt_dist(pt1, pt2):
    """Return the distance between two points"""
    x, y = pt1
    u, v = pt2
    return math.sqrt((x-u)**2 + (y-v)**2)

def find_fillet_pts(r, commonpt, end1, end2):
    """Return ctr of fillet (radius r) and tangent pts for corner
    defined by a common pt, and two adjacent corner pts."""
    line1 = cnvrt_2pts_to_coef(commonpt, end1)
    line2 = cnvrt_2pts_to_coef(commonpt, end2)
    # find 'interior' clines
    cl1a, cl1b = para_lines(line1, r)
    p2a = proj_pt_on_line(cl1a, end2)
    p2b = proj_pt_on_line(cl1b, end2)
    da = pt2pt_dist(p2a, end2)
    db = pt2pt_dist(p2b, end2)
    if da <= db: cl1 = cl1a
    else: cl1 = cl1b
    cl2a, cl2b = para_lines(line2, r)
    p1a = proj_pt_on_line(cl2a, end1)
    p1b = proj_pt_on_line(cl2b, end1)
    da = pt2pt_dist(p1a, end1)
    db = pt2pt_dist(p1b, end1)
    if da <= db: cl2 = cl2a
    else: cl2 = cl2b
    pc = intersection(cl1, cl2)
    p1 = proj_pt_on_line(line1, pc)
    p2 = proj_pt_on_line(line2, pc)
    return (pc, p1, p2)

def cline_box_intrsctn(line, box):
    """Return tuple of pts where line intersects edges of box."""
    x0, y0, x1, y1 = box
    pts = []
    segments = [((x0, y0), (x1, y0)),
                ((x1, y0), (x1, y1)),
                ((x1, y1), (x0, y1)),
                ((x0, y1), (x0, y0))]
    for seg in segments:
        pt = intersection(line, cnvrt_2pts_to_coef(seg[0], seg[1]))
        if pt:
            if pt2pt_dist(pt, seg[0]) <= pt2pt_dist(seg[0], seg[1]) and \
               pt2pt_dist(pt, seg[1]) <= pt2pt_dist(seg[0], seg[1]):
                if pt not in pts:
                    pts.append(pt)
    return tuple(pts)

class cadWindow:
    def __init__(self):
        self.botm = None # bottommost group (child of root)
        self.cons = None # base group for construction lines
        self.tmpg = None # base group for new geometry
        self.geom = None # base group for geometry lines
        self.pnts = None # base group for intersection points
        self.dims = None # base group for dimensions
        self.canvas_size = (-4200, -4200, 4200, 4200)
        self.width = 800    # screen size
        self.height = 600   # screen size
        self.pt_size = 10.0
        self.viewscale = 1.0
        self.viewscale_prev = 1.0
        self.scroll_prev = (0, 0)
        self.linewidth = 2
        self.dimgap = 10
        self.dim_line_color = 'yellow'
        self.dim_text_color = 'white'
        self.hcl = []
        self.vcl = []
        self.acl = []
        self.tmp_geom_list = []
        self.hdims = []
        self.vdims = []
        self.cdict = {}
        self.pdict = {}
        self.ldict = {}
        self.cirdict = {}
        self.arcdict = {}
        self.dimdict = {}
        self.dim_parts_dict = {}
        self.float_stack = []   # float values
        self.pt_stack = []      # (x, y) values
        self.line_stack = []    # (pt1, pt2) values
        self.cline_stack = []   # (a, b, c) values
        self.item_stack = []    # for selected canvas widgets
        self.op = ''            # operation in progress
        self.op_stack = []
        self.flip_flag = True   # If True, +Y values go up
        self.pt_sel_group = ('line', 'polyline', 'rect', 'arc', 'circle',
                             'cline', 'hcline', 'vcline', 'hvcline',
                             'par_cl', 'perp_cl', 'lin_bsctr', 'hdim',
                             'vdim', 'dim_move', 'fillet', 'window')
        self.line_sel_group = ('par_cl', 'perp_cl', 'ang_bsctr')
        self.unit_dict = {'mm': 1.0, 'inches': 25.4, 'feet': 304.8}
        self.units = 'mm'
        self.unitscale = 1.0

    def set_units(self, units):
        if units in self.unit_dict.keys():
            self.units = units
            self.unitscale = self.unit_dict.get(units)
            self.units_display.set_text("Units: %s" % self.units)
            self.remake_dims()
            
    def printer(self):
        '''Send postscript drawing to printer. Scale to fit geometry.'''
        l, t, r, b = self.geom.get_bounds() # geom bounding box (world)
        if not r-l:
            return
        sf = 72 * 8.5/(r-l) # scale factor points/mm
        x = (r+l)/2 - 4.25*72/sf # corresponds to botm edge of paper
        y = (t+b)/2 - 5.50*72/sf # corresponds to left edge of paper
        nsf = sf * 7.5 / 8.5 # leave 1/2 inch margins on sides
        s = '%!\n'
        s = s + 'newpath\n'
        for line in self.ldict.values():
            p1, p2 = line
            s = s + '%f %f moveto\n' % ((p1[0]-x)*nsf+36, (p1[1]-y)*nsf)
            s = s + '%f %f lineto\n' % ((p2[0]-x)*nsf+36, (p2[1]-y)*nsf)
        for arc in self.arcdict.values():
            xc, yc, r, a0, a1 = arc
            x0 = r * math.cos(a0)
            y0 = r * math.sin(a0)
            a0 = a0 * 180 / math.pi # convert to degrees
            a1 = a1 * 180 / math.pi # convert to degrees
            s = s + '%f %f moveto\n' % ((xc-x)*nsf+36, (yc-y)*nsf)
            s = s + '%f %f rmoveto\n' % (x0*nsf, y0*nsf)
            s = s + '%f %f %f %f %f arc\n' % ((xc-x)*nsf+36, (yc-y)*nsf,
                                              r*nsf, a0, a1)
        for circle in self.cirdict.values():
            xc, yc, r = circle
            s = s + '%f %f moveto\n' % ((xc-x)*nsf+36, (yc-y)*nsf)
            s = s + '%f 0 rmoveto\n' % (r*nsf)  #TypeError w/out parens ??
            s = s + '%f %f %f 0 360 arc\n' % ((xc-x)*nsf+36,
                                              (yc-y)*nsf, r*nsf)        
        for part in self.dim_parts_dict.values():
            for item in part[1:]:
                x1, y1, x2, y2 = item.get_property('points')
                s = s + '%f %f moveto\n' % ((x1-x)*nsf+36, (y1-y)*nsf)
                s = s + '%f %f lineto\n' % ((x2-x)*nsf+36, (y2-y)*nsf)
        s = s + 'stroke\n'
        s = s + '/Times-Roman findfont\n'
        s = s + '15 scalefont\n'
        s = s + 'setfont\n'
        for part in self.dim_parts_dict.values():
            textitem = part[0]
            text = textitem.get_property('text')
            xc = textitem.get_property('x')
            yc = textitem.get_property('y')
            width = textitem.get_property('text_width')
            height = textitem.get_property('text_height')
            s = s + '%s %s moveto\n' % ((xc-x)*nsf+36, (yc-y)*nsf)
            s = s + '%s %s rmoveto\n' % ((-width/2)*nsf, (-height/2)*nsf)
            s = s + '(%s) show\n' % (text)
        s = s + 'showpage\n'
        
        # print to file
        f = open('printfoo.ps', 'w')
        f.write(s)
        f.close()
        '''
        # send to printer
        p = os.popen('lpr', 'w')
        p.write(s)
        p.close
        '''

    def save(self, file):
        """Save drawing contents to file."""
        save_dict = {"constr": self.cdict.values(),
                     "line": self.ldict.values(),
                     "circ": self.cirdict.values(),
                     "arc": self.arcdict.values(),
                     "hdim": [self.dimdict.get(key) for key in self.hdims],
                     "vdim": [self.dimdict.get(key) for key in self.vdims]}
        f = open(file, 'w')
        pickle.dump(save_dict, f)
        f.close

    def load(self, file):
        """Load drawing contents from file."""
        f = open(file, 'r')
        load_dict = pickle.load(f)
        f.close
        for line in load_dict.get("line"):
            self.line_gen(line)
        for circ in load_dict.get("circ"):
            self.circle_gen(circ)
        for arc in load_dict.get("arc"):
            self.arc_gen(arc)
        self.reparent_tmp_geom()
        for cline in load_dict.get("constr"):
            self.cline2(cline)
        self.remake_pnts()
        for hdim in load_dict.get("hdim"):
            p1, p2, p3 = hdim
            self.dim_gen(p1, p2, p3, 'h')
        for vdim in load_dict.get("vdim"):
            p1, p2, p3 = vdim
            self.dim_gen(p1, p2, p3, 'v')

    def mainquit(self, *args):
        """Delete all references to canvas items."""
        print 'cdict length was %s' % (len(self.cdict))
        print 'pdict length was %s' % (len(self.pdict))
        print 'ldict length was %s' % (len(self.ldict))
        print 'cirdict length was %s' % (len(self.cirdict))
        print 'arcdict length was %s' % (len(self.arcdict))
        print 'dimdict length was %s' % (len(self.dimdict))
        print 'dim_parts_dict length was %s' % (len(self.dim_parts_dict))
        print 'item_stack list length was %s' % (len(self.item_stack))
        print 'hcl list length was %s' % (len(self.hcl))
        print 'vcl list length was %s' % (len(self.vcl))
        print 'acl list length was %s' % (len(self.acl))
        print 'tmp_geom_list length was', len(self.tmp_geom_list)
        print 'hdims list length was %s' % (len(self.hdims))
        print 'vdims list length was %s' % (len(self.vdims))
        
        self.cdict.clear()
        self.pdict.clear()
        self.ldict.clear()
        self.cirdict.clear()
        self.arcdict.clear()
        self.dimdict.clear()
        self.dim_parts_dict.clear()
        self.item_stack = []
        self.hcl = []
        self.vcl = []
        self.acl = []
        self.tmp_geom_list = []
        self.hdims = []
        self.vdims = []
        self.dims = None
        self.pnts = None
        self.geom = None
        self.tmpg = None
        self.cons = None
        self.botm = None
        
    def zoom(self, widget, scale=1.0):
        """ if scale > 0: Change view scale to scale.
        
        elif scale == 0: scale to fit geometry.
        else: (scale < 0)show previous view.
        """
        if scale > 0:
            self.viewscale_prev = self.viewscale
            self.scroll_prev = self.acanvas.get_scroll_offsets()
            self.viewscale = self.viewscale * scale
            if self.viewscale > 15.0:
                self.viewscale = 15.0
            self.acanvas.set_pixels_per_unit(self.viewscale)
            self.remake_pnts()
            self.remake_dims()
        elif scale == 0:   # zoom to fit
            self.viewscale_prev = self.viewscale
            self.scroll_prev = self.acanvas.get_scroll_offsets()
            l, t, r, b = self.geom.get_bounds() # geom bounding box (world)
            cx, cy = ((l+r)/2, (t+b)/2) # location of cntr of bb
            sx, sy = (r-l, b-t) # size of bb
            if sx == 0 or sy == 0:
                return gtk.TRUE
            else:
                pass
            pad = 1.3   # make view a little bigger than bb
            vx, vy = (pad*sx, pad*sy)   # size of new view
            rsx = float(self.width)/vx     # required scale in x
            rsy = float(self.height)/vy    # required scale in y
            if rsx < rsy:
                self.viewscale = rsx
            else:
                self.viewscale = rsy
            if self.viewscale > 15.0:
                self.viewscale = 15.0
	    self.acanvas.set_pixels_per_unit(self.viewscale)
	    # map bb ctr to canvas ctr
	    if self.flip_flag:
                cxc, cyc = self.acanvas.w2c(cx, -cy)
            else:
                cxc, cyc = self.acanvas.w2c(cx, cy)
            self.acanvas.scroll_to(cxc-self.width/2, cyc-self.height/2)
	    self.remake_pnts()
	    self.remake_dims()
	else:   # display previous view
            # previous settings
            temp_scale = self.viewscale_prev
            temp_scroll = self.scroll_prev
            # store current settings
            self.viewscale_prev = self.viewscale
            self.scroll_prev = self.acanvas.get_scroll_offsets()
            # set new settings
            self.viewscale = temp_scale
            self.acanvas.set_pixels_per_unit(self.viewscale)
            self.acanvas.scroll_to(temp_scroll[0], temp_scroll[1])
            self.remake_pnts()
            self.remake_dims()

    def window(self, widget=None, data= None):
        """Set view window to box defined by 2 clicked points.

        If the selected window becomes too small, the
        self.acanvas.set_pixels_per_unit() call takes "forever",
        consumes all available memory and swap, and may result
        in a segmentation fault. It's just not healthy."""
        if self.op != 'window':
            self.op = 'window'
            self.pt_stack = []
            self.status_line.set_text(
                "Pick first corner of new view window.")
            return
        elif len(self.pt_stack) == 1:
            self.status_line.set_text(
                "Pick opposite corner of new view window.")
            return
        elif len(self.pt_stack) == 2:
            self.op = ''
            self.status_line.set_text('')
            x2, y2 = p2 = self.pt_stack.pop()
            x1, y1 = p1 = self.pt_stack.pop()
            xc, yc = (x1+x2)/2, (y1+y2)/2   # cntr of new view
            sx, sy = abs(x2-x1), abs(y2-y1) # size of new view
            if sx == 0 and sy == 0:
                return
            else:
                pass
            self.viewscale_prev = self.viewscale
            self.scroll_prev = self.acanvas.get_scroll_offsets()
            rsx = float(self.width)/sx     # required scale in x
            rsy = float(self.height)/sy    # required scale in y
            if rsx < rsy:
                self.viewscale = rsx
            else:
                self.viewscale = rsy
            if self.viewscale > 15.0:
                self.viewscale = 15.0
            # map new view ctr to canvas ctr
	    if self.flip_flag:
                cxc, cyc = self.acanvas.w2c(xc, -yc)
            else:
                cxc, cyc = self.acanvas.w2c(xc, yc)
            self.acanvas.scroll_to(cxc-self.width/2, cyc-self.height/2)
            self.acanvas.set_pixels_per_unit(self.viewscale)
	    self.remake_pnts()
	    self.remake_dims()
	    self.end()

####   Event callback functions    ####
    
    def mouse_event(self, widget, event=None, data=None):
        if event.type == gtk.gdk.BUTTON_PRESS:
            if event.button == 1:
                if (event.state & gtk.gdk.CONTROL_MASK):
                    self.remember_x = event.x
                    self.remember_y = event.y
                    self.xbuf = [0 for i in range(8)] # smooth jittery motion
                    self.ybuf = [0 for i in range(8)]
                    return gtk.TRUE
                elif self.op in self.pt_sel_group or \
                     self.op == 'del_box':
                    x = event.x
                    y = event.y
                    if self.flip_flag:
                        y = -y
                    try:
                        pt = (float(x), float(y))
                        self.pt_stack.append(pt)
                        func = "self.%s()" % self.op
                        eval(func)
                    except:
                        pass
                    return gtk.TRUE
                elif self.op_stack:
                    if self.op_stack[-1] == 'dim_sel':
                        self.item_stack.append(widget)
                        self.op_stack.pop()
                        self.op = self.op_stack[-1]
                        func = "self.%s()" % self.op
                        eval(func)
                        return gtk.TRUE
            elif event.button == 2:
                self.end(widget)
            elif event.button == 3:
                menu = gtk.Menu()
                entries = [(gtk.STOCK_ZOOM_IN, self.zoom, 1.35),
                           (gtk.STOCK_ZOOM_OUT, self.zoom, 0.65),
                           (gtk.STOCK_ZOOM_FIT, self.zoom, 0),
                           ("Zoom Previous", self.zoom, -1),
                           ("Zoom Window", self.window, None)]
                for stock_id,callback,scale in entries:
                    item = gtk.ImageMenuItem(stock_id)
                    item.connect("activate", callback, scale)
                    item.show()
                    menu.append(item)
                menu.popup(None,None,None,event.button,event.time)
                return gtk.TRUE
        elif event.type == gtk.gdk.MOTION_NOTIFY:
            self.reparent_tmp_geom()
            x = event.x
            y = event.y
            if self.flip_flag:
                y = -y
            self.mouse_xy.set_text("X= %.2f  Y= %.2f" % (x/self.unitscale,
                                                         y/self.unitscale))
            if (event.state & gtk.gdk.BUTTON1_MASK) and \
               (event.state & gtk.gdk.CONTROL_MASK):
                new_x = event.x
                new_y = event.y
                delta_x = self.remember_x - new_x 
                delta_y = self.remember_y - new_y
                self.xbuf.insert(0, delta_x)
                self.xbuf.pop()
                sum = 0
                for bi in self.xbuf:
                    sum = sum + bi
                delta_x = sum / float(len(self.xbuf))
                
                self.ybuf.insert(0, delta_y)
                self.ybuf.pop()
                sum = 0
                for bi in self.ybuf:
                    sum = sum + bi
                delta_y = sum / len(self.ybuf)
                
                x, y = self.acanvas.get_scroll_offsets()
                self.acanvas.scroll_to(x + delta_x, y + delta_y)
                self.remember_x = new_x
                self.remember_y = new_y
            return gtk.TRUE
        else:
            return gtk.FALSE

    def cline_sel_event(self, widget, event=None, data=None):
        """Place coefficients of selected construction line on stack."""
        if event.type == gtk.gdk.BUTTON_PRESS and \
           event.button == 1:
            if self.op == 'del':
                self.delete_cgeom(widget)
                return gtk.TRUE
            elif self.op in self.line_sel_group:
                self.cline_stack.append(self.cdict.get(widget))
                func = "self.%s()" % self.op
                eval(func)
                return gtk.TRUE
        elif event.type == gtk.gdk.ENTER_NOTIFY:
            if self.op in self.line_sel_group or \
               self.op == 'del':
                widget.set(fill_color='cyan')
                return gtk.TRUE
        elif event.type == gtk.gdk.LEAVE_NOTIFY:
            if self.op in self.line_sel_group or \
               self.op == 'del':
                widget.set(fill_color='magenta')
                return gtk.FALSE
        else:
            return gtk.FALSE

    def pt_sel_event(self, widget, event=None, data=None):
        if self.op in self.pt_sel_group:
            if event.type == gtk.gdk.BUTTON_PRESS:
                if event.button == 1:
                    p = self.pdict.get(widget)
                    self.pt_stack.append(p)
                    func = "self.%s()" % self.op
                    eval(func)
                return gtk.TRUE
            elif event.type == gtk.gdk.ENTER_NOTIFY:
                widget.set(fill_color='yellow')
                return gtk.TRUE
            elif event.type == gtk.gdk.LEAVE_NOTIFY:
                widget.set(fill_color_rgba=0)
                return gtk.FALSE
        return gtk.FALSE

    def geom_sel_event(self, widget, event=None, data=None):
        if event.type == gtk.gdk.BUTTON_PRESS:
            if event.button == 1:
                if self.op == "del":
                    self.delete_geom(widget)
                    return gtk.TRUE
                elif self.op in self.pt_sel_group:
                    # Place ctr or end pt of selected geom elem on pt_stack.
                    # CTRL-SHIFT-LMB selects mid point of element.
                    # LMB selects closest end point.
                    if self.flip_flag:
                        clicked_pt = (event.x, -event.y)
                    else:
                        clicked_pt = (event.x, event.y)
                    if widget in self.ldict.keys():
                        p1, p2 = self.ldict.get(widget)
                        pc = midpoint(p1, p2)
                    elif widget in self.arcdict.keys():
                        param = self.arcdict.get(widget)
                        p1, p2 = arc_end_pts(param)
                        pc = (param[0], param[1])
                    elif widget in self.cirdict.keys():
                        param = self.cirdict.get(widget)
                        pc = (param[0], param[1])
                    p = None
                    if (event.state & gtk.gdk.CONTROL_MASK) and \
                       (event.state & gtk.gdk.SHIFT_MASK):
                        p = pc
                    else:
                        p = closer(clicked_pt, p1, p2)
                    if p:
                        self.pt_stack.append(p)
                    func = "self.%s()" % self.op
                    eval(func)
                    return gtk.TRUE
        elif event.type == gtk.gdk.ENTER_NOTIFY:
            if self.op in self.pt_sel_group or \
               self.op == "del":
                if widget in self.cirdict.keys():
                    widget.set(outline_color='yellow')
                else:
                    widget.set(fill_color='yellow')
                return gtk.TRUE
        elif event.type == gtk.gdk.LEAVE_NOTIFY:
            if self.op in self.pt_sel_group or \
               self.op == "del":
                if widget in self.cirdict.keys():
                    widget.set(outline_color='white')
                else:
                    widget.set(fill_color='white')
                return gtk.FALSE
        else:
            return gtk.FALSE

    def dim_sel_event(self, widget, event=None, data=None):
        if event.type == gtk.gdk.BUTTON_PRESS:
            if event.button == 1:
                if self.op == "del":
                    self.delete_dim(widget)
                    return gtk.TRUE
                elif self.op_stack:
                    if self.op_stack[-1] == 'dim_sel':
                        self.item_stack.append(widget)
                        self.op_stack.pop()
                        self.op = self.op_stack[-1]
                        func = "self.%s()" % self.op
                        eval(func)
                        return gtk.TRUE
                else:
                    return gtk.FALSE

    def enter_data_event(self, widget, entry):
        """Take float and point values entered by user and place on stack.

        Act as a transparent surrogate for screen picks."""
        entry_text = entry.get_text()
        if self.op in ('line', 'rect', 'arc', 'hvcline', 'perp_cl',
                       'hdim', 'vdim'):
            try:
                x, y = entry_text.split(',')
                pt = (float(x)*self.unitscale, float(y)*self.unitscale)
                self.pt_stack.append(pt)
                func = "self.%s()" % self.op
                eval(func)
            except:
                pass
        elif self.op in ('hcline', 'vcline', 'par_cl', 'ang_bsctr', 'fillet'):
            try:
                self.float_stack.append(float(entry_text)*self.unitscale)
                func = "self.%s()" % self.op
                eval(func)
            except:
                pass
        elif self.op in ('circle', 'cline', 'lin_bsctr'):
            try:
                lst = entry_text.split(',')
                if len(lst) == 2:
                    x, y = lst
                    pt = (float(x)*self.unitscale, float(y)*self.unitscale)
                    self.pt_stack.append(pt)
                elif len(lst) == 1:
                    self.float_stack.append(float(lst[0])*self.unitscale)
                func = "self.%s()" % self.op
                eval(func)
            except:
                pass
        entry.set_text("")
        return gtk.TRUE
        
####    Construction Lines  ####
        
    def add_cline(self, widget, event=None, data=None):
        """Initialize operation for creating a c-line defined by 2 points

        or pt & angle. Lines may be horiz, vert or angled."""
        self.status_line.set_text(
            "Pick first point or enter x,y coordinates for construction line.")
        self.op = 'cline'
        self.pt_stack = []
        self.float_stack = []
        self.entry.grab_focus()

    def add_hvcline(self, widget, event=None, data=None):
        """Initialize operation for creating a pair of h & v c-lines."""
        self.status_line.set_text(
            "Pick point for horz & vert construction lines or enter x,y values")
        self.op = "hvcline"
        self.pt_stack = []
        self.entry.grab_focus()

    def hvcline(self):
        """Create a horizontal and vertical pair of c-lines through a pt."""
        x, y = self.pt_stack.pop()
        self.float_stack.append(y)
        self.hcline()
        self.float_stack.append(x)
        self.vcline()

    def cline(self):
        """Create a c-line from 2 pts or from 1 pt and angle on stack.

        Dispatch to appropriate h, v or a function."""
        if len(self.pt_stack) == 1 and not self.float_stack:
            self.status_line.set_text(
                "Specify second point or angle (deg) for construction line.")
            return
        else:
            if self.float_stack:
                ang = self.float_stack.pop()
                p1 = self.pt_stack[-1]
                ang = ang * math.pi / 180
                dx = math.cos(ang)
                dy = math.sin(ang)
                p2 = (p1[0]+dx, p1[1]+dy)
                self.pt_stack.append(p2)
            else:
                p1 = self.pt_stack[-2]
                p2 = self.pt_stack[-1]
            a, b, c = cnvrt_2pts_to_coef(p1, p2)
            if a == 0:
                throw_away = self.pt_stack.pop()
                self.hcline()
            elif b == 0:
                throw_away = self.pt_stack.pop()
                self.vcline()
            else:
                self.acline()
            self.status_line.set_text(
                "Pick first point or enter x,y coordinates for construction line.")

    def cline2(self, cline=None):
        """Create a c-line from c-line passed as argument or on stack.

        Dispatch to appropriate h, v, or a function."""
        if not cline:
            cline = self.cline_stack.pop()
        a, b, c = cline
        if a == 0:
            self.float_stack.append(-c/b)
            self.hcline()
        elif b == 0:
            self.float_stack.append(-c/a)
            self.vcline()
        else:
            self.cline_stack.append(cline)
            self.acline(frm_pts=0)
        

    def acline(self, frm_pts=1):
        """Create an angled c-line from values on stack.

        if frm_pts=1: create from 2 point values (default)
        if frm_pts=0: create from (a,b,c) coefficients
        Destroy and remake all intersection points."""
        if frm_pts:
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            cline = cnvrt_2pts_to_coef(p1, p2)
        else:
            cline = self.cline_stack.pop()
        self.acline_gen(cline)
        self.remake_pnts()

    def acline_gen(self, cline):
        """Generate an angled c-line from cline coefficients (a, b, c)"""
        pts = cline_box_intrsctn(cline, self.canvas_size)
        if len(pts) == 2:
            p1, p2 = pts
            points = (p1[0], p1[1], p2[0], p2[1])
            w = self.cons.add('GnomeCanvasLine',
                              points=points,
                              fill_color='magenta',
                              line_style='GDK_LINE_ON_OFF_DASH',
                              width_pixels=self.linewidth/2)
            w.connect("event", self.cline_sel_event)
            self.cdict[w] = cline
            self.acl.append(w)

    def add_hcline(self, widget, event=None, data=None):
        """Initialize horizontal construction line creation operation."""
        self.status_line.set_text(
            "Pick point for horizontal construction line or enter y value")
        self.op = "hcline"
        self.float_stack = []
        self.entry.grab_focus()

    def hcline(self):
        """Create a horizontal c line from pt or float value on stack.

        Destroy and remake all intersection points."""
        if self.pt_stack:
            x, y = self.pt_stack.pop()
        elif self.float_stack:
            y = self.float_stack.pop()
        self.hcline_gen(y)
        self.remake_pnts()

    def hcline_gen(self, y):
        """Create a horizontal c line from float passed as arg."""
        w = self.cons.add('GnomeCanvasLine',
                          points=(self.canvas_size[0], y,
                                  self.canvas_size[2], y),
                          fill_color='magenta',
                          line_style='GDK_LINE_ON_OFF_DASH',
                          width_pixels=self.linewidth/2)
        line_coef = cnvrt_2pts_to_coef((0.0, y), (600.0, y))
        w.connect("event", self.cline_sel_event)
        self.cdict[w] = line_coef
        self.hcl.append(w)

    def add_vcline(self, widget, event=None, data=None):
        """Initialize vertical construction line creation operation."""
        self.status_line.set_text(
            "Pick point for vertical construction line or enter x value")
        self.op = "vcline"
        self.float_stack = []
        self.entry.grab_focus()
        
    def vcline(self):
        """Create a vertical c line from pt or float value on stack.

        Destroy and remake all intersection points."""
        if self.pt_stack:
            x, y = self.pt_stack.pop()
        elif self.float_stack:
            x = self.float_stack.pop()
        self.vcline_gen(x)
        self.remake_pnts()
        
    def vcline_gen(self, x):
        """Create a vertical c line from float passed as arg."""
        w = self.cons.add('GnomeCanvasLine',
                          points=(x, self.canvas_size[1],
                                  x, self.canvas_size[3]),
                          fill_color='magenta',
                          line_style='GDK_LINE_ON_OFF_DASH',
                          width_pixels=self.linewidth/2)
        line_coef = cnvrt_2pts_to_coef((x, 0.0), (x, 500.0))
        w.connect("event", self.cline_sel_event)
        self.cdict[w] = line_coef
        self.vcl.append(w)

    def add_parcl(self, widget, event=None, data=None):
        """Initialize parallel construction line creation operation."""
        self.status_line.set_text("Pick line for parallel construction")
        self.op = "par_cl"
        self.pt_stack = []
        self.cline_stack = []
        self.float_stack = []
        
    def par_cl(self):
        """Create 1 parallel c line through a point, or 2 straddling lines

        offset by a distance, from values on stack."""
        if (not self.float_stack) and (not self.pt_stack):
            self.status_line.set_text(
                "Pick point for parallel line or enter offset distance")
            return
        elif self.cline_stack:
            baseline = self.cline_stack.pop()
            d = 0
            if self.float_stack:
                d = self.float_stack.pop()
                cline1, cline2 = para_lines(baseline, d)
                self.cline_stack.extend([cline1, cline2])
                self.cline2()
                self.cline2()
            else:
                p = self.pt_stack.pop()
                newline = para_line(baseline, p)
                self.cline_stack.append(newline)
                self.cline2()
            self.remake_pnts()
            self.status_line.set_text("Pick line for parallel construction")

    def add_perpcl(self, widget, event=None, data=None):
        """Initialize perpendicular construction line creation operation."""
        self.status_line.set_text("Pick line for perpendicular construction")
        self.op = "perp_cl"
        self.cline_stack = []
        self.pt_stack = []

    def perp_cl(self):
        """Create a perpendicular c line through a selected point."""
        if not self.pt_stack:
            self.status_line.set_text(
                "Pick point for perpendicular c-line or enter point coordinates x,y")
            return
        if self.cline_stack:
            baseline = self.cline_stack.pop()
            pt = self.pt_stack.pop()
            cline = perp_line(baseline, pt)
            self.cline_stack.append(cline)
            self.cline2()
            self.remake_pnts()
            self.status_line.set_text(
                "Pick line for perpendicular construction")

    def add_ang_bsctr(self, widget, event=None, data=None):
        """Initialize angular bisector construction line creation operation."""
        self.status_line.set_text(
            "Enter bisector factor (def=.5) or pick first construction line")
        self.op = "ang_bsctr"
        self.cline_stack = []
        self.float_stack = []

    def ang_bsctr(self):
        """Create an angular bisector construction line."""
        if not self.cline_stack:
            self.status_line.set_text("Pick first construction line")
            return
        elif len(self.cline_stack) == 1:
            self.status_line.set_text("Pick second construction line")
            return
        f = .5
        if self.float_stack:
            f = self.float_stack.pop()
        a2, b2, c2 = line2 = self.cline_stack.pop()
        a1, b1, c1 = line1 = self.cline_stack.pop()
        x, y = p0 = intersection(line1, line2)
        ang1 = math.atan2(-a1, b1)
        ang2 = math.atan2(-a2, b2)
        deltang = ang2 - ang1
        ang3 = f * deltang + ang1
        dy3 = math.sin(ang3)
        dx3 = math.cos(ang3)
        p3 = (p0[0] + dx3, p0[1] + dy3)
        line3 = cnvrt_2pts_to_coef(p0, p3)
        self.cline_stack.append(line3)
        self.cline2()
        self.remake_pnts()
        self.status_line.set_text(
            "Enter bisector factor (def=.5) or pick first construction line")

    def add_lin_bsctr(self, widget, event=None, data=None):
        """Initialize linear bisector construction line creation operation."""
        self.status_line.set_text(
            "Enter bisector factor (def=.5) or pick first point")
        self.op = "lin_bsctr"
        self.pt_stack = []
        self.float_stack = []

    def lin_bsctr(self):
        """Create a linear bisector construction line."""
        if not self.pt_stack:
            self.status_line.set_text("Pick first point")
            return
        elif len(self.pt_stack) == 1:
            self.status_line.set_text("Pick second point")
            return
        f = .5
        if self.float_stack:
            f = self.float_stack.pop()
        p2 = self.pt_stack.pop()
        p1 = self.pt_stack.pop()
        p0 = midpoint(p1, p2, f)
        baseline = cnvrt_2pts_to_coef(p1, p2)
        newline = perp_line(baseline, p0)
        self.cline_stack.append(newline)
        self.cline2()
        self.remake_pnts()
        self.status_line.set_text(
            "Enter bisector factor (def=.5) or pick first point")

    def add_ccirc(self, widget, event=None, data=None):
        """Initialize construction circle creation operation."""
        self.status_line.set_text("Not yet implemented")

    def remake_pnts(self):
        """Delete all the existing intersection pts and make all new ones.

        These are "clickable points" created at the intersections of all
        construction lines. In order to maintain their correct apparent
        size, they need to be recreated whenever the scale changes.
        This may happen fairly frequently. They also need to be created
        and deleted as construction lines are created and deleted.
        Rather than keep track of which points go with which lines, just
        do a wholesale purge and remake whenever construction lines are
        created or deleted or whenever the zoom changes."""
        
        for i in xrange(len(self.pdict)):
            k, v = self.pdict.popitem()
            k.destroy()

        for vcl in self.vcl:
            for hcl in self.hcl:
                pt = intersection(self.cdict.get(vcl), self.cdict.get(hcl))
                if pt:
                    self.point_gen(pt)
        for acl in self.acl:
            for hcl in self.hcl:
                pt = intersection(self.cdict.get(acl), self.cdict.get(hcl))
                if pt:
                    self.point_gen(pt)
            for vcl in self.vcl:
                pt = intersection(self.cdict.get(acl), self.cdict.get(vcl))
                if pt:
                    self.point_gen(pt)
        for i in range(len(self.acl)):
            ac0 = self.acl[i]
            for acl in self.acl[i:]:
                pt = intersection(self.cdict.get(ac0), self.cdict.get(acl))
                if pt:
                    self.point_gen(pt)

    def point_gen(self, pt):
        """Create an intersection point."""
        pt_rad = self.pt_size/2/self.viewscale
        if pt not in self.pdict.values():
            p = self.pnts.add('GnomeCanvasEllipse',
                              x1 = pt[0] - pt_rad,
                              y1 = pt[1] - pt_rad,
                              x2 = pt[0] + pt_rad,
                              y2 = pt[1] + pt_rad,
                              fill_color=None,
                              outline_color_rgba=0,
                              width_pixels=self.linewidth)
            p.connect("event", self.pt_sel_event)
            self.pdict[p] = pt

    def delete_cgeom(self, widget):
        """ Delete selected construction element."""
        if widget in self.hcl:
            self.hcl.remove(widget)
            del self.cdict[widget]
            widget.destroy()
            self.remake_pnts()
        elif widget in self.vcl:
            self.vcl.remove(widget)
            del self.cdict[widget]
            widget.destroy()
            self.remake_pnts()
        elif widget in self.acl:
            self.acl.remove(widget)
            del self.cdict[widget]
            widget.destroy()
            self.remake_pnts()

    def del_all_c(self, widget, event=None, data=None):
        """Delete all construction elements."""
        self.op = ""
        self.hcl = []
        self.vcl = []
        self.acl = []
        for i in xrange(len(self.cdict)):
            k, v = self.cdict.popitem()
            k.destroy()
        self.remake_pnts()
        
####    Geometry Lines     ####
        
    def add_line(self, widget, event=None, data=None):
        """Initialize geometry line segment creation operation."""
        self.pt_stack = []
        self.op = "line"
        self.status_line.set_text(
            "Click first point for line or enter coordinates x,y")
        self.entry.grab_focus()

    def line(self):
        """Create line from values on point stack."""
        if len(self.pt_stack) == 1:
            self.status_line.set_text(
                "Click second point for line or enter coordinates x,y")
        else:
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            line = (p1, p2)
            self.line_gen(line)
            self.status_line.set_text(
                "Click first point for line or enter coordinates x,y")

    def line_gen(self, line):
        """Generate straight line segment from line = (p1, p2)."""
        p1, p2 = line
        points=(p1[0], p1[1], p2[0], p2[1])
        w = self.tmpg.add('GnomeCanvasLine',
                          points=points,
                          fill_color='white',
                          width_pixels=self.linewidth)
        w.connect("event", self.geom_sel_event)
        self.ldict[w] = line
        self.tmp_geom_list.append(w)

    def add_polyline(self, widget, event=None, data=None):
        """Initialize geometry poly-line creation operation."""
        self.pt_stack = []
        self.op = "polyline"
        self.status_line.set_text(
            "Click first point for polyline or enter coordinates x,y")
        self.entry.grab_focus()

    def polyline(self):
        """Create a chain of independent line segments from values

        on point stack. If chain closes on starting point, close chain
        and start a new one."""
        if len(self.pt_stack) == 1:
            self.status_line.set_text(
                "Click second point for line or enter coordinates x,y")
        elif len(self.pt_stack) == 2:
            p2 = self.pt_stack[-1]
            p1 = self.pt_stack[0]
            self.line_gen((p1, p2))
            self.status_line.set_text(
                "Click next point for line or enter coordinates x,y")
        elif len(self.pt_stack) == 3:
            p2 = self.pt_stack[-1]
            p1 = self.pt_stack[-2]
            p0 = self.pt_stack[0]
            if p2 == p0:
                self.line_gen((p1, p2))
                self.pt_stack = []
                self.status_line.set_text(
                    "Click first point for polyline or enter coordinates x,y")
            else:
                self.pt_stack.pop(1)
                self.line_gen((p1, p2))
                self.status_line.set_text(
                    "Click next point for line or enter coordinates x,y")
                

    def add_rect(self, widget, event=None, data=None):
        """Initialize geometry rectangle creation operation."""
        self.pt_stack = []
        self.op = "rect"
        self.status_line.set_text(
            "Click First Corner of Rectangle or enter coordinates x,y")
        self.entry.grab_focus()

    def rect(self):
        """Create rectangle of 4 independent line segments from  2 pts

        on stack representing diagonally opposite corners."""
        if len(self.pt_stack) == 1:
            self.status_line.set_text(
                "Click Opposite Corner of Rectangle or enter coordinates x,y")
        else:
            x2, y2 = self.pt_stack.pop()
            x1, y1 = self.pt_stack.pop()
            pt_list = [((x1, y1), (x1, y2)),
                       ((x1, y2), (x2, y2)),
                       ((x2, y2), (x2, y1)),
                       ((x2, y1), (x1, y1))]
            for pts in pt_list:
                p1, p2 = pts
                self.line_gen((p1, p2))
            self.status_line.set_text(
                "Click First Corner of Rectangle or enter coordinates x,y")

    def add_circ(self, widget, event=None, data=None):
        """Initialize geometry circle creation operation."""
        self.op = "circle"
        self.status_line.set_text(
            "Click pnt for cntr of circle or enter coords")
        self.entry.grab_focus()
        self.pt_stack = []
        self.float_stack = []

    def circle(self):
        """Create a circle from a point and a float value on stack."""
        if len(self.pt_stack) == 1 and not self.float_stack:
            self.status_line.set_text("Click pnt on circle or enter radius")
            return
        elif self.float_stack:
            x, y = self.pt_stack.pop()
            r = self.float_stack.pop()
        else:
            p1 = self.pt_stack.pop()
            p0 = self.pt_stack.pop()
            x, y = p0
            r = pt2pt_dist(p1, p0)
        param = x, y, r
        self.circle_gen(param)
        self.status_line.set_text(
            "Click pnt for cntr of circle or enter coords")

    def circle_gen(self, param):
        """Generate a circle from param = (x, y, r) """
        x, y, r = param
        w = self.tmpg.add('GnomeCanvasEllipse',
                          x1=x-r, y1=y-r, x2=x+r, y2=y+r,
                          outline_color='white',
                          width_pixels=self.linewidth)
        w.connect("event", self.geom_sel_event)
        self.cirdict[w] = param
        self.tmp_geom_list.append(w)

    def add_arc(self, widget, event=None, data=None):
        """Initialize geometry arc creation operation."""
        self.pt_stack = []
        self.op = "arc"
        self.status_line.set_text("Pick center of arc or enter coordinates.")
        self.entry.grab_focus()

    def arc(self):
        """Create arc from 3 pt values on stack: cntr pt, start pt, end pt."""
        if len(self.pt_stack) == 1:
            self.status_line.set_text("Pick starting point for arc.")
            return
        elif len(self.pt_stack) == 2:
            self.status_line.set_text("Pick ending point for arc.")
            return
        else:
            x1, y1 = self.pt_stack.pop()    # end pt
            x0, y0 = self.pt_stack.pop()    # start pt
            x, y = self.pt_stack.pop()      # cntr
            r = math.sqrt((x0-x)**2 + (y0-y)**2)    # radius
            a0 = math.atan2((y0-y), (x0-x))   # start angle
            a1 = math.atan2((y1-y), (x1-x))   # end angle
            param = (x, y, r, a0, a1)
            self.arc_gen(param)
            self.status_line.set_text(
                "Pick center of arc or enter coordinates.")

    def arc_gen(self, param):
        """Generate arc from param = (xc, yc, r, a0, a1)

        Simulate arc with a portion of a regular polygon."""
        xc, yc, r, a0, a1 = param   # a0 & a1 expressed in radians
        da = (a1-a0)
        if da < 0:
            da = 2 * math.pi + da
        g = 0.15 # angle (radians) subtended by each facet
        n = int(da/g) # number of facets in arc
        points = []
        a = a0
        for i in range(n+1):
            u = xc + r*math.cos(a)
            points.append(u)
            v = yc + r*math.sin(a)
            points.append(v)
            a = a + da/n
        w = self.tmpg.add('GnomeCanvasLine',
                          points=points,
                          smooth=gtk.TRUE,
                          fill_color='white',
                          width_pixels=self.linewidth)
        w.connect("event", self.geom_sel_event)
        self.arcdict[w] = param
        self.tmp_geom_list.append(w)

    def add_fillet(self, widget, event=None, data=None):
        """Initialize operation to create fillet."""
        self.pt_stack = []
        self.float_stack = []
        self.op = "fillet"
        self.status_line.set_text("Click line near corner to fillet.")

    def fillet(self):
        """Generate fillet on corner by shortening two adjacent line
        segments and adding arc."""
        if self.pt_stack and not self.float_stack:
            self.status_line.set_text("Enter fillet radius.")
            return
        elif self.pt_stack and self.float_stack:
            p = self.pt_stack.pop()
            r = self.float_stack[-1]
            lines = []
            for k, v in self.ldict.items():
                if p in v:
                    lines.append(k)
            if len(lines) != 2:
                print "Couldn't find exactly 2 adjacent lines."
                return
            endpts = []
            for k in lines:
                p1, p2 = self.ldict[k]
                if p1 == p:
                    endpts.append(p2)
                else:
                    endpts.append(p1)
            pc, p1, p2 = pts = find_fillet_pts(r, p, endpts[0], endpts[1])
            xc, yc = pc
            x1, y1 = p1
            x2, y2 = p2
            for k in lines:
                del self.ldict[k]
                k.destroy()
            self.line_gen((p1, endpts[0]))
            self.line_gen((p2, endpts[1]))
            a0 = math.atan2(y1-yc, x1-xc)
            a1 = math.atan2(y2-yc, x2-xc)
            if (-math.pi < a1-a0 < 0) or \
               (math.pi < a1-a0 < 2*math.pi):
                a0, a1 = a1, a0
            arcparam = (xc, yc, r, a0, a1)
            self.arc_gen(arcparam)
            self.status_line.set_text(
                "Click line near corner to fillet or enter new radius.")

    def reparent_tmp_geom(self):
        """Reparent canvas items in tmpg group to geom group."""
        if len(self.tmp_geom_list):
            for i in range(len(self.tmp_geom_list)):
                item = self.tmp_geom_list.pop()
                item.reparent(self.geom)
            self.acanvas.set_pixels_per_unit(self.viewscale)

    def delete_geom(self, widget):
        if widget in self.ldict.keys():
            del self.ldict[widget]
            widget.destroy()
        elif widget in self.cirdict.keys():
            del self.cirdict[widget]
            widget.destroy()
        elif widget in self.arcdict.keys():
            del self.arcdict[widget]
            widget.destroy()

    def del_all_geom(self, widget, event=None, data=None):
        """Delete all geometry elements."""
        self.op = ""
        for i in xrange(len(self.ldict)):
            k, v = self.ldict.popitem()
            k.destroy()
        for i in xrange(len(self.cirdict)):
            k, v = self.cirdict.popitem()
            k.destroy()
        for i in xrange(len(self.arcdict)):
            k, v = self.arcdict.popitem()
            k.destroy()

####    Dimensions   ####

    def add_hdim(self, widget, event=None, data=None):
        '''Initialize horizontal dimensioning operation.'''
        self.op = 'hdim'
        self.status_line.set_text(
            "Select first point for horizontal dimension")
        self.pt_stack = []

    def hdim(self):
        """Place horizontal dimension between two points."""
        if len(self.pt_stack) == 1:
            self.status_line.set_text(
                "Select second point for horizontal dimension")
            return
        if len(self.pt_stack) == 2:
            self.status_line.set_text(
                "Click location for dimension text")
            return
        if len(self.pt_stack) == 3:
            p3 = self.pt_stack.pop()
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            self.dim_gen(p1, p2, p3, 'h')
            self.status_line.set_text(
                "Select first point for horizontal dimension")

    def add_vdim(self, widget, event=None, data=None):
        '''Initialize vertical dimensioning operation.'''
        self.op = 'vdim'
        self.status_line.set_text(
            "Select first point for vertical dimension")
        self.pt_stack = []

    def vdim(self):
        """Place vertical dimension between two points."""
        if len(self.pt_stack) == 1:
            self.status_line.set_text(
                "Select second point for vertical dimension")
            return
        if len(self.pt_stack) == 2:
            self.status_line.set_text(
                "Click location for dimension text")
            return
        if len(self.pt_stack) == 3:
            p3 = self.pt_stack.pop()
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            self.dim_gen(p1, p2, p3, 'v')
            self.status_line.set_text(
                "Select first point for vertical dimension")

    def dim_gen(self, p1, p2, p3, type):
        """Generate dim between p1 & p2, dim loc at p3, type = h or v"""
        gap = self.dimgap/self.viewscale
        dimg = self.dims.add('GnomeCanvasGroup', x=0, y=0)
        dimg.connect('event', self.dim_sel_event)
        dimg.connect('event', self.mouse_event)
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3
        if type == 'h':
            value = abs(x2-x1)/self.unitscale
            if not value:
                return
            if y3 > y1:
                epts = (x1, y1+gap, x1, y3+gap)
            else:
                epts = (x1, y1-gap, x1, y3-gap)
            if y3 > y2:
                fpts = (x2, y2+gap, x2, y3+gap)
            else:
                fpts = (x2, y2-gap, x2, y3-gap)
        elif type == 'v':
            value = abs(y2-y1)/self.unitscale
            if not value:
                return
            if x3 > x1:
                epts = (x1+gap, y1, x3+gap, y1)
            else:
                epts = (x1-gap, y1, x3-gap, y1)
            if x3 > x2:
                fpts = (x2+gap, y2, x3+gap, y2)
            else:
                fpts = (x2-gap, y2, x3-gap, y2)
        # extension lines
        e = dimg.add('GnomeCanvasLine',
                     points=epts,
                     fill_color=self.dim_line_color,
                     width_pixels=self.linewidth)
        f = dimg.add('GnomeCanvasLine',
                     points=fpts,
                     fill_color=self.dim_line_color,
                     width_pixels=self.linewidth)
        # text
        text = '%.3f' % value
        t = dimg.add('GnomeCanvasText',
                     text=text,
                     x=x3,
                     y=y3,
                     anchor='center',
                     size=15360,
                     fill_color=self.dim_text_color)
        if type == 'h':
            w = t.get_property('text_width')
            if x1 < x3 < x2:
                q1 = (x3-w/2-gap, y3, x1, y3)
                q2 = (x3+w/2+gap, y3, x2, y3)
            elif x2 < x3 < x1:
                q1 = (x3-w/2-gap, y3, x2, y3)
                q2 = (x3+w/2+gap, y3, x1, y3)

        elif type == 'v':
            h = t.get_property('text_height')
            if y1 < y3 < y2:
                q1 = (x3, y3-h/2-gap, x3, y1)
                q2 = (x3, y3+h/2+gap, x3, y2)
            elif y2 < y3 < y1:
                q1 = (x3, y3-h/2-gap, x3, y2)
                q2 = (x3, y3+h/2+gap, x3, y1)
        # dimension lines
        l1 = dimg.add('GnomeCanvasLine',
                      points=q1,
                      fill_color=self.dim_line_color,
                      first_arrowhead=False,
                      last_arrowhead=True,
                      arrow_shape_a=5,
                      arrow_shape_b=10,
                      arrow_shape_c=5,
                      width_pixels=self.linewidth)
        l2 = dimg.add('GnomeCanvasLine',
                      points=q2,
                      fill_color=self.dim_line_color,
                      first_arrowhead=False,
                      last_arrowhead=True,
                      arrow_shape_a=5,
                      arrow_shape_b=10,
                      arrow_shape_c=5,
                      width_pixels=self.linewidth)
        self.dimdict[dimg] = (p1, p2, p3)
        self.dim_parts_dict[dimg] = (t, e, f, l1, l2)
        if type == 'h':
            self.hdims.append(dimg)
        elif type == 'v':
            self.vdims.append(dimg)

    def remake_dims(self):
        hdimvals = [self.dimdict.get(key) for key in self.hdims]
        vdimvals = [self.dimdict.get(key) for key in self.vdims]
        self.del_all_dims()
        for each in hdimvals:
            p1, p2, p3 = each
            self.dim_gen(p1, p2, p3, 'h')
        for each in vdimvals:
            p1, p2, p3 = each
            self.dim_gen(p1, p2, p3, 'v')

    def move_dim(self, widget, event=None, data=None):
        """Initialize dimension move operation."""
        self.op_stack = ['dim_move', 'dim_sel']
        self.op = []
        self.status_line.set_text("Click on dimension to move")
        self.pt_stack = []
        self.item_stack = []

    def dim_move(self):
        """Move selected dimension to new location."""
        if self.item_stack and not self.pt_stack:
            self.status_line.set_text("Click new location for dimension")
            return
        if self.item_stack and self.pt_stack:
            new_pt = self.pt_stack.pop()
            dim = self.item_stack.pop()
            dimvals = self.dimdict.get(dim)
            p1, p2, p3 = dimvals[0], dimvals[1], new_pt
            if dim in self.hdims:
                type = 'h'
            elif dim in self.vdims:
                type = 'v'
            self.delete_dim(dim)
            self.dim_gen(p1, p2, p3, type)
            self.pt_stack = []
            self.op_stack = ['dim_move', 'dim_sel']
            self.op = []
            self.status_line.set_text("Click on dimension to move")

    def delete_dim(self, widget):
        """Delete specified dimension."""
        if widget in self.hdims:
            self.hdims.remove(widget)
        if widget in self.vdims:
            self.vdims.remove(widget)
        for part in self.dim_parts_dict[widget]:
            part.destroy()
        del self.dimdict[widget]
        del self.dim_parts_dict[widget]
        widget.destroy()

    def del_all_dims(self, widget=None, event=None, data=None):
        """Delete all dimensions."""
        self.hdims = []
        self.vdims = []
        for i in xrange(len(self.dimdict)):
            grp, val = self.dimdict.popitem()
            for part in self.dim_parts_dict[grp]:
                part.destroy()
            del self.dim_parts_dict[grp]
            grp.destroy()

####    Miscellaneous   ####

    def end(self, widget=None, event=None, data=None):
        """End current operation and clear stacks."""
        self.float_stack = []
        self.pt_stack = []
        self.line_stack = []
        self.cline_stack = []
        self.op_stack = []
        self.op = ""
        self.status_line.set_text("Choose a command from the menubar.")

    def del_elem(self, widget, event=None, data=None):
        """Callback for deleting individual elements."""
        self.op = "del"
        self.status_line.set_text("Click on element to delete")

    def del_inbox(self, type):
        """Initialize operation of deleting elements in box."""
        self.pt_stack = []
        self.op = "del_box"
        self.op_stack = [type]
        self.status_line.set_text("Pick first corner of box.")

    def del_box(self):
        """Delete elements of type=op_stack[-1] in box defined by pt_stack[-2:]."""
        if len(self.pt_stack) == 1:
            self.status_line.set_text("Pick opposite corner of box.")
            return
        elif len(self.pt_stack) == 2:
            type = self.op_stack.pop()
            x2, y2 = p2 = self.pt_stack.pop()
            x1, y1 = p1 = self.pt_stack.pop()
            box = (x1, y1, x2, y2)
            if type in ('all', 'cgeo'):
                for k, v in self.cdict.items():
                    if cline_box_intrsctn(v, box):
                        self.delete_cgeom(k)
            else:
                print "Deleting %s in box not yet implemented" % type
            self.end()

####    CAD screen  ####

    def main(self):

        vbox = gtk.VBox()
        vbox.show()

        # Create canvas
        self.acanvas = gnome.canvas.Canvas(aa=gtk.TRUE)
        self.acanvas.set_size_request(self.width, self.height)
        x1, y1, x2, y2 = self.canvas_size
        self.acanvas.set_scroll_region(x1, y1, x2, y2)
        xc, yc = self.acanvas.w2c(0,0)
        # set initial screen location for 0,0 coordinates
        #self.acanvas.scroll_to(xc-self.width/2, yc-self.height/2)  # ctr
        self.acanvas.scroll_to(xc-self.width/10, yc-self.height*.9) # bot-lft
        vbox.pack_start(self.acanvas)
        self.acanvas.show()
        self.botm = self.acanvas.root().add('GnomeCanvasGroup', x=0, y=0)
        w = self.botm.add('GnomeCanvasRect',
                          x1=x1, y1=y1, x2=x2, y2=y2,
                          fill_color='black',
                          width_pixels=self.linewidth)
        w.connect("event", self.mouse_event)
        self.cons = self.botm.add('GnomeCanvasGroup', x=0, y=0)
        self.tmpg = self.botm.add('GnomeCanvasGroup', x=0, y=0)
        self.geom = self.botm.add('GnomeCanvasGroup', x=0, y=0)
        self.pnts = self.botm.add('GnomeCanvasGroup', x=0, y=0)
        self.dims = self.geom.add('GnomeCanvasGroup', x=0, y=0)
        # Flip botm group so +y goes up
        if self.flip_flag:
            self.botm.affine_absolute((1, 0, 0, -1, 0, 0))

        hbox = gtk.HBox(spacing=10)
        vbox.pack_start(hbox)
        hbox.show()

        # Create entry widget
        self.entry = gtk.Entry()
        self.entry.connect("activate", self.enter_data_event, self.entry)
        hbox.pack_start(self.entry, expand=gtk.FALSE)
        self.entry.show()
        
	# Create status line
	self.status_line = gtk.Label("Start by creating some Consruction Lines")
	hbox.pack_start(self.status_line, expand=gtk.FALSE)
	self.status_line.show()

	# Create display for mouse coordinates
	self.mouse_xy = gtk.Label("X= 0.00  Y= 0.00")
	hbox.pack_end(self.mouse_xy, expand=gtk.FALSE, padding=10)
	self.mouse_xy.show()
        
	# Create label for display of current units
	self.units_display = gtk.Label("Units: %s" % self.units)
	hbox.pack_end(self.units_display, expand=gtk.FALSE)
	self.units_display.show()

        return vbox

####    Application main window     ####

class fileLoadSave:
    def __init__(self, action):
        '''action: 0 -> load; 1-> save; 2 -> save_as.'''
        self.action = action
        self.filew = gtk.FileSelection("File selection")
        self.filew.connect("destroy", self.destroy)
        self.filew.ok_button.connect("clicked", self.file_ok_sel)
        self.filew.cancel_button.connect("clicked", self.destroy)
        self.filew.show()
    
    def file_ok_sel(self, w):
        file = "%s" % self.filew.get_filename()
        if self.action:
            cw.save(file)
        else:
            cw.load(file)
        self.filew.destroy()

    def destroy(self, w):
        self.action = None
        self.filew.destroy()

def units_cb(window, action, widget):
    if action == 1:
        cw.set_units('mm')
    elif action == 2:
        cw.set_units('inches')
    elif action == 3:
        cw.set_units('feet')

def del_inbox_cb(window, action, widget):
    if action == 0:
        cw.del_inbox('all')
    elif action == 1:
        cw.del_inbox('cgeo')
    elif action == 2:
        cw.del_inbox('geom')
    elif action == 3:
        cw.del_inbox('dims')

def fileload_cb(window, action, widget):
    fs = fileLoadSave(action)

def filesave_cb(window, action, widget):
    fs = fileLoadSave(action)

def print_cb(window, action, widget):
    cw.printer()

def quit_cb(window, action, widget):
    mainquit(widget)
    
def main():
    global cw
    cw = cadWindow()
    
    # Create the toplevel window
    win = gtk.Window()
    win.set_title('CADvas')
    win.connect('destroy', mainquit)

    table = gtk.Table(1, 3, gtk.FALSE)
    win.add(table)

    # Create the menubar

    menu_items = (
        ('/_File',                      None,   None,               0,  '<Branch>' ),
        ('/File/_Open',         '<control>O',   fileload_cb,        0,  '<StockItem>', gtk.STOCK_OPEN),
        ('/File/_Save',         '<control>S',   filesave_cb,        1,  '<StockItem>', gtk.STOCK_SAVE),
        ('/File/Save _As...',           None,   filesave_cb,        2,  '<StockItem>', gtk.STOCK_SAVE),
        ('/File/_Print',        '<control>P',   print_cb,           0,  '<StockItem>', gtk.STOCK_PRINT),
        ('/File/sep1',                  None,   None,               0,  '<Separator>'),
        ('/File/_Quit',         '<control>Q',   quit_cb,            0,  '<StockItem>', gtk.STOCK_QUIT),
        
        ('/_Preferences',               None,   None,               0,  '<Branch>'),
        ('/Preferences/_Units',         None,   None,               0,  '<Branch>'),
        ('/Preferences/Units/_mm',      None,   units_cb,           1,  '<RadioItem>'),
        ('/Preferences/Units/_inches',  None,   units_cb,           2,  '/Preferences/Units/mm'),
        ('/Preferences/Units/_feet',    None,   units_cb,           3,  '/Preferences/Units/mm'),
        
        ('/C-Geo',                      None,   None,               0,  '<Branch>' ),
        ('/C-Geo/H+V',                  None,   cw.add_hvcline,     0,  '<Item>' ),
        ('/C-Geo/Horizontal',           None,   cw.add_hcline,      0,  '<Item>' ),
        ('/C-Geo/Vertical',             None,   cw.add_vcline,      0,  '<Item>' ),
        ('/C-Geo/Angled',               None,   cw.add_cline,       0,  '<Item>' ),
        ('/C-Geo/Parallel',             None,   cw.add_parcl,       0,  '<Item>' ),
        ('/C-Geo/Perpendicular',        None,   cw.add_perpcl,      0,  '<Item>' ),
        ('/C-Geo/Ang Bisector',         None,   cw.add_ang_bsctr,   0,  '<Item>' ),
        ('/C-Geo/Lin Bisector',         None,   cw.add_lin_bsctr,   0,  '<Item>' ),
        ('/C-Geo/Circle',               None,   cw.add_ccirc,       0,  '<Item>' ),

        ('/Geometry',                   None,   None,               0,  '<Branch>' ),
        ('/Geometry/Line',              None,   cw.add_line,        0,  '<Item>' ),
        ('/Geometry/Poly Line',         None,   cw.add_polyline,    0,  '<Item>' ),
        ('/Geometry/Rectangle',         None,   cw.add_rect,        0,  '<Item>' ),
        ('/Geometry/Circle',            None,   cw.add_circ,        0,  '<Item>' ),
        ('/Geometry/Arc',               None,   cw.add_arc,         0,  '<Item>' ),
        ('/Geometry/Fillet',            None,   cw.add_fillet,      0,  '<Item>' ),
        
        ('/Dimension',                  None,   None,               0,  '<Branch>' ),
        ('/Dimension/Horiz',            None,   cw.add_hdim,        0,  '<Item>' ),
        ('/Dimension/Vert',             None,   cw.add_vdim,        0,  '<Item>' ),
        ('/Dimension/Move',             None,   cw.move_dim,        0,  '<Item>' ),
        
        ('/End',                        None,   cw.end,             0,  '<Item>' ),

        ('/Delete',                     None,   None,               0,  '<Branch>' ),
        ('/Delete/Element',             None,   cw.del_elem,        0,  '<Item>' ),
        ('/Delete/All Constr',          None,   cw.del_all_c,       0,  '<Item>' ),
        ('/Delete/All Geom',            None,   cw.del_all_geom,    0,  '<Item>' ),
        ('/Delete/All Dims',            None,   cw.del_all_dims,    0,  '<Item>' ),
        ('/Delete/In Box',              None,   None,               0,  '<Branch>' ),
        ('/Delete/In Box/All',          None,   del_inbox_cb,       0,  '<Item>' ),
        ('/Delete/In Box/Construct',    None,   del_inbox_cb,       1,  '<Item>' ),
        ('/Delete/In Box/Geometry',     None,   del_inbox_cb,       2,  '<Item>' ),
        ('/Delete/In Box/Dimensions',   None,   del_inbox_cb,       3,  '<Item>' ),

        ('/_Help',                      None,   None,               0,  '<Branch>'),
        ('/Help/_About',                None,   None,               0,  ''),
        )
    
    accel_group = gtk.AccelGroup()
    win.add_accel_group(accel_group)
    
    item_factory = gtk.ItemFactory(gtk.MenuBar, '<main>', accel_group)
    
    # create menu items

    item_factory.create_items(menu_items, win)
    
    table.attach(item_factory.get_widget('<main>'),
                 # X direction           Y direction
                 0, 1,                      0, 1,
                 gtk.EXPAND | gtk.FILL,     0,
                 0,                         0)

    table.attach(cw.main(),
                 # X direction           Y direction
                 0, 1,                   2, 3,
                 gtk.EXPAND | gtk.FILL,  gtk.EXPAND | gtk.FILL,
                 0,                      0)

    win.show_all()
    gtk.main()

def mainquit(widget):
    cw.mainquit()
    gtk.main_quit()

if __name__ == '__main__':
    main()
