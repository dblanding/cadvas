#!/usr/bin/env python
#
# CADvas 
# A 2D CAD application written in Python3 and using the Tkinter canvas.
# The latest  version of this file can be found at:
# https://github.com/dblanding/cadvas
#
# Author: Doug Blanding   <dblanding at gmail dot com>
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

import math
import pickle
import os
from   tkinter import *
from   tkinter.filedialog import *
import Pmw
from   zooming import Zooming
import AppShell
from   toolbarbutton import ToolBarButton
import tkrpncalc
import entities
import pprint

version = '0.5.0'
date = 'Sept 2, 2019'

geomcolor = 'white'     # color of geometry entities
constrcolor = 'magenta' # color of construction entities
textcolor = 'cyan'      # text color
dimcolor = 'red'        # dimension color
rubbercolor = 'yellow'  # color of (temporary) rubber elements

#===========================================================================
# 
# Math & geometry utility functions
# 
#===========================================================================

def intersection(cline1, cline2):
    """Return intersection (x,y) of 2 clines expressed as (a,b,c) coeff."""
    a,b,c = cline1
    d,e,f = cline2
    i = b*f-c*e
    j = c*d-a*f
    k = a*e-b*d
    if k:
        return (i/k, j/k)
    else:
        return None

def cnvrt_2pts_to_coef(pt1, pt2):
    """Return (a,b,c) coefficients of cline defined by 2 (x,y) pts."""
    x1, y1 = pt1
    x2, y2 = pt2
    a = y2 - y1
    b = x1 - x2
    c = x2*y1-x1*y2
    return (a, b, c)

def proj_pt_on_line(cline, pt):
    """Return point which is the projection of pt on cline."""
    a, b, c = cline
    x, y = pt
    denom = a**2 + b**2
    if not denom:
        return pt
    xp = (b**2*x - a*b*y -a*c)/denom
    yp = (a**2*y - a*b*x -b*c)/denom
    return (xp, yp)

def pnt_in_box_p(pnt, box):
    '''Point in box predicate: Return True if pnt is in box.'''
    x, y = pnt
    x1, y1, x2, y2 = box
    if x1<x<x2 and y1<y<y2: return True
    else: return False

def midpoint(p1, p2, f=.5):
    """Return point part way (f=.5 by def) between points p1 and p2."""
    return (((p2[0]-p1[0])*f)+p1[0], ((p2[1]-p1[1])*f)+p1[1])

def p2p_dist(p1, p2):
    """Return the distance between two points"""
    x, y = p1
    u, v = p2
    return math.sqrt((x-u)**2 + (y-v)**2)

def p2p_angle(p0, p1):
    """Return angle (degrees) from p0 to p1."""
    return math.atan2(p1[1]-p0[1], p1[0]-p0[0])*180/math.pi

def add_pt(p0, p1):
    return (p0[0]+p1[0], p0[1]+p1[1])

def sub_pt(p0, p1):
    return (p0[0]-p1[0], p0[1]-p1[1])

def line_circ_inters(x1, y1, x2, y2, xc, yc, r):
    '''Return list of intersection pts of line defined by pts x1,y1 and x2,y2
    and circle (cntr xc,yc and radius r).
    Uses algorithm from Paul Bourke's web page.'''
    intpnts = []
    num = (xc - x1)*(x2 - x1) + (yc - y1)*(y2 - y1)
    denom = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)
    if denom == 0:
        return
    u = num / denom
    xp = x1 + u*(x2-x1)
    yp = y1 + u*(y2-y1)

    a = (x2 - x1)**2 + (y2 - y1)**2
    b = 2*((x2-x1)*(x1-xc) + (y2-y1)*(y1-yc))
    c = xc**2+yc**2+x1**2+y1**2-2*(xc*x1+yc*y1)-r**2
    q = b**2 - 4*a*c
    if q == 0:
        intpnts.append((xp, yp))
    elif q:
        u1 = (-b+math.sqrt(abs(q)))/(2*a)
        u2 = (-b-math.sqrt(abs(q)))/(2*a)
        intpnts.append(((x1 + u1*(x2-x1)), (y1 + u1*(y2-y1))))
        intpnts.append(((x1 + u2*(x2-x1)), (y1 + u2*(y2-y1))))
    return intpnts

def circ_circ_inters(x1, y1, r1, x2, y2, r2):
    '''Return list of intersection pts of 2 circles.
    Uses algorithm from Robert S. Wilson's web page.'''
    pts = []
    D = (x2-x1)**2 + (y2-y1)**2
    if not D:
        return pts  # circles have same cntr; no intersection
    try:
        q = math.sqrt(abs(((r1+r2)**2-D)*(D-(r2-r1)**2)))
    except:
        return pts  # circles don't interect
    pts = [((x2+x1)/2+(x2-x1)*(r1**2-r2**2)/(2*D)+(y2-y1)*q/(2*D),
            (y2+y1)/2+(y2-y1)*(r1**2-r2**2)/(2*D)-(x2-x1)*q/(2*D)),
           ((x2+x1)/2+(x2-x1)*(r1**2-r2**2)/(2*D)-(y2-y1)*q/(2*D),
            (y2+y1)/2+(y2-y1)*(r1**2-r2**2)/(2*D)+(x2-x1)*q/(2*D))]
    if same_pt_p(pts[0], pts[1]):
        pts.pop()   # circles are tangent
    return pts

def same_pt_p(p1, p2):
    '''Return True if p1 and p2 are within 1e-10 of each other.'''
    if p2p_dist(p1, p2) < 1e-6:
        return True
    else:
        return False

def cline_box_intrsctn(cline, box):
    """Return tuple of pts where line intersects edges of box."""
    x0, y0, x1, y1 = box
    pts = []
    segments = [((x0, y0), (x1, y0)),
                ((x1, y0), (x1, y1)),
                ((x1, y1), (x0, y1)),
                ((x0, y1), (x0, y0))]
    for seg in segments:
        pt = intersection(cline, cnvrt_2pts_to_coef(seg[0], seg[1]))
        if pt:
            if p2p_dist(pt, seg[0]) <= p2p_dist(seg[0], seg[1]) and \
               p2p_dist(pt, seg[1]) <= p2p_dist(seg[0], seg[1]):
                if pt not in pts:
                    pts.append(pt)
    return tuple(pts)

def para_line(cline, pt):
    """Return coeff of newline thru pt and parallel to cline."""
    a, b, c = cline
    x, y = pt
    cnew = -(a*x + b*y)
    return (a, b, cnew)

def para_lines(cline, d):
    """Return 2 parallel lines straddling line, offset d."""
    a, b, c = cline
    c1 = math.sqrt(a**2 + b**2)*d
    cline1 = (a, b, c + c1)
    cline2 = (a, b, c - c1)
    return (cline1, cline2)

def perp_line(cline, pt):
    """Return coeff of newline thru pt and perpend to cline."""
    a, b, c = cline
    x, y = pt
    cnew = a*y - b*x
    return (b, -a, cnew)

def closer(p0, p1, p2):
    """Return closer of p1 or p2 to point p0."""
    d1 = (p1[0] - p0[0])**2 + (p1[1] - p0[1])**2
    d2 = (p2[0] - p0[0])**2 + (p2[1] - p0[1])**2
    if d1 < d2: return p1
    else: return p2

def farther(p0, p1, p2):
    """Return farther of p1 or p2 from point p0."""
    d1 = (p1[0] - p0[0])**2 + (p1[1] - p0[1])**2
    d2 = (p2[0] - p0[0])**2 + (p2[1] - p0[1])**2
    if d1 > d2: return p1
    else: return p2

def find_fillet_pts(r, commonpt, end1, end2):
    """Return ctr of fillet (radius r) and tangent pts for corner
    defined by a common pt, and two adjacent corner pts."""
    line1 = cnvrt_2pts_to_coef(commonpt, end1)
    line2 = cnvrt_2pts_to_coef(commonpt, end2)
    # find 'interior' clines
    cl1a, cl1b = para_lines(line1, r)
    p2a = proj_pt_on_line(cl1a, end2)
    p2b = proj_pt_on_line(cl1b, end2)
    da = p2p_dist(p2a, end2)
    db = p2p_dist(p2b, end2)
    if da <= db: cl1 = cl1a
    else: cl1 = cl1b
    cl2a, cl2b = para_lines(line2, r)
    p1a = proj_pt_on_line(cl2a, end1)
    p1b = proj_pt_on_line(cl2b, end1)
    da = p2p_dist(p1a, end1)
    db = p2p_dist(p1b, end1)
    if da <= db: cl2 = cl2a
    else: cl2 = cl2b
    pc = intersection(cl1, cl2)
    p1 = proj_pt_on_line(line1, pc)
    p2 = proj_pt_on_line(line2, pc)
    return (pc, p1, p2)

def find_common_pt(apair, bpair):
    """Return (common pt, other pt from a, other pt from b), where a and b
    are coordinate pt pairs in (p1, p2) format."""
    p0, p1 = apair
    p2, p3 = bpair
    if same_pt_p(p0, p2):
        cp = p0     # common pt
        opa = p1    # other pt a
        opb = p3    # other pt b
    elif same_pt_p(p0, p3):
        cp = p0
        opa = p1
        opb = p2
    elif same_pt_p(p1, p2):
        cp = p1
        opa = p0
        opb = p3
    elif same_pt_p(p1, p3):
        cp = p1
        opa = p0
        opb = p2
    else:
        return
    return (cp, opa, opb)

def cr_from_3p(p1, p2, p3):
    """Return ctr pt and radius of circle on which 3 pts reside.
    From Paul Bourke's web page."""
    chord1 = cnvrt_2pts_to_coef(p1, p2)
    chord2 = cnvrt_2pts_to_coef(p2, p3)
    radial_line1 = perp_line(chord1, midpoint(p1, p2))
    radial_line2 = perp_line(chord2, midpoint(p2, p3))
    ctr = intersection(radial_line1, radial_line2)
    if ctr:
        radius =  p2p_dist(p1, ctr)
        return (ctr, radius)

def extendline(p0, p1, d):
    """Return point which lies on extension of line segment p0-p1,
    beyond p1 by distance d."""
    pts = line_circ_inters(p0[0], p0[1], p1[0], p1[1], p1[0], p1[1], d)
    if pts:
        return farther(p0, pts[0], pts[1])
    else:
        return

def shortenline(p0, p1, d):
    """Return point which lies on line segment p0-p1,
    short of p1 by distance d."""
    pts = line_circ_inters(p0[0], p0[1], p1[0], p1[1], p1[0], p1[1], d)
    if pts:
        return closer(p0, pts[0], pts[1])
    else:
        return

def line_tan_to_circ(circ, p):
    """Return tan pts on circ of line through p."""
    c, r = circ
    d = p2p_dist(c, p)
    ang0 = p2p_angle(c, p)*math.pi/180
    theta = math.asin(r/d)
    ang1 = ang0+math.pi/2-theta
    ang2 = ang0-math.pi/2+theta
    p1 = (c[0]+(r*math.cos(ang1)), c[1]+(r*math.sin(ang1)))
    p2 = (c[0]+(r*math.cos(ang2)), c[1]+(r*math.sin(ang2)))
    return (p1, p2)

def line_tan_to_2circs(circ1, circ2):
    """Return tangent pts on line tangent to 2 circles.
    Order of circle picks determines which tangent line."""
    c1, r1 = circ1
    c2, r2 = circ2
    d = p2p_dist(c1, c2)    # distance between centers
    ang_loc = p2p_angle(c2, c1)*math.pi/180  # angle of line of centers
    f = (r2/r1-1)/d # reciprocal dist from c1 to intersection of loc & tan line
    theta = math.asin(r1*f)    # angle between loc and tangent line
    ang1 = (ang_loc + math.pi/2 - theta)
    ang2 = (ang_loc - math.pi/2 + theta)
    p1 = (c1[0]+(r1*math.cos(ang1)), c1[1]+(r1*math.sin(ang1)))
    p2 = (c2[0]+(r2*math.cos(ang1)), c2[1]+(r2*math.sin(ang1)))
    return (p1, p2)

def angled_cline(pt, angle):
    """Return cline through pt at angle (degrees)"""
    ang = angle * math.pi / 180
    dx = math.cos(ang)
    dy = math.sin(ang)
    p2 = (pt[0]+dx, pt[1]+dy)
    cline = cnvrt_2pts_to_coef(pt, p2)
    return cline

def ang_bisector(p0, p1, p2, f=0.5):
    """Return cline coefficients of line through vertex p0, factor=f
    between p1 and p2."""
    ang1 = math.atan2(p1[1]-p0[1], p1[0]-p0[0])
    ang2 = math.atan2(p2[1]-p0[1], p2[0]-p0[0])
    deltang = ang2 - ang1
    ang3 = (f * deltang + ang1)*180/math.pi
    return angled_cline(p0, ang3)


def pt_on_RHS_p(pt, p0, p1):
    """Return True if pt is on right hand side going from p0 to p1."""
    angline = p2p_angle(p0, p1)
    angpt = p2p_angle(p0, pt)
    if angline >= 0:
        if angline > angpt > angline-180:
            return True
    else:
        angline += 360
        if angpt < 0:
            angpt += 360
        if angline > angpt > angline-180:
            return True

def rotate_pt(pt, ang, ctr):
    """Return coordinates of pt rotated ang (deg) CCW about ctr.

    This is a 3-step process:
    1. translate to place ctr at origin.
    2. rotate about origin (CCW version of Paul Bourke's algorithm.
    3. apply inverse translation. """
    x, y = sub_pt(pt, ctr)
    A = ang * math.pi / 180
    u = x * math.cos(A) - y * math.sin(A)
    v = y * math.cos(A) + x * math.sin(A)
    return add_pt((u, v), ctr)
    

#===========================================================================

class Draw(AppShell.AppShell):
    """A 2D CAD application using the Tkinter canvas. The canvas is wrapped
    by 'Zooming', (slightly modified) which adds a 'world' coordinate system
    and smooth, mouse controlled zoom (ctrl-RMB) and pan (ctrl-LMB).
    The framework for the application inherits from John Grayson's AppShell
    PMW megawidget, modified slightly. 

    Here's how it works:
    All CAD operations are initiated through a dispatch method, which after
    first initalizing things, and saving the name of the operation as self.op,
    then calls the method (whose name is saved in self.op). Within the
    operation method, the selection mode is set, determining what types of
    data (points or canvas items) are needed from the user and an appropriate
    message prompt is displayed at the bottom of the application window.
    The user then follows the instructions of the message prompt, and clicks
    the mouse on the screen, or enters data using the keyboard or calculator.
    Event handlers detect user input, save the data onto the appropriate
    stack (point_stack, float_stack, or object_stack) and then call the
    'self.op' method again. Some operations, may allow items to be "box
    selected" or assembled into a list. If an operation allows a list of
    items to be selected, the RMB popup menu will include 2 additional
    buttons: "Start list" and "End list".
    If the operation wants to show a hypothetical result, such as a 'rubber
    line', the mouse_motion event passes a screen coordinate to the method as
    an argument. When all the needed data have been entered and stored in
    the appropriate stack, the method pops the data off the stacks and
    completes the operation. When the user wants to quit the current
    operation, he can click the MMB (which calls the end() method), or he
    can just click on another operation, (which causes the dispatch method
    to run again, which in turn calls the end() method as part of the
    initialization sequence).

    3 Coordinate systems:
    The tk Canvas has its own Coordinate System (CCS) with 0,0 in the top left
    corner and increasing Y values going down. In order to facilitate zooming
    and panning, the canvas is wrapped by 'Zooming', which introduces a World
    Coordinate System (WCS) which has the benefit of remaining invariant in
    size, but is still inverted, like the CCS. For CAD, it is conventional
    to work in an environment where positive Y values go up. Therefore,
    an Engineering Coordinate System (ECS) is introduced, which is an X-axis
    reflection of the WCS. The ECS has a 1:1 relation to the CAD  model
    (in millimeters). Wherever possible, calculations are done in ECS,
    converting to or from canvas units as needed. Working in the WCS is
    discouraged, because the negative Y-values and negative angles can cause
    a lot of confusion, especially when calculating angles.
    
    Keeping track of items on the canvas:
    Drawing entities of various types, such as construction lines, geometry
    lines, circles, etc are encapsulated in entity objects which save their
    types, coordinates, etc as attributes. These objects are stored as values
    in a dictionary, accessible by a unique key (integer) which has been
    assigned by the tk canvas. A single dictionary contains all the drawing
    entities in the current display.
    
    File save/load:
    For the purpose of being able to save and load the current display to file,
    these entity objects are disassembled and saved as individual dictionaries
    of exactly one key/value pair. The key is the entity type and the value is
    a tuple containing all the other attributes of the entity.
    When it comes time to reload the data, the original objects are first
    reassembled (the key determines the type of entity object to create and the
    tuple of attributes is supplied as the parameter), then the new objects are
    submitted to type-specific generator methods which display them on the
    canvas and use the canvas generated handle to rebuild the dictionary of
    displayed entities. 

    Calculator:
    An RPN calculator can be launched by running many of the "Measure"
    functions. When measurements are made, values are sent to the x-register
    of the calculator. Also, if an operation is prompting the user to enter a
    float value, the buttons to the left of the calculator registers will send
    the associated value to the CAD operation."""
    
    usecommandarea  = 0
    appversion      = version
    appname         = 'CADvas'
    copyright       = 'Copyright GPL %s' % date
    contactname     = 'Doug Blanding'
    contactemail    = 'dblanding%sgmail%scom' % ('@', '.')
    frameWidth      = 840
    frameHeight     = 600
    catchCntr = False
    catch_pnt = None    # ID of (temporary) catch point
    catch_radius = 5    # radius of catch region
    catch_pnt_size = 5  # size of displayed catch point
    rubber = None       # ID of (temporary) rubber element
    rtext = None        # ID of (temporary) rubber text
    sel_boxID = None    # ID of (temporary) selection box
    op = ''             # current CAD operation (create or modify)
    op_stack = []
    text_entry_enable = 0
    text = ''
    curr = {}           # all entities in curr dwg {k=handle: v=entity}
    prev = {}
    allow_list = 0      # enable/disable item selection in list mode
    sel_mode = ''       # selection mode for screen picks
    float_stack = []    # float values (unitless)
    pt_stack = []       # points, in ECS (mm) units
    obj_stack = []      # canvas items picked from the screen
    sel_box_crnr = None # first corner of selection box, if any
    undo_stack = []     # list of dicts of sets of entities
    redo_stack = []     # data popped off undo_stack
    filename = None     # name of file currently loaded (or saved as)
    dimgap = 10         # extension line gap (in canvas units) 
    textsize = 10       # default text size
    textstyle = 'Calibri'   # default text style

    shift_key_advice = ' (Use SHIFT key to select center of element)'
    unit_dict = {'mm': 1.0,
                 'inches': 25.4,
                 'feet': 304.8}
    units = 'mm'
    unitscale = unit_dict[units]
    calculator = None
    popup = None
    
    #=======================================================================
    # Functions for converting between canvas CS and engineering CS
    #=======================================================================

    def ep2cp(self, pt):
        """Convert pt from ECS to CCS."""
        return self.canvas.world2canvas(pt[0], -pt[1])

    def cp2ep(self, pt):
        """Convert pt from CCS to ECS."""
        x, y = self.canvas.canvas2world(pt[0], pt[1])
        return (x, -y)

    #=======================================================================
    # File, View, Units and Measure commands
    #=======================================================================

    def printps(self):
        openfile = None
        ftypes = [('PostScript file', '*.ps'),
                  ('All files', '*')]
        openfile = asksaveasfilename(filetypes=ftypes)
        if openfile:
            outfile = os.path.abspath(openfile)
            self.ipostscript(outfile)

    def ipostscript(self, file='drawing.ps'):
        ps = self.canvas.postscript()
        ps = ps.replace('1.000 1.000 1.000 setrgbcolor',
                        '0.000 0.000 0.000 setrgbcolor')
        fd = open(file, 'w')
        fd.write(ps)
        fd.close()

    def fileOpen(self):
        openfile = None
        ftypes = [('CADvas dwg', '*.pkl'),
                  ('All files', '*')]
        openfile = askopenfilename(filetypes=ftypes,
                                   defaultextension='.pkl')
        if openfile:
            infile = os.path.abspath(openfile)
            self.load(infile)

    def fileImport(self):
        openfile = None
        ftypes = [('DXF format', '*.dxf'),
                  ('All files', '*')]
        openfile = askopenfilename(filetypes=ftypes,
                                   defaultextension='.dxf')
        if openfile:
            infile = os.path.abspath(openfile)
            self.load(infile)

    def fileSave(self):
        openfile = self.filename
        if openfile:
            outfile = os.path.abspath(openfile)
            self.save(outfile)
        else:
            self.fileSaveas()

    def fileSaveas(self):
        ftypes = [('CADvas dwg', '*.pkl'),
                  ('All files', '*')]
        openfile = asksaveasfilename(filetypes=ftypes,
                                     defaultextension='.pkl')
        if openfile:
            self.filename = openfile
            outfile = os.path.abspath(openfile)
            self.save(outfile)

    def fileExport(self):
        ftypes = [('DXF format', '*.dxf'),
                  ('All files', '*')]
        openfile = asksaveasfilename(filetypes=ftypes,
                                     defaultextension='.dxf')
        if openfile:
            outfile = os.path.abspath(openfile)
            self.save(outfile)

    def save(self, file):

        drawlist = []
        for entity in self.curr.values():
            drawlist.append({entity.type: entity.get_attribs()})

        fext = os.path.splitext(file)[-1]
        if fext == '.dxf':
            import dxf
            dxf.native2dxf(drawlist, file)
        elif fext == '.pkl':
            with open(file, 'wb') as f:
                pickle.dump(drawlist, f)
            self.filename = file
        elif not fext:
            print("Please type entire filename, including extension.")
        else:
            print("Save files of type {fext} not supported.")

    def load(self, file):
        """Load CAD data from file.

        Data is saved/loaded as a list of dicts, one dict for each
        drawing entity, {key=entity_type: val=entity_attribs} """
        
        fext = os.path.splitext(file)[-1]
        if fext == '.dxf':
            import dxf
            drawlist = dxf.dxf2native(file)
        elif fext == '.pkl':
            with open(file, 'rb') as f:
                drawlist = pickle.load(f)
            self.filename = file
        else:
            print("Load files of type {fext} not supported.")
        for ent_dict in drawlist:
            if 'cl' in ent_dict:
                attribs = ent_dict['cl']
                e = entities.CL(attribs)
                self.cline_gen(e.coords)  # This method takes coords
            elif 'cc' in ent_dict:
                attribs = ent_dict['cc']
                e = entities.CC(attribs)
                self.cline_gen(e)
            elif 'gl' in ent_dict:
                attribs = ent_dict['gl']
                e = entities.GL(attribs)
                self.gline_gen(e)
            elif 'gc' in ent_dict:
                attribs = ent_dict['gc']
                e = entities.GC(attribs)
                self.gcirc_gen(e)
            elif 'ga' in ent_dict:
                attribs = ent_dict['ga']
                e = entities.GA(attribs)
                self.garc_gen(e)
            elif 'dl' in ent_dict:
                attribs = ent_dict['dl']
                e = entities.DL(attribs)
                self.dim_gen(e)
            elif 'tx' in ent_dict:
                attribs = ent_dict['tx']
                e = entities.TX(attribs)
                self.text_gen(e)
        self.view_fit()
        self.save_delta()  # undo/redo thing

    def close(self):
        self.quit()

    def view_fit(self):
        bbox = self.canvas.bbox('g', 'd', 't')
        if bbox:
            xsize, ysize = bbox[2]-bbox[0], bbox[3]-bbox[1]
            xc, yc = (bbox[2]+bbox[0])/2, (bbox[3]+bbox[1])/2
            w, h = self.canvas.winfo_width(), self.canvas.winfo_height()
            self.canvas.move_can(w/2-xc, h/2-yc)
            wm, hm = .9 * w, .9 * h
            xscale, yscale = wm/float(xsize), hm/float(ysize)
            if xscale > yscale:
                scale = yscale
            else:
                scale = xscale
            self.canvas.scale(w/2, h/2, scale, scale)
            self.regen()

    def regen(self, event=None):
        #self.regen_all_cl()
        #self.regen_all_dims()
        self.regen_all_text()

    def set_units(self, units):
        if units in self.unit_dict.keys():
            self.units = units
            self.unitscale = self.unit_dict.get(units)
            self.unitsDisplay.configure(text="Units: %s" % self.units)
            self.regen_all_dims()

    def meas_dist(self, obj=None):
        """Measure distance between 2 points."""
        self.op = 'meas_dist'
        if not self.pt_stack:
            self.updateMessageBar('Pick 1st point for distance measurement.')
            self.set_sel_mode('pnt')
        elif len(self.pt_stack) == 1:
            self.updateMessageBar('Pick 2nd point for distance measurement.')
        elif len(self.pt_stack) > 1:
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            dist = p2p_dist(p1, p2)/self.unitscale
            self.updateMessageBar('%s %s'%(dist, self.units))
            self.launch_calc()
            self.calculator.putx(dist)

    def itemcoords(self, obj=None):
        """Print coordinates (in ECS) of selected element."""
        if not self.obj_stack:
            self.updateMessageBar('Pick element from drawing.')
            self.set_sel_mode('items')
        elif self.obj_stack:
            elem = self.obj_stack.pop()[0]
            x1, y1, x2, y2 = self.canvas.coords(elem)
            print(self.cp2ep((x1, y1)), self.cp2ep((x2, y2)))

    def itemlength(self, obj=None):
        """Print length (in current units) of selected line, circle, or arc."""
        if not self.obj_stack:
            self.updateMessageBar('Pick element from drawing.')
            self.set_sel_mode('items')
        elif self.obj_stack:
            elem = None
            for item in self.obj_stack.pop():
                if 'g' in self.canvas.gettags(item):
                    elem = self.curr[item]
            length = 0
            if elem:
                if elem.type is 'gl':
                    p1, p2 = elem.coords
                    length = p2p_dist(p1, p2) / self.unitscale
                elif elem.type is 'gc':
                    length = math.pi*2*elem.coords[1]/self.unitscale
                elif elem.type is 'cc':
                    length = math.pi*2*elem.coords[1]/self.unitscale
                elif elem.type is 'ga':
                    pc, r, a0, a1 = elem.coords
                    ang = float(self.canvas.itemcget(item, 'extent'))
                    length = math.pi*r*ang/180/self.unitscale
                if length:
                    self.launch_calc()
                    self.calculator.putx(length)

    def launch_calc(self):
        if not self.calculator:
            self.calculator = tkrpncalc.Calculator(self)
            self.calculator.geometry('+800+50')

    #=======================================================================
    # GUI configuration
    #=======================================================================
 
    def createBase(self):
        self.toolbar = self.createcomponent('toolbar', (), None,
                  Frame, (self.interior(),), background="gray80")
        self.toolbar.pack(fill=X)

        self.canvas = self.createcomponent('canvas', (), None,
                  Zooming, (self.interior(),), background="black")
        self.canvas.pack(side=LEFT, expand=YES, fill=BOTH)
        self.canvas.panbindings()
        self.canvas.zoombindings()
        Widget.bind(self.canvas, "<Motion>", self.mouseMove)
        Widget.bind(self.canvas, "<Button-1>", self.lftClick)
        Widget.bind(self.canvas, "<Button-2>", self.midClick)
        Widget.bind(self.canvas, "<Button-3>", self.rgtClick)
        self.root.bind("<Key>", self.setCC)
        self.root.bind("<KeyRelease>", self.setCC)
        self.root.bind("<Control-B1-ButtonRelease>", self.regen_all_cl)
        self.root.bind("<Control-B3-ButtonRelease>", self.regen)

    def createMenus(self):
        self.menuBar.deletemenuitems('File', 0)
        self.menuBar.addmenuitem('File', 'command', 'Print drawing',
                                 label='Print', command=self.printps)
        self.menuBar.addmenuitem('File', 'command', 'Open drawing',
                                 label='Open...', command=self.fileOpen)
        self.menuBar.addmenuitem('File', 'command', 'Save drawing',
                                 label='Save', command=self.fileSave)
        self.menuBar.addmenuitem('File', 'command', 'Save drawing',
                                 label='SaveAs...', command=self.fileSaveas)
        self.menuBar.addmenuitem('File', 'command', 'Import DXF',
                                 label='Import DXF', command=self.fileImport)
        self.menuBar.addmenuitem('File', 'command', 'Export DXF',
                                 label='Export DXF', command=self.fileExport)
        self.menuBar.addmenuitem('File', 'separator')
        self.menuBar.addmenuitem('File', 'command', 'Exit program',
                                 label='Exit', command=self.quit)
        self.menuBar.addmenu('Edit', 'Undo / Redo')
        self.menuBar.addmenuitem('Edit', 'command', 'Undo',
                                 label='Undo', command=self.undo)
        self.menuBar.addmenuitem('Edit', 'command', 'Redo',
                                 label='Redo', command=self.redo)
        self.menuBar.addmenuitem('Edit', 'command', 'Clear Redo',
                                 label='Clr Redo', command=self.clear_redo)
        self.menuBar.addmenuitem('Edit', 'command', 'Clear Undo',
                                 label='Clr Undo', command=self.clear_undo)
        self.menuBar.addmenu('View', 'View commands')
        self.menuBar.addmenuitem('View', 'command', 'Fit geometry to screen',
                                 label='Fit', command=self.view_fit)
        self.menuBar.addmenu('Units', 'Switch units')
        self.menuBar.addmenuitem('Units', 'command', 'Set units=mm',
                                 label='mm',
                                 command=lambda k='mm': self.set_units(k))
        self.menuBar.addmenuitem('Units', 'command', 'Set units=inches',
                                 label='inches',
                                 command=lambda k='inches': self.set_units(k))
        self.menuBar.addmenuitem('Units', 'command', 'Set units=feet',
                                 label='feet',
                                 command=lambda k='feet': self.set_units(k))
        self.menuBar.addmenu('Measure', 'Measure')
        self.menuBar.addmenuitem('Measure', 'command', 'measure distance',
                                 label='pt-pt distance', command=self.meas_dist)
        self.menuBar.addmenuitem('Measure', 'command', 'print item coords',
                                 label='item coords',
                                 command=lambda k='itemcoords':self.dispatch(k))
        self.menuBar.addmenuitem('Measure', 'command', 'print item length',
                                 label='item length',
                                 command=lambda k='itemlength':self.dispatch(k))
        self.menuBar.addmenuitem('Measure', 'command', 'launch calculator',
                                 label='calculator',
                                 command=self.launch_calc)
        self.menuBar.addmenu('Dimension', 'Dimensions')
        self.menuBar.addmenuitem('Dimension', 'command', 'Horizontal dimension',
                                 label='Dim Horizontal',
                                 command=lambda k='dim_h':self.dispatch(k))
        self.menuBar.addmenuitem('Dimension', 'command', 'Vertical dimension',
                                 label='Dim Vertical',
                                 command=lambda k='dim_v':self.dispatch(k))
        self.menuBar.addmenuitem('Dimension', 'command', 'Parallel dimension',
                                 label='Dim Parallel',
                                 command=lambda k='dim_par':self.dispatch(k))
        self.menuBar.addmenu('Text', 'Text')
        self.menuBar.addmenuitem('Text', 'command', 'Enter text',
                                 label='Create text',
                                 command=lambda k='text_enter':self.dispatch(k))
        self.menuBar.addmenuitem('Text', 'command', 'Move text',
                                 label='Move text',
                                 command=lambda k='text_move':self.dispatch(k))
        self.menuBar.addmenu('Delete', 'Delete drawing elements')
        self.menuBar.addmenuitem('Delete', 'command',
                                 'Delete individual element',
                                 label='Del Element',
                                 command=lambda k='del_el':self.dispatch(k))
        self.menuBar.addmenuitem('Delete', 'separator')
        self.menuBar.addmenuitem('Delete', 'command', 'Delete all construct',
                                 label='All Cons', command=self.del_all_c)
        self.menuBar.addmenuitem('Delete', 'command', 'Delete all geometry',
                                 label='All Geom', command=self.del_all_g)
        self.menuBar.addmenuitem('Delete', 'command', 'Delete all dimensions',
                                 label='All Dims', command=self.del_all_d)
        self.menuBar.addmenuitem('Delete', 'command', 'Delete all text',
                                 label='All Text', command=self.del_all_t)
        self.menuBar.addmenuitem('Delete', 'separator')
        self.menuBar.addmenuitem('Delete', 'command', 'Delete all',
                                 label='Delete All', command=self.del_all)
        self.menuBar.addmenu('Debug', 'Debug')
        self.menuBar.addmenuitem('Debug', 'command', 'Show Curr',
                                 label='Show Curr',
                                 command=lambda k='show_curr':self.dispatch(k))
        self.menuBar.addmenuitem('Debug', 'command', 'Show Prev',
                                 label='Show Prev',
                                 command=lambda k='show_prev':self.dispatch(k))
        self.menuBar.addmenuitem('Debug', 'command', 'Show Undo',
                                 label='Show Undo',
                                 command=lambda k='show_undo':self.dispatch(k))
        self.menuBar.addmenuitem('Debug', 'command', 'Show Redo',
                                 label='Show Redo',
                                 command=lambda k='show_redo':self.dispatch(k))
        

    def createTools(self):
        self.func      = {}
        self.transFunc = {}
        for key, balloon in [
            ('hcl',     'horizontal construction line'),
            ('vcl',     'vertical construction line'),
            ('hvcl',    'horz & vert construction line'),
            ('acl',     'angled construction line'),
            ('clrefang','construction line angled wrt reference'),
            ('abcl',    'angular bisector construction line'),
            ('lbcl',    'linear bisector construction line'),
            ('parcl',   'parallel construction line'),
            ('perpcl',  'perpendicular construction line'),
            ('cltan1',  'construction line tangent to circle'),
            ('cltan2',  'construction line tangent to 2 circles'),
            ('ccirc',   'construction circle'),
            ('cc3p',    'construction circle by 3 pts'),
            ('cccirc',  'concentric construction circle'),
            #('cctan2',  'construction circle with 2 tangent pts'),
            #('cctan3',  'construction circle with 3 tangent pts'),
            ('sep',     ''),
            ('line',    'line'),
            ('poly',    'polyline'),
            ('rect',    'rectangle'),
            ('circ',    'circle'),
            ('arcc2p',  'arc by cntr & 2 points'),
            ('arc3p',   'arc by 3 points'),
            ('slot',    'slot'),
            ('sep',     ''),
            ('split',   'split line'),
            ('join',    'join 2 adjacent lines'),
            ('fillet',  'fillet corner'),
            ('translate', 'translate geometry element by 2 pts'),
            ('rotate',  'rotate geometry element by angle'),
            #('array',   'copy element(s) into an array'),
            #('stretch', 'stretch elements')
            ]:
            if key == 'sep':
                ToolBarButton(self, self.toolbar, 'sep', 'sep.gif',
                              width=10, state='disabled')
            else:
                ToolBarButton(self, self.toolbar, key, '%s.gif' % key,
                              command=self.dispatch,
                              balloonhelp=balloon)

    def dispatch(self, key):
        """Dispatch commands initiated by menubar & toolbar buttons."""
        self.end()
        self.set_sel_mode('pnt')
        self.op = key
        func = 'self.%s()' % self.op
        eval(func)
        self.entry.focus()

    def set_sel_mode(self, mode=''):
        '''Set selection mode and cursor style.
        Selection mode should be controlled by current operation
        in order to determine what is returned from screen picks.'''
        cursordict = {''    :   'top_left_arrow',
                      'pnt' :   'crosshair',
                      'items':  'right_ptr',
                      'list':   'right_ptr'}
        if mode in cursordict.keys():
            self.sel_mode = mode
            self.canvas.config(cursor=cursordict[mode])

    def createInterface(self):
        AppShell.AppShell.createInterface(self)
        self.createMenus()
        self.createBase()
        self.createTools()
        self.canvas.move_can(60,420)    # Put 0,0 near lower left corner

    #=======================================================================
    # Debug Tools
    #=======================================================================

    def show_curr(self):
        pprint.pprint(self.curr)
        self.end()

    def show_prev(self):
        pprint.pprint(self.prev)
        self.end()

    def show_undo(self):
        pprint.pprint(self.undo_stack)
        self.end()
        
    def show_redo(self):
        pprint.pprint(self.redo_stack)
        self.end()

    #=======================================================================
    # Construction
    # construction lines (clines) are "infinite" length lines
    # described by the equation:            ax + by + c = 0
    # they are defined by coefficients:     (a, b, c)
    #
    # circles are defined by coordinates:   (pc, r)
    #=======================================================================

    def cline_gen(self, cline, rubber=0):
        '''Generate infinite length clines .

        cline coords (a,b,c) are in ECS (mm) values.'''
        
        # extend clines 500 canvas units beyond edge of canvas
        w, h = self.canvas.winfo_width(), self.canvas.winfo_height()
        toplft = self.cp2ep((-500, -500))
        botrgt = self.cp2ep((w+500, h+500))
        trimbox = (toplft[0], toplft[1], botrgt[0], botrgt[1])
        endpts = cline_box_intrsctn(cline, trimbox)
        if len(endpts) == 2:
            p1 = self.ep2cp(endpts[0])
            p2 = self.ep2cp(endpts[1])
            if rubber:
                if self.rubber:
                    self.canvas.coords(self.rubber, p1[0], p1[1], p2[0], p2[1])
                else:
                    self.rubber = self.canvas.create_line(p1[0], p1[1],
                                                          p2[0], p2[1],
                                                          fill=constrcolor,
                                                          tags='r')
            else:
                if self.rubber:
                    self.canvas.delete(self.rubber)
                    self.rubber = None
                handle = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1],
                                                 fill=constrcolor, tags='c')
                self.canvas.tag_lower(handle)
                attribs = (cline, constrcolor)
                e = entities.CL(attribs)
                self.curr[handle] = e

    def regen_all_cl(self, event=None):
        """Delete existing clines, remove them from self.curr, and regenerate

        This needs to be done after pan or zoom because the "infinite" length
        clines are not really infinite, they just hang off the edge a bit. So
        when zooming out, new clines need to be generated so they extend over
        the full canvas. Also, when zooming in, some clines are completely off
        the canvas, so we need a way to keep them from getting lost."""
        
        cl_keylist = [k for k, v in self.curr.items() if v.type is 'cl']
        cl_list = [v.coords for v in self.curr.values() if v.type is 'cl']
        for handle in cl_keylist:
            self.canvas.delete(handle)
            del self.curr[handle]
        for cline in cl_list:
            self.cline_gen(cline)

    def ccirc_gen(self, cc, tag='c'):
        """Create constr circle from a CC object. Save to self.curr."""

        coords, color = cc.get_attribs()
        handle = self.circ_gen(coords, color, tag=tag)
        self.curr[handle] = cc
        self.canvas.tag_lower(handle)

    def hcl(self, pnt=None):
        """Create horizontal construction line from one point or y value."""

        message = 'Pick pt or enter value for horizontal constr line'
        message += self.shift_key_advice
        self.updateMessageBar(message)
        proceed = 0
        if self.pt_stack:
            p = self.pt_stack.pop()
            proceed = 1
        elif self.float_stack:
            y = self.float_stack.pop()*self.unitscale
            p = (0, y)
            proceed = 1
        elif pnt:
            p = self.cp2ep(pnt)
            cline = angled_cline(p, 0)
            self.cline_gen(cline, rubber=1)
        if proceed:
            cline = angled_cline(p, 0)
            self.cline_gen(cline)

    def vcl(self, pnt=None):
        """Create vertical construction line from one point or x value."""

        message = 'Pick pt or enter value for vertical constr line'
        message += self.shift_key_advice
        self.updateMessageBar(message)
        proceed = 0
        if self.pt_stack:
            p = self.pt_stack.pop()
            proceed = 1
        elif self.float_stack:
            x = self.float_stack.pop()*self.unitscale
            p = (x, 0)
            proceed = 1
        elif pnt:
            p = self.cp2ep(pnt)
            cline = angled_cline(p, 90)
            self.cline_gen(cline, rubber=1)
        if proceed:
            cline = angled_cline(p, 90)
            self.cline_gen(cline)

    def hvcl(self, pnt=None):
        """Create a horizontal & vertical construction line pair at a point."""

        message = 'Pick pnt or enter coords for vertical & horizontal constr lines'
        message += self.shift_key_advice
        self.updateMessageBar(message)
        if self.pt_stack:
            p = self.pt_stack.pop()
            self.cline_gen(angled_cline(p, 0))
            self.cline_gen(angled_cline(p, 90))

    def acl(self, pnt=None):
        """Create construction line thru a point, at a specified angle."""
        
        if not self.pt_stack:
            message = 'Pick pnt for angled construction line or enter coordinates'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif self.pt_stack and self.float_stack:
            p0 = self.pt_stack[0]
            ang = self.float_stack.pop()
            cline = angled_cline(p0, ang)
            self.cline_gen(cline)
        elif len(self.pt_stack) > 1:
            p0 = self.pt_stack[0]
            p1 = self.pt_stack.pop()
            cline = cnvrt_2pts_to_coef(p0, p1)
            self.cline_gen(cline)
        elif self.pt_stack and not self.float_stack:
            message = 'Specify 2nd point or enter angle in degrees'
            message += self.shift_key_advice
            self.updateMessageBar(message)
            if pnt:
                p0 = self.pt_stack[0]
                p1 = self.cp2ep(pnt)
                ang = p2p_angle(p0, p1)
                cline = angled_cline(p0, ang)
                self.cline_gen(cline, rubber=1)

    def clrefang(self, p3=None):
        """Create a construction line at an angle relative to a reference."""
        
        if not self.pt_stack:
            message = 'Specify point for new construction line'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif not self.float_stack:
            self.updateMessageBar('Enter offset angle in degrees')
        elif len(self.pt_stack) == 1:
            message = 'Pick first point on reference line'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif len(self.pt_stack) == 2:
            message = 'Pick second point on reference line'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif len(self.pt_stack) == 3:
            p3 = self.pt_stack.pop()
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            baseangle = p2p_angle(p2, p3)
            angoffset = self.float_stack.pop()
            ang = baseangle + angoffset
            cline = angled_cline(p1, ang)
            self.cline_gen(cline)

    def abcl(self, pnt=None):
        """Create an angular bisector construction line."""
        
        if not self.float_stack and not self.pt_stack:
            message = 'Enter bisector factor (Default=.5) or specify vertex'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif not self.pt_stack:
            message = 'Specify vertex point'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif len(self.pt_stack) == 1:
            self.updateMessageBar('Specify point on base line')
        elif len(self.pt_stack) == 2:
            self.updateMessageBar('Specify second point')
            if pnt:
                f = .5
                if self.float_stack:
                    f = self.float_stack[-1]
                p2 = self.cp2ep(pnt)
                p1 = self.pt_stack[-1]
                p0 = self.pt_stack[-2]
                cline = ang_bisector(p0, p1, p2, f)
                self.cline_gen(cline, rubber=1)
        elif len(self.pt_stack) == 3:
            f = .5
            if self.float_stack:
                f = self.float_stack[-1]
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            p0 = self.pt_stack.pop()
            cline = ang_bisector(p0, p1, p2, f)
            self.cline_gen(cline)

    def lbcl(self, pnt=None):
        """Create a linear bisector construction line."""
        
        if not self.pt_stack and not self.float_stack:
            message = 'Enter bisector factor (Default=.5) or specify first point'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif not self.pt_stack:
            message = 'Specify first point'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif len(self.pt_stack) == 1:
            message = 'Specify second point'
            message += self.shift_key_advice
            self.updateMessageBar(message)
            if pnt:
                f = .5
                if self.float_stack:
                    f = self.float_stack[-1]
                p2 = self.cp2ep(pnt)
                p1 = self.pt_stack[-1]
                p0 = midpoint(p1, p2, f)
                baseline = cnvrt_2pts_to_coef(p1, p2)
                newline = perp_line(baseline, p0)
                self.cline_gen(newline, rubber=1)
        elif len(self.pt_stack) == 2:
            f = .5
            if self.float_stack:
                f = self.float_stack[-1]
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            p0 = midpoint(p1, p2, f)
            baseline = cnvrt_2pts_to_coef(p1, p2)
            newline = perp_line(baseline, p0)
            self.cline_gen(newline)

    def parcl(self, pnt=None):
        """Create parallel clines in one of two modes:

        1) At a specified offset distance from selected straight element, or
        2) Parallel to a selected straight element through a selected point."""
        
        if not self.obj_stack and not self.float_stack:
            self.updateMessageBar(
                'Pick a straight element or enter an offset distance')
            self.set_sel_mode('items')
        elif self.float_stack:      # mode 1
            if not self.obj_stack:
                self.set_sel_mode('items')
                self.updateMessageBar(
                    'Pick a straight element to be parallel to')
            elif not self.pt_stack:
                self.set_sel_mode('pnt')
                self.updateMessageBar('Pick on (+) side of line')
            else:
                obj = self.obj_stack.pop()
                p = self.pt_stack.pop()
                item = obj[0]
                baseline = (0,0,0)
                if self.canvas.type(item) == 'line':
                    if 'c' in self.canvas.gettags(item):
                        baseline = self.curr[item].coords
                    elif 'g' in self.canvas.gettags(item):
                        p1, p2 = self.curr[item].coords
                        baseline = cnvrt_2pts_to_coef(p1, p2)
                d = self.float_stack[-1]*self.unitscale
                cline1, cline2 = para_lines(baseline, d)
                p1 = proj_pt_on_line(cline1, p)
                p2 = proj_pt_on_line(cline2, p)
                d1 = p2p_dist(p1, p)
                d2 = p2p_dist(p2, p)
                if d1<d2:
                    self.cline_gen(cline1)
                else:
                    self.cline_gen(cline2)
        elif self.obj_stack:        # mode 2
            obj = self.obj_stack[-1]
            if not obj:
                return
            item = obj[0]
            baseline = (0,0,0)
            if self.canvas.type(item) == 'line':
                if 'c' in self.canvas.gettags(item):
                    baseline = self.curr[item].coords
                elif 'g' in self.canvas.gettags(item):
                    p1, p2 = self.curr[item].coords
                    baseline = cnvrt_2pts_to_coef(p1, p2)
            if not self.pt_stack:
                self.set_sel_mode('pnt')
                message = 'Select point for new parallel line'
                message += self.shift_key_advice
                self.updateMessageBar(message)
                if pnt:
                    p = self.cp2ep(pnt)
                    parline = para_line(baseline, p)
                    self.cline_gen(parline, rubber=1) 
            else:
                p = self.pt_stack.pop()
                newline = para_line(baseline, p)
                self.cline_gen(newline)

    def perpcl(self, pnt=None):
        """Create a perpendicular cline through a selected point."""
        
        if not self.obj_stack:
            self.updateMessageBar('Pick line to be perpendicular to')
            self.set_sel_mode('items')
        else:
            message = 'Select point for perpendicular construction'
            message += self.shift_key_advice
            self.updateMessageBar(message)
            self.set_sel_mode('pnt')
            obj = self.obj_stack[0]
            if not obj:
                return
            item = obj[0]
            baseline = (0,0,0)
            if self.canvas.type(item) == 'line':
                if 'c' in self.canvas.gettags(item):
                    baseline = self.curr[item].coords
                elif 'g' in self.canvas.gettags(item):
                    p1, p2 = self.curr[item].coords
                    baseline = cnvrt_2pts_to_coef(p1, p2)
            if self.pt_stack:
                p = self.pt_stack.pop()
                newline = perp_line(baseline, p)
                self.cline_gen(newline)
                self.obj_stack.pop()
            elif pnt:
                p = self.cp2ep(pnt)
                newline = perp_line(baseline, p)
                self.cline_gen(newline, rubber=1)

    def cltan1(self, p1=None):
        '''Create a construction line through a point, tangent to a circle.'''
        
        if not self.obj_stack:
            self.updateMessageBar('Pick circle')
            self.set_sel_mode('items')
        elif self.obj_stack and not self.pt_stack:
            self.updateMessageBar('specify point')
            self.set_sel_mode('pnt')
        elif self.obj_stack and self.pt_stack:
            item = self.obj_stack.pop()[0]
            p = self.pt_stack.pop()
            circ = None
            if self.curr[item].type in ('gc', 'cc'):
                circ = self.curr[item].coords
            if circ:
                p1, p2 = line_tan_to_circ(circ, p)
                cline1 = cnvrt_2pts_to_coef(p1, p)
                cline2 = cnvrt_2pts_to_coef(p2, p)
                self.cline_gen(cline1)
                self.cline_gen(cline2)

    def cltan2(self, p1=None):
        '''Create a construction line tangent to 2 circles.'''
        
        if not self.obj_stack:
            self.updateMessageBar('Pick first circle')
            self.set_sel_mode('items')
        elif len(self.obj_stack) == 1:
            self.updateMessageBar('Pick 2nd circle')
        elif len(self.obj_stack) == 2:
            item1 = self.obj_stack.pop()[0]
            item2 = self.obj_stack.pop()[0]
            circ1 = circ2 = None
            if self.curr[item1].type in ('gc', 'cc'):
                circ1 = self.curr[item1].coords
            if self.curr[item2].type in ('gc', 'cc'):
                circ2 = self.curr[item2].coords
            if circ1 and circ2:
                p1, p2 = line_tan_to_2circs(circ1, circ2)
                cline = cnvrt_2pts_to_coef(p1, p2)
                self.cline_gen(cline)

    def ccirc_gen(self, cc, tag='c'):
        """Create constr circle from a CC object. Save to self.curr."""

        coords, color = cc.get_attribs()
        handle = self.circ_draw(coords, color, tag=tag)
        self.curr[handle] = cc
        self.canvas.tag_lower(handle)

    def ccirc(self, p1=None):
        '''Create a construction circle from center point and
        perimeter point or radius.'''
        
        self.circ(p1=p1, constr=1)

    def cccirc(self, p1=None):
        '''Create a construction circle concentric to an existing circle,
        at a "relative" radius.'''
        
        if not self.obj_stack:
            self.set_sel_mode('items')
            self.updateMessageBar('Select existing circle')
        elif self.obj_stack and not (self.float_stack or self.pt_stack):
            item = self.obj_stack[0][0]
            self.coords = None
            if self.curr[item].type in ('cc', 'gc'):
                self.coords = self.curr[item].coords
            self.set_sel_mode('pnt')
            self.updateMessageBar(
                'Enter relative radius or specify point on new circle')
            if self.coords and p1:
                pc, r0 = self.coords
                ep = self.cp2ep(p1)
                r = p2p_dist(pc, ep)
                self.circ_builder((pc, r), rubber=1)
        elif self.coords and self.float_stack:
            pc, r0 = self.coords
            self.obj_stack.pop()
            r = self.float_stack.pop()*self.unitscale + r0
            self.circ_builder((pc, r), constr=1)
        elif self.coords and self.pt_stack:
            pc, r0 = self.coords
            self.obj_stack.pop()
            p = self.pt_stack.pop()
            r = p2p_dist(pc, p)
            self.circ_builder((pc, r), constr=1)

    def cc3p(self, p3=None):
        """Create a constr circle from 3 pts on circle."""
        
        if not self.pt_stack:
            self.updateMessageBar('Pick first point on circle')
        elif len(self.pt_stack) == 1:
            self.updateMessageBar('Pick second point on circle')
        elif len(self.pt_stack) == 2:
            self.updateMessageBar('Pick third point on circle')
            if p3:
                p3 = self.cp2ep(p3)
                p2 = self.pt_stack[1]
                p1 = self.pt_stack[0]
                tuple = cr_from_3p(p1, p2, p3)
                if tuple:
                    pc, r = tuple
                    self.circ_builder((pc, r,), rubber=1)
        elif len(self.pt_stack) == 3:
            p3 = self.pt_stack.pop()
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            pc, r = cr_from_3p(p1, p2, p3)
            self.circ_builder((pc, r), constr=1)

    #=======================================================================
    # Geometry
    # geometry line parameters are stored in GL objects.
    # geometry lines are finite length segments between 2 pts: p1, p2
    # lines are defined by coordinates:         (p1, p2)
    #
    #=======================================================================

    def line_draw(self, coords, color, arrow=None, tag='g'):
        """Create and display line segment between two pts. Return ID.

        This is a low level method that accesses the canvas directly &
        returns tkid. The caller can save to self.curr if needed."""
        
        p1, p2 = coords
        xa, ya = self.ep2cp(p1)
        xb, yb = self.ep2cp(p2)
        tkid = self.canvas.create_line(xa, ya, xb, yb,
                                       fill=color, tags=tag, arrow=arrow)
        return tkid

    def gline_gen(self, gl):
        """Create line segment from gl object. Store {id: obj} in self.curr.

        This provides access to line_gen using a gl object."""

        coords, color = gl.get_attribs()
        tkid = self.line_draw(coords, color)
        self.curr[tkid] = gl
        
    def line(self, p1=None):
        '''Create line segment between 2 points. Enable 'rubber line' mode'''
        
        rc = rubbercolor
        if not self.pt_stack:
            message = 'Pick start point of line or enter coords'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif self.pt_stack and p1:
            p0 = self.pt_stack[-1]
            x, y = self.ep2cp(p0)   # fixed first point (canvas coords)
            xr, yr = p1             # rubber point (canvas coords)
            x0, y0 = p0             # fixed first point (ECS)
            x1, y1 = self.cp2ep(p1) # rubber point (ECS)
            strcoords = "(%1.3f, %1.3f)" % ((x1-x0)/self.unitscale,
                                            (y1-y0)/self.unitscale)
            if self.rubber:
                self.canvas.coords(self.rubber, x, y, xr, yr)
            else:
                self.rubber = self.canvas.create_line(x, y, xr, yr,
                                                      fill=rc, tags='r')
            if self.rtext:
                self.canvas.delete(self.rtext)
            self.rtext = self.canvas.create_text(xr+20, yr-20,
                                                 text=strcoords,
                                                 fill=textcolor)
            self.updateMessageBar('Specify end point of line')
        elif len(self.pt_stack) > 1:
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            coords = (p1, p2)
            attribs = (coords, geomcolor)
            e = entities.GL(attribs)
            self.gline_gen(e)
            if self.rubber:
                self.canvas.delete(self.rubber)
                self.rubber = None
            if self.rtext:
                self.canvas.delete(self.rtext)
                self.rtext = None

    def poly(self, p1=None):
        '''Create chain of line segments, enabling 'rubber line' mode.'''
        
        if not self.pt_stack:
            self.poly_start_pt = None
            message = 'Pick start point or enter coords'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        elif self.pt_stack and p1:
            if not self.poly_start_pt:
                self.poly_start_pt = self.pt_stack[-1]
            self.line(p1)   # This will generate rubber line
            self.updateMessageBar('Pick next point or enter coords')
        elif len(self.pt_stack) > 1:
            lastpt = self.pt_stack[-1]
            self.line()     # This will pop 2 points off stack
            if not same_pt_p(self.poly_start_pt, lastpt):
                self.pt_stack.append(lastpt)
    
    def rect(self, p2=None):
        '''Generate a rectangle from 2 diagonally opposite corners.'''
        
        rc = rubbercolor
        if not self.pt_stack:
            self.updateMessageBar(
                'Pick first corner of rectangle or enter coords')
        elif len(self.pt_stack) == 1 and p2:
            self.updateMessageBar(
                'Pick opposite corner of rectangle or enter coords')
            p1 = self.pt_stack[0]
            x1, y1 = self.ep2cp(p1)
            x2, y2 = p2
            if self.rubber:
                self.canvas.coords(self.rubber, x1, y1, x2, y2)
            else:
                self.rubber = self.canvas.create_rectangle(x1, y1, x2, y2,
                                                           outline=rc,
                                                           tags='r')
        elif len(self.pt_stack) > 1:
            x2, y2 = self.pt_stack.pop()
            x1, y1 = self.pt_stack.pop()
            a = (x1, y1)
            b = (x2, y1)
            c = (x2, y2)
            d = (x1, y2)
            sides = ((a, b), (b, c), (c, d), (d, a))
            for p in sides:
                coords = (p[0], p[1])
                attribs = (coords, geomcolor)
                e = entities.GL(attribs)
                self.gline_gen(e)
            if self.rubber:
                self.canvas.delete(self.rubber)
                self.rubber = None

    #=======================================================================
    # geometry circle parameters are stored in GC objects.
    # circles are defined by coordinates:       (pc, r)
    #=======================================================================

    def circ_draw(self, coords, color, tag):
        """Draw a circle on the canvas and return the tkid handle.

        This low level method accesses the canvas directly & returns tkid.
        The caller should save handle & entity_obj to self.curr if needed."""

        ctr, rad = coords
        x, y = self.ep2cp(ctr)
        r = self.canvas.w2c_dx(rad)
        handle = self.canvas.create_oval(x-r, y-r, x+r, y+r,
                                         outline=color,
                                         tags=tag)
        return handle

    def gcirc_gen(self, gc, tag='g'):
        """Create geometry circle from a GC object. Save to self.curr."""

        coords, color = gc.get_attribs()
        handle = self.circ_draw(coords, color, tag=tag)
        self.curr[handle] = gc

    def circ_builder(self, coords, rubber=0, constr=0):
        """Create circle at center pc, radius r in engineering (mm) coords.

        Handle rubber circles, construction, and geom circles."""
        
        ctr, rad = coords       # ECS
        x, y = self.ep2cp(ctr)
        r = self.canvas.w2c_dx(rad)
        if rubber:
            color = rubbercolor
            tag = 'r'
            if self.rubber:
                self.canvas.coords(self.rubber, x-r, y-r, x+r, y+r)
            else:
                self.rubber = self.canvas.create_oval(x-r, y-r, x+r, y+r,
                                                      outline=color,
                                                      tags=tag)
        else:
            if constr:  # Constr circle
                attribs = (coords, constrcolor)
                e = entities.CC(attribs)
                self.ccirc_gen(e)
            else:  # geom circle
                attribs = (coords, geomcolor)
                e = entities.GC(attribs)
                self.gcirc_gen(e)
            if self.rubber:
                self.canvas.delete(self.rubber)
                self.rubber = None
            
    def circ(self, p1=None, constr=0):
        '''Create a circle from center pnt and perimeter pnt or radius.'''
        
        finish = 0
        if not self.pt_stack:
            self.updateMessageBar('Pick center of circle or enter coords')
        elif len(self.pt_stack) == 1 and p1 and not self.float_stack:
            self.updateMessageBar('Specify point on circle or radius')
            pc = self.pt_stack[0]
            p1 = self.cp2ep(p1)
            r = p2p_dist(pc, p1)
            coords = (pc, r)
            self.circ_builder(coords, rubber=1)
        elif len(self.pt_stack) > 1:
            p1 = self.pt_stack.pop()
            pc = self.pt_stack.pop()
            r = p2p_dist(pc, p1)
            finish = 1
        elif self.pt_stack and self.float_stack:
            pc = self.pt_stack.pop()
            r = self.float_stack.pop()*self.unitscale
            finish = 1
        if finish:
            self.circ_builder((pc, r), constr=constr)

    #=======================================================================
    # geometry arc parameters are stored in GA objects
    # arcs are defined by coordinates:  (pc, r, a0, a1)
    # where:    pc = (x, y) coords of center point
    #           r = radius
    #           a0 = start angle in degrees
    #           a1 = end angle in degrees
    #=======================================================================

    def garc_gen(self, ga, tag='g'):
        """Create geometry arc from GA object (coords in ECS)

        pc  = arc center pt
        rad = radius of arc center in mm
        a0  = start angle in degrees measured CCW from 3 o'clock position
        a1  = end angle in degrees measured CCW from 3 o'clock position
        """
        coords, color = ga.get_attribs()
        pc, rad, a0, a1 = coords
        ext = a1-a0
        if ext<0:
            ext += 360
        x, y = self.ep2cp(pc)
        r = self.canvas.w2c_dx(rad)
        x1 = x-r
        y1 = y-r
        x2 = x+r
        y2 = y+r
        if tag is 'r':
            if self.rubber:
                self.canvas.coords(self.rubber, x1, y1, x2, y2,)
                self.canvas.itemconfig(self.rubber, start=a0, extent=ext)
            else:
                self.rubber = self.canvas.create_arc(x1, y1, x2, y2,
                                                     start=a0, extent=ext,
                                                     style='arc', tags=tag,
                                                     outline=color)
        else:
            handle = self.canvas.create_arc(x1, y1, x2, y2,
                                            start=a0, extent=ext, style='arc',
                                            outline=color, tags=tag)
            self.curr[handle] = ga

    def arcc2p(self, p2=None):
        """Create an arc from center pt, start pt and end pt."""
        
        if not self.pt_stack:
            self.updateMessageBar('Specify center of arc')
        elif len(self.pt_stack) == 1:
            self.updateMessageBar('Specify start point of arc')
        elif len(self.pt_stack) == 2:
            self.updateMessageBar('Specify end point of arc')
            if p2:
                p2 = self.cp2ep(p2)
                p1 = self.pt_stack[1]
                p0 = self.pt_stack[0]
                r = p2p_dist(p0, p1)
                ang1 = p2p_angle(p0, p1)
                ang2 = p2p_angle(p0, p2)
                coords = (p0, r, ang1, ang2)
                attribs = (coords, rubbercolor)
                e = entities.GA(attribs)
                self.garc_gen(e, tag='r')
        elif len(self.pt_stack) == 3:
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            p0 = self.pt_stack.pop()
            r = p2p_dist(p0, p1)
            ang1 = p2p_angle(p0, p1)
            ang2 = p2p_angle(p0, p2)
            coords = (p0, r, ang1, ang2)
            attribs = (coords, geomcolor)
            e = entities.GA(attribs)
            self.garc_gen(e)

    def arc3p(self, p3=None):
        """Create an arc from start pt, end pt, and 3rd pt on the arc."""
        
        if not self.pt_stack:
            self.updateMessageBar('Specify start of arc')
        elif len(self.pt_stack) == 1:
            self.updateMessageBar('Specify end of arc')
        elif len(self.pt_stack) == 2:
            self.updateMessageBar('Specify point on arc')
            if p3:
                p3 = self.cp2ep(p3)
                p2 = self.pt_stack[1]
                p1 = self.pt_stack[0]
                tuple = cr_from_3p(p1, p2, p3)
                if tuple:   # tuple=None if p1, p2, p3 are colinear
                    pc, r = tuple
                    ang1 = p2p_angle(pc, p1)
                    ang2 = p2p_angle(pc, p2)
                    if not pt_on_RHS_p(p3, p1, p2):
                        ang2, ang1 = ang1, ang2
                    coords = (pc, r, ang1, ang2)
                    attribs = (coords, rubbercolor)
                    e = entities.GA(attribs)
                    self.garc_gen(e, tag='r')
        elif len(self.pt_stack) == 3:
            p3 = self.pt_stack.pop()
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            pc, r = cr_from_3p(p1, p2, p3)
            ang1 = p2p_angle(pc, p1)
            ang2 = p2p_angle(pc, p2)
            if not pt_on_RHS_p(p3, p1, p2):
                ang2, ang1 = ang1, ang2
            coords = (pc, r, ang1, ang2)
            attribs = (coords, geomcolor)
            e = entities.GA(attribs)
            self.garc_gen(e)
            if self.rubber:
                self.canvas.delete(self.rubber)
                self.rubber = None 

    def slot(self, p1=None):
        if not self.pt_stack:
            self.updateMessageBar('Specify first point for slot')
        elif len(self.pt_stack) == 1:
            self.updateMessageBar('Specify second point for slot')
        elif len(self.pt_stack) == 2 and not self.float_stack:
            self.updateMessageBar('Enter slot width')
        elif len(self.pt_stack) == 2 and self.float_stack:
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            w = self.float_stack.pop()*self.unitscale
            baseline = cnvrt_2pts_to_coef(p1, p2)
            crossline1 = perp_line(baseline, p1)
            crossline2 = perp_line(baseline, p2)
            circ1 = (p1, w/2)
            circ2 = (p2, w/2)
            p1e = extendline(p2, p1, w/2)
            paraline1, paraline2 = para_lines(baseline, w/2)
            p1a = intersection(paraline1, crossline1)
            p1b = intersection(paraline2, crossline1)
            self.pt_stack.extend([p1a, p1b, p1e])
            self.arc3p()
            p2e = extendline(p1, p2, w/2)
            p2a = intersection(paraline1, crossline2)
            p2b = intersection(paraline2, crossline2)
            self.pt_stack.extend([p2a, p2b, p2e])
            self.arc3p()
            self.gline_gen(entities.GL(((p1a, p2a), geomcolor)))
            self.gline_gen(entities.GL(((p1b, p2b), geomcolor)))

    #=======================================================================
    # Modify geometry
    #=======================================================================

    def split(self, p1=None):
        """Split 1 line segment into 2 (at a selected point.)"""
        
        if not self.obj_stack:
            self.set_sel_mode('items')
            self.updateMessageBar('Pick straight line to split')
        elif self.obj_stack and not self.pt_stack:
            self.set_sel_mode('pnt')
            message = 'Pick point for split'
            message += self.shift_key_advice
            self.updateMessageBar(message)
        else:
            # When picking a geometry line that overlays a
            # construction line, need to ignore the c-line
            item_tuple = self.obj_stack.pop()
            for item in item_tuple:
                entity = self.curr[item]
                if entity.type is 'gl':
                    line = item
                    p0 = self.pt_stack.pop()
                    (p1, p2), clr = self.curr[line].get_attribs()
                    del self.curr[line]
                    self.canvas.delete(line)
                    self.gline_gen(entities.GL(((p0, p1), geomcolor)))
                    self.gline_gen(entities.GL(((p0, p2), geomcolor)))

    def join(self, p1=None):
        """Join 2 adjacent line segments into 1. """
        
        if not self.obj_stack:
            self.set_sel_mode('items')
            self.updateMessageBar('Pick first line to join')
        elif len(self.obj_stack) == 1:
            self.updateMessageBar('Pick second line to join')
        elif len(self.obj_stack) == 2:
            item2 = self.obj_stack.pop()[0]
            item1 = self.obj_stack.pop()[0]
            for item in (item1, item2):
                if not (self.canvas.type(item) == 'line' and \
                        'g' in self.canvas.gettags(item)):
                    print('Incorrect types of items picked for join')
                    return
            gl1 = self.curr[item1]
            gl2 = self.curr[item2]
            coords1, clr = gl1.get_attribs()
            coords2, clr = gl2.get_attribs()
            pts = find_common_pt(coords1, coords2)
            if pts:
                cp, ep1, ep2 = pts
            else:
                print('No common pt found')
                return
            for item in (item1, item2):
                del self.curr[item]
                self.canvas.delete(item)
            self.gline_gen(entities.GL(((ep1, ep2), geomcolor)))

    def fillet(self, p1=None):
        """Create a fillet of radius r at the common corner of 2 lines."""
        
        if not self.obj_stack and not self.float_stack:
            self.updateMessageBar('Enter radius for fillet')
        elif not self.obj_stack:
            self.set_sel_mode('items')
            self.updateMessageBar('Pick corner to apply fillet')
        elif self.obj_stack and self.float_stack:
            rw = self.float_stack[-1]*self.unitscale
            rc = self.canvas.w2c_dx(rw)
            found = self.obj_stack.pop()
            items = []
            for item in found:
                if self.canvas.type(item) == 'line' and \
                   'g' in self.canvas.gettags(item):
                    items.append(item)
            if len(items) == 2:
                line1coords, color = self.curr[items[0]].get_attribs()
                line2coords, color = self.curr[items[1]].get_attribs()
                pts = find_common_pt(line1coords,
                                     line2coords)
                if pts:
                    # common pt, other end pt1, other end pt2
                    cp, ep1, ep2 = pts
                else:
                    print('No common point found')
                    return
                # find arc center and tangent points
                ctr, tp1, tp2 = find_fillet_pts(rw, cp, ep1, ep2)
                # shorten adjacent sides
                for item in items:
                    del self.curr[item]
                    self.canvas.delete(item)
                self.gline_gen(entities.GL(((ep1, tp1), geomcolor)))
                self.gline_gen(entities.GL(((ep2, tp2), geomcolor)))

                # make arc, but first, get the order of tp1 and tp2 right
                a1 = math.atan2(tp1[1]-ctr[1], tp1[0]-ctr[0])
                a2 = math.atan2(tp2[1]-ctr[1], tp2[0]-ctr[0])
                if (a2-a1) > math.pi or -math.pi < (a2-a1) < 0:
                    tp1, tp2 = tp2, tp1
                self.pt_stack = [ctr, tp1, tp2]
                self.arcc2p()

    def translate(self, p=None):
        """Move (or copy) selected geometry item(s) by two points.

        To copy items, enter number of copies.
        Otherwise, item(s) will be moved (not copied)."""
        
        if not self.obj_stack and not self.pt_stack and \
           not self.float_stack:
            self.set_sel_mode('items')
            self.allow_list = 1
            msg = 'Specify number of copies or select geometry item(s) to move'
            self.updateMessageBar(msg)
        elif not self.obj_stack and not self.pt_stack:
            self.updateMessageBar('Select geometry item(s) to move')
        elif self.obj_stack and not self.pt_stack:
            self.set_sel_mode('pnt')
            self.allow_list = 0
            self.updateMessageBar('Select "FROM" point')
        elif self.obj_stack and len(self.pt_stack) == 1:
            self.updateMessageBar('Select "TO" point')
        elif self.obj_stack and len(self.pt_stack) == 2:
            if self.float_stack:
                repeat = int(self.float_stack.pop())
            else:
                repeat = 0
            p1 = self.pt_stack.pop()
            p0 = self.pt_stack.pop()
            handles = self.obj_stack.pop()
            dp = sub_pt(p1, p0)
            cx, cy = sub_pt(self.ep2cp(p1), self.ep2cp(p0))
            items = [self.curr[handle] for handle in handles]
            delete_original = False
            if not repeat:  # move, (not copy)
                delete_original = True
                repeat = 1
            for item in items:
                if item.type is 'gl':
                    pnts, _ = item.get_attribs()
                    for x in range(repeat):
                        pnts = (add_pt(pnts[0], dp),
                                add_pt(pnts[1], dp))
                        gl = entities.GL((pnts, geomcolor))
                        self.gline_gen(gl)
                elif item.type is 'gc':
                    pnts, _ = item.get_attribs()
                    for x in range(repeat):
                        pnts = (add_pt(pnts[0], dp), pnts[1])
                        gc = entities.GC((pnts, geomcolor))
                        self.gcirc_gen(gc)
                elif item.type is 'ga':
                    pnts, _ = item.get_attribs()
                    for x in range(repeat):
                        pnts = (add_pt(pnts[0], dp),
                                pnts[1], pnts[2], pnts[3])
                        ga = entities.GA((pnts, geomcolor))
                        self.garc_gen(ga)
                else:
                    print('Only geometry type items can be moved with this command.')
            if delete_original:
                for handle in handles:
                    self.canvas.delete(handle)
                    del self.curr[handle]

    def rotate(self, p=None):
        """Move (or copy) selected geometry item(s) by rotating about a point.

        To copy items, enter number of copies.
        Otherwise, item(s) will be moved (not copied)."""
        
        if not self.obj_stack and not self.pt_stack and not self.float_stack:
            self.repeat = 0   # No copies. "move" mode is intended.
            self.set_sel_mode('items')
            self.allow_list = 1
            msg = 'Specify number of copies or select geometry item(s) to move'
            self.updateMessageBar(msg)
        elif not self.obj_stack and not self.pt_stack:
            self.updateMessageBar('Select geometry item(s) to move')
        elif self.obj_stack and not self.pt_stack:
            if self.float_stack:
                self.repeat = int(self.float_stack.pop())   # number of copies
            self.set_sel_mode('pnt')
            self.allow_list = 0
            self.updateMessageBar('Select center of rotation')
        elif self.obj_stack and self.pt_stack and not self.float_stack:
            self.updateMessageBar('Specify angle of rotation in degrees')
        elif self.obj_stack and self.pt_stack and self.float_stack:
            ctr = self.pt_stack.pop()
            handles = self.obj_stack.pop()
            A = self.float_stack.pop()
            items = [self.curr[handle] for handle in handles]
            delete_original = False
            if not self.repeat:  # move, (not copy)
                delete_original = True
                self.repeat = 1
            for item in items:
                if item.type is 'gl':
                    pnts, _ = item.get_attribs()
                    for x in range(self.repeat):
                        pnts = (rotate_pt(pnts[0], A, ctr),
                                rotate_pt(pnts[1], A, ctr))
                        gl = entities.GL((pnts, geomcolor))
                        self.gline_gen(gl)
                elif item.type is 'gc':
                    pnts, _ = item.get_attribs()
                    for x in range(self.repeat):
                        pnts = (rotate_pt(pnts[0], A, ctr), pnts[1])
                        gc = entities.GC((pnts, geomcolor))
                        self.gcirc_gen(gc)
                elif item.type is 'ga':
                    pnts, _ = item.get_attribs()
                    for x in range(self.repeat):
                        pnts = (rotate_pt(pnts[0], A, ctr),
                                pnts[1], pnts[2] + A, pnts[3] + A)
                        ga = entities.GA((pnts, geomcolor))
                        self.garc_gen(ga)
                else:
                    print('Only geometry type items can be moved with this command.')
            if delete_original:
                for handle in handles:
                    self.canvas.delete(handle)
                    del self.curr[handle]
                
    #=======================================================================
    # Dimensions
    # linear dimensions have coords:    (p1, p2, p3, d)
    # where p1 and p2 are the points being dimensioned,
    # d is the direction along which the dimension is being measured,
    # represented by the coefficients of a cline: d = (a, b, c)
    # and p3 is the location of the center of the dimension text.
    #=======================================================================

    def dim_draw(self, dim_obj):
        """Create a linear dimension from dim_obj and return handle.
        
        There are 5 individual components that make up a linear dimension:
        The text, 2 dimension lines, and 2 extension lines. Each component
        shares a tag which is unique to this 'group' of 5 components. This
        permits all components to be found when any component is selected
        on the canvas. It is intended to treat dimensions as 'disposable'.
        For example, to move a dimension, just delete all 5 components,
        then regenerate them in the new position."""

        (p1, p2, p3, c), color = dim_obj.get_attribs()
        dimdir = para_line(c, p3)
        p1b = proj_pt_on_line(dimdir, p1)
        p2b = proj_pt_on_line(dimdir, p2)
        d = p2p_dist(p1b, p2b) / self.unitscale
        text = '%.3f' % d
        x3, y3 = self.ep2cp(p3)
        tkid = self.canvas.create_text(x3, y3, fill=color, text=text)
        dgidtag = 'd%s' % tkid  # unique dimension group ID tag
        self.canvas.itemconfig(tkid, tags=('d', dgidtag))
        # create dimension lines
        xa, ya, xb, yb = self.canvas.bbox(tkid)
        xa, ya = self.cp2ep((xa, ya))
        xb, yb = self.cp2ep((xb, yb))
        innerpts = cline_box_intrsctn(dimdir, (xa, ya, xb, yb))
        ip1 = closer(p1b, innerpts[0], innerpts[1])
        ip2 = closer(p2b, innerpts[0], innerpts[1])
        self.line_draw((ip1, p1b), color=color, tag=('d', dgidtag), arrow=LAST)
        self.line_draw((ip2, p2b), color=color, tag=('d', dgidtag), arrow=LAST)
        # create extension lines
        # make ext line gap appear same size irrespective of zoom
        gap = self.canvas.c2w_dx(self.dimgap)
        p1a = shortenline(p1b, p1, gap)
        p2a = shortenline(p2b, p2, gap)
        p1c = extendline(p1, p1b, gap)
        p2c = extendline(p2, p2b, gap)
        if p1a and p2a and p1c and p2c:
            self.line_draw((p1a, p1c), color=color, tag=('d', dgidtag))
            self.line_draw((p2a, p2c), color=color, tag=('d', dgidtag))
        return dgidtag
        

    def dim_gen(self, dim_obj):
        """Generate dimension from dim_obj and save to self.curr."""
        
        dgid = self.dim_draw(dim_obj)
        self.curr[dgid] = dim_obj

    def regen_all_dims(self, event=None):
        """Delete all existing dimensions, and regenerate.

        This needs to be done after zoom because the dimension text does
        not change size with zoom."""
        
        dimlist = [v.coords for v in self.curr.values() if v.type is 'dl']
        self.del_all_d()
        for coords in dimlist:
            self.dim_gen(coords)

    def dim_lin(self, p=None, d=(0,1,0)):
        """Manually create a linear dimension obj. Add to self.curr."""

        rc = rubbercolor
        if not self.pt_stack:
            self.updateMessageBar('Pick 1st point.')
        elif len(self.pt_stack) == 1:
            self.updateMessageBar('Pick 2nd point.')
        elif len(self.pt_stack) == 2 and p:
            self.updateMessageBar('Pick location for dimension text.')
            p3 = self.cp2ep(p)
            p2 = self.pt_stack[1]
            p1 = self.pt_stack[0]
            if not same_pt_p(p3, p2):
                if self.rubber:
                    for each in self.canvas.find_withtag(self.rubber):
                        self.canvas.delete(each)
                att = ((p1, p2, p3, d), rc)
                rubber_ent = entities.DL(att)
                self.rubber = self.dim_draw(rubber_ent)
        elif len(self.pt_stack) == 3:
            if self.rubber:
                for each in self.canvas.find_withtag(self.rubber):
                    self.canvas.delete(each)
            p3 = self.pt_stack.pop()
            p2 = self.pt_stack.pop()
            p1 = self.pt_stack.pop()
            coords = (p1, p2, p3, d)
            attribs = (coords, dimcolor)
            dl = entities.DL(attribs)
            dgid = self.dim_draw(dl)
            self.curr[dgid] = dl

    def dim_h(self, p=None):
        """Create a horizontal dimension"""
        
        self.dim_lin(p)

    def dim_v(self, p=None):
        """Create a vertical dimension"""
        
        self.dim_lin(p, d=(1,0,0))

    def dim_par(self, p=None):
        """Create a dimension parallel to a selected line element."""
        
        if not self.obj_stack:
            self.set_sel_mode('items')
            self.updateMessageBar(
                'Pick linear element to define direction of dimension.')
        elif self.obj_stack:
            self.set_sel_mode('pnt')
            item = self.obj_stack[-1][0]
            if self.canvas.type(item) == 'line':
                tags = self.canvas.gettags(item)
                d = None
                if 'c' in tags:
                    d = self.curr[item].coords
                elif 'g' in tags:
                    p1, p2 = self.curr[item].coords
                    d = cnvrt_2pts_to_coef(p1, p2)
                if d:
                    self.dim_lin(p, d)

    #=======================================================================
    # Text
    # Text parameters are stored as attributes of a TX object.
    # attribs = (x,y), text, style, size, color
    # where (x, y) are the coordinates of the center of the text.
    # style, size, color define the font.
    #=======================================================================

    def text_draw(self, tx, tag='t'):
        """Draw text on canvas and return handle."""

        x, y = tx.coords
        text = tx.text
        style = tx.style
        size = tx.size
        color =tx. color
        u, v = self.ep2cp((x, y))
        zoom_scale = self.canvas.scl.x
        zoomed_font_size = int(size * zoom_scale)
        font = (style, zoomed_font_size)
        handle = self.canvas.create_text(u, v, text=text, tags=tag,
                                         fill=color, font=font)
        return handle

    def text_gen(self, tx, tag='t'):
        """Generate text from a TX object and save to self.curr."""

        handle = self.text_draw(tx, tag=tag)
        self.curr[handle] = tx

    def regen_all_text(self, event=None):
        """Delete all existing text, clear tx_dict, and regenerate.

        This needs to be done after zoom because text size is defined
        in terms of canvas pixels and doesn't change size with zoom."""
        
        tx_list = [tx for tx in self.curr.values() if tx.type == 'tx']
        attribs_list = [tx.get_attribs() for tx in tx_list]
        self.del_all_t()
        for attribs in attribs_list:
            tx = entities.TX(attribs)
            self.text_gen(tx)

    def text_enter(self, p=None):
        """Place new text on drawing."""
        
        rc = rubbercolor
        if not self.text:
            self.text_entry_enable = 1
            self.updateMessageBar('Enter text')
        elif not self.pt_stack:
            self.updateMessageBar('Pick location for center of text')
            if p:
                x, y = p
                if self.rubber:
                    self.canvas.delete(self.rubber)
                self.rubber = self.canvas.create_text(x, y, text=self.text,
                                                      fill=rc, tags='r')
        elif self.pt_stack:
            p = self.pt_stack.pop()
            attribs = (p, self.text, self.textstyle,
                       self.textsize, textcolor)
            tx = entities.TX(attribs)
            self.text_gen(tx)
            self.text = None
            if self.rubber:
                self.canvas.delete(self.rubber)
            if self.op_stack:
                self.op = self.op_stack.pop()

    def text_move(self, p=None):
        """Move existing text to new point."""

        if not self.obj_stack:
            self.set_sel_mode('items')
            self.updateMessageBar('Select text to move.')
        elif not self.pt_stack:
            if not self.rubber:
                for item_tuple in self.obj_stack:
                    for item in item_tuple:
                        if 't' in self.canvas.gettags(item):
                            if item in self.curr:
                                old_tx = self.curr[item]
                                old_attribs = old_tx.get_attribs()
                                self.rubber_tx = entities.TX(old_attribs)
                                self.rubber = self.text_draw(self.rubber_tx,
                                                             tag='r')
            if self.rubber:
                self.canvas.delete(self.rubber)
            if p:  # cursor coordinates supplied by mouseMove
                p = self.cp2ep(p)  # coords of p are in CCS
                self.rubber_tx.coords = p
            self.rubber = self.text_draw(self.rubber_tx, tag='r')
            self.updateMessageBar('Pick new location for center of text')
            self.set_sel_mode('pnt')
        elif self.pt_stack:
            newpoint = self.pt_stack.pop()
            handle = self.obj_stack.pop()[0]
            if handle in self.curr:
                tx = self.curr[handle]
                attribs = list(tx.get_attribs())
                attribs[0] = newpoint
                attribs = tuple(attribs)
                new_tx = entities.TX(attribs)
                self.text_gen(new_tx)
                del self.curr[handle]
                self.canvas.delete(handle)
            if self.rubber:
                self.canvas.delete(self.rubber)
                self.rubber = None
                del self.rubber_tx
            self.regen_all_text()
            
    #=======================================================================
    # Delete
    #=======================================================================

    def del_el(self, item_tuple=None):
        '''Delete individual elements.'''
        
        self.set_sel_mode('items')
        self.allow_list = 1
        self.updateMessageBar('Pick element(s) to delete.')
        if self.obj_stack:
            item_tuple = self.obj_stack.pop()
            for item in item_tuple:
                if item in self.curr:
                    del self.curr[item]
                    self.canvas.delete(item)
                else:
                    tags = self.canvas.gettags(item)
                    if 'd' in tags:
                        dgid = tags[1]
                        dim_items = self.canvas.find_withtag(dgid)
                        for each in dim_items:
                            self.canvas.delete(each)
                        del self.curr[dgid]
             
    def del_all_c(self):
        '''Delete All construction.'''
        
        delete = [k for k, v in self.curr.items() if v.type in ('cl', 'cc')]
        for k in delete:
            del self.curr[k]
        for item in self.canvas.find_withtag('c'):
            self.canvas.delete(item)
             
    def del_all_g(self):
        '''Delete all geometry.'''
        
        delete = [k for k, v in self.curr.items() if v.type is 'gl']
        for k in delete:
            del self.curr[k]
        delete = [k for k, v in self.curr.items() if v.type is 'gc']
        for k in delete:
            del self.curr[k]
        delete = [k for k, v in self.curr.items() if v.type is 'ga']
        for k in delete:
            del self.curr[k]
        for item in self.canvas.find_withtag('g'):
            self.canvas.delete(item)

    def del_all_d(self):
        '''Delete all dimensions.'''
        
        delete = [k for k, v in self.curr.items() if v.type is 'dl']
        for k in delete:
            del self.curr[k]
        for item in self.canvas.find_withtag('d'):
            self.canvas.delete(item)

    def del_all_t(self):
        '''Delete all text.'''
        
        delete = [k for k, v in self.curr.items() if v.type is 'tx']
        for k in delete:
            del self.curr[k]
        for item in self.canvas.find_withtag('t'):
            self.canvas.delete(item)

    def del_all(self):
        '''Delete all.'''
        
        self.curr.clear()
        self.canvas.delete(ALL)
        

    #=======================================================================
    # Undo / Redo
    #=======================================================================

    """
    When drawing entities are generated and displayed, their parameters are
    stored in objects that are specific to their 'type'. The objects which
    encapsulate them each have a .type attribute mirroring the type of the
    entity being encapsulated. The types are as follows:
    
    'cl'    construction line
    'cc'    construction circle
    'gl'    geometry line
    'gc'    geometry circle
    'ga'    geometry arc
    'dl'    linear dimension
    'tx'    text

    Information about all the entities currently in the drawing is kept in a
    dictionary named self.curr, whose values are the objects encapsulating each
    entity and whose keys are the canvas generated handles associated with each
    entity.
    In order to implement undo and redo, it is neccesary to detect whenever
    there is a change in self.curr. To do this, a copy of self.curr (named
    self.prev) is maintained. Whenever a CAD operation ends, the save_delta()
    method is called. This method first compares self.curr with self.prev to
    see if they are equal. If not, a set containing the values in self.curr is
    compared with a set containing the values in self.prev. The difference is
    loaded onto the undo_stack. The curr config is then copied to self.prev.
                             __________
                            |  Change  |
                            |_detected_|
                                 ||
                                 ||1
                                 \/          2
     ____________            __________    diff    ______________
    | redo stack |          |   Curr   |    -->   |  Undo stack  |
    |____________|          |__________|          |______________|
                                 ||
                                 ||3
                                 \/
                             __________
                            |   Prev   |
                            |__________|

    1. difference detected between curr and prev.
    2. diff (delta) pushed onto undo_stack.
    3. copy of curr saved to prev.
    
    
    The undo & redo buttons work as shown in the diagram below.

     ____________     2      __________ 3       1  ______________
    | redo stack |   <--    |   Curr   |    <--   |  Undo stack  |
    |____________|          |__________|          |______________|
                                 ||
                                 ||4
                                 \/
                             __________
                            |   Prev   |
                            |__________|

    For example, when the Undo button is clicked:
    1. undo_data is popped off the undo_stack.
    2. undo data is pushed onto the redo_stack.
    3. curr is updated with undo_data.
    4. copy of curr is save to prev.


     ____________ 1       3  __________      2     ______________
    | redo stack |   -->    |   Curr   |    -->   |  Undo stack  |
    |____________|          |__________|          |______________|
                                 ||
                                 ||4
                                 \/
                             __________
                            |   Prev   |
                            |__________|

    Similarly, if the Redo button is clicked:
    1. redo_data is popped off the redo_stack.
    2. redo data is pushed onto the undo_stack.
    3. curr is updated with redo_data.
    4. copy of curr is saved to prev.

    Typically, after clicking undo / redo buttons one or more times,
    the user will resume running operations that create, modify or
    delete CAD data. Once these operations resume, the data on the redo
    stack will no longer be relevant and needs to be discarded.
    For now, there is a button under the edit menu that allows the user
    to manually clear the redo stack.

    ToDo: Figure out how best to clear the redo_stack automatically.
    """

    def save_delta(self):
        """After a drawing change, save deltas on undo stack."""
        
        if self.curr != self.prev:
            plus = set(self.curr.values()) - set(self.prev.values())
            minus = set(self.prev.values()) - set(self.curr.values())
            delta = {'+': plus, '-': minus}
            self.undo_stack.append(delta)
            self.prev = self.curr.copy()

    def undo(self):
        """Pop data off undo, push onto redo, and update curr."""
        
        if self.undo_stack:
            undo_data = self.undo_stack.pop()
            self.redo_stack.append(undo_data)
            for item in undo_data['+']:
                self.rem_draw(item)
            for item in undo_data['-']:
                self.add_draw(item)
            self.prev = self.curr.copy()
        else:
            print("No more Undo steps available.")

    def redo(self):
        """Pop data off redo, push onto undo, and update curr."""

        if self.redo_stack:
            redo_data = self.redo_stack.pop()
            self.undo_stack.append(redo_data)
            for item in redo_data['+']:
                self.add_draw(item)
            for item in redo_data['-']:
                self.rem_draw(item)
            self.prev = self.curr.copy()
        else:
            print("No more Redo steps available.")

    def add_draw(self, entity):
        """Add entity to current drawing."""

        if entity.type is 'cl':
            self.cline_gen(entity.coords)  # This one takes coords
        elif entity.type is 'cc':
            self.ccirc_gen(entity)
        elif entity.type is 'gl':
            self.gline_gen(entity)
        elif entity.type is 'gc':
            self.gcirc_gen(entity)
        elif entity.type is 'ga':
            self.garc_gen(entity)
        elif entity.type is 'dl':
            self.dim_gen(entity)
        elif entity.type is 'tx':
            self.text_gen(entity)
        

    def rem_draw(self, entity):
        """Remove entity from current drawing."""

        kvlist = list(self.curr.items())
        for k, v in kvlist:
            if v == entity:
                self.canvas.delete(k)
                del self.curr[k]

    def clear_redo(self):
        self.redo_stack.clear()

    def clear_undo(self):
        self.undo_stack.clear()

    #=======================================================================
    # Event handling
    #=======================================================================

    def end(self):
        '''End current operation'''
        
        if self.rubber:
            self.canvas.delete(self.rubber)
            self.rubber = None
        if self.rtext:
            self.canvas.delete(self.rtext)
            self.rtext = None
        if self.catch_pnt:
            self.canvas.delete(self.catch_pnt)
            self.catch_pnt = None
        if self.op:
            self.op = ''
        self.sel_box_crnr = None
        self.canvas.delete(self.sel_boxID)
        self.sel_boxID = None
        self.text = ''
        self.pt_stack = []
        self.float_stack = []
        self.obj_stack = []
        self.text_entry_enable = 0
        self.set_sel_mode('')
        self.allow_list = 0
        self.quitpopup()
        self.save_delta()
        self.updateMessageBar('CTRL-LMB to pan.  CTRL-RMB to zoom.')

    def enterfloat(self, str_value):
        """Receive string value (from calculator) and do the right thing."""
        
        if str_value:
            val = float(str_value)
            self.float_stack.append(val)
            func = 'self.%s()' % self.op
            eval(func)

    def keybrdEntry(self, event):
        """Store user entered values on stack.

        POINTS:
        points are stored in mm units in ECS on self.pt_stack.
        This is one of the places where unit scale is applied.

        FLOATS:
        floats are stored as unitless numbers on self.float_stack. Because a
        float value may be used for anything: radius, angle, x value, y value,
        whatever; it is not possible to know here how a float value will
        be used. It remains the responsibility of the using function to
        condition the float value appropriately by applying unitscale for
        distances, etc."""
        
        if self.op:
            text = self.entry.get()
            self.entry.delete(0, len(text))
            if self.text_entry_enable:
                self.text = text
            else:
                list = text.split(',')
                if len(list) == 1:
                    val = list[0]
                    self.float_stack.append(float(val))
                elif len(list) == 2 and self.sel_mode == 'pnt':
                    # user entered points are already in ECS units
                    x, y = list
                    x = float(x) * self.unitscale
                    y = float(y) * self.unitscale
                    self.pt_stack.append((x, y))
            func = 'self.%s()' % self.op
            eval(func)

    def lftClick(self, event):
        '''Place screen picks on appropriate stack, call method named by self.op.

        In "point" mode, put x,y coords of "catch point", if any, on point
        stack, otherwise put pointer x,y coords on stack.
        In "items" mode, put a tuple of selected items on "object stack".
        If first click does not find one or more items within its
        "catch radius", enter "box select mode" and look for objects that
        lie completely inside box defined by 1st and 2nd clicks.
        '''
        
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        cr = self.catch_radius
        if self.sel_mode == 'pnt':
            # convert screen coords to ECS units and put on pt_stack
            if self.catch_pnt:
                l, t, r, b = self.canvas.coords(self.catch_pnt)
                x = (r + l)/2
                y = (t + b)/2
            p = self.cp2ep((x, y))
            self.pt_stack.append(p)
            func = 'self.%s()' % self.op
            eval(func)
        elif self.sel_mode in ('items', 'list'):
            items = self.canvas.find_overlapping(x-cr, y-cr, x+cr, y+cr)
            if not items and not self.sel_box_crnr:
                self.sel_box_crnr = (x, y)
                return
            elif self.sel_box_crnr:
                x1, y1 = self.sel_box_crnr
                items = self.canvas.find_enclosed(x1, y1, x, y)
                self.sel_box_crnr = None
                self.canvas.delete(self.sel_boxID)
                self.sel_boxID = None
            if self.sel_mode == 'items':
                self.obj_stack.append(items)
                func = 'self.%s()' % self.op
                eval(func)
            elif self.sel_mode == 'list':
                if not self.obj_stack:
                    self.obj_stack.append([])
                for item in items:
                    if item not in self.obj_stack[-1]:
                        self.obj_stack[-1].append(item)

    def midClick(self, event):
        self.end()

    def rgtClick(self, event):
        '''Popup menu for view options.'''
        
        if self.popup:
            self.popup.destroy()
        self.popup = Toplevel()
        self.popup.overrideredirect(1)
        frame = Frame(self.popup)
        Button(frame, text='View Fit',
               command=lambda:(self.view_fit(), self.quitpopup())).pack()
        if self.allow_list:
            Button(frame, text='Start list',
                   command=lambda:(self.set_sel_mode('list'), self.quitpopup())).pack()
            Button(frame, text='End list',
                   command=lambda:(self.set_sel_mode('items'), eval('self.%s()' % self.op),
                                   self.quitpopup())).pack()
        frame.pack()
        size, x, y = self.winfo_toplevel().winfo_geometry().split('+')
        x = int(x)
        y = int(y)
        if self.allow_list:
            self.popup.geometry('60x90+%s+%s' % (x+event.x, y+event.y+30))
        else:
            self.popup.geometry('60x30+%s+%s' % (x+event.x, y+event.y+30))

    def quitpopup(self):
        if self.popup:
            self.popup.destroy()
            self.popup = None

    def genCatchPnt(self, x, y, color='yellow', regen=0):
        '''Generate (or regenerate) a catch point at coordinates x, y.'''
        
        ps = self.catch_pnt_size
        if regen:
            self.canvas.coords(self.catch_pnt, x-ps, y-ps, x+ps, y+ps)
        else:
            self.catch_pnt = self.canvas.create_rectangle(x-ps, y-ps,
                                                          x+ps, y+ps,
                                                          outline=color)

    def setCC(self, event):
        '''Set center catch flag'''
        
        if event.type == '2' and event.keysym == 'Shift_L':
            self.catchCntr = True
        else:
            self.catchCntr = False

    def mouseMove(self, event):
        '''Display a catch point (ID=self.catch_pnt) on a line within
        self.catch_radius of the cursor. Catch point should be "sticky"
        at midpoints, ends and intersections.'''
        
        x = self.canvas.canvasx(event.x)
        y = self.canvas.canvasy(event.y)
        
        if self.sel_mode == 'pnt':
            cr = self.catch_radius
            found = self.canvas.find_overlapping(x-cr, y-cr, x+cr, y+cr)
            items = []
            for each in found:
                if self.canvas.type(each) in ('line', 'oval', 'arc') and\
                   'r' not in self.canvas.gettags(each):
                    items.append(each)
            cp = self.find_catch_pt(items, x, y)
            if cp:
                x, y = cp
                if self.catch_pnt:
                    self.genCatchPnt(x, y, regen=1)
                else:
                    self.genCatchPnt(x, y)
            else:
                if self.catch_pnt:
                    self.canvas.delete(self.catch_pnt)
                    self.catch_pnt = 0
            p1 = (x, y)  # func wants canvas coords to make rubber element 
            func = 'self.%s(%s)' % (self.op, p1)
            eval(func)
        elif self.sel_box_crnr:
            x1, y1 = self.sel_box_crnr
            if self.sel_boxID:
                self.canvas.coords(self.sel_boxID, x1, y1, x, y)
            else:
                self.sel_boxID = self.canvas.create_rectangle(x1, y1, x, y,
                                                              outline='cyan',
                                                              tags='sb')
        elif self.sel_mode == 'items':
            func = 'self.%s()' % self.op
            eval(func)

    def find_catch_pt(self, items, x, y):
        cr = self.catch_radius
        if len(items) == 1:
            item = items[0]
            if self.canvas.type(item) == 'arc':
                x0, y0, x1, y1 = self.canvas.coords(item)
                xc = (x0+x1)/2
                yc = (y0+y1)/2
                r = (x1-x0)/2
                a0 = float(self.canvas.itemcget(item, 'start'))
                a1 = a0 + float(self.canvas.itemcget(item, 'extent'))
                a0 = -a0*math.pi/180
                a1 = -a1*math.pi/180
                p0 = (xc+r*math.cos(a0), yc+r*math.sin(a0))
                p1 = (xc+r*math.cos(a1), yc+r*math.sin(a1))
                arc_end_pts = (p0, p1)
                if self.catchCntr:
                    return (xc, yc)
                caught = None
                for pt in arc_end_pts:
                    if pnt_in_box_p((pt[0], pt[1]),
                                    (x-cr, y-cr, x+cr, y+cr)):
                        caught = pt
                if caught:
                    return caught
                else:
                    ip = line_circ_inters(xc, yc, x, y, xc, yc, r)
                    for pt in ip:
                        if p2p_dist(pt, (x,y)) < cr:
                            return pt
            elif self.canvas.type(item) == 'oval':
                x0, y0, x1, y1 = self.canvas.coords(item)
                xc, yc = ctr = midpoint((x0, y0), (x1, y1))
                r = (x1-x0)/2
                if self.catchCntr:
                    return (xc, yc)
                else:
                    inters_pts = line_circ_inters(xc, yc, x, y, xc, yc, r)
                    for pt in inters_pts:
                        if p2p_dist(pt, (x,y)) < cr:
                            return (pt[0], pt[1])
            elif self.canvas.type(item) == 'line':
                x0, y0, x1, y1 = self.canvas.coords(item)  # end points
                xm, ym = midpoint((x0, y0), (x1, y1))   # mid point
                pts = ((x0, y0), (x1, y1), (xm, ym))
                caught = None
                for pt in pts:
                    if 'g' in self.canvas.gettags(item) and \
                       pnt_in_box_p((pt[0], pt[1]), (x-cr, y-cr, x+cr, y+cr)):
                        caught = pt
                if caught:
                    return caught
                else:
                    line = cnvrt_2pts_to_coef((x0, y0), (x1, y1))
                    u, v = proj_pt_on_line(line, (x, y))
                    if x0<u<x1 or x0>u>x1 or y0<v<y1 or y0>v>y1:
                        return (u, v)
        
        elif len(items) > 1:  # intersection found                   
            if self.canvas.type(items[0]) == 'line' and\
               self.canvas.type(items[1]) == 'line':
                a,b,c,d = self.canvas.coords(items[0])
                e,f,g,h = self.canvas.coords(items[1])
                line1 = cnvrt_2pts_to_coef((a,b), (c,d))
                line2 = cnvrt_2pts_to_coef((e,f), (g,h))
                if line1 == line2:  # colinear; toss one and try again
                    items.pop()
                    return self.find_catch_pt(items, x, y)
                ip = intersection(line1, line2)
                if not ip:
                    items.pop(0)
                    return self.find_catch_pt(items, x, y)
                elif ip:
                    return ip
            elif self.canvas.type(items[0]) in ('oval', 'arc') and\
                 self.canvas.type(items[1]) in ('oval', 'arc'):
                a,b,c,d = self.canvas.coords(items[0])
                x1, y1 = midpoint((a,b), (c,d))
                r1 = (c-a)/2
                e,f,g,h = self.canvas.coords(items[1])
                x2, y2 = midpoint((e,f), (g,h))
                r2 = (g-e)/2
                ip = circ_circ_inters(x1, y1, r1, x2, y2, r2)
                if ip:
                    for pt in ip:
                        if p2p_dist(pt, (x,y)) < cr:
                            return pt
            elif self.canvas.type(items[0]) in ('oval', 'arc') and\
                 self.canvas.type(items[1]) == 'line':
                items[0], items[1] = items[1], items[0]
            if self.canvas.type(items[0]) == 'line' and\
               self.canvas.type(items[1]) in ('oval', 'arc'):
                x1,y1,x2,y2 = self.canvas.coords(items[0])
                line = cnvrt_2pts_to_coef((x1,y1), (x2,y2))
                e,f,g,h = self.canvas.coords(items[1])
                xc, yc = cntr = midpoint((e,f), (g,h))
                r = (g-e)/2
                ip = line_circ_inters(x1, y1, x2, y2, xc, yc, r)
                for pt in ip:
                    if p2p_dist(pt, (x,y)) < cr:
                        return pt
        
if __name__ == '__main__':
    draw = Draw()
    draw.run()
