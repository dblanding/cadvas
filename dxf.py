#!/usr/bin/env python
"""Utilities for translating between dxf and native cadvas (.pkl) format"""
import math
import ezdxf

def pnt_n_vctr_to_coef(pnt, vector):
    (u, v, w) = pnt
    (x, y, z) = vector
    pt1 = (u, v)
    pt2 = (u + x, v + y)
    return cnvrt_2pts_to_coef(pt1, pt2)

def normalize_vector(vctr):
    x, y, z = vctr
    diag = math.sqrt(x * x + y * y)
    return (x / diag, y / diag, 0)

def coef_to_pnt_n_vctr(coords):
    a, b, c = coords
    if abs(b) >= abs(a):
        y_intercept = -(c / b)
        p0 = (0, y_intercept, 0)
    else:
        x_intercept = -(c / a)
        p0 = (x_intercept, 0, 0)
    vector = normalize_vector((b, -a, 0))
    return (p0, vector)

def cnvrt_2pts_to_coef(pt1, pt2):
    """Return (a,b,c) coefficients of cline defined by 2 (x,y) pts."""
    x1, y1 = pt1
    x2, y2 = pt2
    a = y2 - y1
    b = x1 - x2
    c = x2*y1-x1*y2
    return (a, b, c)

def dxf2native(filename):
    """Generate cadvas native CAD data from dxf entities."""

    cllist = []
    gldict = {}
    gcdict = {}
    gadict = {}
    txdict = {}
    k = 0
    dwg = ezdxf.readfile(filename)
    for e in dwg.modelspace():  # e = dxf entity
        if e.dxftype() == 'XLINE':
            print(e.dxfattribs())
            coords = pnt_n_vctr_to_coef(e.dxf.start, e.dxf.unit_vector)
            k+=1
            cllist.append(coords)
        if e.dxftype() == 'LINE':
            # print(e.dxfattribs())
            coords = (e.dxf.start, e.dxf.end)
            k+=1
            gldict[k] = coords
        elif e.dxftype() == 'CIRCLE':
            # print(e.dxfattribs())
            coords = (e.dxf.center, e.dxf.radius)
            k+=1
            gcdict[k] = coords
        elif e.dxftype() == 'ARC':
            # print(e.dxfattribs())
            coords = (e.dxf.center, e.dxf.radius,
                      e.dxf.start_angle, e.dxf.end_angle)
            k+=1
            gadict[k] = coords
        elif e.dxftype() == 'TEXT':
            # print(e.dxfattribs())
            x, y = e.dxfattribs()['align_point']
            text = e.dxfattribs()['text']
            k += 1
            txdict[k] = (x, y, text)

    drawdict = {}
    drawdict['cl'] = cllist
    drawdict['gl'] = gldict
    drawdict['gc'] = gcdict
    drawdict['ga'] = gadict
    drawdict['tx'] = txdict
    return drawdict


def native2dxf(drawdict, dxf_filename):
    """Generate .dxf file format from native CADvas drawing."""
    # Create a new DXF R2010 drawing
    dwg = ezdxf.new('R2010')  # Official DXF version name: 'AC1024'
    msp = dwg.modelspace()  # Create new model space
    # Add new entities to the model space
    for coef in drawdict['cl']:
        pnt, vctr = coef_to_pnt_n_vctr(coef)
        msp.add_xline(pnt, vctr)
    for p0, p1 in drawdict['gl'].values():
        msp.add_line(p0, p1)
    for center, radius in drawdict['gc'].values():
        msp.add_circle(center, radius) 
    for center, radius, start, end in drawdict['ga'].values():
        msp.add_arc(center, radius, start, end)
    for x, y, text in drawdict['tx'].values():
        dxfattribs = dict((('align_point',(x, y)),
                           ('halign', 2),
                           ('height', 1.0),
                           ('insert', (x, y)),
                           ('layer', '0'),
                           ('oblique', 0.0),
                           ('paperspace', 0),
                           ('rotation', 0.0),
                           ('style', 'STANDARD'),
                           ('text', text),
                           ('text_generation_flag', 0),
                           ('valign', 2),
                           ('width', 1.0)))
        msp.add_text(text, dxfattribs)
    dwg.saveas(dxf_filename)
