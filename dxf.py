#!/usr/bin/env python
"""Utilities for translating between dxf and native cadvas (.pkl) format"""
import math
import ezdxf
import entities
from cadvas import geomcolor
from cadvas import constrcolor

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

    drawlist = []
    dwg = ezdxf.readfile(filename)
    for e in dwg.modelspace():  # e = dxf entity
        if e.dxftype() == 'XLINE':
            # print(e.dxfattribs())
            coords = pnt_n_vctr_to_coef(e.dxf.start, e.dxf.unit_vector)
            cc = entities.CC((coords, constrcolor))
            drawlist.append(cc)
        if e.dxftype() == 'LINE':
            # print(e.dxfattribs())
            coords = (e.dxf.start, e.dxf.end)
            gl = entities.GL((coords, geomcolor))
            drawlist.append(gl)
        elif e.dxftype() == 'CIRCLE':
            # print(e.dxfattribs())
            coords = (e.dxf.center, e.dxf.radius)
            gc = entities.GC((coords, geomcolor))
            drawlist.append(gc)
        elif e.dxftype() == 'ARC':
            # print(e.dxfattribs())
            coords = (e.dxf.center, e.dxf.radius,
                      e.dxf.start_angle, e.dxf.end_angle)
            ga = entities.GA((coords, geomcolor))
            drawlist.append(ga)
        elif e.dxftype() == 'TEXT':
            # print(e.dxfattribs())
            coords = e.dxfattribs()['align_point']
            text = e.dxfattribs()['text']
            style = e.dxfattribs()['style']
            size = e.dxfattribs()['height']
            attribs = (coords, text, style, size, 'cyan')
            tx = entities.TX(attribs)
            drawlist.append(tx)
            
    return drawlist


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
    for attribs in drawdict['tx'].values():
        coords = attribs['coords']
        text = attribs['text']
        style = attribs['style']
        size = attribs['size']
        color = attribs['color']
        dxfattribs = dict((('align_point', coords),
                           ('halign', 2),
                           ('height', size),
                           ('insert', coords),
                           ('layer', '0'),
                           ('oblique', 0.0),
                           ('paperspace', 0),
                           ('rotation', 0.0),
                           ('style', style),
                           ('text', text),
                           ('text_generation_flag', 0),
                           ('valign', 2),
                           ('width', 1.0)))
        msp.add_text(text, dxfattribs)
    dwg.saveas(dxf_filename)
