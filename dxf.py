#!/usr/bin/env python
"""Utilities for translating between dxf and native cadvas (.pkl) format"""

import ezdxf


def dxf2native(filename):
    """Generate cadvas native CAD data from dxf entities."""

    gldict = {}
    gcdict = {}
    gadict = {}
    k = 0
    dwg = ezdxf.readfile(filename)
    for e in dwg.modelspace():  # e = dxf entity
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

    drawdict = {}
    drawdict['gl'] = gldict
    drawdict['gc'] = gcdict
    drawdict['ga'] = gadict
    return drawdict


def native2dxf(drawdict, dxf_filename):
    """Generate .dxf file format from native CADvas drawing."""
    # Create a new DXF R2010 drawing
    dwg = ezdxf.new('R2010')  # Official DXF version name: 'AC1024'
    msp = dwg.modelspace()  # Create new model space
    # Add new entities to the model space
    for p0, p1 in drawdict['gl'].values():
        msp.add_line(p0, p1)
    for center, radius in drawdict['gc'].values():
        msp.add_circle(center, radius) 
    for center, radius, start, end in drawdict['ga'].values():
        msp.add_arc(center, radius, start, end) 
    dwg.saveas(dxf_filename)
