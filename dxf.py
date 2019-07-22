#!/usr/bin/env python
"""Utilities for translating between dxf and native cadvas (.pkl) format.
"""

##########################
#
# dxf to native
#
##########################

def _parsepair(a):
    """Parse a pair of lines containing 'group code' and 'value'."""
    groupcode = a.next().strip()
    value = a.next().strip()
    return (groupcode, value)

def _gotosection(a, secname):
    """Go to secname and stop."""
    while 1:
        gc, val = _parsepair(a)
        if gc == '2' and val == secname:
            return
            
def _get_units(a):
    """Parse through HEADER section and detect whether units
    are set for metric(1) or english(0)."""
    _gotosection(a, 'HEADER')
    units = 0   # assume inches by default
    while 1:
        gc, val = _parsepair(a)
        if gc == '9' and val == '$MEASUREMENT':
            gc, val = _parsepair(a)
            if gc == '70':
                units = int(val)
        elif gc == '0' and val == 'ENDSEC':
            return units

def _process_entity(a):
    """Return a dictionary of groupcodes : values for the next
    entity in the ENTITIES section. Go until groupcode == 0."""
    entitydict = {}
    flag = 1
    while 1:
        gc, val = _parsepair(a)
        if gc == '0':
            if val == 'ENDSEC':
                flag = 0    # Done with ENTITIES section
            return (entitydict, flag, val)
        else:
            entitydict[gc] = val

def _parse_file(f):
    """Parse contents of the dxf file, looking for units and
    all the drawing entities."""
    a = iter(open(f))
    units = _get_units(a)
    _gotosection(a, 'ENTITIES')
    lines = []
    circles = []
    arcs = []
    entities = [lines, circles, arcs]
    gc, val = _parsepair(a)
    while 1:
        if val == 'LINE':
            ed, f, val = _process_entity(a)
            lines.append(ed)
        elif val == 'CIRCLE':
            ed, f, val = _process_entity(a)
            circles.append(ed)
        elif val == 'ARC':
            ed, f, val = _process_entity(a)
            arcs.append(ed)
        else:
            ed, f, val = _process_entity(a)
        if not f:
            return (units, entities)

def _gen_cadvas_entities(f):
    """Generate cadvas native CAD data from dxf entities."""
    units, entities = _parse_file(f)
    lines, circles, arcs = entities
    if units:
        scale = 1.0
    else:
        scale = 25.4
    gldict = {}
    gcdict = {}
    gadict = {}
    k = 0
    for line in lines:
        p1 = (float(line['10'])*scale,
              float(line['20'])*scale)
        p2 = (float(line['11'])*scale,
              float(line['21'])*scale)
        coords = (p1, p2)
        k+=1
        gldict[k] = coords
    k = 0
    for circ in circles:
        cntr = (float(circ['10'])*scale,
                float(circ['20'])*scale)
        radius = float(circ['40'])*scale
        coords = (cntr, radius)
        k+=1
        gcdict[k] = coords
    k = 0
    for arc in arcs:
        cntr = (float(arc['10'])*scale,
                float(arc['20'])*scale)
        radius = float(arc['40'])*scale
        a0 = float(arc['50'])
        a1 = float(arc['51'])
        coords = (cntr, radius, a0, a1)
        k+=1
        gadict[k] = coords
    return [gldict, gcdict, gadict]

def dxf2native(inf):
    drawdict = {}
    gldict, gcdict, gadict = _gen_cadvas_entities(inf)
    drawdict['gl'] = gldict
    drawdict['gc'] = gcdict
    drawdict['ga'] = gadict
    return drawdict

##########################
#
# native to dxf
#
##########################

template1 = """999
DXF file created from CADvas
0
SECTION
2
HEADER
9
$ACADVER
1
AC1006
9
$INSBASE
10
0.0
20
0.0
30
0.0
9
$EXTMIN
10
0.0
20
0.0
9
$EXTMAX
10
1000.0
20
1000.0
0
ENDSEC
0
SECTION
2
TABLES
0
TABLE
2
LTYPE
70
1
0
LTYPE
2
CONTINUOUS
70
64
3
Solid line
72
65
73
0
40
0.000000
0
ENDTAB
0
TABLE
2
LAYER
70
6
0
LAYER
2
1
70
64
62
7
6
CONTINUOUS
0
LAYER
2
2
70
64
62
7
6
CONTINUOUS
0
ENDTAB
0
TABLE
2
STYLE
70
0
0
ENDTAB
0
ENDSEC
0
SECTION
2
BLOCKS
0
ENDSEC
  0
SECTION
  2
ENTITIES"""

template2 = """
  0
ENDSEC
  0
EOF"""

def native2dxf(drawdict, dxf):
    scale = 25.4    # default value (cadvas units are mm, dxf units = inches)
    gldict = drawdict['gl']
    gcdict = drawdict['gc']
    gadict = drawdict['ga']
    entitytext = ""
    for p0, p1 in gldict.values():
        entitytext = entitytext + """
  0
LINE
  8
0
 10
%f
 20
%f
 11
%f
 21
%f""" % (p0[0]/scale, p0[1]/scale, p1[0]/scale, p1[1]/scale)
    for pc, r in gcdict.values():
        entitytext = entitytext + """
  0
CIRCLE
  8
0
 10
%f
 20
%f
 40
%f""" % (pc[0]/scale, pc[1]/scale, r/scale)
    for pc, r, a0, a1 in gadict.values():
        entitytext = entitytext + """
  0
ARC
  8
0
 10
%f
 20
%f
 40
%f
 50
%f
 51
%f""" % (pc[0]/scale, pc[1]/scale, r/scale, a0, a1)
    dxftext = template1 + entitytext + template2
    f = open(dxf, 'w')
    f.write(dxftext)
    f.close()
    
