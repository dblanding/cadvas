import sys
from tkinter import *
import math
import entities

def but(root, text, row, col, com=None, span=2, clr='darkslateblue', pad=1): 
    w = Button(root, text=text, command=com, bg=clr, fg='white', padx=pad)
    w.grid(row=row, column=col, columnspan=span, sticky=E+W)

def ent(root, var, row, col=2, span=10):
    e = Entry(root, textvariable=var, relief=SUNKEN)
    e.grid(row=row, column=col, columnspan=span)

def f2s(f):
    """Convert float to string with 12 significant figures."""
    return '%1.12f' % f

class TxtDialog(Toplevel):
    """Dimension dialog for editing dimension parameters."""
    
    def __init__(self, caller=None):
        Toplevel.__init__(self)
        self.caller = caller    # ref to Draw instance
        self.title('Text Parameters')
        self.protocol("WM_DELETE_WINDOW", self.quit)
        self.resizable(width=0, height=0)
        if caller:
            self.transient(caller)
        
        but(self, 'Style', 0, 0, lambda r='t': self.pr(r), clr='dimgray')
        but(self, 'Size', 1, 0, lambda r='z': self.pr(r), clr='dimgray')
        but(self, 'Color', 2, 0, lambda r='y': self.pr(r), clr='dimgray')
        but(self, 'Text', 3, 0, lambda r='x': self.pr(r), clr='dimgray')

        self.tdisplay = StringVar()
        self.zdisplay = StringVar()
        self.ydisplay = StringVar()
        self.xdisplay = StringVar()
        ent(self, self.tdisplay, 0)
        ent(self, self.zdisplay, 1)
        ent(self, self.ydisplay, 2)
        ent(self, self.xdisplay, 3)

        but(self, 'Change Dim Parameters', 4, 0, self.mm2in, span=12, clr='dimgray')
        

    def quit(self):
        if self.caller:
            self.caller.txtdialog = None
        self.destroy()

    def mm2in(self):  # save modified text object to caller
        attribs = (self.coords,
                   self.xdisplay.get().strip("'"),
                   self.tdisplay.get().strip("'"),
                   int(float(self.zdisplay.get().strip("'"))),
                   self.ydisplay.get().strip("'"))
        print(attribs)
        tx = entities.TX(attribs)
        self.caller.modified_text_object = tx

    def putx(self, value):
        self.xdisplay.set(repr(value))
        self.keip = False
        
    def puty(self, value):
        self.ydisplay.set(repr(value))
        self.keip = False
        
    def putz(self, value):
        self.zdisplay.set(repr(value))
        self.keip = False
        
    def putt(self, value):
        self.tdisplay.set(repr(value))
        self.keip = False
        
if __name__ == '__main__':
    TxtDialog(None).mainloop()
