import sys
from tkinter import *
import math

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

        but(self, 'mm->in', 4, 0, self.mm2in, span=4, clr='dimgray')
        but(self, 'in->mm', 4, 4, self.in2mm, span=4, clr='dimgray')
        but(self, 'Sto', 4, 8, self.storex, clr='darkgreen')
        but(self, 'Rcl', 4, 10, self.recallx, clr='darkgreen')
        

    def quit(self):
        if self.caller:
            self.caller.txtdialog = None
        self.destroy()

    def pr(self, val):
        """Send value in register to CADvas."""
        # There must be a better way to get this value
        str_value = repr(eval('self.'+val+'display.get()')).strip("'")
        self.caller.enterfloat(str_value)
        self.keip = False
        self.needrup = True
        

    def keyin(self, c):
        if self.keip:
            self.xdisplay.set(self.xdisplay.get()+c)
        else:
            self.keip = True
            if self.needrup:
                self.rotateup(loop=0)
            self.clearx()
            self.keyin(c)

    def enter(self):
        self.tdisplay.set(self.zdisplay.get())
        self.zdisplay.set(self.ydisplay.get())
        self.ydisplay.set(self.xdisplay.get())
        self.keip = False
        self.needrup = False

    def calc(self, op):
        pass
        """Arithmetic calculations between x and y registers, then rotate down."""
        try:
            if op == '+/-':
                self.xdisplay.set(repr(eval('-'+self.xdisplay.get())))
            else:
                x = repr(eval(self.ydisplay.get()+op+self.xdisplay.get()))
                self.xdisplay.set(x)
                self.ydisplay.set(self.zdisplay.get())
                self.zdisplay.set(self.tdisplay.get())
            self.keip = False
            self.needrup = True
        except:
            self.xdisplay.set("ERROR")
        

    def func(self, op, in_cnvrt=0, out_cnvrt=0):
        """Evaluate function op then put result in x-register, don't rotate stack.
        if in_cnvrt: convert input value of x-register from degrees to radians.
        if out_cnvrt: convert output value of x-register from radians to degrees."""
        try:
            x = float(self.xdisplay.get())
        except:
            x = 0
        try:
            y = float(self.ydisplay.get())
        except:
            y = 0
        if in_cnvrt:
            x = x * math.pi / 180
        result = eval(op)
        if out_cnvrt:
            result = result * 180 / math.pi
        self.xdisplay.set(f2s(result))
        self.keip = False
        self.needrup = True

    def mm2in(self):
        if self.xdisplay.get():
            self.xdisplay.set(repr(eval(self.xdisplay.get()+'/25.4')))
            self.keip = False
            self.needrup = True

    def in2mm(self):
        if self.xdisplay.get():
            self.xdisplay.set(repr(eval(self.xdisplay.get()+'*25.4')))
            self.keip = False
            self.needrup = True

    def storex(self):
        self.mem = self.xdisplay.get()
        self.keip = False
        self.needrup = True

    def recallx(self):
        self.rotateup()
        self.xdisplay.set(self.mem)
        self.keip = False
        self.needrup = True

    def rotateup(self, loop=1):
        x = self.tdisplay.get()
        self.tdisplay.set(self.zdisplay.get())
        self.zdisplay.set(self.ydisplay.get())
        self.ydisplay.set(self.xdisplay.get())
        if loop:
            self.xdisplay.set(x)

    def rotatedn(self):
        x = self.xdisplay.get()
        self.xdisplay.set(self.ydisplay.get())
        self.ydisplay.set(self.zdisplay.get())
        self.zdisplay.set(self.tdisplay.get())
        self.tdisplay.set(x)

    def trimx(self):
        self.xdisplay.set(self.xdisplay.get()[:-1])

    def swapxy(self):
        x = self.xdisplay.get()
        y = self.ydisplay.get()
        self.xdisplay.set(y)
        self.ydisplay.set(x)

    def clearx(self):
        self.xdisplay.set('')
        
    def clearall(self):
        self.xdisplay.set('')
        self.ydisplay.set('')
        self.zdisplay.set('')
        self.tdisplay.set('')

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
    Dialog(None).mainloop()
