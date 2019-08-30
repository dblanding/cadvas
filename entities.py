    
class TxObject:
    """Text object initialized with a tuple of attributes.

    attribs = (coords, text, style, size, color)
    """

    def __init__(self, attribs):
        self.coords, self.text, self.style, self.size, self.color = attribs
        self.type = 'tx'  # text entity
        self.show = True  # control show/hide visability

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and
                self.coords == other.coords and
                self.text == other.text and
                self.style == other.style and
                self.size == other.size and
                self.color == other.color)

    def __repr__(self):
        return f"{self.type} object with coordinates {self,coords}"

    def get_attribs(self):
        return [self.coords, self.text, self.style, self.size, self.color]

if __name__ == "__main__":
    attribs = ((50,50), "this is some text", 'Verdana', 10, 'cyan',)
    t1 = TxObject(attribs)
    print(t1)
    print(t1.coords)
    print(t1.get_attribs())
    t1.coords = (100, 100)
    
