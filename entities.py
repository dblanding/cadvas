class Entity:
    """A drawing entity such as text, containing all params except handle."""

    def __init__(self, attribs, show=True):
        self.attribs = attribs
        self.show = show

    def __hash__(self):
        return hash(self.attribs)

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and
                self.attribs == other.attribs)

    def __repr__(self):
        return f"Entity of type {etype} with attributes {attribs}"

class TxEntity(Entity):
    """Text entity"""

    def __init__(self, attribs):
        super().__init__(attribs)
        self.coords, self.text, self.style, self.size, self.color = self.attribs
        self.type = 'tx'

if __name__ == "__main__":
    attribs = ((100, 100), "rain is coming", 'Courier', 10, 'cyan')
    txt = TxEntity(attribs)
    txt1 = TxEntity(attribs)
    attribs = ((100, 100), "rain isn't coming", 'Courier', 10, 'cyan')
    txt2 = TxEntity(attribs)
    print(txt == txt1)
    print(txt1 == txt2)
    print(txt == txt2)
    print(txt.coords, txt.text, txt.type)
