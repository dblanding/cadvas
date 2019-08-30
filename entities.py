class Entity:
    """A drawing entity such as text, containing all params except handle."""

    def __init__(self, attribs):
        self.attribs = attribs  # needs to be a list to allow editing
        self.show = True  # to control show/hide visability

    def __hash__(self):
        return hash(self.attribs)

    def __eq__(self, other):
        return (self.__class__ == other.__class__ and
                self.attribs == other.attribs)

    def __repr__(self):
        return f"Entity of type {etype} with attributes {attribs}"

class TxEntity(Entity):
    """Text entity object initialized with a list of attributes.

    attribs = [coords, text, style, size, color]
    """

    def __init__(self, attribs):
        super().__init__(attribs)
        self.coords, self.text, self.style, self.size, self.color = self.attribs
        self.type = 'tx'

