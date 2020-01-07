import sys
import os

import mbuild as mb
import foyer


class Ethane(mb.Compound):
    def __init__(self):
        super(Ethane, self).__init__()
        # PDB contains coordinates and some atom names
        # But we need to touch up the atom names 
        # and add bonds
        mb.load('TraPPE_UA_2_ethane_dimer.pdb', compound=self,
                relative_to_module=self.__module__)
        self.name = 'Eth'
        particles = [a for a in self.particles()]
        particles[0].name = particles[0].name[1:]
        particles[1].name = particles[1].name[1:]
        self.add_bond((particles[0], particles[1]))

        # For convenience, each mbcompound carries its xml file
        xml_path = os.path.realpath(
                sys.modules[self.__module__].__file__)
        file_dir = os.path.dirname(xml_path)
        self.xml = os.path.join(file_dir, "TraPPE_UA_2_fully_flexible_ethane.xml" )
