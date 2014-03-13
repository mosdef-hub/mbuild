from scimath.units import example_units
from scimath.units.unit_db import UnitDB
from scimath.units.unit_system import UnitSystem
import os
from scimath.units.unit_manager import unit_manager

__author__ = 'sallai'


print "Loading additional units into unit manager " + str(unit_manager)
# instantiate a UnitDB object using our text files:
udb = UnitDB()
udb.get_family_members_from_file(filename=os.path.join( 'units',
                   'unit_family_membership.txt' ))

udb.get_unit_families_from_file(filename=os.path.join( 'units',
                   'unit_families.txt' ))


# # Add each unit system from unit_db
# for udb_sys in udb.unit_systems:
#
#     us = UnitSystem(udb_sys)
#
#     column_name = '%s_%s' % (udb_sys.upper(), 'UNITS')
#
#     for fam in udb.unit_names:
#         us.add_family(fam, \
#             udb.unit_names[fam][udb.column_names[column_name]],
#             description=udb.unit_names[fam][udb.column_names['DESCRIPTION']],
#             inverse=udb.unit_names[fam][udb.column_names['INVERSE']])
#
#     unit_manager.add_unit_system(us)

# Add unit members/preferred names from unit_db
unit_manager.unit_members.update(udb.member_names)
unit_manager.preferred_names.update(udb.preferred_names)

for name in udb.member_names:
    if name.find('*') != -1:
        unit_manager._wildcards.append(name)

# # Load unit converters from default_unit_converters file
# self.unit_converters = default_unit_converters
# self.default_system = self.unit_systems[0]

