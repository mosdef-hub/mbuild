import os

import numpy as np
import parmed as pmd
import pytest

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.utils.io import get_fn, has_intermol, has_openbabel
from mbuild.tests.base_test import BaseTest

class TestCompound(BaseTest):

    def test_load_and_create(self):
        mb.load(get_fn('methyl.pdb'))

    def test_update_from_file(self, ch3):
        ch3.update_coordinates(get_fn("methyl.pdb"))

    def test_save_simple(self, ch3):
        extensions = ['.xyz', '.pdb', '.mol2']
        for ext in extensions:
            outfile = 'methyl_out' + ext
            ch3.save(filename=outfile)
            assert os.path.exists(outfile)

    def test_save_box(self, ch3):
        extensions = ['.mol2', '.pdb', '.hoomdxml', '.gro']
        box_attributes = ['mins', 'maxs', 'lengths']
        custom_box = mb.Box([.8, .8, .8])
        for ext in extensions:
            outfile_padded = 'padded_methyl' + ext
            outfile_custom = 'custom_methyl' + ext
            ch3.save(filename=outfile_padded, box=None, overwrite=True)
            ch3.save(filename=outfile_custom, box=custom_box, overwrite=True)
            padded_ch3 = mb.load(outfile_padded)
            custom_ch3 = mb.load(outfile_custom)
            for attr in box_attributes:
                pad_attr = getattr(padded_ch3.boundingbox, attr)
                custom_attr = getattr(custom_ch3.boundingbox, attr)
                assert np.array_equal(pad_attr, custom_attr)

    def test_save_overwrite(self, ch3):
        extensions = ['.gsd', '.hoomdxml', '.lammps', '.lmp', '.top', '.gro']
        for ext in extensions:
            outfile = 'lyhtem' + ext
            ch3.save(filename=outfile)
            ch3.save(filename=outfile, overwrite=True)
            with pytest.raises(IOError):
                ch3.save(filename=outfile, overwrite=False)

    def test_save_forcefield(self, methane):
        exts = ['.gsd', '.hoomdxml', '.lammps', '.lmp', '.top', '.gro',
                '.mol2', '.pdb', '.xyz']
        for ext in exts:
            methane.save('lythem' + ext,
                         forcefield_name='oplsaa',
                         overwrite=True)

    def test_save_resnames(self, ch3, h2o):
        system = mb.Compound([ch3, h2o])
        system.save('resnames.gro', residues=['CH3', 'H2O'])
        struct = pmd.load_file('resnames.gro')

        assert struct.residues[0].name == 'CH3'
        assert struct.residues[1].name == 'H2O'

    def test_save_resnames_single(self, c3, n4):
        system = mb.Compound([c3, n4])
        system.save('resnames_single.gro', residues=['C3', 'N4'])
        struct = pmd.load_file('resnames_single.gro')
        assert struct.residues[0].number ==  1
        assert struct.residues[1].number ==  2

    def test_save_references(self, methane):
        methane.save('methyl.mol2', forcefield_name='oplsaa',
                     references_file='methane.bib')
        assert os.path.isfile('methane.bib')

    def test_save_combining_rule(self, methane):
        combining_rules = ['lorentz', 'geometric']
        gmx_rules = {'lorentz': 2, 'geometric': 3}
        for combining_rule in combining_rules:
            methane.save('methane.top', forcefield_name='oplsaa',
                         combining_rule=combining_rule, overwrite=True)
            with open('methane.top') as fp:
                for i, line in enumerate(fp):
                    if i == 18:
                        gmx_rule = int(line.split()[1])
                        assert gmx_rule == gmx_rules[combining_rule]

    def test_batch_add(self, ethane, h2o):
        compound = mb.Compound()
        compound.add([ethane, h2o])
        assert compound.n_particles == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds1(self, ethane):
        compound = mb.Compound(ethane)
        assert compound.n_particles == 8
        assert compound.n_bonds == 7

    def test_init_with_subcompounds2(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])
        assert compound.n_particles == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds3(self, ethane, h2o):
        compound = mb.Compound([ethane, [h2o, mb.clone(h2o)]])
        assert compound.n_particles == 8 + 2 * 3
        assert compound.n_bonds == 7 + 2 * 2

    def test_init_with_bad_name(self):
        with pytest.raises(ValueError):
            mb.Compound(name=1)

    def test_add_wrong_input(self, ethane):
        with pytest.raises(ValueError):
            ethane.add('water')

    def test_add_existing_parent(self, ethane, h2o):
        water_in_water = mb.clone(h2o)
        h2o.add(water_in_water)
        with pytest.raises(MBuildError):
            ethane.add(water_in_water)

    def test_add_label_exists(self, ethane, h2o):
        ethane.add(h2o, label='water')
        with pytest.raises(MBuildError):
            ethane.add(mb.clone(h2o), label='water')

    def test_set_pos(self, ethane):
        with pytest.raises(MBuildError):
            ethane.pos = [0, 0, 0]

    def test_xyz(self, ethane):
        xyz = ethane.xyz
        assert xyz.shape == (8, 3)

        xyz = ethane.xyz_with_ports
        assert xyz.shape == (24, 3)

    def test_mirror_std_axis(self, labeled_tetrahedral):
       origdict = {piece.my_label: piece.pos for piece in labeled_tetrahedral}
       labeled_tetrahedral.mirror(about_vectors=[(-1,6,0),(1,1,0)])
       for k, v in origdict.items():
           v[2]*=-1
           assert np.allclose(labeled_tetrahedral[k].pos, v, atol=1e-8)

    def test_mirror_colinearity(self, methane):
       with pytest.raises(ValueError):
           methane.mirror(about_vectors=[[-1,-1,-1],[1,1,1]])

    def test_mirror_overdefined_4points(self, labeled_tetrahedral):
       with pytest.raises(ValueError):
           labeled_tetrahedral.mirror(about_vectors= None,
                                      mirror_plane_points= [["H"],
                                                            ["F"],
                                                            ["N"],
                                                            ["O"]])

    def test_mirror_overdefined_3vecs(self, labeled_tetrahedral):
       with pytest.raises(ValueError):
           labeled_tetrahedral.mirror(about_vectors= [(0,1,0), (1,0,0), (-1,0,0)])

    def test_mirror_overdefined_2vecs_2points(self, labeled_tetrahedral):
       with pytest.raises(ValueError):
           labeled_tetrahedral.mirror(about_vectors= [(1,0,0), (0,0,1)],
                                      mirror_plane_points=[["H[0]"],["F[0]"]])

    def test_mirror_overdefined_1vec_3points(self, labeled_tetrahedral):
       with pytest.raises(ValueError):
           labeled_tetrahedral.mirror(about_vectors= [(1,0,0)],
                                      mirror_plane_points=[["H[0]"],["F[0]"], ["N"]])

    def test_mirror_underdefined_1vec(self, labeled_tetrahedral):
       with pytest.raises(ValueError):
           labeled_tetrahedral.mirror(about_vectors= [(1,0,0)])

    def test_mirror_underdefined_2points(self, labeled_tetrahedral):
       with pytest.raises(ValueError):
           labeled_tetrahedral.mirror(mirror_plane_points=[["H[0]"],["F[0]"]])

    def test_mirror_3points_and_anchor(self,labeled_tetrahedral):
       l = labeled_tetrahedral
       anchor = np.array([1,1,1])
       l.translate(anchor)
       f = l["F[0]"].pos
       n = l["N[0]"].pos
       o = l["O[0]"].pos
       hfuture = 2*anchor-l["H[0]"].pos
       dist = 2*(np.mean([o,n,f], axis=0) - anchor)
       l.mirror(mirror_plane_points= [["F"], ["O"], ["N"]], anchor_point=anchor)
       assert np.allclose(anchor, l["C[0]"].pos, atol=1e-15)
       assert np.allclose(hfuture, l["H[0]"].pos, atol=2e-3)
       assert np.allclose((l["F[0]"].pos + dist), f, atol=7e-4)
       assert np.allclose(l["N[0]"].pos+dist, n, atol = 7e-4)
       assert np.allclose(l["O[0]"].pos+dist, o, atol = 7e-4)


    def test_mirror_2points_1vec(self, benzene):
       b=benzene
       c0, c1, c2, c3, c4, c5 = (b["C[{}]".format(ii)].pos for ii in range(6))
       b.mirror(about_vectors=[[0,0,1]], mirror_plane_points=[["C[0]"], ["C[3]"]])
       assert np.allclose(b["C[0]"].pos, c0, atol=1e-6)
       assert np.allclose(b["C[3]"].pos, c3, atol=5e-4)
       assert np.allclose(b["C[1]"].pos, c5, atol=5e-4)
       assert np.allclose(b["C[2]"].pos, c4, atol=5e-4)
       assert np.allclose(b["C[4]"].pos, c2, atol=5e-4)
       assert np.allclose(b["C[5]"].pos, c1, atol=5e-4)

    def test_mirror_2vec(self, labeled_tetrahedral):
       l = labeled_tetrahedral
       l.translate([1,1,1])
       old_N = l["N[0]"].pos
       old_O = l["O[0]"].pos
       old_H = l["H[0]"].pos
       old_F = l["F[0]"].pos
       men = np.mean([old_O, old_H, old_F], axis=0)
       dist = 2*(men-old_N)
       vecs = [old_O-old_H,old_O-old_F]
       l.mirror(about_vectors=vecs, anchor_point= men)
       print(l["N[0]"].pos-dist)
       print(l["N[0]"].pos+dist)
       print(old_N)
       # print(dist)
       # print(l["N[0]"].pos)
       # print(old_N-(l["N[0]"].pos-dist))
       # print(old_N+l["N[0]"].pos)
       assert np.allclose(l["N[0]"].pos-dist, old_N, atol=5e-3)
       print(l["F[0]"].pos)
       print(old_F)
       assert np.allclose(l["F[0]"].pos, old_F, atol= 1e-10)
       assert np.allclose(l["O[0]"].pos, old_O, atol= 1e-10)
       assert np.allclose(l["H[0]"].pos, old_H, atol= 1e-10)

    def test_mirror_3pts(self, labeled_tetrahedral):
       l = labeled_tetrahedral
       l.translate([1,1,1])
       old_N = l["N[0]"].pos
       old_O = l["O[0]"].pos
       old_H = l["H[0]"].pos
       old_F = l["F[0]"].pos
       men = np.mean([old_O, old_H, old_N], axis=0)
       dist = 2*(men-old_F)
       l.mirror(mirror_plane_points=[["O"], ["H"], ["N"]])
       assert np.allclose(l["F[0]"].pos-dist, old_F, atol=5e-3)
       assert np.allclose(l["N[0]"].pos, old_N, atol= 1e-10)
       assert np.allclose(l["O[0]"].pos, old_O, atol= 1e-10)
       assert np.allclose(l["H[0]"].pos, old_H, atol= 1e-10)


    def test_mirror_pos_neg_vectors(self, labeled_tetrahedral):
       l = deepcopy(labeled_tetrahedral)
       l.mirror(about_vectors=((6,8,9),(-90,-5,-10)))
       labeled_tetrahedral.mirror(about_vectors=((6,8,9),(90,5,10)))
       for ii, jj in zip(l, labeled_tetrahedral):
           assert np.allclose(ii.pos, jj.pos, atol=1e-8)

    def test_mirror_non_unique_anchor(self, methane):
       with pytest.raises(MBuildError):
           methane.mirror(about_vectors=None, mirror_plane_points=None,
                          anchor_point=["H"])

    def test_mirror_non_unique_plane_point1(self, methane):
       with pytest.raises(MBuildError):
           methane.mirror(mirror_plane_points= [["C"], ["H"]],
                          about_vectors=[(1,1,1)])

    def test_mirror_non_unique_plane_point2(self, methane):
       with pytest.raises(MBuildError):
           methane.mirror(mirror_plane_points= [["H"],["H[0]"]],
                          about_vectors=[(1,1,1)])

    def test_mirror_anchor_exterior_coordinates(self, labeled_tetrahedral):
       l = labeled_tetrahedral
       anchor = (0,.1,0)
       l.mirror(anchor_point=anchor)
       assert np.allclose(l["C[0]"].pos, (0,.2,0), atol=1e-10)
       assert np.allclose(l["H[0]"].pos, (.1, .2, -.07), atol=1e-10)
       assert np.allclose(l["O[0]"].pos, (-.1, .2, -.07), atol=1e-10)
       assert np.allclose(l["F[0]"].pos, (0,.1,.07), atol=1e-10)
       assert np.allclose(l["N[0]"].pos, (0, .3, .07), atol=1e-10)

    def test_align_vectors_2D_default_anchor(self,benzene):
       c1 = []
       h1 = []
       c2 = []
       h2 = []
       for b in benzene:
           if b.name == 'H':
               h1.append(b.pos)
           else:
               c1.append(b.pos)
       benzene.align_vectors(align_these=(benzene["C[0]"].pos - benzene.center, (-.0117,.23333,1)-benzene.center),
                       with_these=(benzene["C[1]"].pos - benzene.center, (-.1325, .1637, 1)- benzene.center))
       for b in benzene:
           if b.name == "H":
               h2.append(b.pos)
           else:
               c2.append(b.pos)
       h1.append(h1.pop(0))
       c1.append(c1.pop(0))
       assert np.allclose(h1,h2, atol= 1e-5)
       assert np.allclose(c1,c2, atol =1e-5)

    def test_align_vectors_2D_defined_anchor(self, benzene):
       pre = benzene["C[2]"].pos
       anchor = ["C[4]"]
       benzene.align_vectors(align_these=(benzene["C[0]"].pos - benzene["C[4]"].pos, (0,0,1)),
                             with_these=(benzene["C[2]"].pos- benzene["C[4]"].pos, (0,0,3)),
                             anchor_pt=anchor)
       assert np.allclose(benzene["C[0]"].pos, pre, atol = 1e-4)

    def test_align_vectors_3D_default_anchor(self, labeled_tetrahedral):
       l = labeled_tetrahedral
       l.translate([1,1,1])
       old_O = l["O"][0].pos
       old_H = l["H"][0].pos
       old_N = l["N[0]"].pos
       old_F = l["F[0]"].pos
       l.align_vectors(align_these=(l["C[0]"].pos - l["F[0]"].pos, l["C[0]"].pos - l["H[0]"].pos),
                       with_these=(l["C[0]"].pos - l["O[0]"].pos, l["C[0]"].pos - l["N[0]"].pos),
                       anchor_pt=[1,1,1])
       assert np.allclose(l["C[0]"].pos, (1,1,1), atol=1e-10)
       assert np.allclose(l["O[0]"].pos, old_F, atol= 1e-10)
       assert np.allclose(l["N[0]"].pos, old_H, atol= 1e-10)
       assert np.allclose(l["F[0]"].pos, old_O, atol= 1e-10)
       assert np.allclose(l["H[0]"].pos, old_N, atol= 1e-10)


    def test_align_vectors_3D_defined_anchor1(self, labeled_tetrahedral):
       l = labeled_tetrahedral
       l.translate([1,1,1])
       anchor = l["N[0]"].pos
       old_O = l["O[0]"].pos
       old_H = l["H[0]"].pos
       old_F = l["F[0]"].pos
       l.align_vectors(align_these=(anchor - l["O[0]"].pos, anchor - l["H[0]"].pos),
                       with_these=(anchor - l["F[0]"].pos, anchor - l["O[0]"].pos),
                       anchor_pt=anchor)
       assert np.allclose(l["C[0]"].pos, (1,1,1), atol=8e-4)
       assert np.allclose(l["N[0]"].pos, anchor, atol= 1e-3)
       assert np.allclose(l["O[0]"].pos, old_F, atol= 1e-3)
       assert np.allclose(l["F[0]"].pos, old_H, atol= .00213)
       assert np.allclose(l["H[0]"].pos, old_O, atol= .0015)

    def test_align_vectors_3D_defined_anchor2(self, labeled_tetrahedral):
       l = labeled_tetrahedral
       l.translate([1,1,1])
       anchor = l["N[0]"].pos
       old_O = l["O[0]"].pos
       old_H = l["H[0]"].pos
       old_F = l["F[0]"].pos
       men = np.mean([old_O, old_H, anchor], axis=0)
       dist = 2*(men-old_F)
       l.align_vectors(align_these=(anchor - l["O[0]"].pos, anchor - l["H[0]"].pos),
                       with_these=(anchor - l["H[0]"].pos, anchor - l["O[0]"].pos),
                       anchor_pt=anchor)
       assert np.allclose(l["F[0]"].pos-dist, old_F, atol=5e-3)
       assert np.allclose(l["N[0]"].pos, anchor, atol= 1e-10)
       assert np.allclose(l["O[0]"].pos, old_H, atol= 1e-10)
       assert np.allclose(l["H[0]"].pos, old_O, atol= 1e-10)

    def test_subcompounds_by_name_or_label_particle_case1(self, benzene_from_parts):
       with pytest.raises(ValueError):
           list(benzene_from_parts.subcompounds_by_name_or_label(looking_for= "H"))

    def test_subcompounds_by_name_or_label_particle_case2(self, benzene_from_parts):
       with pytest.raises(ValueError):
           list(benzene_from_parts.subcompounds_by_name_or_label(looking_for= "C[0]"))

    def test_subcompounds_by_name_or_label_label_case1(self, benzene_from_parts):
       bet=0
       Bet=False
       for l in benzene_from_parts.subcompounds_by_name_or_label(looking_for="CH[2]"):
           if l:
               bet+=1
               if np.allclose([  1.73417655e-01,  -2.38361811e-01,   1.70108418e-17], l.pos):
                   Bet = True
       assert bet ==1
       assert Bet is True


    def test_subcompounds_by_name_or_label_hierarchy1(self, mixed_bilayer):
       tot = 0
       d = 0
       a = 0
       for money in list( t for t in mixed_bilayer.subcompounds_by_name_or_label(looking_for="AlkylMonomer") if t):
           tot+=1
           if "DSPC" in list(y.name for y in money.ancestors()):
               d+=1
           elif "ALC" in list(y.name for y in money.ancestors()):
               a+=1
       assert tot == 994
       assert a == 270
       assert d == 444

    def test_subcompounds_by_name_or_label_hierarchy2(self, mixed_bilayer):
       assert 18 == len(list(t for t in mixed_bilayer.subcompounds_by_name_or_label(looking_for="ALC") if t))
       assert 12 == len(list(t for t in mixed_bilayer.subcompounds_by_name_or_label(looking_for="DSPC") if t))

    def test_subcompounds_by_name_or_label_hierarchy3(self, mixed_bilayer):
       assert 2 == len(list(t for t in mixed_bilayer.subcompounds_by_name_or_label(looking_for="ALC[0]") if t))
       assert 2 == len(list(t for t in mixed_bilayer.subcompounds_by_name_or_label(looking_for="ALC[3]") if t))
       assert 2 == len(list(t for t in mixed_bilayer.subcompounds_by_name_or_label(looking_for="DSPC[0]") if t))
       assert 0 == len(list(t for t in mixed_bilayer.subcompounds_by_name_or_label(looking_for="DSPC[6]") if t))


    def test_find_subcompound_in_path_str1(self, mixed_bilayer):
       path = ["AlkylMonomer","DSPC[0]"]
       assert 74 == len(list(t for t in mixed_bilayer.find_subcompounds_in_path(path) if t))

    def test_find_subcompound_in_path_empty(self, mixed_bilayer):
       path = ["COOH", "ALC"]
       assert 0 == len(list(t for t in mixed_bilayer.find_subcompounds_in_path(path) if t))

    def test_find_subcompound_in_path_str2(self, mixed_bilayer):
       path = ("AlkylMonomer","DSPC")
       assert 444 == len(list(t for t in mixed_bilayer.find_subcompounds_in_path(path) if t))

    def test_find_subcompound_in_path_toError(self, mixed_bilayer):
       pass

    def test_find_subcompound_in_path_cmpnd(self, mixed_bilayer):
       parti = list(t for t in mixed_bilayer.find_subcompounds_in_path(["ALC[0]","top_leaflet"]) if t)[0]
       path1 = [ "OH", parti, "lipid_bilayer"]
       path2 = ["OH", parti, "ip"]
       path1parti = list(t for t in mixed_bilayer.find_subcompounds_in_path(pathway=path1) if t)
       path2parti = list(t for t in mixed_bilayer.find_subcompounds_in_path(pathway=path2) if t)
       assert 1 == len(path1parti)
       assert path1parti[0] is path2parti[0]

    def test_find_subcompound_in_path_list(self, mixed_bilayer):
       path = ["AlkylMonomer", ["DSPC[{}]".format(t) for t in range(4)], "top_leaflet"]
       assert len(list(t for t in mixed_bilayer.find_subcompounds_in_path(pathway=path) if t)) == 37*4

    def test_find_subcompound_in_path_tuple(self, mixed_bilayer):
       path = ["AlkylMonomer", tuple("DSPC[{}]".format(t) for t in range(2,6)), "top_leaflet"]
       assert len(list(t for t in mixed_bilayer.find_subcompounds_in_path(pathway=path) if t)) == 37*4

    def test_find_particles_in_path_hierarchy(self, mixed_bilayer):
       path = ["P", tuple("DSPC[{}]".format(t) for t in range(2,6)), "top_leaflet"]
       assert 4 == len(list(t for t in mixed_bilayer.find_particles_in_path(path) if t))

    def test_find_particles_in_path_yield_self1(self, mixed_bilayer):
       for p in mixed_bilayer.particles_by_name("C"):
           i = p
           break
       assert i is list(mixed_bilayer.find_particles_in_path([i]))[0]

    def test_find_particles_in_path_yield_self2(self, mixed_bilayer):
       for p in mixed_bilayer.particles_by_name("C"):
           i = p
           break
       assert i is list(mixed_bilayer.find_particles_in_path(i))[0]

    def test_find_particles_in_path_yield_self3(self, alc):
       assert np.allclose(list(alc.find_particles_in_path(["O"]))[0].pos, np.array([0,0,0]))

    def test_find_particles_in_path_Not_Found(self, mixed_bilayer):
       with pytest.raises(ValueError):
           path = ["AlkylMonomer","DSPC[0]"]
           list(mixed_bilayer.find_particles_in_path(within_path=path))

    def test_find_particles_in_path_bad_input(self, mixed_bilayer):
       with pytest.raises(ValueError):
           path = ["P", "P", "DSPC"]
           list(mixed_bilayer.find_particles_in_path(within_path=path))

    def test_find_particles_in_path_specific_location(self, benzene_from_parts):
       path = ['C', "CH[0]"]
       assert np.allclose(np.array([0,0,0]), list(benzene_from_parts.find_particles_in_path(path))[0].pos)

    def test_particles_in_range(self, ethane):
        group = ethane.particles_in_range(ethane[0], 0.141)
        assert sum([1 for x in group if x.name == 'H']) == 3
        assert sum([1 for x in group if x.name == 'C']) == 2

        group = ethane.particles_in_range(ethane[0], 0.141, max_particles=4)
        assert sum([1 for x in group if x.name == 'H']) == 3
        assert sum([1 for x in group if x.name == 'C']) == 1

    def test_generate_bonds(self, ch3):
        ch3.generate_bonds('H', 'H', dmin=0.01, dmax=2.0)
        assert ch3.n_bonds == 3 + 3

    def test_remove_from_box(self, ethane):
        n_ethanes = 5
        box = mb.fill_box(ethane, n_ethanes, [3, 3, 3])
        box.remove(box.children[3])

        n_ethanes -= 1
        assert box.n_particles == n_ethanes * ethane.n_particles
        assert len(box.children) == n_ethanes
        assert box.n_bonds == n_ethanes * ethane.n_bonds
        assert len([meth.referenced_ports()
                    for eth in box.children
                    for meth in eth.children]) == 2 * n_ethanes

    def test_remove(self, ethane):
        hydrogens = ethane.particles_by_name('H')
        ethane.remove(hydrogens)

        assert ethane.n_particles == 2
        assert ethane.n_bonds == 1
        for part in ethane.children:
            assert part.n_bonds == 0

        carbons = ethane.particles_by_name('C')
        ethane.remove(carbons)
        assert ethane.n_particles == 0
        assert ethane.n_bonds == 0
        assert len(ethane.children) == 2
        assert len(ethane.children[0].children) == 1  # Still contains a port

    def test_remove_many(self, ethane):
        ethane.remove([ethane.children[0], ethane.children[1]])

        assert ethane.n_particles == 1
        assert ethane._n_particles() == 0
        assert ethane.n_bonds == 0
        for part in ethane.children:
            assert isinstance(part, mb.Port)

    def test_remove_subcompound(self, ethane):
        methyl = ethane.children[0]
        ethane.remove(methyl)

        assert ethane.n_particles == 4
        assert ethane.n_bonds == 3
        assert len(ethane.children) == 1
        assert len(ethane.children[0].children) == 5  # Still contains a port

        methyl = ethane.children[0]
        ethane.remove(methyl)

        assert ethane.n_particles == 1
        assert ethane._n_particles() == 0
        assert ethane.n_bonds == 0
        assert len(ethane.children) == 0

    def test_remove_no_bond_graph(self):
        compound = mb.Compound()
        particle = mb.Compound(name='C', pos=[0, 0, 0])
        compound.add(particle, 'test-particle')
        compound.remove(particle)
        assert particle not in compound.particles()

    def test_remove_bond(self, ch3):
        ch_bond = list(ch3.bonds())[0]
        ch3.remove_bond(ch_bond)
        assert ch3.n_bonds == 2

        with pytest.warns(UserWarning):
            ch3.remove_bond(ch_bond)

    def test_center(self, methane):
        assert np.array_equal(methane.center, np.array([0, 0, 0]))
        for orientation in np.identity(3):
            separation = 0.2
            port = mb.Port(anchor=methane[0], orientation=orientation)
            assert np.allclose(port.center, np.array([0.0, 0.0, 0.0]), atol=1e-15)
            port = mb.Port(anchor=methane[0], orientation=orientation,
                           separation=separation)
            assert np.allclose(port.center, separation*orientation, atol=1e-15)
        np.random.seed(0)
        for orientation in np.random.rand(5, 3):
            port = mb.Port(anchor=methane[0], orientation=orientation)
            assert np.allclose(port.center, np.array([0.0, 0.0, 0.0]), atol=1e-15)
            port = mb.Port(anchor=methane[0], orientation=orientation,
                           separation=separation)
            assert np.allclose(port.center,
                               separation*orientation/np.linalg.norm(orientation),
                               atol=1e-15)

    def test_single_particle(self):
        part = mb.Particle(name='A')
        assert part.n_particles == 1
        assert len(list(part.particles())) == 1
        assert part.xyz.shape == (1, 3)
        assert part.root == part
        assert len(list(part.ancestors())) == 0
        assert next(part.particles_by_name('A')) == part

    def test_name(self):
        with pytest.raises(ValueError):
            mb.Compound(name=1)

    def test_my_label_custom(self, ch3):
       ch3["H[0]"].pos = [-.1, 0, -.07]
       ch3["H[1]"].pos = [0,.1,.07]
       ch3["H[2]"].pos = [0,-.1, .07]
       ch3.add(mb.Particle(name="H", pos=[.1,0,-.07]), label="odd_one")
       assert ch3[1].name == ch3[4].name
       assert ch3[4].my_label == "odd_one"
       assert ch3["odd_one"] is ch3[4]
       assert len(ch3["H"])  == 3

    def test_my_label_default(self,ch3):
       ch3["H[0]"].pos = [-.1, 0, -.07]
       ch3["H[1]"].pos = [0,.1,.07]
       ch3["H[2]"].pos = [0,-.1, .07]
       ch3.add(mb.Particle(name="H", pos=[.1,0,-.07]))
       assert len(ch3["H"]) == 4
       assert ch3["H[3]"] is ch3[4]
       assert ch3[4].my_label == "H[3]"
       assert ch3[1].name == ch3[4].name

    def test_my_label_parmed(self, alc):
       # dumalc=alc.to_parmed()
       # a = mb.Compound()
       # a.from_parmed(dumalc)
       #assert np.allclose(d["H[0]"].pos, (.1,0,-.07))
       # i did a test case with methane. it worked. When i tried to
       # do it with alc from parmed it did not preserve the hierarchy.
       # this should have been expected. the labeling convention still
       # is relevant for particles.
       pass


    def test_particle_in_particle(self):
        part = mb.Particle(name='A')
        parent = mb.Compound(part)

        assert part.n_particles == 1
        assert len(list(part.particles())) == 1
        assert part.xyz.shape == (1, 3)
        assert part.root == parent
        assert len(list(part.ancestors())) == 1
        assert next(part.particles_by_name('A')) == part

        assert parent.n_particles == 1
        assert len(list(parent.particles())) == 1
        assert parent.xyz.shape == (1, 3)
        assert parent.root == parent
        assert len(list(parent.ancestors())) == 0
        assert next(parent.particles_by_name('A')) == part

    def test_reload(self):
        from mbuild.examples.pmpc.brush import Brush
        from numpy import pi
        # Create a compound and write it to file.
        brush1 = Brush()
        brush1.save("brush1.pdb")

        # Create another compound, rotate it and write it to file.
        brush2 = Brush()
        mb.rotate(brush2, pi/2, [0, 0, 1])
        brush2.save("brush2.pdb")

        # Load brush2.pdb into brush1, modifying the atom positions of brush1.
        brush1.update_coordinates("brush2.pdb")
        brush1.save("modified_brush1.pdb")

        assert brush1['pmpc'].n_particles == 164
        assert brush1['pmpc'].n_bonds == 163
        assert len(brush1['pmpc']['monomer']) == 4
        assert brush1['pmpc']['monomer'][0].n_particles == 41
        assert brush1['pmpc']['monomer'][0].n_bonds == 40

    def test_to_trajectory(self, ethane, c3, n4):
        traj = ethane.to_trajectory()
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 1

        traj = ethane.to_trajectory(residues='CH3')
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 2
        assert 'CH3' in [res.name for res in traj.top.residues]
        assert all(res.n_atoms == 4 for res in traj.top.residues)

        traj = ethane.to_trajectory(chains='CH3')
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 2
        assert traj.n_residues == 2
        assert all(chain.n_atoms == 4 for chain in traj.top.chains)
        assert all(chain.n_residues == 1 for chain in traj.top.chains)

        system = mb.Compound([c3, n4])
        traj = system.to_trajectory(residues=['C', 'N'])
        assert traj.n_atoms == 2
        assert traj.top.n_bonds == 0
        assert traj.n_chains == 1
        assert traj.n_residues == 2

        traj = system.to_trajectory(chains=['C', 'N'])
        assert traj.n_atoms == 2
        assert traj.top.n_bonds == 0
        assert traj.n_chains == 2
        assert traj.n_residues == 2

        methyl = next(iter(ethane.children))
        traj = methyl.to_trajectory()
        assert traj.n_atoms == 4
        assert traj.top.n_bonds == 3
        assert traj.n_chains == 1
        assert traj.n_residues == 1

    def test_resnames_mdtraj(self, h2o, ethane):
        system = mb.Compound([h2o, mb.clone(h2o), ethane])
        traj = system.to_trajectory(residues=['Ethane', 'H2O'])
        residues = list(traj.top.residues)
        assert traj.n_residues == 3
        assert residues[0].name == 'H2O'
        assert residues[1].name == 'H2O'
        assert residues[2].name == 'Ethane'

        traj = system.to_trajectory(residues='Ethane')
        residues = list(traj.top.residues)
        assert traj.n_residues == 2
        assert residues[0].name == 'RES'
        assert residues[1].name == 'Ethane'

        traj = system.to_trajectory(residues=['Ethane'])
        residues = list(traj.top.residues)
        assert traj.n_residues == 2
        assert residues[0].name == 'RES'
        assert residues[1].name == 'Ethane'

        traj = system.to_trajectory()
        residues = list(traj.top.residues)
        assert traj.n_residues == 1
        assert residues[0].name == 'RES'

    def test_chainnames_mdtraj(self, h2o, ethane):
        system = mb.Compound([h2o, mb.clone(h2o), ethane])
        traj = system.to_trajectory(chains=['Ethane', 'H2O'])
        assert traj.n_chains == 3

        traj = system.to_trajectory(chains='Ethane')
        assert traj.n_chains == 2

        traj = system.to_trajectory(chains=['Ethane'])
        assert traj.n_chains == 2

        traj = system.to_trajectory()
        assert traj.n_chains == 1

    @pytest.mark.skipif(not has_intermol, reason="InterMol is not installed")
    def test_intermol_conversion1(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])

        intermol_system = compound.to_intermol()
        assert len(intermol_system.molecule_types) == 1
        assert 'Compound' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Compound'].bonds) == 9

        assert len(intermol_system.molecule_types['Compound'].molecules) == 1
        molecules = list(intermol_system.molecule_types['Compound'].molecules)
        assert len(molecules[0].atoms) == 11

    @pytest.mark.skipif(not has_intermol, reason="InterMol is not installed")
    def test_intermol_conversion2(self, ethane, h2o):
        # 2 distinct Ethane objects.
        compound = mb.Compound([ethane, mb.clone(ethane), h2o])

        molecule_types = [type(ethane), type(h2o)]
        intermol_system = compound.to_intermol(molecule_types=molecule_types)
        assert len(intermol_system.molecule_types) == 2
        assert 'Ethane' in intermol_system.molecule_types
        assert 'H2O' in intermol_system.molecule_types
        assert len(intermol_system.molecule_types['Ethane'].bonds) == 7
        assert len(intermol_system.molecule_types['H2O'].bonds) == 2

        assert len(intermol_system.molecule_types['Ethane'].molecules) == 2
        ethanes = list(intermol_system.molecule_types['Ethane'].molecules)
        assert len(ethanes[0].atoms) == len(ethanes[1].atoms) == 8

        assert len(intermol_system.molecule_types['H2O'].molecules) == 1
        h2os = list(intermol_system.molecule_types['H2O'].molecules)
        assert len(h2os[0].atoms) == 3

    def test_parmed_conversion(self, ethane, h2o):
        compound = mb.Compound([ethane, h2o])

        structure = compound.to_parmed()
        assert structure.title == 'Compound'

        structure = compound.to_parmed(title='eth_h2o')
        assert structure.title == 'eth_h2o'

        assert len(structure.atoms) == 11
        assert len([at for at in structure.atoms if at.element == 6]) == 2
        assert len([at for at in structure.atoms if at.element == 1]) == 8
        assert len([at for at in structure.atoms if at.element == 8]) == 1

        assert len(structure.bonds) == 9

        assert (sum(len(res.atoms) for res in structure.residues) ==
                len(structure.atoms))

        compound2 = mb.Compound()
        compound2.from_parmed(structure)

        assert compound2.n_particles == 11
        assert len([at for at in compound2.particles() if at.name == 'C']) == 2
        assert len([at for at in compound2.particles() if at.name == 'H']) == 8
        assert len([at for at in compound2.particles() if at.name == 'O']) == 1

        assert compound2.n_bonds == 9

    def test_resnames_parmed(self, h2o, ethane):
        system = mb.Compound([h2o, mb.clone(h2o), ethane])
        struct = system.to_parmed(residues=['Ethane', 'H2O'])
        assert len(struct.residues) == 3
        assert struct.residues[0].name == 'H2O'
        assert struct.residues[1].name == 'H2O'
        assert struct.residues[2].name == 'Ethane'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

        struct = system.to_parmed(residues=['Ethane', 'H2O'])
        assert len(struct.residues) == 3
        assert struct.residues[0].name == 'H2O'
        assert struct.residues[1].name == 'H2O'
        assert struct.residues[2].name == 'Ethane'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

        struct = system.to_parmed(residues='Ethane')
        assert len(struct.residues) == 2
        assert struct.residues[0].name == 'RES'
        assert struct.residues[1].name == 'Ethane'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

        struct = system.to_parmed()
        assert len(struct.residues) == 1
        assert struct.residues[0].name == 'RES'
        assert sum(len(res.atoms) for res in struct.residues) == len(struct.atoms)

    def test_parmed_element_guess(self):
        compound = mb.Particle(name='foobar')
        with pytest.warns(UserWarning):
            _ = compound.to_parmed()

        compound = mb.Particle(name='XXXXXX')
        with pytest.warns(UserWarning):
            _ = compound.to_parmed()

    def test_min_periodic_dist(self, ethane):
        compound = mb.Compound(ethane)
        C_pos = np.array([atom.pos for atom in list(compound.particles_by_name('C'))])
        assert round(compound.min_periodic_distance(C_pos[0], C_pos[1]), 2) == 0.14
        compound.periodicity = np.array([0.2, 0.2, 0.2])
        assert round(compound.min_periodic_distance(C_pos[0], C_pos[1]), 2) == 0.06

    def test_bond_graph(self, ch3):
        compound = mb.Compound()
        compound.add(ch3)
        assert compound.n_bonds == 3
        assert all(compound.bond_graph.has_node(particle)
                   for particle in ch3.particles())

        ch3_nobonds = mb.clone(ch3)
        for bond in ch3_nobonds.bonds():
            ch3_nobonds.remove_bond(bond)
        compound.add(ch3_nobonds)
        assert compound.n_bonds == 3
        assert not any(compound.bond_graph.has_node(particle)
                       for particle in ch3_nobonds.particles())

        carbons = list(compound.particles_by_name('C'))
        compound.add_bond((carbons[0], carbons[1]))
        assert compound.n_bonds == 4
        assert all(compound.bond_graph.has_node(particle)
                   for particle in carbons)
        assert any(compound.bond_graph.has_node(particle)
                   for particle in ch3_nobonds.particles())

        compound.remove_bond((carbons[0], carbons[1]))
        assert not any(compound.bond_graph.has_node(particle)
                       for particle in ch3_nobonds.particles())

    def test_update_coords_update_ports(self, ch2):
        distances = np.round([ch2.min_periodic_distance(port.pos, ch2[0].pos)
                              for port in ch2.referenced_ports()], 5)
        orientations = np.round([port.pos - port.anchor.pos
                                 for port in ch2.referenced_ports()], 5)

        ch2_clone = mb.clone(ch2)
        ch2_clone[0].pos += [1, 1, 1]
        ch2_clone.save('ch2-shift.pdb')

        ch2.update_coordinates('ch2-shift.pdb')
        updated_distances = np.round([ch2.min_periodic_distance(port.pos, ch2[0].pos)
                                      for port in ch2.referenced_ports()], 5)
        updated_orientations = np.round([port.pos - port.anchor.pos
                                         for port in ch2.referenced_ports()], 5)

        assert np.array_equal(distances, updated_distances)
        assert np.array_equal(orientations, updated_orientations)

    def test_charge(self, ch2, ch3):
        compound = mb.Compound(charge=2.0)
        assert compound.charge == 2.0
        compound2 = mb.Compound()
        assert compound2.charge == 0.0

        ch2[0].charge = 0.5
        ch2[1].charge = -0.25
        ch3[0].charge = 1.0
        compound.add([ch2, ch3])
        assert compound.charge == 1.25
        assert ch2.charge == 0.25
        assert compound[0].charge == 0.5

        with pytest.raises(AttributeError):
            compound.charge = 2.0

    def test_charge_subcompounds(self, ch2, ch3):
        ch2[0].charge = 0.5
        ch2[1].charge = -0.25
        compound = mb.Compound(subcompounds=ch2)
        assert compound.charge == 0.25

        with pytest.raises(MBuildError):
            compound = mb.Compound(subcompounds=ch3, charge=1.0)

    def test_charge_neutrality_warn(self, benzene):
        benzene[0].charge = 0.25
        with pytest.warns(UserWarning):
            benzene.save('charge-test.mol2')

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization(self, octane):
        octane.energy_minimization()

    @pytest.mark.skipif(has_openbabel, reason="Open Babel package is installed")
    def test_energy_minimization_openbabel_warn(self, octane):
        with pytest.raises(MBuildError):
            octane.energy_minimization()

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_ff(self, octane):
        for ff in ['UFF', 'GAFF', 'MMFF94', 'MMFF94s', 'Ghemical']:
            octane.energy_minimization(forcefield=ff)
        with pytest.raises(MBuildError):
            octane.energy_minimization(forcefield='fakeFF')

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_algorithm(self, octane):
        for algorithm in ['cg', 'steep', 'md']:
            octane.energy_minimization(algorithm=algorithm)
        with pytest.raises(MBuildError):
            octane.energy_minimization(algorithm='fakeAlg')

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_non_element(self, octane):
        for particle in octane.particles():
            particle.name = 'Q'
        with pytest.raises(MBuildError):
            octane.energy_minimization()

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel package not installed")
    def test_energy_minimization_ports(self, octane):
        distances = np.round([octane.min_periodic_distance(port.pos, port.anchor.pos)
                              for port in octane.all_ports()], 5)
        orientations = np.round([port.pos - port.anchor.pos
                                 for port in octane.all_ports()], 5)

        octane.energy_minimization()

        updated_distances = np.round([octane.min_periodic_distance(port.pos,
                                                                   port.anchor.pos)
                                      for port in octane.all_ports()], 5)
        updated_orientations = np.round([port.pos - port.anchor.pos
                                         for port in octane.all_ports()], 5)

        assert np.array_equal(distances, updated_distances)
        assert np.array_equal(orientations, updated_orientations)

    def test_clone_outside_containment(self, ch2, ch3):
        compound = mb.Compound()
        compound.add(ch2)
        mb.force_overlap(ch3, ch3['up'], ch2['up'])
        with pytest.raises(MBuildError):
            ch3_clone = mb.clone(ch3)

    def test_load_mol2_mdtraj(self):
        with pytest.raises(KeyError):
            mb.load(get_fn('benzene-nonelement.mol2'))
        mb.load(get_fn('benzene-nonelement.mol2'), use_parmed=True)

    def test_siliane_bond_number(self, silane):
        assert silane.n_bonds == 4
