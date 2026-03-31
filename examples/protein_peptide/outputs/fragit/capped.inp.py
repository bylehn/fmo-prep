from pymol import stored

fragments_data="1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28:29,30,31,32,33,34,35:36,37,38,39,40,41,42:43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64:65,66,67,68,69,70,71,72,73:172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,197,198,199,200,201,202,203,204,205,206,207,208:209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225:258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285:286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305:306,307,308,309,310,311,312,313,316,317,318,319,320,321,322,323,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468:314,315,324,325,326,327,328,329,330,331,332,333,334,335:336,337,338,339,340,341,342:343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359:360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379:380,381,382,383,384,385,386:387,388,389,390,391,392,393,394,395,396:397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418:419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439:11:74:75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95:96:256:257:97:98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,118,119,120,121,122,123,124,125,126,127,128,129,130,131:132:117:133,134,135,136,137,138,139,140,141,142,143,144,145,146,147:148:149:150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169:170:171:196:226:227,228,229,230,231,232,233,234,235,236:237:238:239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255"

buffer_data=""

active_data=""

backbone_data="1,2,3,4,9,10,12,13,14,15,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,65,66,67,68,69,70,71,72,73,75,76,77,78,94,95,98,99,100,101,115,116,118,119,120,121,130,131,133,134,135,136,146,147,150,151,152,153,168,169,172,173,174,175,194,195,197,198,207,208,209,210,211,212,213,214,224,225,227,228,229,230,235,236,239,240,241,242,254,255,262,264,265,266,267,286,287,288,289,290,291,306,307,312,313,314,315,322,323,324,325,326,327,336,337,338,339,340,341,342,343,344,345,346,347,348,360,361,362,363,364,365,380,381,382,383,384,385,386,387,388,389,390,394,395,396,397,398,399,400,404,418,419,420,421,422,423,424,440,441,442,443,452,453,454,455,469,473,477,481,486,490,495,499"
fragment_charges="0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"


def atom_data_to_lists(l):
    lists = list()
    if len(l) == 0: return []
    fraglists = l.split(":")
    for fraglist in fraglists:
        atmlists = fraglist.split(",")
        lists.append(list(map(int, atmlists)))
    return lists


def fragment_data_to_list(l):
    if len(l) == 0: return []
    fragment_properties = map(int, l.split(","))
    return fragment_properties


def select_fragment_by_id(id):
    fragments = atom_data_to_lists(fragments_data)
    idx = int(id) - 1
    selection = "select fragment-%03i, none" % (idx + 1)
    string = "".join([" or id %i" % atom for atom in fragments[idx]])
    cmd.do(selection + string)
    return "fragment-%03i" % (idx + 1)


def select_backbone():
    select_region("backbone", backbone_data)


def flatten_list(list_to_flatten):
    flat_list = []
    for items in list_to_flatten:
        flat_list.extend(items)
    return flat_list


def select_active_region():
    select_region("active", active_data)


def select_buffer_region():
    select_region("buffer", buffer_data)


def select_frozen_region():
    all_atoms = flatten_list(atom_data_to_lists(fragments_data))
    buffer_atoms = flatten_list(atom_data_to_lists(buffer_data))
    if len(buffer_atoms) == 0: return
    frozen_atoms = []
    for atom in all_atoms:
        if atom not in buffer_atoms:
            frozen_atoms.append(atom)
    frozen_data = ",".join(map(str, frozen_atoms))
    select_region("frozen", frozen_data)


def select_region(type, data, cut=40):
    atoms = flatten_list(atom_data_to_lists(data))
    if len(atoms) == 0: return
    selection_string = "select sele-%s, %s"
    sel_name = "none"
    while len(atoms) > 0:
        t_atoms = atoms[:cut]
        atoms = atoms[cut:]
        selection = selection_string % (type, sel_name)
        string = "".join([" or id %i" % atom for atom in t_atoms])
        cmd.do(selection + string)
        sel_name = "sele-%s" % type


def make_selection(type="fragment", id="1"):
    if type == "fragment" or type == "fid":
        select_fragment_by_id(id)
    elif type == "backbone" or type == "bb":
        select_backbone()
    elif type == "active":
        select_active_region()
    elif type == "buffer":
        select_buffer_region()
    elif type == "frozen":
        select_frozen_region()
    else:
        cmd.do("select sele, none")


def make_fragment_selections():
    fragments = atom_data_to_lists(fragments_data)
    for i, fragments in enumerate(fragments):
        selection_name = select_fragment_by_id("%i" % (i + 1))
        cmd.do("pseudoatom lbl-frag%i, selection='%s', label='Frag-%03i'" % (i + 1, selection_name, i + 1))
    cmd.do("group fragments, fragment-*")
    cmd.do("group labels, lbl-*")
    cmd.do("disable labels")


def make_selections():
    make_selection("active")
    make_selection("buffer")
    make_selection("frozen")
    select_backbone()
    make_fragment_selections()
    cmd.do("group selections, sele-*")


def color_selection(sel="all", color="green"):
    cmd.color(color, sel)


def color_atoms(data, color="green"):
    selection = "all"
    for d in data:
        sel = selection + " and id %i" % (d)
        color_selection(sel, color)


def hex_to_float(value):
    (ri, gi, bi) = hex_to_rgb(value)
    r = ri / 255.0
    g = gi / 255.0
    b = bi / 255.0
    return (r, g, b)


def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + 2], 16) for i in [0, 2, 4])


def get_colors_for_fragments(list_of_fragments):
    color_names = ["color1", "color2", "color3", "color4", "color5", "color6", "color7", "color8", "color9", "color10"]
    colors_hex = ["#006837", "#1A9850", "#66BD63", "#A6D96A", "#D9EF8B", "#F46D43", "#ED5D3C", "#E0422F", "#D73027",
                  "#A50026"]
    for c, n in zip(colors_hex, color_names):
        (r, g, b) = hex_to_float(c)
        cmd.do("set_color %s, [%3.1f, %3.1f, %3.1f]" % (n, r, g, b))
    col = color_names[:]
    while len(col) < len(list_of_fragments):
        col.extend(color_names)
    return col


def color_all_fragments():
    frags = atom_data_to_lists(fragments_data)
    colors = get_colors_for_fragments(frags)
    for i, atomlist in enumerate(frags):
        color_atoms(atomlist, colors[i])


def color_fragments_by_charge():
    """ Colors fragments by charges """
    charges = fragment_data_to_list(fragment_charges)
    fragments = atom_data_to_lists(fragments_data)
    for i, (fragment, charge) in enumerate(zip(fragments, charges)):
        if charge == -1:
            color_atoms(fragment, "red")
        elif charge == +1:
            color_atoms(fragment, "blue")
        else:
            color_atoms(fragment, "white")


def color_fragments(sele="fragments"):
    """ PyMOL GUI coloring function

        This function is invoked by the 'ColorFragments' option
        in the PyMOL user interface.
    """
    cmd.bg_color("white")
    if sele == "fragments":
        color_all_fragments()

    elif sele == "buffer" or sele == "layers":
        cmd.do("color gray, sele-frozen")
        cmd.do("color marine, sele-buffer")

    elif sele == "active":
        cmd.do("color gray, sele-frozen")
        cmd.do("color marine, sele-buffer")
        cmd.do("color raspberry, sele-active")

    elif "charge" in sele:
        color_fragments_by_charge()


# iterate over atoms in a fragment
def name_all_fragments():
    cmd.do("enable labels")


def name_fragments(action="show"):
    if action == "show":
        name_all_fragments()


def setup_rendering_defaults():
    cmd.do("hide spheres")
    cmd.do("set antialias, 2")
    cmd.do("set ray_trace_mode, 1")
    cmd.do("set ray_shadow, off")

# default commands we need to execute
# to set up the environment correctly
cmd.do("load /home/fabian/Documents/kuano/fmo-prep/examples/protein_peptide/outputs/capped.pdb")


color_fragments("fragments")
make_selections()

cmd.do("show sticks, all")
cmd.do("select sele, none")
cmd.extend("NameFragments", name_fragments)
cmd.extend("ColorFragments", color_fragments)
cmd.extend("SelectFragments", make_selection)

setup_rendering_defaults()
