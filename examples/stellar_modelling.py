"""
Example for modifying ATLAS, DETAIL and SURFACE command files.
Not all parameters can be edited this way, but the calculation of a grid, everything relevant is editable.
"""
from pkg_resources import resource_filename
from astro_scripts_uibk import stellar_modelling


def main():
    # define your stellar parameters
    t_eff = 9000
    log_g = 1.4
    x = 2
    he = 0.1
    z = -0.1
    abun = {'C': -4.2, 'MG': -3.4}

    # A surface file can have multiple detail pops files as input.
    # Make sure that the keywords (e.g. 'cpops' and 'mgpops') are already in the surface file.
    detail_pops_dict = {'cpops': 'path_to_detail_c_file', 'mgpops': 'path_to_detail_c_file'}

    # define paths for output files
    atlas_struct = 'path_to_file.struct'
    atlas_control = 'path_to_file.atout'
    detail_out = 'path_to_file.dout'
    detail_pops = 'path_to_file.dpops'
    surface_out = 'path_to_file.sout'
    surface_flux = 'path_to_file.flux'

    # ATLAS
    # define paths of command files
    a_com_path = resource_filename("astro_scripts_uibk", "test_data/asu_test.acom")
    d_com_path = resource_filename("astro_scripts_uibk", "test_data/asu_test.dcom")
    s_com_path = resource_filename("astro_scripts_uibk", "test_data/asu_test.scom")

    # load the atlas command file
    with open(a_com_path) as f:
        lines = f.readlines()

    # modify lines
    lines = stellar_modelling.modify_atlas_com(lines, t_eff=t_eff, log_g=log_g, microturbulence=x, z=z, he=he,
                                               out_control_path=atlas_control, out_struct_path=atlas_struct)

    # save the modified atlas command file
    with open(a_com_path + '.tmp', 'w') as f:
        f.writelines(lines)

    # DETAIL
    # load the detail command file
    with open(d_com_path) as f:
        lines = f.readlines()

    # modify lines
    lines = stellar_modelling.modify_detail_com(lines, microturbulence=x, abundance_dict=abun,
                                                atlas_struct_path=atlas_struct, out_control_path=detail_out,
                                                out_pops_path=detail_pops)

    # save the modified detail command file
    with open(d_com_path + '.tmp', 'w') as f:
        f.writelines(lines)

    # SURFACE
    # load the surface command file
    with open(s_com_path) as f:
        lines = f.readlines()

    # modify lines
    lines = stellar_modelling.modify_surface_com(lines, microturbulence=x, abundance_dict=abun,
                                                 atlas_struct_path=atlas_struct, detail_pops_dict=detail_pops_dict,
                                                 out_control_path=surface_out, out_flux_path=surface_flux)

    # save the modified detail command file
    with open(s_com_path + '.tmp', 'w') as f:
        f.writelines(lines)


if __name__ == '__main__':
    main()
