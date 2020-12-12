import pymol

def generate_3structure_image(pdb_id,pdbs_to_process,output_path,logger):
    pymol.cmd.fetch(' '.join(pdbs_to_process))    

    pymol.cmd.alignto(pdb_id,object="all_to_" + pdb_id)

    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)

    pymol.cmd.save(output_path + "/"+pdb_id + ".pse")
    pymol.cmd.png(output_path + "/" + "%s.png"%(pdb_id))
