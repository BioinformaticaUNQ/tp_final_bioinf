import pymol

def generate_3structure_image(pdb_id,pdbs_to_process,output_path):
    pymol.cmd.fetch(' '.join(pdbs_to_process))    
    # for pdb in pdbs_to_process:
    #     pymol.cmd.enable(pdb)
    pymol.cmd.alignto(pdb_id,object="all_to_" + pdb_id)
    print(pymol.cmd.get_names())
    #pymol.cmd.disable("all")

    pymol.cmd.hide('all')
    #pymol.cmd.enable("all_to_" + pdb_id)
    pymol.cmd.show('cartoon')
    pymol.cmd.set('ray_opaque_background', 0)

    # for pdb in pdbs_to_process:
    #     if pdb_id != pdb:
    #         pymol.cmd.align(pdb_id,pdb)
    pymol.cmd.save(output_path + "/"+pdb_id + ".pse")
    pymol.cmd.png(output_path + "/" + "%s.png"%(pdb_id))
    #pymol.cmd.quit()
