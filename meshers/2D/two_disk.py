import gmsh
import math
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Creates two circle meshes positioned on `ypos` and `xpos1`/`xpos2`")
    parser.add_argument('--radius_1', type=float, default=1.0, help='radius of ball 1')
    parser.add_argument('--radius_2', type=float, default=1.0, help='radius of ball 2')
    parser.add_argument('--xpos_1', type=float, default=1.25, help='x position of ball 1')
    parser.add_argument('--xpos_2', type=float, default=4.25, help='x position of ball 2')
    parser.add_argument('--ypos', type=float, default=1.25, help='y position of both balls')
    parser.add_argument('-m', '--mesh_density', type=int, default=4, help='the \"density\" of the mesh - defines the number of elements on each side of the square') 
    parser.add_argument('-o', '--output', type=str, default="meshes/2D/", help='the output directory of the mesh file')
    parser.add_argument('-g', '--graphics', action='store_true', help='if --graphics is given, display the GMSH GUI of the built mesh')
    parser.add_argument('--mesh_algorithm', type=int, default=8, help='the GMSH meshing algorithm number (see gmsh.info for details)')
    parser.add_argument('--recomb_algorithm', type=int, default=1, help='the GMSH recombination algorithm number (see gmsh.info for details)')
    parser.add_argument('--subdiv_algorithm', type=int, default=0, help='the GMSH subdivision algorithm number (see gmsh.info for details)')   

    args = parser.parse_args()
    return args


def create_mesh(args):
    print("MSHPY | Initializing GMSH")
    gmsh.initialize()
    gmsh.model.add("TwoDisk")

    yposition = args.ypos
    xposition_1 = args.xpos_1
    xposition_2 = args.xpos_2
    radius_1 = args.radius_1
    radius_2 = args.radius_2

    file_name = args.output + f"twodisk_r{radius_1}_r{radius_2}_xa{xposition_1}_xb{xposition_2}_msh{args.mesh_density}.msh"

    print("MSHPY | Creating the circle geometries")

    gmsh.model.occ.addDisk(xposition_1, yposition, 0.0, radius_1, radius_1, tag=1)
    gmsh.model.occ.addDisk(xposition_2, yposition, 0.0, radius_2, radius_2, tag=2)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 1.0 / args.mesh_density)

    # gmsh.model.mesh.setTransfiniteCurve(1, args.mesh_density + 1)
    # gmsh.model.mesh.setTransfiniteCurve(2, args.mesh_density + 1)
    # gmsh.model.mesh.setTransfiniteCurve(3, args.mesh_density + 1)
    # gmsh.model.mesh.setTransfiniteCurve(4, args.mesh_density + 1)

    bulk_objects = [1, 2]
    surf_objects = [1, 2]
    print("MSHPY | Defining GMSH Physical Objects")
    bulktag = gmsh.model.addPhysicalGroup(2, bulk_objects, name="Bulk")
    surftag = gmsh.model.addPhysicalGroup(1, surf_objects, name="Surface")
    print("MSHPY |     Bulk Physical Tag:", str(bulktag))
    print("MSHPY |     Surface Physical Tag:", str(surftag))

    print("MSHPY | Creating 2D Mesh with Algorithms:")
    print("MSHPY |     Meshing Algorithm", args.mesh_algorithm)
    print("MSHPY |     Recombination Algorithm", args.recomb_algorithm)
    print("MSHPY |     Subdivision Algorithm", args.subdiv_algorithm)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", args.mesh_algorithm)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", args.recomb_algorithm)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", args.subdiv_algorithm)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()
    # gmsh.model.mesh.refine()
    
    print("MSHPY | Written mesh to", file_name)
    gmsh.write(file_name)

    if args.graphics:
        gmsh.fltk.run()
    
    print("MSHPY | Closing GMSH instance")
    gmsh.finalize()

    return

if __name__ == "__main__":
    args = parse_arguments()
    print("MSHPY | Creating .msh File With These Parameters:")
    print("MSHPY |     Mesh Density =", args.mesh_density)
    print("MSHPY |     Output Directory:", args.output)
    print("MSHPY |     Meshing Algorithm =", args.mesh_algorithm)
    print("MSHPY |     Recombination Algorithm =", args.recomb_algorithm)
    print("MSHPY |     Subdivision Algorithm =", args.subdiv_algorithm)
    create_mesh(args)
    print("MSHPY | Finished Mesh Creation")
