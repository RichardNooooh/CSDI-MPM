import gmsh
import math
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Creates a 2D mesh a box.")
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
    gmsh.model.add("rectangle")

    file_name = args.output + "rectangle" + "_msh" + str(args.mesh_density) + ".msh"

    length = 4
    height = 1

    print("MSHPY | Creating the rectangle geometry")
    gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, length, height, tag=1)
    gmsh.model.occ.synchronize()

    # this is the modified mesh density following refinement so that args.mesh_density matches the actual mesh
    # refined_density = int(args.mesh_density / 2)

    # gmsh.model.mesh.setTransfiniteCurve(1, math.floor(math.pi * refined_density) + 1)
    # gmsh.model.mesh.setTransfiniteCurve(2, 2*args.mesh_density + 1)
    # gmsh.model.mesh.setTransfiniteSurface(1)
    mesh_density = args.mesh_density
    gmsh.model.mesh.setTransfiniteCurve(1, length * mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(2, height * mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(3, length * mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(4, height * mesh_density + 1)
    gmsh.model.mesh.setTransfiniteSurface(1)

    # gmsh.model.mesh.setSize(gmsh.model.getEntities(0), 1 / refined_density)

    bulk_objects = [1, 2]
    surf_objects = [1, 2, 3, 4]
    print("MSHPY | Defining GMSH Physical Objects")
    bulktag = gmsh.model.addPhysicalGroup(2, bulk_objects, name="Bulk")
    surftag = gmsh.model.addPhysicalGroup(1, surf_objects, name="Surface")
    print("MSHPY |     Bulk Physical Tag:", str(bulktag))
    print("MSHPY |     Surface Physical Tag:", str(surftag))

    print("MSHPY | Creating 2D Mesh with Algorithms:")
    print("MSHPY |     Meshing Algorithm", args.mesh_algorithm)
    print("MSHPY |     Recombination Algorithm", args.recomb_algorithm)
    print("MSHPY |     Subdivision Algorithm", args.subdiv_algorithm)

    gmsh.option.setNumber("Mesh.Algorithm", args.mesh_algorithm)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", args.recomb_algorithm)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.recombine()

    # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", args.subdiv_algorithm)
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
    # print("MSHPY |     Meshing Algorithm =", args.mesh_algorithm)
    # print("MSHPY |     Recombination Algorithm =", args.recomb_algorithm)
    # print("MSHPY |     Subdivision Algorithm =", args.subdiv_algorithm)
    create_mesh(args)
    print("MSHPY | Finished Mesh Creation")
