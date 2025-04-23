import gmsh
import math
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Creates a 3D mesh of a cube.")
    parser.add_argument('-m', '--mesh_density', type=int, default=2, help='the \"density\" of the mesh - defines the number of elements on each side of the square') 
    parser.add_argument('-o', '--output', type=str, default="meshes/3D/", help='the output directory of the mesh file')
    parser.add_argument('-g', '--graphics', action='store_true', help='if --graphics is given, display the GMSH GUI of the built mesh')

    args = parser.parse_args()
    return args


def create_mesh(args):
    print("MSHPY | Initializing GMSH")
    gmsh.initialize()
    gmsh.model.add("CubeDroplet")

    file_name = args.output + "cube" +"_msh" + str(args.mesh_density) + ".msh"

    print("MSHPY | Creating the cube geometry")
    gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1)
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteCurve(1, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(2, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(3, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(4, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(5, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(6, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(7, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(8, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(9, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(10, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(11, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteCurve(12, args.mesh_density + 1)
    gmsh.model.mesh.setTransfiniteSurface(1)
    gmsh.model.mesh.setTransfiniteSurface(2)
    gmsh.model.mesh.setTransfiniteSurface(3)
    gmsh.model.mesh.setTransfiniteSurface(4)
    gmsh.model.mesh.setTransfiniteSurface(5)
    gmsh.model.mesh.setTransfiniteSurface(6)
    gmsh.model.mesh.setTransfiniteVolume(1)

    bulk_objects = [1]
    surf_objects = [1, 2, 3, 4, 5, 6]
    print("MSHPY | Defining GMSH Physical Objects")
    bulktag = gmsh.model.addPhysicalGroup(3, bulk_objects, name="Bulk")
    surftag = gmsh.model.addPhysicalGroup(2, surf_objects, name="Surface")
    print("MSHPY |     Bulk Physical Tag:", str(bulktag))
    print("MSHPY |     Surface Physical Tag:", str(surftag))

    # gmsh.option.setNumber("Mesh.RecombineAll", 1)
    # gmsh.option.setNumber("Mesh.Algorithm", args.mesh_algorithm)
    # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", args.recomb_algorithm)
    # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", args.subdiv_algorithm)

    gmsh.model.mesh.generate(3)
    # gmsh.model.mesh.recombine()
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
    create_mesh(args)
    print("MSHPY | Finished Mesh Creation")
