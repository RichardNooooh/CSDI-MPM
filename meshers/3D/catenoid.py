import gmsh
from math import cos, sin, pi, floor
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Creates a 3D mesh of 1/4th of a tube")
    parser.add_argument('-m', '--mesh_density', type=int, default=8, help='the \"density\" of the mesh - defines the number of elements on each side of the square') 
    parser.add_argument('--height', type=float, default=2.0, help='the height of the cylinder')
    parser.add_argument('--radius', type=float, default=2.0, help='the outer radius of the tube')
    parser.add_argument('-o', '--output', type=str, default="meshes/3D/", help='the output directory of the mesh file')
    parser.add_argument('-g', '--graphics', action='store_true', help='if --graphics is given, display the GMSH GUI of the built mesh')

    args = parser.parse_args()
    return args


def create_mesh(args):
    print("MSHPY | Initializing GMSH")
    gmsh.initialize()
    gmsh.model.add("Tube")

    file_name = f"{args.output}catenoid_msh{args.mesh_density}.msh"
    mesh_density = args.mesh_density

    print("MSHPY | Creating the catenoid geometry")
    gmsh.model.occ.addCylinder(0.0, 0.0, 0.0,
                               0.0, 0.0, args.height,
                               args.radius, angle=pi/2)
    gmsh.model.occ.remove([(3, 1), 
                           (2, 2), (2, 3), (2, 4), (2, 5), 
                           (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), 
                           (0, 5), (0, 6)])
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteCurve(1, floor(pi / 2 * mesh_density * args.radius) + 1)
    gmsh.model.mesh.setTransfiniteCurve(3, floor(pi / 2 * mesh_density * args.radius) + 1)
    gmsh.model.mesh.setTransfiniteCurve(2, int(args.height * mesh_density) + 1)
    gmsh.model.mesh.setTransfiniteCurve(4, int(args.height * mesh_density) + 1)
    gmsh.model.mesh.setTransfiniteSurface(1)

    # mesh_size = 1 / args.mesh_density
    # gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

    # bulk_objects = [1]
    surf_objects = [1]
    print("MSHPY | Defining GMSH Physical Objects")
    # bulktag = gmsh.model.addPhysicalGroup(3, bulk_objects, name="Bulk")
    surftag = gmsh.model.addPhysicalGroup(2, surf_objects, name="Surface")
    # print("MSHPY |     Bulk Physical Tag:", str(bulktag))
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
    print("MSHPY |     Mesh Density = ", f"{args.mesh_density}")
    print("MSHPY |     Output Directory:", args.output)
    create_mesh(args)
    print("MSHPY | Finished Mesh Creation")
