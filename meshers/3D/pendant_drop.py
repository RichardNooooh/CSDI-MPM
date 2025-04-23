import gmsh
from math import cos, sin, pi, floor
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Creates a 3D mesh of 1/4th of a spherical pendant")
    parser.add_argument('-m', '--mesh_density', type=int, default=32, help='the \"density\" of the mesh - defines the number of elements on each side of the square') 
    parser.add_argument('--ceiling_diff', type=float, default=0.65, help='the distance of the center of the sphere to z=1.0')
    parser.add_argument('--radius', type=float, default=0.75, help='the radius of the pendant')
    parser.add_argument('-o', '--output', type=str, default="meshes/3D/", help='the output directory of the mesh file')
    parser.add_argument('-g', '--graphics', action='store_true', help='if --graphics is given, display the GMSH GUI of the built mesh')

    args = parser.parse_args()
    return args


def create_mesh(args):
    print("MSHPY | Initializing GMSH")
    gmsh.initialize()
    gmsh.model.add("PendantDrop")

    file_name = args.output + "pendant_msh" + str(args.mesh_density) + ".msh"

    print("MSHPY | Creating the catenoid geometry")
    gmsh.model.occ.addSphere(0.0, 0.0, 3.0-args.ceiling_diff,
                             args.radius, angle1=-pi/2, angle2=pi/2, angle3=pi/2, tag=1)
    gmsh.model.occ.addBox(-1.0, -1.0, 3.0, 2.0, 2.0, 2.0, tag=2)
    gmsh.model.occ.cut([(3, 1)], [(3, 2)], tag=3)

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    mesh_size = 1 / args.mesh_density
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

    bulk_objects = [3]
    surf_objects = [1, 2, 3, 4]
    print("MSHPY | Defining GMSH Physical Objects")
    bulktag = gmsh.model.addPhysicalGroup(3, bulk_objects, name="Bulk")
    surftag = gmsh.model.addPhysicalGroup(2, surf_objects, name="Surface")
    print("MSHPY |     Bulk Physical Tag:", str(bulktag))
    print("MSHPY |     Surface Physical Tag:", str(surftag))


    gmsh.model.mesh.generate(3)
    
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
