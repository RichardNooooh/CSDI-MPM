import gmsh
from math import cos, sin, pi, floor
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Creates a 3D mesh of weird fork thing")
    parser.add_argument('-m', '--mesh_density', type=int, default=8, help='the \"density\" of the mesh - defines the number of elements on each side of the square') 
    parser.add_argument('-o', '--output', type=str, default="meshes/3D/", help='the output directory of the mesh file')
    parser.add_argument('-g', '--graphics', action='store_true', help='if --graphics is given, display the GMSH GUI of the built mesh')

    args = parser.parse_args()
    return args


def create_mesh(args):
    print("MSHPY | Initializing GMSH")
    gmsh.initialize()
    gmsh.model.add("Cylinder")

    file_name = f"{args.output}forktubes_msh{args.mesh_density}.msh"

    print("MSHPY | Creating the geometry")

    translateX = 2.25
    translateY = 1.25
    
    gmsh.model.occ.addCylinder(0 + translateX, 0 + translateY, 0.0,
                               0.0, 0.0, 0.50,
                               1.0)
    gmsh.model.occ.remove([(3, 1), (2, 3)])

    gmsh.model.occ.addCircle(0.5 + translateX, 0.0 + translateY, 0.5, 0.40, tag=100)
    gmsh.model.occ.addDisk(0.5 + translateX, 0.0 + translateY, 0.5, 0.40, 0.40, tag=500)
    gmsh.model.occ.addCurveLoop([100], 100)
    gmsh.model.occ.addCircle(1.1 + translateX, 0.0 + translateY, 1.0, 0.8, 200)
    gmsh.model.occ.addCurveLoop([200], 200)
    gmsh.model.occ.addThruSections([100, 200], 1000, makeSolid=False, makeRuled=True)

    gmsh.model.occ.addCircle(-0.5 + translateX, 0.0 + translateY, 0.5, 0.40, tag=110)
    gmsh.model.occ.addDisk(-0.5 + translateX, 0.0 + translateY, 0.5, 0.40, 0.40, tag=510)
    gmsh.model.occ.addCurveLoop([110], 110)
    gmsh.model.occ.addCircle(-1.1 + translateX, 0.0 + translateY, 1.0, 0.8, 210)
    gmsh.model.occ.addCurveLoop([210], 210)
    gmsh.model.occ.addThruSections([110, 210], 1010, makeSolid=False, makeRuled=True)

    # gmsh.model.occ.addBox(0.15, 0.15, 0.40, 4.2, 2.2, 0.2)
    # gmsh.model.occ.remove([(3, 1)])

    gmsh.model.occ.cut([(2, 2)], [(2, 500)])
    gmsh.model.occ.cut([(2, 2)], [(2, 510)])

    # gmsh.model.occ.addCylinder(1.15, 1.1, 0.6,
    #                            0.0, 0.0, 0.4,
    #                            0.75)
    # gmsh.model.occ.remove([(3, 1), (2, 10)])
    # gmsh.model.occ.cut([(2, 8)], [(2, 11)])



    # gmsh.model.occ.addCylinder(3.35, 1.4, 0.6,
    #                            0.0, 0.0, 0.4,
    #                            0.75)
    # gmsh.model.occ.remove([(3, 1), (2, 11)])
    # gmsh.model.occ.cut([(2, 8)], [(2, 12)])


    gmsh.model.occ.synchronize()


    mesh_size = 1 / args.mesh_density
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)

    # bulk_objects = [1]
    surf_objects = [1, 2, 1000, 1010]
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
    print("MSHPY |     Mesh Density =", args.mesh_density)
    print("MSHPY |     Output Directory:", args.output)
    create_mesh(args)
    print("MSHPY | Finished Mesh Creation")
