import gmsh
import math
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Creates a 3D mesh of a cube.")
    parser.add_argument('-m', '--mesh_density', type=int, default=64, help='the \"density\" of the mesh') 
    parser.add_argument('-r', '--radius', type=float, default=0.5, help='the radius of the quarter circles') 
    parser.add_argument('-o', '--output', type=str, default="meshes/3D/", help='the output directory of the mesh file')
    parser.add_argument('-g', '--graphics', action='store_true', help='if --graphics is given, display the GMSH GUI of the built mesh')

    args = parser.parse_args()
    return args

def create_boundary(args):
    print("MSHPY | Initializing GMSH for rigid mesh")
    R = args.radius
    mesh_size = 1.0 / args.mesh_density
    file_name = args.output + "schwarz"
    file_name_end = "_msh" + str(args.mesh_density) + ".msh"

    gmsh.initialize()
    gmsh.model.add("Schwarz_P_Cube")

    p001 = gmsh.model.occ.addPoint(0, 0, 1)
    p010 = gmsh.model.occ.addPoint(0, 1, 0)
    p100 = gmsh.model.occ.addPoint(1, 0, 0)

    p01R = gmsh.model.occ.addPoint(0, 1, R)
    p10R = gmsh.model.occ.addPoint(1, 0, R)
    p0R1 = gmsh.model.occ.addPoint(0, R, 1)
    p1R0 = gmsh.model.occ.addPoint(1, R, 0)
    pR01 = gmsh.model.occ.addPoint(R, 0, 1)
    pR10 = gmsh.model.occ.addPoint(R, 1, 0)

    Lx1 = gmsh.model.occ.addLine(p100, p1R0)
    Lx3 = gmsh.model.occ.addLine(p100, p10R)
    Ly1 = gmsh.model.occ.addLine(p010, p01R)
    Ly3 = gmsh.model.occ.addLine(p010, pR10)
    Lz1 = gmsh.model.occ.addLine(p001, pR01)
    Lz3 = gmsh.model.occ.addLine(p001, p0R1)
                                 
    Lcx = gmsh.model.occ.addCircleArc(p1R0, p100, p10R, center=True)
    Lcy = gmsh.model.occ.addCircleArc(p01R, p010, pR10, center=True)
    Lcz = gmsh.model.occ.addCircleArc(pR01, p001, p0R1, center=True)

    sx_loop = gmsh.model.occ.addCurveLoop([Lx1, Lcx, -Lx3])
    sy_loop = gmsh.model.occ.addCurveLoop([Ly1, Lcy, -Ly3])
    sz_loop = gmsh.model.occ.addCurveLoop([Lz1, Lcz, -Lz3])

    sx = gmsh.model.occ.addPlaneSurface([sx_loop])
    sy = gmsh.model.occ.addPlaneSurface([sy_loop])
    sz = gmsh.model.occ.addPlaneSurface([sz_loop])

    rigid_objects = [sx, sy, sz]

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    


    print("MSHPY | Defining GMSH Physical Objects")
    rigidtag = gmsh.model.addPhysicalGroup(2, rigid_objects, name="Rigid")
    rigid_file_name = file_name + "_rigid" + file_name_end
    
    gmsh.model.mesh.generate(2)    
    
    gmsh.write(rigid_file_name)
    print("MSHPY |     Rigid Physical Tag:", str(rigidtag))
    print("MSHPY |         Written rigid mesh to", rigid_file_name)



    if args.graphics:
        gmsh.fltk.run()
    
    print("MSHPY | Closing GMSH instance")
    gmsh.finalize()

def create_mesh(args):
    print("MSHPY | Initializing GMSH for membrane mesh")
    gmsh.initialize()
    gmsh.model.add("Schwarz_P_Cube")

    R = args.radius
    mesh_size = 1.0 / args.mesh_density
    file_name = args.output + "schwarz"
    file_name_end = "_msh" + str(args.mesh_density) + ".msh"

    print("MSHPY | Creating the plane geometry")
    p001 = gmsh.model.occ.addPoint(0, 0, 1)
    p010 = gmsh.model.occ.addPoint(0, 1, 0)
    p011 = gmsh.model.occ.addPoint(0, 1, 1)
    p100 = gmsh.model.occ.addPoint(1, 0, 0)
    p101 = gmsh.model.occ.addPoint(1, 0, 1)
    p110 = gmsh.model.occ.addPoint(1, 1, 0)
    p111 = gmsh.model.occ.addPoint(1, 1, 1)

    p01R = gmsh.model.occ.addPoint(0, 1, R)
    p10R = gmsh.model.occ.addPoint(1, 0, R)
    p0R1 = gmsh.model.occ.addPoint(0, R, 1)
    p1R0 = gmsh.model.occ.addPoint(1, R, 0)
    pR01 = gmsh.model.occ.addPoint(R, 0, 1)
    pR10 = gmsh.model.occ.addPoint(R, 1, 0)

    Lx2 = gmsh.model.occ.addLine(p1R0, p110)
    Lx4 = gmsh.model.occ.addLine(p10R, p101)

    Ly2 = gmsh.model.occ.addLine(p01R, p011)
    Ly4 = gmsh.model.occ.addLine(pR10, p110)

    Lz2 = gmsh.model.occ.addLine(pR01, p101)
    Lz4 = gmsh.model.occ.addLine(p0R1, p011)
    
    Lvx = gmsh.model.occ.addLine(p011, p111)
    Lvy = gmsh.model.occ.addLine(p101, p111)
    Lvz = gmsh.model.occ.addLine(p110, p111)

    Lcx = gmsh.model.occ.addCircleArc(p1R0, p100, p10R, center=True)
    Lcy = gmsh.model.occ.addCircleArc(p01R, p010, pR10, center=True)
    Lcz = gmsh.model.occ.addCircleArc(pR01, p001, p0R1, center=True)

    scx_loop = gmsh.model.occ.addCurveLoop([Lx2, Lvz, -Lvy, -Lx4, -Lcx])
    scy_loop = gmsh.model.occ.addCurveLoop([Ly2, Lvx, -Lvz, -Ly4, -Lcy])
    scz_loop = gmsh.model.occ.addCurveLoop([Lz2, Lvy, -Lvx, -Lz4, -Lcz])

    scx = gmsh.model.occ.addPlaneSurface([scx_loop])
    scy = gmsh.model.occ.addPlaneSurface([scy_loop])
    scz = gmsh.model.occ.addPlaneSurface([scz_loop])

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    
    print("MSHPY | Defining GMSH Physical Objects")
    surf_objects = [scx, scy, scz]
    surftag = gmsh.model.addPhysicalGroup(2, surf_objects, name="Surface")
    surf_file_name = file_name + "_surf" + file_name_end

    gmsh.model.mesh.generate(2)

    gmsh.write(surf_file_name)
    print("MSHPY |     Surface Physical Tag:", str(surftag))
    print("MSHPY |         Written surface mesh to", surf_file_name)

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
    create_boundary(args)
    print("MSHPY | Finished Mesh Creation")

