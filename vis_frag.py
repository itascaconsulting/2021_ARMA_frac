import numpy as np
import scipy.spatial
import networkx as nx
import itasca as it
from itasca import ballarray as ba
from collections import defaultdict
import itertools

def draw_fragment(points,faces,i):
    gname = "tmp{}.geom".format(i)
    f = open(gname, "w")
    print("ITASCA GEOMETRY3D", file=f)
    print("NODES", file=f)
    for j, p in enumerate(points):
        x,y,z = p
        print(j+1, x, y, z, "EXTRA 1 0", file=f)

    print("POLYS", file=f)
    for j, t in enumerate(faces):
        n0, n1, n2 = t
        print(j+1, "NODE", n0+1, n1+1, n2+1, "EXTRA 1 0", file=f)

    f.close()
    it.command("geom import {} format geometry".format(gname))



def mesh_fragments(min_size=10, tol=4.0, explode=0, explode_xyz=False):
    current_bond_list = set(((c.end1().id(), c.end2().id()) for c in \
                             it.contact.list("mechanical") if
                             (c.model()=='linearpbond' and
                              c.prop("pb_state")==3)))
    G=nx.Graph()
    G.add_edges_from(current_bond_list)
    clusters = nx.algorithms.connected_components(G)
    big_fragments = [l for l in clusters if len(l)>min_size]
    max_rad = ba.radius().max()
    for fid, fragment in enumerate(big_fragments):  # ball ids.
        keep = []
        points = np.array([it.ball.find(p).pos() for p in fragment])
        tets = scipy.spatial.Delaunay(points, qhull_options="Qbb Qc Qz Qx")
        edges = ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
        indices = ((0,1,2), (1,2,3), (0,2,3), (0,1,3))

        for i, nl in enumerate(tets.simplices): # for each tet
            edge_lengths = np.array([np.linalg.norm(tets.points[nl[e0]]-tets.points[nl[e1]])
                                     for e0, e1 in edges])
            if (edge_lengths < 4*max_rad).all():
                keep.append(nl)

        face_hash = defaultdict(int)
        for tet in keep: # indices into point list
            for n0,n1,n2 in indices:
                face = [tet[n0], tet[n1],tet[n2]]
                face.sort()
                face_hash[tuple(face)] += 1
        outside_faces = [k for k,v in face_hash.items() if v == 1]
        uniques = list(set(itertools.chain(*outside_faces)))
        pmap = {v : i for i,v in enumerate(uniques)}
        face_points = points[uniques]
        centroid = np.sum(face_points, axis=0)/len(face_points)
        if not explode_xyz:
            centroid[2] = 0
        offset = centroid/np.linalg.norm(centroid) * explode
        face_points += offset
        draw_fragment(face_points, [[pmap[t[j]] for j in range(3)]
                                    for t in outside_faces],fid)


if __name__ == '__main__':
    #sf = "ablast{}.p3sav".format(13)
    #base3d_5.p3sav
    #sf = "basefine103_0.p3sav"
    #it.command("res {}".format(sf))
    it.command("geometry delete")
    mesh_fragments(explode=0.0, explode_xyz=True)
