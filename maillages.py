import matplotlib.pyplot as plt
import math
import numpy as np
import copy
#Classe Maillage permetttant de générer un maillage à partir d'un fichier et effectuer différentes opérations dessus
class Maillage(object):
    def __init__(self, path):
        (self.nbn,
         self.nbe,
         self.nba,
         self.coord,
         self.tri,
         self.ar,
         self.refn,
         self.reft,
         self.refa) = self.lit_fichier_msh(path)


    def lit_fichier_msh(self, path):
        f = open(path, "r")
        coord = []
        tri = []
        ar = []
        refn = []
        reft = []
        refa = []
        (nbn, nbe, nba) = f.readline().split()
        nbn = int(nbn)
        nbe = int(nbe)
        nba = int(nba)
        for i in range(nbn):
            coord.append(f.readline().split())
            refn.append(int(coord[i][-1]))
            del(coord[i][-1])
            for j in range (2):
                coord[i][j] = float(coord[i][j])
        for i in range(nbe):
            tri.append(f.readline().split())
            reft.append(int(tri[i][-1]))
            del(tri[i][-1])
            for j in range (3):
                tri[i][j] = int(tri[i][j])
        for i in range(nba):
            ar.append(f.readline().split())
            refa.append(int(ar[i][-1]))
            del(ar[i][-1])
            for j in range (2):
                ar[i][j] = int(ar[i][j])
        f.close()

        return (nbn, nbe, nba, coord, tri, ar, refn, reft, refa)

    def trace_maillage_ind(self, nbn, nbe, nba, coord, tri, ar):
        xy = np.asarray(coord)
        x = np.degrees(xy[:, 0])
        y = np.degrees(xy[:, 1])
        triangles = np.asarray(tri)

        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.triplot(x, y, triangles-1, 'go-', lw=1.0)
        plt.show()

    def charge_et_affiche_maillage(fichierMaillage):
        m = Maillage(fichierMaillage)
        m.trace_maillage_ind(m.nbn, m.nbe, m.nba, m.coord, m.tri, m.ar)
        return m;

    def pas_et_qualite_maillage(self, nbn, nbe, nba, coord, tri):
        #calcul distance entre chaque point du triangle et prendre le max
        #pour le pas d'un maillage -> prendre le pas le plus haut sur ts les triangles
        #qualité : prendre le pas du triangle + son rayon isncrit et faire des opérations dessus
        listePasTriangles = []
        listeQualiteTriangles = []
        for triangle in tri:
            points = []
            dists = []
            point1 = triangle[0] -1
            points.append(point1)
            point2 = triangle[1] - 1
            points.append(point2)
            point3 = triangle[2] - 1
            points.append(point3)
            a = Maillage.distance(coord, point1, point2)
            dists.append(a)
            b = Maillage.distance(coord, point2, point3)
            dists.append(b)
            c = Maillage.distance(coord, point3, point1)
            dists.append(c)
            aire = Maillage.calcul_aire(coord, points, dists)
            rayon = 2 * aire / (a + b + c)
            pas = max(dists[0], dists[1], dists[2])
            listePasTriangles.append(pas)
            qual = (math.sqrt(3)/6) * (pas/rayon)
            listeQualiteTriangles.append(qual)
        return max(listePasTriangles), max(listeQualiteTriangles)

    def distance(coord, point1, point2):
        return math.sqrt((coord[point1][0] - coord[point2][0])**2 + (coord[point1][1] - coord[point2][1])**2)

    def calcul_aire(coord, points, dists):
        a, b, c = dists[0], dists[1], dists[2]
        s =(a+b+c)/2
        return math.sqrt(s * (s - a) * (s - b) * (s - c))

    def triangle_to_coordonnees(tri, coord):
        triangle = copy.deepcopy(tri)
        for i in range(0,len(triangle)):
            for j in range(0,len(triangle[i])):
                triangle[i][j] = coord[triangle[i][j]-1]
        return triangle


if __name__ == '__main__':
    path = "Maillages/m"
    typeFile = ".msh"
    #for i in range (1,5):
        #m = Maillage.charge_et_affiche_maillage(path+str(i)+typeFile)
        #pas, qualite = m.pas_et_qualite_maillage(m.nbn, m.nbe, m.nba, m.coord, m.tri)
        #print(f'Pas : {pas :.2f}, Qualité : {qualite:.2f}')
    Maillage.charge_et_affiche_maillage("Maillages/m00.msh")
