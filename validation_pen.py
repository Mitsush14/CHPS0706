from maillages import *
from math import *

def aireTriangle(triangle):
    x1, y1, x2, y2, x3, y3 = triangle[0][0], triangle[0][1],triangle[1][0], triangle[1][1],triangle[2][0], triangle[2][1]
    return 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

def fct_u(x,y):
    return 1 + sin((pi/2)* x) + x * (x-4) * cos((pi/2) * y)

def fct_uE():
    return 1

def fct_f(x,y):
    return ((pi**2)/4) * sin((pi/2) * x) + (((pi**2/4)*(x**2) - (pi**2)*x) - 2) * cos((pi/2)*y)

def fct_kappa():
    return 1

def fct_alpha():
    return 10**8

def coeffelem_P1_rigid(triangle):
    kappa = fct_kappa()
    k = np.array([[0.0 for x in range(3)]for y in range(3)])
    coeff = kappa/(4 * aireTriangle(triangle))
    x1,x2,x3,y1,y2,y3 = triangle[0][0], triangle[1][0], triangle[2][0], triangle[0][1], triangle[1][1], triangle[2][1]
    k[0][0] = coeff * (((x2-x3)**2) + (y2-y3)**2)
    k[1][1] = coeff * (((x3-x1)**2) + (y3-y1)**2)
    k[2][2] = coeff * (((x1-x2)**2) + (y1-y2)**2)
    k[0][1] = k[1][0] = coeff *(-(x1-x3)*(x2-x3) - (y1-y3)*(y2-y3))
    k[0][2] = k[2][0] = coeff *(-(x3-x2)*(x1-x2) - (y3-y2)*(y1-y2))
    k[1][2] = k[2][1] = coeff *(-(x2-x1)*(x3-x1) - (y2-y1)*(y3-y1))
    return k

def coeffelem_P1_source(triangle):
    coeff = aireTriangle(triangle)/3
    x,y = (triangle[0][0] + triangle[1][0] + triangle[2][0])/3,(triangle[0][1] + triangle[1][1] + triangle[2][1])/3
    return coeff * fct_f(x,y) * np.array([1.0 for x in range(3)])

def coeffelem_P1_masse(triangle):
    coeff = aireTriangle(triangle)/3
    return np.array([[1,0,0],[0,1,0], [0,0,1]]) *coeff

def coeffelem_P1_transf(arete):
    pointA, pointB = arete[0], arete[1]
    coeff = dist(pointA, pointB)/2
    return coeff * fct_alpha() * fct_uE()  * np.array([1.0 for x in range(2)])

def coeffelem_P1_poids(arete):
    pointA, pointB = arete[0], arete[1]
    alpha = fct_alpha()
    coeff = dist(pointA, pointB)/6
    return coeff * alpha * np.array([[2.0,1.0],[1.0,2.0]])

def assemblage_EF_P1(N, nbe,tri, tri_coord, gammaF, ar, coord):
    #1ère étape : initialisation
    A = np.array([[0. for x in range(N)] for y in range(N)])
    F =  np.array([0. for x in range(N)])

    #2ème étape : addition des termes volumiques
    for l in range(0,nbe):
        kl = coeffelem_P1_rigid(tri_coord[l])
        fl = coeffelem_P1_source(tri_coord[l])
        I1 = tri[l][0] -1
        I2 = tri[l][1] -1
        I3 = tri[l][2] -1

        A[I1][I1] += kl[0][0]
        A[I1][I2] += kl[0][1]
        A[I1][I3] += kl[0][2]
        F[I1] += fl[0]

        A[I2][I1] += kl[1][0]
        A[I2][I2] += kl[1][1]
        A[I2][I3] += kl[1][2]
        F[I2] += fl[1]

        A[I3][I1] += kl[2][0]
        A[I3][I2] += kl[2][1]
        A[I3][I3] += kl[2][2]
        F[I3] += fl[2]

    #Calcul de K ?
    K = np.copy(A)

    #3ème étape : Addition des termes de bords (Fourier/Robin)
    for a in gammaF:
        pa = coeffelem_P1_poids([coord[a[0]-1],coord[a[1]-1]])
        ea = coeffelem_P1_transf([coord[a[0]-1],coord[a[1]-1]])
        I1 = a[0] -1
        I2 = a[1] -1

        A[I1][I1] += pa[0][0]
        A[I1][I2] += pa[0][1]
        F[I1]+=ea[0]

        A[I2][I1] += pa[1][0]
        A[I2][I2] += pa[1][1]
        F[I2]+=ea[1]

    #calcul de Uh
    Uh = np.linalg.solve(A,F)
    U = np.array([0. for x in range(len(coord))])
    for i in range(len(coord)):
        U[i] = fct_u(coord[i][0], coord[i][1])

    return A, F, Uh, K, U

def fct_remplir_gammaF(aretes, refa, codeDir):
    res = []
    for i in range(len(aretes)):
        if(refa[i] == codeDir):
            res.append(aretes[i])
    return res

def validation(path):
    print("================= validation_pas_a_pas =================")
    print(" * Resultats de calculs elementaires ...")
    print("\n * element triangle: xl = [0,1,0] (abscisses), yl = [0,0,1] (ordonnees)")
    kl = coeffelem_P1_rigid([[0,0],[1,0],[0,1]])
    print("kl = ")
    print(kl)
    ml = coeffelem_P1_masse([[0,0],[1,0],[0,1]])
    print("ml = ")
    print(ml)
    fl = coeffelem_P1_source([[0,0],[1,0],[0,1]])
    print("fl = ")
    print(fl)

    print("\n * element arete: xa = [0,0] (abscisses), ya = [0,1] (ordonnees)")
    pa = coeffelem_P1_poids([[0,0], [0,1]])
    print("pa = ")
    print(pa)
    ea = coeffelem_P1_transf([[0,0], [0,1]])
    print("ea =")
    print(ea)
    maillage = Maillage(path)
    print("================= validation_pas_a_pas =================")
    print(" * Resultats sur le mini_maillage "+path+" ...")
    print("nbn = " + str(maillage.nbn))
    print("nbe = " + str(maillage.nbe))
    print("nba = " + str(maillage.nba))

    #ça fonctionne jusque là :p
    gammaF = fct_remplir_gammaF(maillage.ar, maillage.refa, 1)
    A,F, Uh, K, U = assemblage_EF_P1(maillage.nbn, maillage.nbe, maillage.tri, Maillage.triangle_to_coordonnees(maillage.tri, maillage.coord), gammaF, maillage.ar, maillage.coord)
    print("A = ")
    print(A)
    print("F = ")
    print(F)
    print("Uh = ")
    print(Uh)
    erreur1 = sqrt(np.matmul(np.matmul(np.transpose(U-Uh),K),(U-Uh)))
    erreur2 = max(Uh-U)
    print(" ___---===***   RESULTATS:   ***===---___")
    print(" ________________________________________")
    print(" {            min(Uh) : "+str(round(min(Uh),2)))#minUh
    print(" {            max(Uh) :  "+str(round(max(Uh),2)))#maxUh
    print(" {           mean(Uh) :  "+str(round(Uh.mean(),3)))#moyenneUh
    h,Q = maillage.pas_et_qualite_maillage(maillage.nbn, maillage.nbe, maillage.nba, maillage.coord, maillage.tri)
    print(" {            h       : "+str(round(h,3)))#h
    print(" {            Q       : "+str(round(Q,3)))#Q
    print(" {erreur |uh-rh(u)|_H1: "+str(erreur1))#1ère erreur
    print(" {erreur |Uh-U|_inf   : "+str(erreur2))#2ème erreur
    print(" ________________________________________")
    print("================= validation_pas_a_pas =================")
    return erreur1


if __name__ == '__main__':
    path = "Maillages/m"
    typeFile = ".msh"
    erreurs=[]
    for i in range(1, 5):
        erreurs.append(validation(path+str(i)+typeFile))
    for i in range(0,3):
        print("ln(e"+str(i+1)+"/e"+str(i+2)+")/ ln(2) = " + str(log(erreurs[i]/erreurs[i+1])/log(2)))
    #validation(path+"00"+typeFile)





