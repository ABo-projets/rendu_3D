import math
import numpy as np
import matplotlib.pyplot as plt
import random
import time
import csv

class Noeud :
    """Classe correspondant a un arbre, pour implementer l'arbre k-d"""
    def __init__(self,feuille,u,plan,obj,gauche,droit):
        self.feuille = feuille
        self.axe = u
        self.val = plan
        self.obj = obj
        self.gauche = gauche
        self.droit = droit
    
def extrema(coord,maxi,liste_obj):
    """Trouve les coordonnees extremales d'une liste d'objets"""
    tmp = -maxi*math.inf
    for obj in liste_obj :
        if type(obj) == Sphere :
            if maxi*obj.centre[coord] + maxi*obj.rayon > maxi*tmp :
                tmp = obj.centre[coord] + obj.rayon
        if type(obj) == Triangle :
            if maxi*obj.p1[coord] > maxi*tmp :
                tmp = obj.p1[coord]
            if maxi*obj.p2[coord] > maxi*tmp :
                tmp = obj.p2[coord]
            if maxi*obj.p3[coord] > maxi*tmp :
                tmp = obj.p3[coord]
    return tmp 
             
def construit_arbre(liste_obj,p):
    if len(liste_obj) == 0 :
        return None
    if len(liste_obj) == 1 :
        return Noeud(True,None,None,liste_obj[0],None,None)
    coord = p%3
    moy = 0
    for obj in liste_obj :
        moy += obj.centre[coord]
    moy = moy/(len(liste_obj))
    lst_gauche = []
    lst_droite = []
    for obj in liste_obj :
        if obj.centre[coord] < moy :
            lst_gauche.append(obj)
        else :
            lst_droite.append(obj)
    return Noeud(False,coord,moy,None,construit_arbre(lst_gauche,p+1),construit_arbre(lst_droite,p+1))

class Vecteur :
    """Description d'un vecteur avec ses coordonnees x,y et z et les methodes classiques associees
    Un point sera considere comme un vecteur depuis l'origine"""
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    def __add__(self,p2):
        return Vecteur(self.x+p2.x,self.y+p2.y,self.z+p2.z)
    
    def __sub__(self,p2):
        return Vecteur(self.x-p2.x,self.y-p2.y,self.z-p2.z)
    
    def __mul__(self,k):
        return Vecteur(k*self.x,k*self.y,k*self.z)
    
    def __div__(self,k):
        return Vecteur(self.x/k,self.y/k,self.z/k)
    
    def __str__(self):
        c = "x: {}, y: {}, z: {}".format(self.x,self.y,self.z)
        return c
    
    def __getitem__(self,coord):
        if coord == 0 :
            return self.x
        if coord == 1 :
            return self.y
        else :
            return self.z
    
    def norme(self):
        return math.sqrt(self.x**2+self.y**2+self.z**2)
    
    def scalaire(self,p2):
        return self.x*p2.x+self.y*p2.y+self.z*p2.z
    
    def normalise(self):
        n = self.norme()
        if n != 0 : 
            return self*(1/n)
        else :
            return self
    
    def reflechi(self,axe):
        return (self - axe*2*(self.scalaire(axe))).normalise()
    
    def vectoriel(self,v):
        x = self.y*v.z-self.z*v.y
        y = self.z*v.x-self.x*v.z
        z = self.x*v.y-self.y*v.x
        return Vecteur(x,y,z)
    
class Ray :
    def __init__(self,o,v):
        self.origine = o
        self.direction = v.normalise()
        
class Sphere :
    def __init__(self,c,r,ca,cd,cs,s,ref):
        self.centre = c
        self.rayon = r
        self.coul_amb = ca
        self.coul_diff = cd
        self.coul_spec = cs
        self.shin = s
        self.reflection = ref
    
    def intersection(self,ray,scene):
        b = 2*(ray.direction.scalaire(ray.origine-self.centre))
        c = (ray.origine-self.centre).norme()**2 - self.rayon**2
        delta = b*b - 4*c
        if delta > 0 :
            r1 = (-b-math.sqrt(delta))/2
            r2 = (-b+math.sqrt(delta))/2
            if r1 > 0 and r2 > 0 :
                return min(r1,r2)
        return None
    
    def normale(self,point):
        return (point-self.centre).normalise()

class Triangle :
    def __init__(self,p1,p2,p3,ca,cd,cs,s,ref):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.coul_amb = ca
        self.coul_diff = cd
        self.coul_spec = cs
        self.shin = s
        self.reflection = ref
        self.n = self.normale(None)
        self.p1p2 = (self.p2-self.p1).normalise()
        self.p2p3 = (self.p3-self.p2).normalise()
        self.p1p3 = (self.p1-self.p3).normalise()
        self.centre = (p1+p2+p3)*(1/3)
        self.ang1 = self.p1p2.scalaire((self.p3-self.p1).normalise())
        self.ang2 = (self.p1-self.p2).normalise().scalaire((self.p3-self.p2).normalise())
        self.ang3 = (self.p1-self.p3).normalise().scalaire((self.p2-self.p3).normalise())
       
    def intersection(self,ray,scene):
        if (self.n.scalaire(scene.normale) < 0):
            self.n*=-1.0
        if ray.direction.scalaire(self.n) == 0 :
            return None

        k = self.n.x*self.p1.x + self.n.y*self.p1.y + self.n.z*self.p1.z
        somc = k - self.n.x*ray.origine.x - self.n.y*ray.origine.y - self.n.z*ray.origine.z
        somt = self.n.x*ray.direction.x + self.n.y*ray.direction.y + self.n.z*ray.direction.z
        t = somc/somt
        if t<= 0 :
            return None
        p = ray.origine+ray.direction*t
        if self.appartient(p) :
            return t
        return None
   
    def appartient(self,p):
        angp1 = self.p1p2.scalaire((p-self.p1).normalise())
        angp2 = self.p2p3.scalaire((p-self.p2).normalise())
        angp3 = self.p1p3.scalaire((p-self.p3).normalise())
        return angp1>=self.ang1 and angp2>=self.ang2 and angp3>=self.ang3
   
    def normale(self,_):
        return ((self.p3-self.p2).vectoriel(self.p2-self.p1)).normalise()

class Scene :
    """Definition d'une scene comme on l'a definit dnas l'introduction"""
    def __init__(self):
        self.recul = 2
        self.lumiere = Vecteur(-3,7,1)
        self.amb = np.array([1,1,1])
        self.diff = np.array([1,1,1])
        self.spec = np.array([1,1,1])
        self.p1 = Vecteur(-3,4.7,2)
        self.p2 = Vecteur(3,4.7,2)
        self.p3 = Vecteur(3,0.7,2)
        self.p4 = Vecteur(-3,0.7,2)
        centre = self.p1+(self.p2-self.p1)*(1/2)+(self.p3-self.p2)*(1/2)
        self.normale = ((self.p1-self.p2).vectoriel(self.p2-self.p3)).normalise()
        self.camera = centre-self.normale*self.recul
        self.ratio = (self.p3-self.p2).norme()/(self.p2-self.p1).norme()
        self.width = 600
        self.height = int(self.width*self.ratio)
        self.objet = []
        
    def trouve_position(self,i,j):
        pi = (self.p2-self.p1)*(j/self.width)
        pj = (self.p3-self.p2)*(i/self.height)
        return self.p1+pi+pj
    
    def change_largeur(self,nb_pix):
        self.width = nb_pix
        self.height = int(self.width*self.ratio)

def obj_proche(ray,scene):
    """Trouve l'objet le plus proche naivement"""
    obj_proche = None
    t = math.inf
    for obj in scene.objet :
        t_int = obj.intersection(ray,scene)
        if t_int != None and t_int<t:
            t = t_int
            obj_proche = obj
    return obj_proche,t

def touche(ray,lst_conditions):
    """Determine si un rayon touche le pave droit forme par lst_conditions"""
    xmi = lst_conditions[0]
    ymi = lst_conditions[1]
    zmi = lst_conditions[2]
    xma = lst_conditions[3]
    yma = lst_conditions[4]
    zma = lst_conditions[5]
    point_min = ray.origine + ray.direction*((zmi-ray.origine.z)/ray.direction.z)
    point_max = ray.origine + ray.direction*((zma-ray.origine.z)/ray.direction.z)
    condition_y1 = (point_min.y>yma and point_max.y>yma)
    condition_y2 = (point_max.y<ymi and point_min.y<ymi)
    condition_x1 = (point_min.x>xma and point_max.x>xma)
    condition_x2 = (point_max.x<xmi and point_min.x<xmi)
    return (condition_y1 or condition_y2 or condition_x1 or condition_x2)

def intersection(ray,scene,arbre,lst_conditions):
    """intersection avec heuristique sur l'ordre de traitement des objets"""
    if arbre == None :
        return None, math.inf
    if arbre.feuille :
        t = arbre.obj.intersection(ray,scene)
        if t != None :
            return arbre.obj,t
        else :
            return None,math.inf
    if touche(ray,lst_conditions) :
        if ray.origine[arbre.axe] <= arbre.val :
            lst_new = lst_conditions.copy()
            lst_new[3+arbre.axe] = arbre.val
            obj,dist = intersection(ray,scene,arbre.gauche,lst_new)
            if obj != None :
                return obj,dist
            else :
                lst_new = lst_conditions.copy()
                lst_new[arbre.axe] = arbre.val
                obj,dist = intersection(ray,scene,arbre.droit,lst_new)
                return obj,dist
        else :
            lst_new = lst_conditions.copy()
            lst_new[arbre.axe] = arbre.val
            obj,dist = intersection(ray,scene,arbre.droit,lst_new)
            if obj != None :
                return obj,dist
            else :
                lst_new = lst_conditions.copy()
                lst_new[3+arbre.axe] = arbre.val
                obj,dist = intersection(ray,scene,arbre.gauche,lst_new)
                return obj,dist
    else :
        return None, math.inf

def illumination(scene,obj,L,N,V):
    amb = obj.coul_amb*scene.amb
    diff = obj.coul_diff*scene.diff*L.scalaire(N)
    H = (L + V).normalise()
    spec = obj.coul_spec*scene.spec*(N.scalaire(H))**(obj.shin/4)
    return amb+diff+spec

def raytrace(i,j,scene,indice_reflexion_max) :
    o = scene.trouve_position(i,j)
    d = o-scene.camera
    ray = Ray(o,d)
    couleur = np.zeros((3))
    reflection = 1
    while reflection>=indice_reflexion_max :
        n_obj, t = obj_proche(ray,scene)
        if n_obj != None :
            p_int = ray.origine+ray.direction*t
            N = n_obj.normale(p_int)
            _, dist_av_lum = obj_proche(Ray(p_int+N*1e-10,scene.lumiere-p_int),scene)
            if dist_av_lum>=(scene.lumiere-p_int).norme() :
                L = (scene.lumiere-p_int).normalise()
                V = (scene.camera-p_int).normalise()
                ill = illumination(scene,n_obj,L,N,V)
                couleur += (1-n_obj.reflection)*reflection*ill
                reflection *= n_obj.reflection
                o = p_int+N*1e-3
                d = ray.direction.reflechi(N)
                ray = Ray(o,d)
            else :
                break
        else :
            break
    return np.clip(couleur,0,1)

def raytracing(nb_pixel_largeur,nb_scene,scene,indice_reflexion_max):
    scene.change_largeur(nb_pixel_largeur)
    screen = np.zeros((scene.height,scene.width,3))
    for i in range(scene.height):
        for j in range(scene.width):
            screen[i,j] = raytrace(i,j,scene,indice_reflexion_max)
    nom_image = "Test_Sans_amelioration_"+str(nb_scene)+"_"+".png"
    plt.imsave(nom_image,screen)

def new_raytrace(i,j,scene,arbre,conditions,indice_reflexion_max) :
    """Avec amelioration"""
    o = scene.trouve_position(i,j)
    d = o-scene.camera
    ray = Ray(o,d)
    couleur = np.zeros((3))
    reflection = 1
    while reflection>=indice_reflexion_max :
        n_obj, t = intersection(ray,scene,arbre,conditions)
        if n_obj != None :
            p_int = ray.origine+ray.direction*t
            N = n_obj.normale(p_int)
            _, dist_av_lum = obj_proche(Ray(p_int+N*1e-10,scene.lumiere-p_int),scene)
            if dist_av_lum>=(scene.lumiere-p_int).norme() :
                L = (scene.lumiere-p_int).normalise()
                V = (scene.camera-p_int).normalise()
                ill = illumination(scene,n_obj,L,N,V)
                couleur += (1-n_obj.reflection)*reflection*ill
                reflection *= n_obj.reflection
                o = p_int+N*1e-3
                d = ray.direction.reflechi(N)
                ray = Ray(o,d)
            else :
                break
        else :
            break
    return np.clip(couleur,0,1)

def new_raytracing(nb_pixel_largeur,nb_scene,scene,arbre,conditions,indice_reflexion_max):
    """Avec Amelioration"""
    scene.change_largeur(nb_pixel_largeur)
    screen = np.zeros((scene.height,scene.width,3))
    for i in range(scene.height):
        for j in range(scene.width):
            screen[i,j] = new_raytrace(i,j,scene,arbre,conditions,indice_reflexion_max)
    nom_image = "Test_Avec_amelioration_"+str(nb_scene)+"_"+".png"
    plt.imsave(nom_image,screen)
    
def projection(scene,p):
    """fonction de projection pour la rasterisatoin"""
    centre_ecran = scene.p1+(scene.p2-scene.p1)*(1/2)+(scene.p4-scene.p1)*(1/2)
    x_axis = (scene.p2-scene.p1).normalise()
    y_axis = (scene.p2-scene.p3).normalise()
    z_axis = x_axis.vectoriel(y_axis)
    local_to_world = np.array([[x_axis.x,x_axis.y,x_axis.z,0],
                               [y_axis.x,y_axis.y,y_axis.z,0],
                               [z_axis.x,z_axis.y,z_axis.z,0],
                               [centre_ecran.x,centre_ecran.y,centre_ecran.z,1]])
    world_to_camera = np.linalg.inv(local_to_world)
    p_loc = np.array([p.x,p.y,p.z,1]) @ world_to_camera
    p_ecran_x = p_loc[0]/(-p_loc[2])*scene.recul
    p_ecran_y = p_loc[1]/(-p_loc[2])*scene.recul
    return Vecteur(p_ecran_x,p_ecran_y,-p_loc[2])

def extraire_pixel(v,scene):
    width = (scene.p2-scene.p1).norme()
    height = (scene.p3-scene.p2).norme() 
    p_norm0 = (v.x+width/2)/width
    p_norm1 = (v.y+height/2)/height
    pix0 = math.floor(p_norm0*scene.width)
    pix1 = math.floor((1-p_norm1)*scene.height)
    return[pix1,pix0,v.z]

def fonctioncote(v1,v2,p):
    return (p[0]-v1[0])*(v2[1]-v1[1])-(p[1]-v1[1])*(v2[0]-v1[0])
    
def PixelIn(v1,v2,v3,i,j):
    b1 = fonctioncote(v1,v2,[i,j]) >= 0
    b2 = fonctioncote(v2,v3,[i,j]) >= 0
    b3 = fonctioncote(v3,v1,[i,j]) >= 0
    return b1 and b2 and b3

def distance(i,j,v0,v1,v2,scene):
    lambda0 = fonctioncote(v1,v2,[i,j])
    lambda1 = fonctioncote(v2,v0,[i,j])
    lambda2 = fonctioncote(v0,v1,[i,j])
    return (lambda0*v0[2] + lambda1*v1[2] + lambda2*v2[2])

def rasterisation(nb_pixel_largeur,nb_scene,scene):
    scene.change_largeur(nb_pixel_largeur)
    screen = np.zeros((scene.height,scene.width,3))
    z_buffer = [[math.inf for i in range(scene.width)] for j in range(scene.height)]
    for obj in scene.objet :
        v1 = projection(scene,obj.p1)
        v2 = projection(scene,obj.p2)
        v3 = projection(scene,obj.p3)
        pix1 = extraire_pixel(v1,scene)
        pix2 = extraire_pixel(v2,scene)
        pix3 = extraire_pixel(v3,scene)
        if pix1[:2] != pix2[:2] and pix1[:2] != pix3[:2] and pix2[:2] != pix3[:2]  :
            for i in range(max(0,min(pix1[0],pix2[0],pix3[0])),min(max(pix1[0],pix2[0],pix3[0]),scene.height)):
                for j in range(max(0,min(pix1[1],pix2[1],pix3[1])),min(max(pix1[1],pix2[1],pix3[1]),scene.width)):
                    if PixelIn(pix1,pix2,pix3,i,j):
                        dist = distance(i,j,pix1,pix2,pix3,scene)
                        if dist<z_buffer[i][j]:
                            screen[i][j] = obj.coul_amb
                            z_buffer[i][j] = dist
    return screen


noir = np.array([0, 0, 0])
blanc = np.array([1, 1, 1])
bleu = np.array([0, 0.4, 0])
marron = np.array([0.49, 0.2, 0])
violet = np.array([0.52,0,0.54])
vert = np.array([0.38,0.8,0])
jaune = np.array([0.9,0.87,0])
couleurs = [blanc,bleu,marron,violet,vert,jaune]


scene10 = Scene()
scene50 = Scene()
scene150 = Scene()

for i in range(10):
    x = random.uniform(-8,8)
    y = random.uniform(-2,7)
    z = random.uniform(-13,-1)
    r = random.uniform(0.1,0.5)
    ref = random.uniform(0,0.4)
    c1 = couleurs[random.randint(0,len(couleurs)-1)]
    c2 = couleurs[random.randint(0,len(couleurs)-1)]
    scene10.objet.append(Sphere(Vecteur(x,y,z),r,c1,c2,blanc,100,ref))
    
for i in range(50):
    x = random.uniform(-8,8)
    y = random.uniform(-2,7)
    z = random.uniform(-13,-1)
    r = random.uniform(0.1,0.5)
    ref = random.uniform(0,0.4)
    c1 = couleurs[random.randint(0,len(couleurs)-1)]
    c2 = couleurs[random.randint(0,len(couleurs)-1)]
    scene50.objet.append(Sphere(Vecteur(x,y,z),r,c1,c2,blanc,100,ref))
    
for i in range(150):
    x = random.uniform(-8,8)
    y = random.uniform(-2,7)
    z = random.uniform(-13,-1)
    r = random.uniform(0.1,0.5)
    ref = random.uniform(0,0.4)
    c1 = couleurs[random.randint(0,len(couleurs)-1)]
    c2 = couleurs[random.randint(0,len(couleurs)-1)]
    scene150.objet.append(Sphere(Vecteur(x,y,z),r,c1,c2,blanc,100,ref))

scenes = [scene10,scene50,scene150]

arbres = []
for scene in scenes :
    arbres.append(construit_arbre(scene.objet,0))

conditions = []
for scene in scenes :
    xmin = extrema(0,-1,scene.objet)
    ymin = extrema(1,-1,scene.objet)
    zmin = extrema(2,-1,scene.objet)
    xmax = extrema(0,1,scene.objet)
    ymax = extrema(1,1,scene.objet)
    zmax = extrema(2,1,scene.objet)
    conditions.append([xmin,ymin,zmin,xmax,ymax,zmax])

with open('TestAmelioration.csv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter=';')
    writer.writerow(["Numero de scene", "Nombre d'objets","Sans amelioration","Avec amelioration"])
    for taille in [200]:
        for i in range(len(scenes)) :
            deb1 = time.process_time()
            raytracing(taille,i,scenes[i],0.2)
            fin1 = time.process_time()
            deb2 = time.process_time()
            new_raytracing(taille,i,scenes[i],arbres[i],conditions[i],0.2)
            fin2 = time.process_time()
            writer.writerow([i,len(scenes[i].objet),fin1-deb1,fin2-deb2])