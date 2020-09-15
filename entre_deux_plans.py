# Considérons deux plans parallèles et non confudus, P1 et P2, d'équation :

# (P1) : (A)*x + (B)*y + (C)*z + D1
# (P2) : (A)*x + (B)*y + (C)*z + D2

# Un point, de coordonnée [Xa,Ya,Za] appartient à la portion d'espace situé entre les deux plans ssi :

#    Cas n°1 : D1 > D2   et   D2 < A*Xa + B*Ya + C*Za < D1
#
#    Cas n°2 : D2 > D1   et   D1 < A*Xa + B*Ya + C*Za < D2

# Cette fonction prend en argument les quatre constantes des plans concernés ainsi que les coordonnées du point dont on veut savoir s'il appartient à l'espace situé entre les deux plans.

def entre_deux_plans(plan_1, plan_2, point):
    
    resolution = float( (plan_1[0]*point[0]) + (plan_1[1]*point[1]) + (plan_1[2]*point[2]) )

    if plan_1[3] < plan_2[3] and resolution >= plan_1[3] and resolution <= plan_2[3] :
        return True

    elif plan_1[3] > plan_2[3] and resolution <= plan_1[3] and resolution >= plan_2[3] :
        return True

    else :
        return False

# Finalement, cette fonction rapporte une valeur True ou False.
