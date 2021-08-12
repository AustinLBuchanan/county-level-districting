import gurobipy as gp
from gurobipy import GRB 
import math
import networkx as nx
import csv
from gerrychain import Graph
from geopy.distance import geodesic

state_codes = {
    'WA': '53', 'DE': '10', 'WI': '55', 'WV': '54', 'HI': '15',
    'FL': '12', 'WY': '56', 'NJ': '34', 'NM': '35', 'TX': '48',
    'LA': '22', 'NC': '37', 'ND': '38', 'NE': '31', 'TN': '47', 'NY': '36',
    'PA': '42', 'AK': '02', 'NV': '32', 'NH': '33', 'VA': '51', 'CO': '08',
    'CA': '06', 'AL': '01', 'AR': '05', 'VT': '50', 'IL': '17', 'GA': '13',
    'IN': '18', 'IA': '19', 'MA': '25', 'AZ': '04', 'ID': '16', 'CT': '09',
    'ME': '23', 'MD': '24', 'OK': '40', 'OH': '39', 'UT': '49', 'MO': '29',
    'MN': '27', 'MI': '26', 'RI': '44', 'KS': '20', 'MT': '30', 'MS': '28',
    'SC': '45', 'KY': '21', 'OR': '41', 'SD': '46'
}

congressional_districts_2010 = {
    'WA': 10, 'DE': 1, 'WI': 8, 'WV': 3, 'HI': 2,
    'FL': 27, 'WY': 1, 'NJ': 12, 'NM': 3, 'TX': 36,
    'LA': 6, 'NC': 13, 'ND': 1, 'NE': 3, 'TN': 9, 'NY': 27,
    'PA': 18, 'AK': 1, 'NV': 4, 'NH': 2, 'VA': 11, 'CO': 7,
    'CA': 53, 'AL': 7, 'AR': 4, 'VT': 1, 'IL': 18, 'GA': 14,
    'IN': 9, 'IA': 4, 'MA': 9, 'AZ': 9, 'ID': 2, 'CT': 5,
    'ME': 2, 'MD': 8, 'OK': 5, 'OH': 16, 'UT': 4, 'MO': 8,
    'MN': 8, 'MI': 14, 'RI': 2, 'KS': 4, 'MT': 1, 'MS': 4,
    'SC': 7, 'KY': 6, 'OR': 5, 'SD': 1
}


def build_base_model(m, X, DG, L, U, k):
    # create distance dictionary
    dist = dict()
    for i in DG.nodes:
        for j in DG.nodes:
            loc_i = ( DG.nodes[i]['C_Y'],  DG.nodes[i]['C_X'] )
            loc_j = ( DG.nodes[j]['C_Y'],  DG.nodes[j]['C_X'] )
            dist[i,j] = geodesic(loc_i,loc_j).miles
    
    # Set objective function, based on the "moment-of-inertia" and hop-based distances
    m.setObjective(gp.quicksum(dist[i,j]*dist[i,j]*DG.nodes[i]['TOTPOP']*X[i,j] for i in DG.nodes for j in DG.nodes), GRB.MINIMIZE) 
    
    # Each vertex i assigned to one district
    m.addConstrs(sum(X[i,j] for j in DG.nodes) == 1 for i in DG.nodes)
     
    # Pick k centers
    m.addConstr(sum(X[j,j] for j in DG.nodes) == k)
    
    # Population balance: population assigned to vertex j should be in [L,U]
    m.addConstrs(gp.quicksum(DG.nodes[i]['TOTPOP'] * X[i,j] for i in DG.nodes) <= U * X[j,j] for j in DG.nodes)
    m.addConstrs(gp.quicksum(DG.nodes[i]['TOTPOP'] * X[i,j] for i in DG.nodes) >= L * X[j,j] for j in DG.nodes)
    
    # Add coupling inequalities for added model strength
    m.addConstrs(X[i,j] <= X[j,j] for i in DG.nodes for j in DG.nodes)
    
    # Set branch priority on center vars
    for j in DG.nodes:
        X[j,j].BranchPriority=1
        
    # Choose each district's "center" as its most populous county, by fixing x_ij=0 when population[i]>population[j]
    for i in DG.nodes:
        for j in DG.nodes:
            if DG.nodes[i]['TOTPOP']>DG.nodes[j]['TOTPOP']:
                X[i,j].UB=0
    

def add_contiguity_constraints(m,X,DG,U):
    # F[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    F = m.addVars(DG.nodes,DG.edges,vtype=GRB.CONTINUOUS)
    M = DG.number_of_nodes()-1
    for j in DG.nodes:
        m.addConstr( sum(F[j,u,j] for u in DG.neighbors(j)) == 0 )
        for i in DG.nodes:
            if i!=j:
                m.addConstr( sum(F[j,u,i] for u in DG.neighbors(i)) - sum(F[j,i,u] for u in DG.neighbors(i)) == X[i,j] )
                m.addConstr( sum(F[j,u,i] for u in DG.neighbors(i)) <= M*X[i,j] )
                

if __name__ == '__main__':
    
    
    output_file_name = "outputs-2010.csv"
    
    with open(output_file_name,'w',newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        fields = ['State','n','m','k','L','U','Contiguity?','Feasible?']
        csvwriter.writerow(fields)
        
        for state in state_codes.keys():
            for impose_contiguity in {True,False}:
            
                ###########################
                # INPUTS 
                ###########################
                
                k = congressional_districts_2010[state]
                
                # Read county-level graph from "<state>_county.json"
                filepath = "C:\\districting-data-2010\\"
                filename = state + '_county.json'
                
                G = Graph.from_json( filepath + filename )
                
                DG = nx.DiGraph(G) # bidirected version of G
                
                deviation = 0.01 # +/- 0.5%
                
                total_population = sum( G.nodes[i]['TOTPOP'] for i in G.nodes )
                
                L = math.ceil((1-deviation/2)*total_population/k)
                U = math.floor((1+deviation/2)*total_population/k)
                
                print("L =",L,"U =",U,"k =",k)
                
                row = [state,G.number_of_nodes(),G.number_of_edges(),k,L,U,impose_contiguity]
                
                ###########################
                # BUILD MIP 
                ###########################
                
                m = gp.Model()
                
                # X[i,j]=1 if vertex i is assigned to (district centered at) vertex j
                X = m.addVars(DG.nodes, DG.nodes, vtype=GRB.BINARY) 
                
                # add compactness objective; add constraints (i.e., assignment, population balance, coupling)
                build_base_model(m, X, DG, L, U, k)
                if impose_contiguity:
                    add_contiguity_constraints(m, X, DG, U)
                
                #############################
                # SOLVE MIP AND RETRIEVE SOLUTION
                #############################
                
                m.optimize()
                if m.status == GRB.OPTIMAL:
                    centers = [j for j in DG.nodes if X[j,j].x > 0.5]
                    districts = [[i for i in DG.nodes if X[i,j].x > 0.5] for j in centers]
                    print("districts=",districts)
                    row.append("True")
                elif m.status == GRB.INFEASIBLE:
                    print("No feasible solution!")
                    row.append("False")
                else:
                    print("Gurobi terminated with status = ",m.status)
                    row.append("?")
                csvwriter.writerow(row)