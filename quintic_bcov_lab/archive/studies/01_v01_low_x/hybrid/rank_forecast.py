"""Count potentially nonvanishing primitive BPS degrees among MUM slots.

Uses only the established quintic Castelnuovo bound g <= (d^2+5d+10)/10.
This is a count of *missing geometric inputs*, not automatically a rigorous
rank defect: rank also depends on the full boundary matrix and modular data.
"""

def mum_degree(g):return (2*g-2)//5

def unknown_degrees(g,known_gv=()):
    known_gv=set(known_gv)
    return [d for d in range(1,mum_degree(g)+1)
            if 10*g<=d*d+5*d+10 and (g,d) not in known_gv]

def forecast(genera,known_gv=((51,20),)):
    return [{'genus':g,'mum_degree':mum_degree(g),
             'potentially_unknown_GV_degrees':unknown_degrees(g,known_gv),
             'count':len(unknown_degrees(g,known_gv))} for g in genera]
