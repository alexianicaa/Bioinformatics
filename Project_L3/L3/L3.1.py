import math

def tm_1(seq, Na = 0.0001):
        A = seq.count("A")
        T = seq.count("T")
        G = seq.count("G")
        C = seq.count("C")
        
        tm = 4 * (G + C) + 2 * (A + T)
        return tm

def tm_2(seq, Na = 0.0001):
        G = seq.count("G")
        C = seq.count("C")
        length = len(seq)
        
        tm = -(81.5 + 16.6 * math.log10(Na) + 0.41 * ((G + C)/length*100) - 600/length)
        return tm

seq = "ACGTTT"
print("Sequence: " + seq)
print("Tm with first formula: " + str(tm_1(seq)) + "°C" )
print("Tm with second formula: " + str(tm_2(seq))+ "°C")

