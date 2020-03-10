import pandas as pd

# row = data.loc[data['element']=='Cu']
class FormFactor:
    def __init__(self):
        data = pd.read_csv('formfactors.csv')
        self.FF_coeffs = {}
        for index, row in data.iterrows():
            element = str(row["element"])
            a1 = row["a1"]
            a2 = row["a2"]
            a3 = row["a3"]
            a4 = row["a4"]
            b1 = row["b1"]
            b2 = row["b2"]
            b3 = row["b3"]
            b4 = row["b4"]
            c = row["c"]
            self.FF_coeffs[element] = FF_data(a1, a2, a3, a4, b1, b2, b3, b4, c)
    def getCoeffs(self, element_id):
        return self.FF_coeffs[element_id]

class FF_data:
    def __init__(self, a1, a2, a3, a4, b1, b2, b3, b4, c):
        self.a_coeffs = [a1, a2, a3, a4]
        self.b_coeffs = [b1, b2, b3, b4]
        self.c = c
