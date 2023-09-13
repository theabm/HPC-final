import pandas as pd 
import plotly as plt 

mkl_float = pd.read_csv("./float/results_mkl.csv")
oblas_float = pd.read_csv("./float/results_oblas.csv")
blis_float = pd.read_csv("./float/results_blis.csv")

mkl_double = pd.read_csv("./double/results_mkl.csv")
oblas_double = pd.read_csv("./double/results_oblas.csv")
blis_double = pd.read_csv("./double/results_blis.csv")

# group by blocks of 10
mkl_float = mkl_float.groupby(mkl_float.index // 10)
oblas_float = oblas_float.groupby(oblas_float.index // 10)
blis_float = blis_float.groupby(blis_float.index // 10)

mkl_double = mkl_double.groupby(mkl_double.index // 10)
oblas_double = oblas_double.groupby(oblas_double.index // 10)
blis_double = blis_double.groupby(blis_double.index // 10)

mkl_float_avg = mkl_float.mean()
oblas_float_avg = oblas_float.mean()
blis_float_avg = blis_float.mean()

mkl_double_avg = mkl_double.mean()
oblas_double_avg = oblas_double.mean()
blis_double_avg = blis_double.mean()




