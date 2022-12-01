import numpy as np
from metric_learn import ITML
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def learn_itml(input_matr):
  
  itml = ITML(preprocessor = input_matr)
  tuples_indices = np.array([[0,1],[0,2],[0,3],[0,4],[0,5]])
  y = np.array([1,1,1,1,1])
  try:
    itml.fit(tuples_indices, y)
  except ValueError:
    print("ValueError occured!")
    return None
  return(itml.get_mahalanobis_matrix())



