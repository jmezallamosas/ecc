import pickle
from sklearn.preprocessing import MinMaxScaler, StandardScaler

THRESHOLD = 0.75

# preprocess data, putting some features into standard scaling and others into minmax scaling
def preprocess_data(X):
    stan = StandardScaler()
    stan_cols = ["prim_TSS_dist5p", "prim_TSS_dist3p", "sec_TSS_dist5p", "sec_TSS_dist3p", "length", 
                 "GC_content", "GA_content", "CA_content", "frac_A", "frac_C", "frac_T", "frac_G"]

    minmax = MinMaxScaler()
    minmax_cols = ["longest_A", "longest_C", "longest_T", "longest_G", "mfe"]
    minmax_cols.extend([(str(i) + "_MMs") for i in range(4)])
    minmax_cols.extend([("num_MMs_pos" + str(i + 1)) for i in range(25)])

    X_scaled = X.copy()
    X_scaled[stan_cols] = stan.fit_transform(X_scaled[stan_cols])
    X_scaled[minmax_cols] = minmax.fit_transform(X_scaled[minmax_cols])
    
    return X_scaled

# use the fitted model to predict new sgRNAs
def model_data(X, model_dir):
    with open(model_dir, "rb") as file:
        model = pickle.load(file)
        
    y = model.predict(X)
    
    return y

# given a dataframe of sgRNAs and a targeted strand, determine the scores for each of the sgRNAs
def score_sgRNAs(df, strand):
    df["same_strand"] = (df["strand"] == strand).astype(int)
    
    X = df.drop(columns = ["PAMloc", "strand", "sequence", "struct"])
    
    X_scaled = preprocess_data(X)
    y = model_data(X_scaled, "scoring_model/models/LGBM.pk1")
    
    return y