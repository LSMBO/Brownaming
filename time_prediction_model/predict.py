import pandas as pd
import pickle
import os
from sklearn.linear_model import LinearRegression

def predict_diamond_time(nb_query, dbsize, model_path='../diamond_time_model.pkl'):
    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    features = pd.DataFrame([[nb_query, dbsize]], columns=['nb_query', 'dbsize'])    
    predicted_time = model.predict(features)[0]
    return predicted_time

if __name__ == "__main__":
    test_data = [
        (77470, 16710507),
        (76289, 417608),
        (75588, 141714),
        (73962, 19389125),
        (72192, 25825154),
        (70296, 160384335),
        (67816, 104402)
    ]
    for test in test_data:
        nb_query, db_size = test
        time_estimate = predict_diamond_time(nb_query, db_size)
        print(f"{nb_query} queries against a DB with {db_size} entries -> Estimated time: {time_estimate:.2f} minutes")
        