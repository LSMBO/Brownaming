import pandas as pd
import pickle
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score

def train_diamond_time_model(data_path='data_file.tsv'):
    data = pd.read_csv(data_path, sep='\t')
    
    X = data[['nb_query', 'dbsize']]
    y = data['time']
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    print(f"Training data: {len(X_train)} samples")
    print(f"Test data: {len(X_test)} samples")

    print("Training multiple linear regression model...")
    model = LinearRegression()
    model.fit(X_train, y_train)
    
    y_pred = model.predict(X_test)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    
    print("=== RESULTS ===")
    print(f"R² score: {r2:.4f}")
    print(f"MAE: {mae:.4f} minutes")
    print(f"Model formula:")
    print(f"Time (minutes) = {model.intercept_:.4f} + {model.coef_[0]:.4f} * nb_query + {model.coef_[1]:.4f} * dbsize")
    
    # Save the model
    model_path = '../diamond_time_model.pkl'
    with open(model_path, 'wb') as f:
        pickle.dump(model, f)
    print(f"Model saved to {model_path}")
    
    return model, {'r2': r2, 'mae': mae}

if __name__ == "__main__":
    train_diamond_time_model()
    