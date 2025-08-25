#!/usr/bin/env python
# -*- encoding: utf-8 -*-
import pickle
import seaborn as sns
import warnings
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, precision_recall_curve, auc
from sklearn.ensemble import RandomForestClassifier
from joblib import dump, load

warnings.filterwarnings("ignore")

def train_model(model_data, out, random_state=42):
    """Reproducible 10-fold cross-validation training function"""
    np.random.seed(random_state)
    import random
    random.seed(random_state)
    
    rf_accuracy_scores, rf_recall_scores, rf_precision_scores, rf_f1_scores = [], [], [], []
    rf_precision_pr, rf_recall_pr, rf_auc_pr = [], [], []
    
    data_df = pd.read_excel(model_data, index_col=None, header=0, sheet_name="count.filter.signal.model")
    DEL_X = data_df[["Line_SumCount","Line_N","Line_MaxCount","Line_MeanCount","Line_Count_Ratio","Line_N_Ratio"]]
    DEL_y = data_df[["Line"]]
    
    feature_names = DEL_X.columns.tolist()
    feature_names_path = os.path.join(out, "feature_names.pkl")
    with open(feature_names_path, 'wb') as f:
        pickle.dump(feature_names, f)
    print(f"Feature names saved: {feature_names}")
    
    RF_clf = RandomForestClassifier(
        random_state=random_state,
        n_estimators=100,
        max_depth=None,
        min_samples_split=2,
        n_jobs=-1
    )

    skfolds = StratifiedKFold(n_splits=10, shuffle=True, random_state=random_state)
    
    print(f"Starting 10-fold cross-validation (random_state={random_state})")
    for fold, (train_index, test_index) in enumerate(skfolds.split(DEL_X, DEL_y), 1):
        X_train_fold = DEL_X.iloc[train_index]
        y_train_fold = DEL_y.iloc[train_index]
        X_test_fold = DEL_X.iloc[test_index]
        y_test_fold = DEL_y.iloc[test_index]
        
        RF_clf.fit(X_train_fold, y_train_fold.values.ravel())
        rf_y_pred = RF_clf.predict(X_test_fold)
        rf_y_pred_proba = RF_clf.predict_proba(X_test_fold)[:, 1]
        
        rf_accuracy_scores.append(accuracy_score(y_test_fold, rf_y_pred))
        rf_recall_scores.append(recall_score(y_test_fold, rf_y_pred))
        rf_precision_scores.append(precision_score(y_test_fold, rf_y_pred))
        rf_f1_scores.append(f1_score(y_test_fold, rf_y_pred))

        rf_precision, rf_recall, _ = precision_recall_curve(y_test_fold, rf_y_pred_proba)
        rf_auc_score = auc(rf_recall, rf_precision)
        
        rf_precision_pr.append(rf_precision)
        rf_recall_pr.append(rf_recall)
        rf_auc_pr.append(rf_auc_score)
        
        print(f'Fold {fold}: RF - Acc: {rf_accuracy_scores[-1]:.3f}, Pre: {rf_precision_scores[-1]:.3f}, Rec: {rf_recall_scores[-1]:.3f}, F1: {rf_f1_scores[-1]:.3f}, PR-AUC: {rf_auc_score:.3f}')
        
        fold_out_dir = os.path.join(out, f"fold_{fold}")
        os.makedirs(fold_out_dir, exist_ok=True)
        dump(RF_clf, os.path.join(fold_out_dir, "RF.pkl"))
    
    rf_mean_acc = np.average(rf_accuracy_scores)
    rf_mean_recall = np.average(rf_recall_scores)
    rf_mean_f1 = np.average(rf_f1_scores)
    rf_mean_precision = np.average(rf_precision_scores)
    rf_mean_auc_pr = np.average(rf_auc_pr)

    print("RF - Accuracy: %0.3f (+/- %0.2f)" % (rf_mean_acc, np.std(rf_accuracy_scores) * 2))
    print("RF - Recall: %0.3f (+/- %0.2f)" % (rf_mean_recall, np.std(rf_recall_scores) * 2))
    print("RF - F1_score: %0.3f (+/- %0.2f)" % (rf_mean_f1, np.std(rf_f1_scores) * 2))
    print("RF - PR-AUC: %0.3f (+/- %0.2f)" % (rf_mean_auc_pr, np.std(rf_auc_pr) * 2))
    
    plt.figure(figsize=(12, 12))
    folds = range(1, 11)
    plt.plot(folds, rf_precision_scores, "y-.", linewidth=4)
    plt.plot(folds, rf_accuracy_scores, "r-.", linewidth=4)
    plt.plot(folds, rf_recall_scores, "b-.", linewidth=4)
    plt.plot(folds, rf_f1_scores, "g-.", linewidth=4)
    plt.legend(["Precision: %0.3f" % rf_mean_precision, 
                "Accuracy: %0.3f" % rf_mean_acc, 
                "Recall: %0.3f" % rf_mean_recall, 
                "F1: %0.3f" % rf_mean_f1], fontsize=20)
    plt.ylim(0.5, 1.2)
    plt.xlabel("Fold", fontsize=22)
    plt.xticks(fontsize=20)
    plt.ylabel("Score", fontsize=22)
    plt.yticks(fontsize=20)
    rf_png = os.path.join(out, "RF_10fold.png")
    plt.savefig(rf_png, dpi=300, bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(12, 12))
    mean_recall_rf = np.linspace(0, 1, 100)
    mean_precision_rf = np.zeros_like(mean_recall_rf)
    for i in range(10):
        interp_precision = np.interp(mean_recall_rf, rf_recall_pr[i][::-1], rf_precision_pr[i][::-1])
        mean_precision_rf += interp_precision
    mean_precision_rf /= 10
    plt.plot(mean_recall_rf, mean_precision_rf, color=(0/255, 78/255, 162/255), 
             label=f'RF (PR-AUC = {rf_mean_auc_pr:.3f})', linewidth=4)
    positive_ratio = np.mean(DEL_y.values)
    plt.plot([0, 1], [positive_ratio, positive_ratio], 'k--', 
             label=f'Random (AUC = {positive_ratio:.3f})')
    plt.xlabel('Recall', fontsize=22)
    plt.xticks(fontsize=20)
    plt.ylabel('Precision', fontsize=22)
    plt.yticks(fontsize=20)
    plt.legend(loc='lower left', fontsize=20)
    plt.grid(True, alpha=0.3)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    
    pr_curve_png = os.path.join(out, "PR_Curve_10fold.png")
    plt.savefig(pr_curve_png, dpi=300, bbox_inches='tight')
    plt.close()
    
    return RF_clf, feature_names

def feature_visual(model_path, out_dir, fold_number=None):
    """Visualize feature importance for a model"""
    feature_names_file = os.path.join(out_dir, "feature_names.pkl")
    
    try:
        model = load(model_path)
        with open(feature_names_file, 'rb') as f:
            feature_names = pickle.load(f)
    except FileNotFoundError:
        feature_importances = model.feature_importances_
        n_features = len(feature_importances)
        feature_names = [f"Feature_{i}" for i in range(n_features)]
        print("Using default feature indices as feature names")
    
    feature_importances = model.feature_importances_
    
    features_df = pd.DataFrame({
        "Features": feature_names,
        "Importance": feature_importances
    })
    features_df["Importance"] = features_df["Importance"].apply(lambda x: float(f"{x:.2f}"))
    features_df.sort_values("Importance", inplace=True, ascending=False)
    
    plt.figure(figsize=(12, 12))
    ax = sns.barplot(x=features_df["Importance"], y=features_df["Features"], 
                    color=(0/255, 78/255, 162/255))
    
    plt.xlabel("Importance", fontsize=22)
    plt.ylabel("Features", fontsize=22)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    
    title_suffix = f" (Fold {fold_number})" if fold_number is not None else " (Final Model)"
    plt.title(f"Feature Importances{title_suffix}", fontsize=24)
    
    ax.set_facecolor('white')
    
    for i, v in enumerate(features_df["Importance"]):
        ax.text(v + 0.001, i, f'{v:.2f}', va='center', fontsize=16)
    
    png_out = os.path.join(out_dir, f"fold_{fold_number}", f"feature_importance_fold_{fold_number}.png")
    csv_out = os.path.join(out_dir, f"fold_{fold_number}", f"feature_importance_fold_{fold_number}.csv")
    
    plt.savefig(png_out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    features_df.to_csv(csv_out, index=False)
    
    return features_df

def all_feature_visual(main_output_dir):
    """Visualize feature importance for all folds"""
    fold_dirs = [d for d in os.listdir(main_output_dir) if d.startswith("fold_") and os.path.isdir(os.path.join(main_output_dir, d))]
    
    all_importances = []
    for fold_dir in fold_dirs:
        fold_num = fold_dir.split("_")[1]
        model_path = os.path.join(main_output_dir, fold_dir, "RF.pkl")
        
        if os.path.exists(model_path):
            importance_df = feature_visual(model_path, main_output_dir, fold_number=fold_num)
            all_importances.append(importance_df)
    
    if all_importances:
        combined_df = pd.concat(all_importances)
        mean_importance = combined_df.groupby('Features')['Importance'].mean().reset_index()
        mean_importance.sort_values('Importance', ascending=False, inplace=True)
        mean_importance["Importance"] = mean_importance["Importance"].apply(lambda x: float(f"{x:.2f}"))
        plt.figure(figsize=(12, 10))
        ax = sns.barplot(x=mean_importance["Importance"], y=mean_importance["Features"], 
                        color=(0/255, 78/255, 162/255))
        
        plt.xlabel("Mean Importance", fontsize=22)
        plt.ylabel("Features", fontsize=22)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.title("Average Feature Importances (All Folds)", fontsize=24)
        
        ax.set_facecolor('white')
        for i, v in enumerate(mean_importance["Importance"]):
            ax.text(v + 0.001, i, f'{v:.2f}', va='center', fontsize=16)
        
        avg_png = os.path.join(main_output_dir, "feature_importance_average.png")
        avg_csv = os.path.join(main_output_dir, "feature_importance_average.csv")
        
        plt.savefig(avg_png, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        mean_importance.to_csv(avg_csv, index=False)
        print(f"Average feature importance saved to: {avg_png}, {avg_csv}")
    
    return all_importances

if __name__ == "__main__":
    RF_clf, feature_names = train_model("./model/train_dataset.xlsx", "./model/")
    all_importances = all_feature_visual("./model/")