import random
import pandas as pd
import numpy as np
import tensorflow as tf
import keras_tuner as kt
import matplotlib.pyplot as plt
import argparse

from utils_tunning import (
    LocalLinear1D, 
    create_test_dataset, 
    create_dataset, 
    calrelM, 
    fc_model, 
    cnn_model, 
    mul_model, 
    process_and_save_predictions
)

parser = argparse.ArgumentParser(description="Run a model with configurable parameters.")
parser.add_argument("--doubleT", action="store_true", help="Enable doubleT mode")
parser.add_argument("--Th", action="store_true", help="Enable Th mode")
parser.add_argument("--softm", action="store_true", help="Enable softmax mode")
parser.add_argument("--linear", action="store_true", help="Enable linear mode")
args = parser.parse_args()

doubleT = args.doubleT
Th = args.Th
softm = args.softm
linear = args.linear

seed = 42  # 选择一个固定的种子

#=========================loading and merging phenotype and genotype files========================
genotype_df = pd.read_csv('geno.csv')
genotype_df = genotype_df.rename(columns={genotype_df.columns[0]: "ID"})

if not doubleT:
    phenotype_file = 'phenotype_t.csv' if Th else 'phenotype_c.csv'
    phenotype_df = pd.read_csv(phenotype_file)
    mergedDF = pd.merge(genotype_df, phenotype_df, on='ID')
    genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype']]
else:
    phenotype_c_df = pd.read_csv('phenotype_c.csv')
    phenotype_c_df = phenotype_c_df.rename(columns={'phenotype': 'phenotype_c'})
    phenotype_t_df = pd.read_csv('phenotype_t.csv')
    phenotype_t_df = phenotype_t_df.rename(columns={'phenotype': 'phenotype_t'})
    mergedDF = pd.merge(genotype_df, phenotype_c_df, on='ID')
    mergedDF = pd.merge(mergedDF, phenotype_t_df, on=['ID'])
    genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype_c', 'phenotype_t']]
    mergedDF['phenotype_c_mask'] = (mergedDF['phenotype_c'] != -999).astype(float)

    # if softm:
    #     phenotype_c_values = mergedDF.loc[mergedDF['phenotype_t'] == 1, 'phenotype_c']
    #     percentiles = np.percentile(phenotype_c_values, [20, 50, 80])
    #     conditions = [
    #         (mergedDF['phenotype_t'] == 1) & (mergedDF['phenotype_c'] <= percentiles[0]),
    #         (mergedDF['phenotype_t'] == 1) & (percentiles[0] < mergedDF['phenotype_c']) & (mergedDF['phenotype_c'] <= percentiles[1]),
    #         (mergedDF['phenotype_t'] == 1) & (percentiles[1] < mergedDF['phenotype_c']) & (mergedDF['phenotype_c'] <= percentiles[2]),
    #         (mergedDF['phenotype_t'] == 1) & (mergedDF['phenotype_c'] > percentiles[2]),
    #         (mergedDF['phenotype_t'] == 0)
    #     ]
    #     choices = [1, 2, 3, 4, 0]

    #     mergedDF['phenotype_c'] = np.select(conditions, choices)
    
    # elif linear:

    #     endvalue = np.percentile(mergedDF.loc[mergedDF['phenotype_t'] == 1, 'phenotype_c'],80)

    #     meanvalue = np.mean(mergedDF.loc[mergedDF['phenotype_t'] == 1, 'phenotype_c'])
    
    #     mergedDF.loc[mergedDF['phenotype_c'] > endvalue, 'phenotype_c'] = endvalue

    #     mergedDF.loc[mergedDF['phenotype_t'] == 0, 'phenotype_c'] = 999

    #     mergedDF['phenotype_c_mask'] = (mergedDF['phenotype_c'] != 999).astype(float)

    # else:

    #     medvalue = np.median(mergedDF.loc[mergedDF['phenotype_t'] == 1, 'phenotype_c'])
    
    #     # Assign 1 if 'phenotype_c' is greater than the median, otherwise assign 0
    #     mergedDF['phenotype_c'] = (mergedDF['phenotype_c'] > medvalue).astype(int)

    #     mergedDF.loc[mergedDF['phenotype_t'] == 0, 'phenotype_c'] = 0


genotype_data = mergedDF[genotype_columns].values.astype(np.float32)
#===================================================================================================


#==============================Prepare the relationship features==============================
# triangular_data = np.fromfile('Amat_b.PA.bin', dtype=np.float32)
# relationship_ids = pd.read_csv('Amat_b.PA.id', header=None)[0].tolist()
# num_individuals = int((np.sqrt(8 * len(triangular_data) + 1) - 1) / 2)
# # Initialize an empty square matrix
# relationship_matrix = np.zeros((num_individuals, num_individuals), dtype=np.float32)
# relationship_matrix[np.tril_indices(num_individuals)] = triangular_data
# relationship_matrix = relationship_matrix + relationship_matrix.T - np.diag(relationship_matrix.diagonal())
# rel_id_to_index = {id_: idx for idx, id_ in enumerate(relationship_ids)}
# rel_indices = [rel_id_to_index[id_] for id_ in mergedDF['ID'].tolist()]
# Amat = relationship_matrix[np.ix_(rel_indices, rel_indices)]

Gmat = calrelM(X=genotype_data,Amat=None)

L = np.linalg.cholesky(Gmat)
Z = np.eye(L.shape[0])
relationship_features = Z @ L
relationship_features = pd.DataFrame(relationship_features)

mergedDF = pd.concat([mergedDF, relationship_features], axis=1)

if doubleT:
    rel_columns = [col for col in mergedDF.columns if col not in (['ID', 'phenotype_c', 'phenotype_t', 'FamilyID']+genotype_columns)]
else:
    rel_columns = [col for col in mergedDF.columns if col not in (['ID', 'phenotype', 'FamilyID']+genotype_columns)]


#===============================================================================================


#============loading the train, val and test data set========================
val_ids = pd.read_csv('phe_val.csv')['ID'].tolist()
test_ids = pd.read_csv('phe_test.csv')['ID'].tolist()

# Create masks for validation and training data
is_in_val_ids = mergedDF['ID'].isin(val_ids)
is_in_train_ids = ~mergedDF['ID'].isin(val_ids + test_ids)
is_in_test_ids = mergedDF['ID'].isin(test_ids)


val_df = mergedDF[is_in_val_ids]

train_df = mergedDF[is_in_train_ids]

train_dataset_tf = create_dataset(df = train_df, doubleT = doubleT,genotype_columns = genotype_columns,linear=linear)

val_dataset_tf = create_dataset(df = val_df, doubleT = doubleT,genotype_columns = genotype_columns,linear=linear)

rel_train_dataset_tf = create_dataset(df = train_df, doubleT = doubleT,rel_columns = rel_columns,linear=linear)

rel_val_dataset_tf = create_dataset(df = val_df, doubleT = doubleT,rel_columns = rel_columns,linear=linear)

mul_train_dataset_tf = create_dataset(df = train_df, doubleT = doubleT,genotype_columns = genotype_columns,rel_columns = rel_columns,linear=linear)

mul_val_dataset_tf = create_dataset(df = val_df, doubleT = doubleT,genotype_columns = genotype_columns,rel_columns = rel_columns,linear=linear)
#=======================================================================================================================================





#===============================================#===============================================


#========================================training the model========================================
if not doubleT:
     sfx=""
elif doubleT and softm:
    sfx = "_soft"

elif doubleT and linear:
     sfx = "_linear"

elif doubleT and not softm and not linear:
    sfx = "_mt"

np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)

stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5)

#tunning the hyperparameters

def CNN_builder(hp):
    return cnn_model(hp,num_genotype_features = len(genotype_columns),doubleT=doubleT,Th = Th,softm = softm,linear=linear)

def FC_builder(hp):
    return fc_model(hp,rel_features = len(rel_columns),doubleT=doubleT,Th = Th,softm = softm,linear=linear)

def MUL_builder(hp):
    return mul_model(hp,num_genotype_features = len(genotype_columns),rel_features = len(rel_columns),doubleT=doubleT,Th = Th,softm = softm,linear=linear)

class CNNTuner(kt.Hyperband):
        
        def run_trial(self, trial):

            hp = trial.hyperparameters

            # Get the batch_size hyperparameter
            #batch_size = hp.Int('batch_size', min_value=50, max_value=100, step=16)
            batch_size = hp.Choice('batch_size', [30, 50, 100])
            #batch_size = 100

            # Rebuild data loaders with the batch_size
            train_loader = train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)
            val_loader = val_dataset_tf.batch(batch_size).prefetch(tf.data.AUTOTUNE)

            # Build the model
            model = self.hypermodel.build(hp)

            # Train the model
            history = model.fit(
                train_loader,
                epochs=100,
                validation_data=val_loader,
                callbacks=[stop_early],
            )

            # Save the metrics
            return history
        
tuner_cnn = CNNTuner(
            hypermodel=CNN_builder,
            objective='val_loss',
            max_epochs=100,
            factor=3,
            directory=f'my_dir_cnn{sfx}',
            project_name=f'cnn_tuning{sfx}'
        )

tuner_cnn.search()



np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
class FCTuner(kt.Hyperband):
        
        def run_trial(self, trial):

            hp = trial.hyperparameters

            # Get the batch_size hyperparameter
            #batch_size = hp.Int('batch_size', min_value=50, max_value=100, step=16)
            batch_size = hp.Choice('batch_size', [30, 50, 100])
            #batch_size = 100

            # Rebuild data loaders with the batch_size
            rel_train_loader = rel_train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)
            rel_val_loader = rel_val_dataset_tf.batch(batch_size).prefetch(tf.data.AUTOTUNE)

            # Build the model
            model = self.hypermodel.build(hp)

            # Train the model
            history = model.fit(
                rel_train_loader,
                epochs=100,
                validation_data=rel_val_loader,
                callbacks=[stop_early],
            )

            # Save the metrics
            return history
        
tuner_fc = FCTuner(
            hypermodel=FC_builder,
            objective='val_loss',
            max_epochs=100,
            factor=3,
            directory=f'my_dir_rel{sfx}',
            project_name=f'rel_tuning{sfx}'
        )

tuner_fc.search()


np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
class MULTuner(kt.Hyperband):
        
        def run_trial(self, trial):

            hp = trial.hyperparameters

            # Get the batch_size hyperparameter
            #batch_size = hp.Int('batch_size', min_value=50, max_value=100, step=16)
            batch_size = hp.Choice('batch_size', [30, 50, 100])
            #batch_size = 100

            # Rebuild data loaders with the batch_size
            mul_train_loader = mul_train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)
            mul_val_loader = mul_val_dataset_tf.batch(batch_size).prefetch(tf.data.AUTOTUNE)

            # Build the model
            model = self.hypermodel.build(hp)

            # Train the model
            history = model.fit(
                mul_train_loader,
                epochs=100,
                validation_data=mul_val_loader,
                callbacks=[stop_early],
            )

            # Save the metrics
            return history
        
tuner_mul = MULTuner(
            hypermodel=MUL_builder,
            objective='val_loss',
            max_epochs=100,
            factor=3,
            directory=f'my_dir_mul{sfx}',
            project_name=f'mul_tuning{sfx}'
        )

tuner_mul.search()



best_hps_cnn = tuner_cnn.get_best_hyperparameters(num_trials=1)[0]
print(f"cnn{sfx}最佳超参数：")
for param in best_hps_cnn.values:
    print(f"{param}: {best_hps_cnn.get(param)}")
cnn = tuner_cnn.hypermodel.build(best_hps_cnn)

best_hps_fc = tuner_fc.get_best_hyperparameters(num_trials=1)[0]
print(f"fc{sfx}最佳超参数：")
for param in best_hps_fc.values:
    print(f"{param}: {best_hps_fc.get(param)}")
fc = tuner_fc.hypermodel.build(best_hps_fc)


best_hps_mul = tuner_mul.get_best_hyperparameters(num_trials=1)[0]
print(f"mul{sfx}最佳超参数：")
for param in best_hps_mul.values:
    print(f"{param}: {best_hps_mul.get(param)}")
multimodel = tuner_mul.hypermodel.build(best_hps_mul)


histories = []
epochs = 150

np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
batch_size_cnn = best_hps_cnn.get('batch_size')
train_loader = train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size_cnn).prefetch(tf.data.AUTOTUNE)
val_loader = val_dataset_tf.batch(batch_size_cnn).prefetch(tf.data.AUTOTUNE)
history_cnn = cnn.fit(
            train_loader,
            validation_data=val_loader,
            callbacks=[stop_early],
            epochs=epochs, 
            initial_epoch=0
        )
histories.append(history_cnn)

np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
batch_size_fc = best_hps_fc.get('batch_size')
rel_train_loader = rel_train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size_fc).prefetch(tf.data.AUTOTUNE)
rel_val_loader = rel_val_dataset_tf.batch(batch_size_fc).prefetch(tf.data.AUTOTUNE)
history_fc = fc.fit(
            rel_train_loader,
            validation_data=rel_val_loader,
            callbacks=[stop_early],
            epochs=epochs, 
            initial_epoch=0
        )
histories.append(history_fc)



batch_size_mul = best_hps_mul.get('batch_size')
np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
mul_train_loader = mul_train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size_mul).prefetch(tf.data.AUTOTUNE)
mul_val_loader = mul_val_dataset_tf.batch(batch_size_mul).prefetch(tf.data.AUTOTUNE)
history_mul = multimodel.fit(
            mul_train_loader,
            validation_data=mul_val_loader,
            callbacks=[stop_early],
            epochs=epochs, 
            initial_epoch=0
        )
histories.append(history_mul)

for idx, history in enumerate(histories, 1):
    train_loss = history.history['loss']
    val_loss = history.history['val_loss']
    epochs_range = range(1, len(train_loss) + 1)
    plt.figure(figsize=(8, 5))
    plt.plot(epochs_range, train_loss, label='Train Loss', color='blue', linestyle='-')
    plt.plot(epochs_range, val_loss, label='Validation Loss', color='orange', linestyle='--')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title(f'Loss Curve for Model {idx}')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # 保存图像
    plt.savefig(f'loss_curve_{idx}.pdf', format='pdf')
    plt.show()


#prediction for the test set

test_df = mergedDF[mergedDF['ID'].isin(test_ids)].reset_index(drop=True)
predicted_ids = test_df['ID'].tolist()

test_dataset_tf_cnn = create_test_dataset(df=test_df,batch_size=batch_size_cnn,genotype_columns=genotype_columns)
test_dataset_tf_fc = create_test_dataset(df=test_df,batch_size=batch_size_mul,rel_columns=rel_columns)
test_dataset_tf_mul = create_test_dataset(df=test_df,batch_size=batch_size_fc,genotype_columns = genotype_columns,rel_columns=rel_columns)


predictions_cnn = cnn.predict(test_dataset_tf_cnn, verbose=1)
predictions_fc = fc.predict(test_dataset_tf_fc, verbose=1)
predictions_mul = multimodel.predict(test_dataset_tf_mul, verbose=1)

# Define a dictionary mapping model names to their prediction outputs
models_predictions = {
    "cnn": predictions_cnn,
    "fc": predictions_fc,
    "mul": predictions_mul
}

# Iterate over each model and its predictions
for model_name, predictions in models_predictions.items():
    process_and_save_predictions(
        model_name=model_name,
        predictions=predictions,
        predicted_ids=predicted_ids,
        doubleT=doubleT,
        Th=Th, # Ensure Th is defined in your context
        linear=linear,
        softm=softm  # Ensure softm is defined in your context
    )