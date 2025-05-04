import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models, optimizers
from tensorflow.keras.layers import Conv1D, ReLU
from tensorflow.keras.regularizers import l2

def caldiag(Mat):
    diag_Mat = np.mean(np.diag(Mat))
    
    mask = np.eye(Mat.shape[0], dtype=bool)

    # Extract off-diagonal elements using the mask and compute their mean
    off_diag_Mat = Mat[~mask].mean()

    return diag_Mat,off_diag_Mat

def calrelM(X,Amat):
    
    n, k = X.shape
    pi = X.sum(axis=0) / (2 * n)  # cal maf
    P = pi[np.newaxis, :]  # shape (1, k)
    Z = X - 2 * P 

    G = (Z @ Z.T) / (2 * np.sum(pi * (1 - pi)))
    diag_G, off_diag_G = caldiag(G)

    #A = np.eye(n)
    A = Amat
    diag_A, off_diag_A = caldiag(A)

    b = (diag_A-off_diag_A)/(diag_G-off_diag_G)

    a = diag_A - diag_G * b

    G = a + b * G

    # scaling G
    G = 0.98 * G + 0.02 * A

    return G

def create_dataset(df,doubleT,genotype_columns=None,rel_columns=None,linear=None):
    # process the labels
    if not doubleT:
        labels = df['phenotype'].astype(np.float32).values.reshape(-1, 1)
    else:
        labels = {
            'phenotype_c': df['phenotype_c'].astype(np.float32).squeeze(),
            'phenotype_t': df['phenotype_t'].astype(np.float32).squeeze()
        }
        if linear:
            sample_weights = {
                'phenotype_c': df['phenotype_c_mask'].astype(np.float32).squeeze(),
                'phenotype_t': np.ones_like(df['phenotype_t'].astype(np.float32).squeeze())
            }
    
    # relationship features
    if rel_columns is not None and genotype_columns is not None:
        genotype = df[genotype_columns].astype(np.float32).values
        rel_features = df[rel_columns].astype(np.float32).values
        features = (genotype, rel_features)
    elif genotype_columns is not None and rel_columns is None:
        #genotype features
        genotype = df[genotype_columns].astype(np.float32).values
        features = genotype
    elif rel_columns is not None and genotype_columns is None:
        rel_features = df[rel_columns].astype(np.float32).values
        features = rel_features

    
    # making TensorFlow data set
    if not linear:
        dataset = tf.data.Dataset.from_tensor_slices((features, labels))
    else:
        dataset = tf.data.Dataset.from_tensor_slices((features, labels, sample_weights))
    
    return dataset


def create_test_dataset(df,batch_size,genotype_columns=None,rel_columns=None):

     # relationship features
    if rel_columns is not None and genotype_columns is not None:
        genotype = df[genotype_columns].astype(np.float32).values
        rel_features = df[rel_columns].astype(np.float32).values
        features = (genotype, rel_features)
    elif genotype_columns is not None and rel_columns is None:
        #genotype features
        genotype = df[genotype_columns].astype(np.float32).values
        features = genotype
    elif rel_columns is not None and genotype_columns is None:
        rel_features = df[rel_columns].astype(np.float32).values
        features = rel_features

    dataset = tf.data.Dataset.from_tensor_slices((features,))

    return dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)


class LocalLinear1D(layers.Layer):
    def __init__(self, local_features, kernel_size, stride=1, bias=True, kernel_regularizer=None, bias_regularizer=None, **kwargs):
        super(LocalLinear1D, self).__init__(**kwargs)
        self.local_features = local_features
        self.kernel_size = kernel_size
        self.stride = stride
        self.bias = bias
        self.kernel_regularizer = kernel_regularizer
        self.bias_regularizer = bias_regularizer

    def build(self, input_shape):
        in_features = input_shape[-1]
        self.padding = self.kernel_size - 1
        self.fold_num = (in_features + self.padding - self.kernel_size) // self.stride + 1

        # 初始化权重
        self.weight = self.add_weight(
            shape=(self.fold_num, self.kernel_size, self.local_features),
            initializer='glorot_uniform',
            trainable=True,
            name='weight',
            regularizer=self.kernel_regularizer
        )

        if self.bias:
            self.bias_var = self.add_weight(
                shape=(self.fold_num, self.local_features),
                initializer='zeros',
                trainable=True,
                name='bias',
                regularizer=self.bias_regularizer
            )
        else:
            self.bias_var = None

        super(LocalLinear1D, self).build(input_shape)

    def call(self, inputs):
        # 输入形状: (batch_size, in_features)
        x = tf.pad(inputs, [[0,0], [0, self.padding]], mode='CONSTANT')  # (batch_size, in_features + padding)

        # 提取滑动窗口
        x = tf.signal.frame(x, frame_length=self.kernel_size, frame_step=self.stride, axis=-1)  # (batch_size, fold_num, kernel_size)

        #x = tf.expand_dims(x, axis=-2) # (batch_size, fold_num, 1, kernel_size)

        #weight = tf.expand_dims(self.weight, axis=0) # (1, fold_num, kernel_size, local_features)

        #x = tf.matmul(x, weight) # (batch_size, fold_num, 1, local_features)

        #x = tf.squeeze(x, axis=-2) # (batch_size, fold_num, local_features)

        # 执行批量乘法
        x = tf.einsum('bfk,fkl->bfl', x, self.weight)  # (batch_size, fold_num, local_features)

        # 添加偏置
        if self.bias:
            x = x + self.bias_var  # (batch_size, fold_num, local_features)


        if self.local_features == 1:

            x = tf.squeeze(x, axis=-1) # (batch_size, fold_num)
            
        return x

    def get_config(self):
        config = super(LocalLinear1D, self).get_config()
        config.update({
            'local_features': self.local_features,
            'kernel_size': self.kernel_size,
            'stride': self.stride,
            'bias': self.bias
        })
        return config
    
def cnn_model(num_genotype_features,doubleT,Th,softm,linear):

    learning_rate = 1e-5
    # drate = 0.3
    if not doubleT:
        kerl = 3
        con1 = 64
        ker1 = 8
        stride1 = 6
        dense_units = 128
    elif softm:
        kerl = 3
        con1 = 32
        ker1 = 4
        stride1 = 6
        dense_units = 96
    elif linear:
        kerl = 3
        con1 = 96
        ker1 = 8
        stride1 = 6
        dense_units = 64
    else:
        kerl = 3
        con1 = 32
        ker1 = 8
        stride1 = 6
        dense_units = 96


    X_in = layers.Input(shape=(num_genotype_features,), name='genotype_input')

    x = X_in
    
    # x = LocalLinear1D(
    #     local_features=1,
    #     kernel_size=kerl,
    #     stride=1,
    #     bias=True
    # )(x)
    # x = layers.BatchNormalization()(x)
    # x = ReLU()(x)
    # x = layers.Dropout(drate)(x)
    x = tf.expand_dims(x, axis=-1)

    #residual = x

    x = Conv1D(filters=con1, kernel_size=ker1, strides = stride1,padding='valid', activation='relu')(x)

    x = Conv1D(filters=con1, kernel_size=ker1, strides = stride1,padding='valid', activation='relu')(x)

    x = layers.MaxPooling1D(pool_size=5)(x)

    x = layers.Flatten()(x)

    if doubleT:

        task_t = layers.Dense(dense_units,activation='relu')(x)
        # task_t = layers.BatchNormalization()(task_t)
        # task_t = layers.Dropout(drate)(task_t)

        task_c = layers.Dense(dense_units,activation='relu')(x)
        # task_c = layers.BatchNormalization()(task_c)
        # task_c = layers.Dropout(drate)(task_c)

    else:

        x = layers.Dense(dense_units,activation='relu')(x)
        # x = layers.BatchNormalization()(x)
        # x = layers.Dropout(drate)(x)

    z = x


    if not doubleT:
        
        if Th:
            # 二分类
            output = layers.Dense(1, activation='sigmoid')(z)
            loss = tf.keras.losses.BinaryCrossentropy()
            metrics = ['accuracy', tf.keras.metrics.AUC()]
        else:
            # 回归
            output = layers.Dense(1, activation=None)(z)
            loss = tf.keras.losses.MeanSquaredError()
            metrics = ['mse']
    else:

        output_t = layers.Dense(1, activation='sigmoid', name='phenotype_t')(task_t)

        if not softm and not linear:

            output_c = layers.Dense(1, activation='sigmoid', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': tf.keras.losses.BinaryCrossentropy(),
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy', tf.keras.metrics.AUC()],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

        elif softm:
            output_c = layers.Dense(5, activation='softmax', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': 'sparse_categorical_crossentropy',
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }
        
        elif linear:

            output_c = layers.Dense(1, activation='linear', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': 'mean_squared_error',
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['mean_squared_error', 'mean_absolute_error'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

            weighted_metrics = {
                'phenotype_c': ['mean_squared_error', 'mean_absolute_error'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()]
            }

        
        output = [output_c,output_t]

    model = models.Model(inputs=X_in, outputs=output)

    for layer in model.layers:

        if hasattr(layer, 'kernel_regularizer'):

            layer.kernel_regularizer = l2(0.01)
        
        if hasattr(layer, 'bias_regularizer'):
        
            layer.bias_regularizer = l2(0.01)

    optimizer = optimizers.Adam(learning_rate=learning_rate)
    
    if linear:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics,
            weighted_metrics=weighted_metrics
        )
    else:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics
        )

    return model

def fc_model(rel_features,doubleT,Th,softm,linear):

    learning_rate = 1e-5
    # drate = 0.3

    if not doubleT:

        # num_layers =  1
        # dense_units = 32
        dense_sp = 128
    elif softm:
        # num_layers =  1
        # dense_units = 96
        dense_sp = 64
    elif linear:
        # num_layers =  2
        # dense_units = 128
        dense_sp = 128
    else:
        # num_layers =  1
        # dense_units = 96
        dense_sp = 96


    rel_in = layers.Input(shape=(rel_features,), name='rel_input')
    
    x = rel_in

    # for i in range(num_layers):
    #     x = layers.Dense(dense_units, activation='relu')(x)
    #     x = layers.BatchNormalization()(x)
    #     x = layers.Dropout(drate)(x)

    if doubleT:

        task_t = layers.Dense(dense_sp,activation='relu')(x)
        # task_t = layers.BatchNormalization()(task_t)
        # task_t = layers.Dropout(drate)(task_t)

        task_c = layers.Dense(dense_sp,activation='relu')(x)
        # task_c = layers.BatchNormalization()(task_c)
        # task_c = layers.Dropout(drate)(task_c)
    else:
        x = layers.Dense(dense_sp,activation='relu')(x)
        # x = layers.BatchNormalization()(x)
        # x = layers.Dropout(drate)(x)

    z = x


    if not doubleT:
        if Th:
            # 二分类
            output = layers.Dense(1, activation='sigmoid')(z)
            loss = tf.keras.losses.BinaryCrossentropy()
            metrics = ['accuracy', tf.keras.metrics.AUC()]
        else:
            # 回归
            output = layers.Dense(1, activation=None)(z)
            loss = tf.keras.losses.MeanSquaredError()
            metrics = ['mse']
    else:

        output_t = layers.Dense(1, activation='sigmoid', name='phenotype_t')(task_t)

        if not softm and not linear:

            output_c = layers.Dense(1, activation='sigmoid', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': tf.keras.losses.BinaryCrossentropy(),
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy', tf.keras.metrics.AUC()],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

        elif softm:
            output_c = layers.Dense(5, activation='softmax', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': 'sparse_categorical_crossentropy',
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

        elif linear:
            
            output_c = layers.Dense(1, activation='linear', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': 'mean_squared_error',
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['mean_squared_error', 'mean_absolute_error'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

            weighted_metrics = {
                'phenotype_c': ['mean_squared_error', 'mean_absolute_error'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()]
            }

        
        output = [output_c,output_t]

    model = models.Model(inputs=rel_in, outputs=output)

    for layer in model.layers:

        if hasattr(layer, 'kernel_regularizer'):

            layer.kernel_regularizer = l2(0.01)
        
        if hasattr(layer, 'bias_regularizer'):
        
            layer.bias_regularizer = l2(0.01)

    optimizer = optimizers.Adam(learning_rate=learning_rate)
    
    if linear:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics,
            weighted_metrics=weighted_metrics
        )
    else:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics
        )

    return model

def mul_model(num_genotype_features,rel_features,doubleT,Th,softm,linear):

    learning_rate = 1e-5
    # drate = 0.3

    if not doubleT:

        kerl = 3
        con1 = 64
        ker1 = 8
        stride1 = 6
        # dense_units = 128
        # num_layers =  2
        dense_sp = 128
    elif softm:
        kerl = 6
        con1 = 96
        ker1 = 8
        stride1 = 8
        # dense_units = 32
        # num_layers =  1
        dense_sp = 32
    elif linear:
        kerl = 3
        con1 = 32
        ker1 = 4
        stride1 = 8
        # dense_units = 128
        # num_layers =  3
        dense_sp = 128
    else:
        kerl = 6
        con1 = 96
        ker1 = 4
        stride1 = 8
        # dense_units = 32
        # num_layers =  3
        dense_sp = 128

    X_in = layers.Input(shape=(num_genotype_features,), name='genotype_input') 

    rel_in = layers.Input(shape=(rel_features,), name='rel_input')

    #processing genotype


    x_geno = X_in
    
    x_geno = LocalLinear1D(
        local_features=1,
        kernel_size=kerl,
        stride=1,
        bias=True
    )(x_geno)
    x_geno = layers.BatchNormalization()(x_geno)
    x_geno = ReLU()(x_geno)
    # x_geno = layers.Dropout(drate)(x_geno)
    x_geno = tf.expand_dims(x_geno, axis=-1)

    x_geno = Conv1D(filters=con1, kernel_size=ker1, strides = stride1,padding='valid', activation='relu')(x_geno)

    # x_geno = Conv1D(filters=con2, kernel_size=ker2, strides = stride2,padding='valid', activation='relu')(x_geno)

    x_geno = layers.MaxPooling1D(pool_size=5)(x_geno)

    x_geno = layers.Flatten()(x_geno)

    #processing relationship


    x_rel = rel_in

    # for i in range(num_layers):
        # x_rel = layers.Dense(dense_units, activation='relu')(x_rel)
        # x_rel = layers.BatchNormalization()(x_rel)
        # x_rel = layers.Dropout(drate)(x_rel)

    #processing connected features

    
    x = layers.Concatenate()([x_geno, x_rel])

    if doubleT:

        task_t = layers.Dense(dense_sp,activation='relu')(x)
        # task_t = layers.BatchNormalization()(task_t)
        # task_t = layers.Dropout(drate)(task_t)

        task_c = layers.Dense(dense_sp,activation='relu')(x)
        # task_c = layers.BatchNormalization()(task_c)
        # task_c = layers.Dropout(drate)(task_c)
    else:
        x = layers.Dense(dense_sp,activation='relu')(x)
        # x = layers.BatchNormalization()(x)
        # x = layers.Dropout(drate)(x)


    z = x


    if not doubleT:

        if Th:
            # 二分类
            output = layers.Dense(1, activation='sigmoid')(z)
            loss = tf.keras.losses.BinaryCrossentropy()
            metrics = ['accuracy', tf.keras.metrics.AUC()]
        else:
            # 回归
            output = layers.Dense(1, activation=None)(z)
            loss = tf.keras.losses.MeanSquaredError()
            metrics = ['mse']
    else:

        output_t = layers.Dense(1, activation='sigmoid', name='phenotype_t')(task_t)

        if not softm and not linear:

            output_c = layers.Dense(1, activation='sigmoid', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': tf.keras.losses.BinaryCrossentropy(),
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy', tf.keras.metrics.AUC()],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

        elif softm:
            output_c = layers.Dense(5, activation='softmax', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': 'sparse_categorical_crossentropy',
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

        elif linear:
            
            output_c = layers.Dense(1, activation='linear', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': 'mean_squared_error',
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['mean_squared_error', 'mean_absolute_error'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

            weighted_metrics = {
                'phenotype_c': ['mean_squared_error', 'mean_absolute_error'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()]
            }

        
        output = [output_c,output_t]

    
    model = models.Model(inputs=[X_in, rel_in], outputs=output)
    

    for layer in model.layers:

        if hasattr(layer, 'kernel_regularizer'):

            layer.kernel_regularizer = l2(0.01)

        if hasattr(layer, 'bias_regularizer'):

            layer.bias_regularizer = l2(0.01)


    optimizer = optimizers.Adam(learning_rate=learning_rate)
    
    if linear:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics,
            weighted_metrics=weighted_metrics
        )
    else:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics
        )

    return model


def process_and_save_predictions(model_name, predictions, predicted_ids, doubleT, Th=None, softm=None,linear=None):

    if doubleT:
        if not softm:
            pred_c, pred_t = predictions[0].squeeze(-1), predictions[1].squeeze(-1)
        else:
            pred_c, pred_t = np.dot(predictions[0], np.arange(5)), predictions[1].squeeze(-1)

        df = pd.DataFrame({
            "ID": predicted_ids,
            "Prediction_c": pred_c,
            "Prediction_t": pred_t
        })

        if not softm and not linear:

            filename = f"predictions_c_t_{model_name}.csv"
        
        elif softm:
            filename = f"predictions_c_t_{model_name}_soft.csv"

        elif linear:
            
            filename = f"predictions_c_t_{model_name}_linear.csv"
    else:
        pred = predictions.squeeze(-1)
        df = pd.DataFrame({
            "ID": predicted_ids,
            "Prediction": pred
        })
        suffix = "t" if Th else "c"
        filename = f"predictions_{suffix}_{model_name}.csv"
    
    df.to_csv(filename, index=False)

    print(f"Saved predictions to {filename}")

