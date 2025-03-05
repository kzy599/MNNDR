import random
import pandas as pd
import numpy as np
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from tensorflow.keras import layers, models, optimizers
from tensorflow.keras.layers import Conv1D, ReLU, MultiHeadAttention, TimeDistributed, LayerNormalization
import keras_tuner as kt
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from tensorflow.keras.regularizers import l2

seed = 42  # 选择一个固定的种子

# 设置任务类型：True 表示二分类，False 表示回归
Th = True  # 根据需要调整
new = False
rel = False
nonlinear = False
rasieLDim = False
pureSnp = False
attention = False
uspos = False
addDense = False
search = False
usingPca = False
newcnn = True
usingall = False
doubleT = True
connect = False
softm = True
# ============================
# 步骤1：加载和准备数据
# ============================

# 加载基因型数据
genotype_df = pd.read_csv('geno.csv')
genotype_df = genotype_df.rename(columns={genotype_df.columns[0]: "ID"})

if not doubleT:
    # 加载表型数据
    phenotype_file = 'phenotype_t.csv' if Th else 'phenotype_c.csv'
    phenotype_df = pd.read_csv(phenotype_file)

    # 合并基因型和表型数据
    mergedDF = pd.merge(genotype_df, phenotype_df, on='ID')

    # 提取基因型特征列
    genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype', 'FamilyID']]
else:
    # 加载表型数据
    phenotype_c_df = pd.read_csv('phenotype_c.csv')
    phenotype_c_df = phenotype_c_df.rename(columns={'phenotype': 'phenotype_c'})

    phenotype_t_df = pd.read_csv('phenotype_t.csv')
    phenotype_t_df = phenotype_t_df.rename(columns={'phenotype': 'phenotype_t'})

    # 合并基因型和表型数据
    mergedDF = pd.merge(genotype_df, phenotype_c_df, on='ID')
    mergedDF = pd.merge(mergedDF, phenotype_t_df, on=['ID','FamilyID'])

    #mergedDF = mergedDF[mergedDF['FamilyID'] == mergedDF['FamilyID'].iloc[0]]


    genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype_c', 'phenotype_t', 'FamilyID']]

    if softm:
        # 计算20%、50%和80%的百分位数
        phenotype_c_values = mergedDF.loc[mergedDF['phenotype_t'] == 1, 'phenotype_c']
        percentiles = np.percentile(phenotype_c_values, [20, 50, 80])

        # 根据百分位数和条件直接更新 'phenotype_c'
        conditions = [
            (mergedDF['phenotype_t'] == 1) & (mergedDF['phenotype_c'] <= percentiles[0]),
            (mergedDF['phenotype_t'] == 1) & (percentiles[0] < mergedDF['phenotype_c']) & (mergedDF['phenotype_c'] <= percentiles[1]),
            (mergedDF['phenotype_t'] == 1) & (percentiles[1] < mergedDF['phenotype_c']) & (mergedDF['phenotype_c'] <= percentiles[2]),
            (mergedDF['phenotype_t'] == 1) & (mergedDF['phenotype_c'] > percentiles[2]),
            (mergedDF['phenotype_t'] == 0)
        ]
        choices = [1, 2, 3, 4, 0]

        mergedDF['phenotype_c'] = np.select(conditions, choices)

    
    else:

        medvalue = np.median(mergedDF.loc[mergedDF['phenotype_t'] == 1, 'phenotype_c'])
    
        # Assign 1 if 'phenotype_c' is greater than the median, otherwise assign 0
        mergedDF['phenotype_c'] = (mergedDF['phenotype_c'] > medvalue).astype(int)

        mergedDF.loc[mergedDF['phenotype_t'] == 0, 'phenotype_c'] = 0

genotype_data = mergedDF[genotype_columns].values.astype(np.float32)

if usingPca:
    pca = PCA(n_components=0.95)
    genotype_data = pca.fit_transform(genotype_data)
    genotype_data = pd.DataFrame(genotype_data)
    mergedDF = mergedDF.drop(columns=genotype_columns)
    mergedDF = pd.concat([mergedDF, genotype_data], axis=1)
    if doubleT:
        genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype_c', 'phenotype_t', 'FamilyID']]
    else:
        genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype', 'FamilyID']]

if connect:
    triangular_data = np.fromfile('Amat_b.PA.bin', dtype=np.float32)
    relationship_ids = pd.read_csv('Amat_b.PA.id', header=None)[0].tolist()
    num_individuals = int((np.sqrt(8 * len(triangular_data) + 1) - 1) / 2)
    # Initialize an empty square matrix
    relationship_matrix = np.zeros((num_individuals, num_individuals), dtype=np.float32)
    relationship_matrix[np.tril_indices(num_individuals)] = triangular_data
    relationship_matrix = relationship_matrix + relationship_matrix.T - np.diag(relationship_matrix.diagonal())
    rel_id_to_index = {id_: idx for idx, id_ in enumerate(relationship_ids)}
    rel_indices = [rel_id_to_index[id_] for id_ in mergedDF['ID'].tolist()]
    relationship_matrix_aligned = relationship_matrix[np.ix_(rel_indices, rel_indices)]
    L = np.linalg.cholesky(relationship_matrix_aligned)
    Z = np.eye(L.shape[0])
    relationship_features = Z @ L
    relationship_features = pd.DataFrame(relationship_features)
    mergedDF = pd.concat([mergedDF, relationship_features], axis=1)
    if doubleT:
        rel_columns = [col for col in mergedDF.columns if col not in (['ID', 'phenotype_c', 'phenotype_t', 'FamilyID']+genotype_columns)]
    else:
        rel_columns = [col for col in mergedDF.columns if col not in (['ID', 'phenotype', 'FamilyID']+genotype_columns)]


def calfold(in_features,stride):
            return (in_features  - 1 ) // stride + 1

def clnm(r, k):
    return f"r{r}_k{k}"

if new:
    #triangular_data = np.fromfile('Amat_b.PA.bin', dtype=np.float32)
    #relationship_ids = pd.read_csv('Amat_b.PA.id', header=None)[0].tolist()
    #num_individuals = int((np.sqrt(8 * len(triangular_data) + 1) - 1) / 2)
    # Initialize an empty square matrix
    #relationship_matrix = np.zeros((num_individuals, num_individuals), dtype=np.float32)
    #relationship_matrix[np.tril_indices(num_individuals)] = triangular_data
    #relationship_matrix = relationship_matrix + relationship_matrix.T - np.diag(relationship_matrix.diagonal())
    #rel_id_to_index = {id_: idx for idx, id_ in enumerate(relationship_ids)}
    #rel_indices = [rel_id_to_index[id_] for id_ in mergedDF['ID'].tolist()]
    #Amat = relationship_matrix[np.ix_(rel_indices, rel_indices)]
    def calrelM(X):
        """
        计算基因型相关矩阵 G。
        
        参数：
        X (numpy.ndarray): 基因型矩阵（n 个体 x k 标记）。
        
        返回：
        G (numpy.ndarray): 规范化后的基因型相关矩阵（n 个体 x n 个体）。
        """
        n, k = X.shape
        pi = X.sum(axis=0) / (2 * n)  # 计算等位基因频率
        P = pi[np.newaxis, :]  # 形状为 (1, k)
        Z = X - 2 * P  # 中心化基因型矩阵

        # 计算 G
        G = (Z @ Z.T) / (2 * np.sum(pi * (1 - pi)))

        # 定义 A 作为虚拟家系（身份矩阵）
        A = np.eye(n)
        #A = Amat

        # 使用虚拟家系对 G 进行重新缩放
        G = 0.99 * G + 0.01 * A

        return G

    # 计算全局关系矩阵
    globel_matrix = calrelM(genotype_data)

    if not rel:

        L = np.linalg.cholesky(globel_matrix)
        Z = np.eye(L.shape[0])
        relationship_features = Z @ L
        relationship_features = pd.DataFrame(relationship_features)
        mergedDF.drop(columns=genotype_columns)
        mergedDF = pd.concat([mergedDF, relationship_features], axis=1)
        if doubleT:
            genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype_c', 'phenotype_t', 'FamilyID']]
        else:
            genotype_columns = [col for col in mergedDF.columns if col not in ['ID', 'phenotype', 'FamilyID']]
        

        #new_features = globel_matrix @ genotype_data

        #genotype_scaler = StandardScaler()

        #features_scaler = StandardScaler()

        #genotype_data = genotype_scaler.fit_transform(genotype_data) + features_scaler.fit_transform(new_features)

        #genotype_data = genotype_data + new_features/(2*genotype_data.shape[0])

        #mergedDF[genotype_columns] = genotype_data
    else:
        relationship_vectors = globel_matrix  # 形状为 (num_individuals, num_individuals)

        # 将关系向量转换为 DataFrame，以便与mergedDF合并
        relationship_df = pd.DataFrame(relationship_vectors, index=mergedDF['ID'], columns=mergedDF['ID'])

        # 将关系向量添加到 mergedDF 中，每个个体有一个与所有其他个体的关系向量
        # 这里我们将关系向量存储为嵌套列表（列表中的每个元素是一个列表）
        mergedDF['relationship_vector'] = mergedDF['ID'].apply(lambda id_: relationship_df.loc[id_].values)

val_ids = pd.read_csv('phe_val.csv')['ID'].tolist()
test_ids = pd.read_csv('phe_test.csv')['ID'].tolist()

# Create masks for validation and training data
is_in_val_ids = mergedDF['ID'].isin(val_ids)
is_in_train_ids = ~mergedDF['ID'].isin(val_ids + test_ids)
is_in_test_ids = mergedDF['ID'].isin(test_ids)

if usingall:
    all_train_ids = ~mergedDF['ID'].isin(test_ids)
    all_train_df = mergedDF[all_train_ids]
    labels_train_all = all_train_df['phenotype'].astype(np.float32).values.reshape(-1, 1)

# Create validation and training DataFrames
val_df = mergedDF[is_in_val_ids]
train_df = mergedDF[is_in_train_ids]
#data_not_in_test_ids = mergedDF[~is_in_test_ids]

# Now, split the data_not_in_test_ids into train and validation sets
#train_df, val_df = train_test_split(
#    data_not_in_test_ids,
#    test_size=0.1,
#    random_state=42
#)

print(f"训练集包含 {len(train_df)} 个个体。")
print(f"验证集包含 {len(val_df)} 个个体。")

if not doubleT:
    labels_train = train_df['phenotype'].astype(np.float32).values.reshape(-1, 1)

    labels_val = val_df['phenotype'].astype(np.float32).values.reshape(-1, 1)
else:
    labels_train = {
        'phenotype_c': train_df['phenotype_c'].astype(np.float32).values.reshape(-1, 1),
        'phenotype_t': train_df['phenotype_t'].astype(np.float32).values.reshape(-1, 1)
    }

    labels_val = {
        'phenotype_c': val_df['phenotype_c'].astype(np.float32).values.reshape(-1, 1),
        'phenotype_t': val_df['phenotype_t'].astype(np.float32).values.reshape(-1, 1)
    }

# ============================
# 步骤4：定义自定义的 Local CNN 层
# ============================

class LocalLinear2D(layers.Layer):
    def __init__(self, local_features, kernel_size, stride=1, bias=True, kernel_regularizer=None, bias_regularizer=None, **kwargs):
        super(LocalLinear2D, self).__init__(**kwargs)
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

        super(LocalLinear2D, self).build(input_shape)

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
        config = super(LocalLinear2D, self).get_config()
        config.update({
            'local_features': self.local_features,
            'kernel_size': self.kernel_size,
            'stride': self.stride,
            'bias': self.bias
        })
        return config

class PositionalEncoding(tf.keras.layers.Layer):
    def __init__(self, max_len, embedding_dim):
        super(PositionalEncoding, self).__init__()
        self.max_len = max_len
        self.embedding_dim = embedding_dim
        
    def call(self, inputs):
        # Generate position indices (from 0 to max_len-1)

        pos_enc = np.zeros((self.max_len, self.embedding_dim))

        for j in range(self.max_len):
            pos_enc[j, :] = np.sin(np.arange(self.embedding_dim) / 10000 ** (2 * j / self.max_len))
        
        # Convert to TensorFlow tensor and add a batch dimension
        pos_enc = tf.cast(pos_enc, dtype=tf.float32)
        
        # Add a batch dimension to the positional encodings
        pos_enc = tf.expand_dims(pos_enc, axis=0)
        
        # Add positional embeddings to input tensor (broadcasting)
        return inputs + pos_enc
# ============================
# 步骤5：定义数据集类
# ============================

def create_dataset(df,labels):
    """
    创建 TensorFlow 数据集。
    
    参数：
    df (pd.DataFrame): 数据框，包含基因型特征、亲缘关系向量和标签。
    
    返回：
    tf.data.Dataset: TensorFlow 数据集。
    """
    
    if doubleT:
        labels_c = labels['phenotype_c'].squeeze()  # Shape: (num_individuals,)
        labels_t = labels['phenotype_t'].squeeze()  # Shape: (num_individuals,)

        labels = {
            'phenotype_c': labels_c,
            'phenotype_t': labels_t
        }
    
    if rel:
        # 提取亲缘关系向量
        relationship = np.stack(df['relationship_vector'].values, axis=0)  # 形状: (num_individuals, num_individuals)
        dataset = tf.data.Dataset.from_tensor_slices((relationship, labels))
    elif connect:
        rel_features = np.stack(df[rel_columns].values.astype('float32'), axis=0)
        genotype = np.stack(df[genotype_columns].values.astype('float32'), axis=0)  # 形状: (num_individuals, num_genotype_features)
        dataset = tf.data.Dataset.from_tensor_slices(((genotype,rel_features), labels))
    else:
        # 创建 TensorFlow 数据集
        # 提取基因型特征
        genotype = np.stack(df[genotype_columns].values, axis=0)  # 形状: (num_individuals, num_genotype_features)
        dataset = tf.data.Dataset.from_tensor_slices((genotype, labels))
    
    return dataset

# 创建训练集和验证集
# 创建训练集和验证集
if usingall:
    train_dataset_tf_all = create_dataset(all_train_df, labels_train_all)

train_dataset_tf = create_dataset(train_df, labels_train)
val_dataset_tf = create_dataset(val_df, labels_val)
# 创建数据加载器
#train_loader = train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)
#val_loader = val_dataset_tf.batch(batch_size).prefetch(tf.data.AUTOTUNE)

# ============================
# 步骤7：模型构建（Local CNN + GNN + Laplacian Regularization）
# ============================
def build_model(num_genotype_features):
    """
    构建一个结合 Local CNN 和 GNN 的模型，用于个体级别的预测，并包含拉普拉斯正则化。
    
    参数：
    hp (keras_tuner.HyperParameters): 超参数对象。
    num_features (int): 基因型特征数量。
    
    返回：
    model (tf.keras.Model): 构建好的模型。
    """
    # 输入层
    X_in = layers.Input(shape=(num_genotype_features,), name='genotype_input')  # 节点特征
    
    if connect:
        rel_in = layers.Input(shape=(rel_features,), name='rel_input')  # 节点特征
        x_rel = layers.Dense(128,activation='relu')(rel_in)
        x_rel = layers.BatchNormalization()(x_rel)

     # Local CNN 部分

    if rasieLDim:
        local_features = hp.Int('local_features', min_value=32, max_value=128, step=32)
    elif pureSnp:
        local_features = np.nan
    else:
        local_features = 1

    stride = 1  # 固定步长
    # 超参数：指定网络层数
    #num_layers = hp.Int('num_layers', min_value=1, max_value=5, step=1)

    x = layers.BatchNormalization()(X_in)
    
    
    if np.isnan(local_features):

        if search:

            dense_units = hp.Int('dense_units', min_value=32, max_value=128, step=32)
            
            if doubleT:
                dense_t = hp.Int('dense_units_t', min_value=32, max_value=128, step=32)
                dense_c = hp.Int('dense_units_c', min_value=32, max_value=128, step=32)
        else:
            dense_units = 96

            if doubleT:
                dense_t = 32
                dense_c = 32

        x = layers.Dense(dense_units,activation='relu')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.3)(x)

        x = layers.Dense(64,activation='relu')(x)
        x = layers.BatchNormalization()(x)
        x = layers.Dropout(0.3)(x)

        if doubleT:
            task_t = layers.Dense(dense_t,activation='relu')(x)
            task_c = layers.Dense(dense_c,activation='relu')(x)


    elif local_features > 1:

        kernel_size = hp.Int('kernel_size', min_value=3, max_value=7, step=2)
        residual = x
        x = LocalLinear2D(
                local_features=local_features,
                kernel_size=kernel_size,
                stride=stride,
                bias=True
            )(x)
        x = layers.BatchNormalization()(x)

        x = ReLU()(x)

        if nonlinear:

            x = tf.keras.layers.TimeDistributed(tf.keras.layers.Dense(units=local_features, activation='relu'))(x)
            x = layers.BatchNormalization()(x)

        x = layers.GlobalAveragePooling1D()(x)

        residual = layers.Dense(local_features,activation='relu')(residual)
        residual = layers.BatchNormalization()(residual)

        x = x + residual
        x = ReLU()(x)

    else:
        drate = 0.3
        
        residual  = x

        if newcnn:
                drate = 0.3
                if usingPca:
                    x = tf.expand_dims(x, axis=-1)
                else:
                    if search:
                        kerl = hp.Int('kernel_size_lcl', min_value=3, max_value=6, step=3)
                    else:
                        kerl = 3
                    x = LocalLinear2D(
                        local_features=1,
                        kernel_size=kerl,
                        stride=1,
                        bias=True
                    )(x)
                    x = layers.BatchNormalization()(x)
                    x = ReLU()(x)
                    x = layers.Dropout(drate)(x)
                    x = tf.expand_dims(x, axis=-1)

                residual = x

                if search:
                    con1 = hp.Int('filters1', min_value=32, max_value=96, step=32)
                    con2 = hp.Int('filters2', min_value=32, max_value=96, step=32)
                    ker1 = hp.Int('kernel_size1', min_value=4, max_value=8, step=2)
                    ker2 = hp.Int('kernel_size2', min_value=4, max_value=8, step=2)
                    stride1 = hp.Int('stride1', min_value=4, max_value=8, step=2)
                    stride2 = hp.Int('stride2', min_value=4, max_value=8, step=2)
                    #ker1 = 4
                    #ker2 = 4
                    #stride1 = ker1
                    #stride2 = ker2
                    #strider = stride1*stride2
                    #kerr = ker1*ker2
                else:
                    con1 = 32
                    con2 = 96
                    ker1 = 8
                    ker2 = 4
                    stride1 = 8
                    stride2 = 8
                    #strider = stride1*stride2
                    #kerr = ker1*ker2
                
                x = Conv1D(filters=con1, kernel_size=ker1, strides = stride1,padding='valid', activation='relu')(x)

                x = Conv1D(filters=con2, kernel_size=ker2, strides = stride2,padding='valid', activation='relu')(x)

                #residual = Conv1D(filters=con2, kernel_size=kerr,strides = strider, padding='valid', activation='relu')(residual)

                #x = x + residual
                
                #x = ReLU()(x)

                x = layers.MaxPooling1D(pool_size=5)(x)

                x = layers.Flatten()(x)

                if search:

                    dense_units = hp.Int('dense_units', min_value=64, max_value=128, step=32)
                    
                else:
                    dense_units = 96

                if connect:
                    x= layers.Concatenate()([x, x_rel])

                if doubleT:
                    task_t = layers.Dense(dense_units,activation='relu')(x)

                    task_c = layers.Dense(dense_units,activation='relu')(x)
                else:
                    x = layers.Dense(dense_units,activation='relu')(x)
                
        else:
            x = LocalLinear2D(
                local_features=local_features,
                kernel_size=7,
                stride=stride,
                bias=True
            )(x)

            x = layers.BatchNormalization()(x)
            x = ReLU()(x)
            x = layers.Dropout(drate)(x)

            fold = calfold(genotype_data.shape[1],stride=1)

        if attention:
            drate = 0.3
            x = x + residual
            x = ReLU()(x)
            residual = x
            x = LocalLinear2D(
                local_features=local_features,
                kernel_size=5,
                stride=5,
                bias=True
            )(x)

            x = layers.BatchNormalization()(x)
            x = ReLU()(x)
            x = layers.Dropout(drate)(x)

            fold = calfold(in_features=fold,stride=5)

            residual = layers.Dense(fold,activation='relu')(residual)
            residual = layers.BatchNormalization()(residual)

            x = x + residual
            x = ReLU()(x)
            residual = x

            x = LocalLinear2D(
                local_features=local_features,
                kernel_size=5,
                stride=5,
                bias=True
            )(x)

            x = layers.BatchNormalization()(x)
            x = ReLU()(x)
            x = layers.Dropout(drate)(x)

            fold = calfold(in_features=fold,stride=5)

            residual = layers.Dense(fold,activation='relu')(residual)
            residual = layers.BatchNormalization()(residual)

            x = x + residual
            x = ReLU()(x)

            x = LocalLinear2D(
                local_features=12,
                kernel_size=3,
                stride=3,
                bias=True
            )(x)

            x = layers.BatchNormalization()(x)
            x = ReLU()(x)
            x = layers.Dropout(drate)(x)

            fold = calfold(in_features=fold,stride=3)

            if uspos:
                    x = PositionalEncoding(max_len=fold,embedding_dim=12)(x)

            residual = x

            if search:
                num_attention_heads = hp.Int('num_attention_heads', min_value=2, max_value=8, step=2)
                key_dim = hp.Int('key_dim', min_value=2, max_value=8, step=2)
                dropout=hp.Float('attention_dropout', min_value=0.1, max_value=0.5, step=0.1)
            else:
                num_attention_heads = 2
                key_dim = 8
                dropout = 0.3

                #num_attention_heads = 4
                #key_dim = 2
                #dropout = 0.2

            attention_output = MultiHeadAttention(
                                    num_heads=num_attention_heads,
                                    key_dim=key_dim,
                                    dropout = dropout
                                    )(x, x)

            x = layers.BatchNormalization()(attention_output + residual)

            # Optionally, add a feed-forward layer
            x = layers.GlobalAveragePooling1D()(x)


            
            #x = tf.reduce_sum(attention_output, axis=1)  # Aggregate into a global feature

        elif not newcnn:
            drate = 0.3
            x = x + residual
            x = ReLU()(x)
            residual = x
            x = LocalLinear2D(
                local_features=local_features,
                kernel_size=3,
                stride=stride,
                bias=True
            )(x)
            x = layers.BatchNormalization()(x)
            x = ReLU()(x)
            x = layers.Dropout(drate)(x)

            x = x + residual
            x = ReLU()(x)
            
            if addDense:

                x = layers.Dense(64,activation='relu')(x)

    genotype_features = x

    #z = layers.LayerNormalization()(genotype_features)
    z = genotype_features
    # 输出层：针对每个节点进行预测

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


        if not softm:

            output_c = layers.Dense(1, activation='sigmoid', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': tf.keras.losses.BinaryCrossentropy(),
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy', tf.keras.metrics.AUC()],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

        else:
            output_c = layers.Dense(5, activation='softmax', name='phenotype_c')(task_c)

            loss = {
                'phenotype_c': 'sparse_categorical_crossentropy',
                'phenotype_t': tf.keras.losses.BinaryCrossentropy(),
            }

            metrics = {
                'phenotype_c': ['accuracy'],
                'phenotype_t': ['accuracy', tf.keras.metrics.AUC()],
            }

        
        output = [output_c,output_t]
        
    
    # 构建模型
    if connect:
        model = models.Model(inputs=[X_in, rel_in], outputs=output)
    else:
        model = models.Model(inputs=X_in, outputs=output)

    for layer in model.layers:
        if hasattr(layer, 'kernel_regularizer'):
            layer.kernel_regularizer = l2(0.01)
        if hasattr(layer, 'bias_regularizer'):
            layer.bias_regularizer = l2(0.01)

    
    # 编译模型
    if search:
        #learning_rate = hp.Float('learning_rate', min_value=1e-7, max_value=1e-5, sampling='log')
        learning_rate = 1e-5
    else:
        #learning_rate = 0.00408
        #learning_rate = 0.005307
        learning_rate = 1e-5#fc
        #learning_rate = 4.44853401907757e-07#attention

    optimizer = optimizers.Adam(learning_rate=learning_rate)
    
    model.compile(
        optimizer=optimizer,
        loss=loss,
        metrics=metrics
    )
    
    return model

# ============================
# 步骤10：超参数调优与训练
# ============================

# 定义早停回调
stop_early = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5)
np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
if search:
    def model_builder(hp):
        if not rel:
            nF = len(genotype_columns)
        else:
            nF = relationship_vectors.shape[1]

        if connect:
            return build_model(hp,num_genotype_features=nF,rel_features=mergedDF.shape[0])
        else:
            return build_model(hp,num_genotype_features=nF)
#kt.GridSearch kt.Hyperband
    class MyTuner(kt.Hyperband):
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


    tuner = MyTuner(
            hypermodel=model_builder,
            objective='val_loss',
            max_epochs=100,
            factor=3,
            directory='my_dir',
            project_name='multimodal_model_tuning'
        )
    tuner.search()

    # 创建调优器
    #tuner = kt.Hyperband(
    #    model_builder,
    #    objective='val_loss',
    #    max_epochs=50,
    #    factor=3,
    #    directory='my_dir',
    #    project_name='multimodal_model_tuning'
    #)

    #tuner = kt.GridSearch(
    #    hypermodel=model_builder,
    #    objective='val_loss',
    #    max_trials=100,  # Total number of combinations to try
    #     directory='grid_search_dir',
    #    project_name='grid_search_example'
    #)

    #tuner = kt.BayesianOptimization(
    #    model_builder,
    #    objective='val_loss',
    #    max_trials=50,  # Number of trials
    #    num_initial_points=20,
    #    directory='my_dir',
    #    project_name='bayesian_tuning'
    #)

    #tuner.search(
    #    train_loader,
    #    epochs=50,
    #    validation_data=val_loader,
    #    callbacks=[stop_early]
    #)

    # 获取最佳超参数
    best_hps = tuner.get_best_hyperparameters(num_trials=1)[0]
    print("最佳超参数：")
    for param in best_hps.values:
        print(f"{param}: {best_hps.get(param)}")

    # 构建并训练最佳模型
    model = tuner.hypermodel.build(best_hps)
else:
    if not rel:
        nF = len(genotype_columns)
    else:
        nF = relationship_vectors.shape[1]
    
    if connect:
        model = build_model(num_genotype_features=nF,rel_features=mergedDF.shape[0])
    else:
        model = build_model(num_genotype_features=nF)

#model.save('best_model.h5')

np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
if search:
    batch_size = best_hps.get('batch_size')
    train_loader = train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)
    val_loader = val_dataset_tf.batch(batch_size).prefetch(tf.data.AUTOTUNE)
    history = model.fit(
        train_loader,
        #epochs=best_hps.get("tuner/epochs"),
        epochs=100,
        initial_epoch = 0,
        validation_data=val_loader,
        callbacks=[stop_early]
    )
else:
    batch_size = 30
    train_loader = train_dataset_tf.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)
    val_loader = val_dataset_tf.batch(batch_size).prefetch(tf.data.AUTOTUNE)
    epochs = 100
    if usingall:
        epochs = 6
        train_loader_all = train_dataset_tf_all.shuffle(buffer_size=1024).batch(batch_size).prefetch(tf.data.AUTOTUNE)
        history = model.fit(
            train_loader_all,
            #validation_data=val_loader,
            #callbacks=[stop_early],
            epochs=epochs, 
            initial_epoch=0
        )
    else:
        history = model.fit(
                    train_loader,
                    validation_data=val_loader,
                    callbacks=[stop_early],
                    epochs=epochs, 
                    initial_epoch=0
                )

# 直接获取训练和验证损失值
if not usingall:
    train_loss = history.history['loss']
    val_loss = history.history['val_loss']

    # 绘制损失曲线
    plt.plot(history.epoch, train_loss, label='Train Loss')
    plt.plot(history.epoch, val_loss, label='Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.legend()
    # 保存图像为 PDF
    plt.savefig('loss_curve.pdf', format='pdf')
    plt.show()

# ============================
# 步骤11：在测试数据上进行预测
# ============================

# 加载测试集 ID

# 筛选测试集数据
test_df = mergedDF[mergedDF['ID'].isin(test_ids)].reset_index(drop=True)

# 准备基因型数据并标准化（如果启用了标准化）
if rel:
    test_relationship = np.stack(test_df['relationship_vector'].values, axis=0)
else:
    test_genotype = test_df[genotype_columns].values.astype(np.float32)
    if connect:
        test_rel = test_df[rel_columns].values.astype(np.float32)

def create_test_dataset(features):
    """
    创建测试 TensorFlow 数据集。
    
    参数：
    genotype (np.ndarray): 基因型特征数组。
    relationship (np.ndarray): 亲缘关系向量数组。
    
    返回：
    tf.data.Dataset: TensorFlow 数据集。
    """
    dataset = tf.data.Dataset.from_tensor_slices((features,))
    return dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

# 创建测试数据集
if rel:
    test_dataset_tf = create_test_dataset(test_relationship)
elif connect:
    test_dataset_tf = create_test_dataset((test_genotype,test_rel))
else:
    test_dataset_tf = create_test_dataset(test_genotype)

predictions = model.predict(test_dataset_tf, verbose=1)
# 获取预测对应的个体ID
predicted_ids = test_df['ID'].tolist()

if doubleT:
    if not softm:
        predictions_c = predictions[0].squeeze(-1)  # 形状: (num_test_individuals,)
    else:
        predictions_c = np.dot(predictions[0], np.arange(5)) 
        
    predictions_t = predictions[1].squeeze(-1)  # 形状: (num_test_individuals,)
    # 构建最终的预测DataFrame
    predictions_df = pd.DataFrame({
        "ID": predicted_ids,
        "Prediction_c": predictions_c,
        "Prediction_t": predictions_t
    })

    # 保存预测结果
    output_file = "predictions_c_t.csv"
else:
    predictions = predictions.squeeze(-1)

    # 构建最终的预测DataFrame
    predictions_df = pd.DataFrame({
        "ID": predicted_ids,
        "Prediction": predictions
    })

    # 保存预测结果
    output_file = "predictions_t.csv" if Th else "predictions_c.csv"

predictions_df.to_csv(output_file, index=False)
print(f"预测结果已保存到 {output_file}。")