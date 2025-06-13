#!/usr/bin/env python
# coding: utf-8



# ===============================# SECTION 0: 0 Instructions# ===============================


# ===============================# SECTION 0.1: 0.1 Feature Selection: 22 features# ===============================
# select 22 features:   
# selected_features = [
#     'num_seg_p2',
#     'num_variant_window',
#     'exon_density', 
#     'mean_introg_anc', 
#     'exon_window', 
#     'introg_anc_window', 
#     'divergence_p3_p1', 
#     'divergence_p3_p2', 
#     'watterson_theta_p3', 
#     'windowed_tajima_d_p3', 
#     'D', 
#     'Het', 
#     'Q95', 
#     'U0', 
#     'U20', 
#     'U50', 
#     'U80', 
#     'num_seg_p1', 
#     'num_seg_p3', 
#     'num_private_seg_p1', 
#     'num_private_seg_p2', 
#     'num_private_seg_p3'
# ]
# 
# Feature importances (sorted): 
# | Feature               | Importance   |
# |-----------------------|--------------| 
# |num_seg_p2| 0.179797947406768| 
# |num_private_seg_p2| 0.14682884514331818|
# |exon_window| 0.10894051194190979|
# |num_variant_window| 0.09169841557741165|  
# |mean_introg_anc| 0.0881299301981926|  
# |watterson_theta_p3| 0.06600210070610046|  
# |divergence_p3_p1| 0.0538877509534359|  
# |exon_density| 0.027985380962491035|  
# |num_private_seg_p3| 0.02782861702144146|
# |divergence_p3_p2| 0.025785868987441063|
# |num_seg_p1| 0.023056915029883385|
# |recrate_window| 0.016651177778840065|
# |mean_recrate| 0.015112540684640408|
# |df_p3_p1| 0.0138973668217659|
# |U50| 0.008871671743690968|
# |U0| 0.00883196946233511|
# |Het| 0.00876428559422493|
# |num_seg_p3| 0.008294460363686085|
# |garud_h1| 0.007584222126752138|
# |RD| 0.007552091032266617|
# |D| 0.007459886837750673|
# |hap_diversity_p3| 0.007423382252454758|
# |df_p3_p2| 0.006998920347541571|
# |Q95| 0.006585646886378527|
# |introg_anc_window| 0.006539533846080303|
# |U20| 0.006470245774835348|
# |fD| 0.006391296163201332|
# |num_private_seg_p1| 0.0059153544716537|
# |windowed_tajima_d_p3| 0.00569456210359931|
# |U80| 0.005019065923988819|
# |garud_h2_h1| 0.0|
# |garud_h12| 0.0|
# 
# Feature importances (sorted):
# | Feature               | Importance            | CNN Model | Other Models |
# |-----------------------|-----------------------|-----------|--------------|
# | num_seg_p2            | 0.179797947406768     | ✔️        |              |
# | num_private_seg_p2    | 0.14682884514331818   | ✔️        | ✔️           |
# | exon_window           | 0.10894051194190979   | ✔️        | ✔️           |
# | num_variant_window    | 0.09169841557741165   | ✔️        |              |
# | mean_introg_anc       | 0.0881299301981926    | ✔️        | ✔️           |
# | watterson_theta_p3    | 0.06600210070610046   | ✔️        | ✔️           |
# | divergence_p3_p1      | 0.0538877509534359    | ✔️        | ✔️           |
# | exon_density          | 0.027985380962491035  | ✔️        | ✔️           |
# | num_private_seg_p3    | 0.02782861702144146   | ✔️        | ✔️           |
# | divergence_p3_p2      | 0.025785868987441063  | ✔️        | ✔️           |
# | num_seg_p1            | 0.023056915029883385  | ✔️        | ✔️           |
# | recrate_window        | 0.016651177778840065  |           |              |
# | mean_recrate          | 0.015112540684640408  |           |              |
# | df_p3_p1              | 0.0138973668217659    |           |              |
# | U50                   | 0.008871671743690968  | ✔️        | ✔️           |
# | U0                    | 0.00883196946233511   | ✔️        | ✔️           |
# | Het                   | 0.00876428559422493   | ✔️        | ✔️           |
# | num_seg_p3            | 0.008294460363686085  | ✔️        | ✔️           |
# | garud_h1              | 0.007584222126752138  |           |              |
# | RD                    | 0.007552091032266617  |           |              |
# | D                     | 0.007459886837750673  | ✔️        | ✔️           |
# | hap_diversity_p3      | 0.007423382252454758  |           |              |
# | df_p3_p2              | 0.006998920347541571  |           |              |
# | Q95                   | 0.006585646886378527  | ✔️        | ✔️           |
# | introg_anc_window     | 0.006539533846080303  | ✔️        | ✔️           |
# | U20                   | 0.006470245774835348  | ✔️        | ✔️           |
# | fD                    | 0.006391296163201332  |           |              |
# | num_private_seg_p1    | 0.0059153544716537    | ✔️        | ✔️           |
# | windowed_tajima_d_p3  | 0.00569456210359931   | ✔️        | ✔️           |
# | U80                   | 0.005019065923988819  | ✔️        | ✔️           |
# | garud_h2_h1           | 0.0                   |           |              |
# | garud_h12             | 0.0                   |           |              |
# 
# 


# ===============================# SECTION 0.2: 0.2 Strategy# ===============================
# 1) split all data into test and training data
# 2) stratify test data into: exon_density >200 (all); exon_density >400; exon_density >600; exon_density >800
# 3) normalize training data and test data using the (mean,std) of training data
# 4) normalize real data using the (mean,std) of real data



# ===============================# SECTION 1: 1 CNN (Sigmoid) with gradient clipping# ===============================


# ===============================# SECTION 1.1: 1.1 A 3-layer CNN (Sigmoid) with gradient clipping# ===============================
# 1DConvolutionalArchitecture：
# * Conv1d(32 3x1 kernels) + ReLU + MaxPool1d(2x1)
# * Conv1d(64 3x1 kernels) + ReLU + MaxPool1d(2x1)
# * Conv1d(32 3x1 kernels) + ReLU + MaxPool1d(2x1)
# * Dense Layer (512 units) + ReLU
# * Dense Layer (number of classes units) + Sigmoid (assuming classification problem)

# In[1]:


import torch 
from torch.utils.data import Dataset, DataLoader, Subset
import pandas as pd
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm
import pickle
import seaborn as sns
from numpy.lib.format import read_magic, read_array_header_1_0, read_array_header_2_0
import time
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import io
import contextlib
from torchvision import transforms

from itertools import combinations
from sklearn.metrics import confusion_matrix
import openpyxl
from matplotlib.colors import LinearSegmentedColormap


# In[2]:


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

set_seed(1)

class WindowDataset(Dataset):
    def __init__(self, csv_file, mean_std_dir, transform=None):
        self.csv_file = csv_file
        self.transform = transform
        self.column_names = [
            'num_seg_p2','num_variant_window','exon_density', 'mean_introg_anc', 'exon_window', 'introg_anc_window', 'divergence_p3_p1', 'divergence_p3_p2', 'watterson_theta_p3', 'windowed_tajima_d_p3', 'D', 
            'Het', 'Q95', 'U0', 'U20', 'U50', 'U80', 'num_seg_p1', 'num_seg_p3', 'num_private_seg_p1', 'num_private_seg_p2', 'num_private_seg_p3'
        ]
        
        # Load entire dataset to fetch column indices
        all_data = np.genfromtxt(csv_file, delimiter=',', names=True, dtype=None)
        self.data = np.column_stack([all_data[name] for name in self.column_names])
        # Check if 'dominance' column exists
        if 'dominance' in all_data.dtype.names:
            self.labels = all_data['dominance'].astype(int)
            self.has_labels = True
        else:
            self.has_labels = False

        # Load means and stds for the selected columns only
        mean_std_data = pd.read_csv(mean_std_dir)
        filtered_mean_std_data = mean_std_data.set_index('Unnamed: 0').reindex(self.column_names).reset_index()
        self.means = filtered_mean_std_data['Mean'].to_numpy()
        self.stds = filtered_mean_std_data['Std'].to_numpy()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        x = self.data[idx].astype(float)
        x = (x - self.means) / self.stds  # Normalize x
        x = np.stack([x, x, x], axis=0)  # Duplicate the single channel to create 3 channels
        x = torch.tensor(x, dtype=torch.float32)
        if self.transform:
            x = self.transform(x)
        
        if self.has_labels:
            y = self.labels[idx]
            return x, torch.tensor(y, dtype=torch.long)
        else:
            return x

transform = transforms.Compose([
    transforms.Lambda(lambda x: x.unsqueeze(0))
])

# Emprical data
data_dir = '/home/zhulx/lab/summary_statistics/data/final/empirical data'
def create_datasets(data_dir):
    datasets = {}
    for file in os.listdir(data_dir): # Loop through all files in the directory
        if file.endswith('_all.csv'):
            base_name = file.replace('_empirical-stats_all.csv', '')
            meanstd_file = f'{base_name}_empirical-stats_all_meanstd.csv'
            if meanstd_file in os.listdir(data_dir):# Check if the meanstd file exists
                file_path = os.path.join(data_dir, file) # Full paths for files
                meanstd_path = os.path.join(data_dir, meanstd_file)
                datasets[base_name] = WindowDataset(file_path, meanstd_path) # Create WindowDataset object (assuming WindowDataset is a predefined class or function)
    return datasets

start_time = time.time()
train_dataset = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_train.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger200 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_200.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger400 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_400.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger600 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_600.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger800 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_800.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
data_dir = '/home/zhulx/lab/summary_statistics/data/final/empirical data'
CEU_dataset_meanstdReal = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_CEU_empirical-stats_all.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_CEU_empirical-stats_all_meanstd.csv')#,transform=transform)
real_datasets = create_datasets(data_dir) # print(datasets.keys())
end_time = time.time()
print(f"Loading data took {end_time - start_time} seconds.") 
print(f"Training Dataset size: {len(train_dataset)}; additive size: {len([i for i, label in enumerate(train_dataset.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(train_dataset.labels) if label == 1])}")
print(f"Testing>200 Dataset size: {len(test_dataset_larger200)}; additive size: {len([i for i, label in enumerate(test_dataset_larger200.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger200.labels) if label == 1])}") 
print(f"Testing>400 Dataset size: {len(test_dataset_larger400)}; additive size: {len([i for i, label in enumerate(test_dataset_larger400.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger400.labels) if label == 1])}") 
print(f"Testing>600 Dataset size: {len(test_dataset_larger600)}; additive size: {len([i for i, label in enumerate(test_dataset_larger600.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger600.labels) if label == 1])}") 
print(f"Testing>800 Dataset size: {len(test_dataset_larger800)}; additive size: {len([i for i, label in enumerate(test_dataset_larger800.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger800.labels) if label == 1])}") 
print(f"CEU Dataset size: {len(CEU_dataset_meanstdReal)}")
print(f"TSI Dataset size: {len(real_datasets['hg19_TSI'])}")
print(f"GBR Dataset size: {len(real_datasets['hg19_GBR'])}")
print(f"IBS Dataset size: {len(real_datasets['hg19_IBS'])}")
print(f"CHB Dataset size: {len(real_datasets['hg19_CHB'])}")
print(f"CHS Dataset size: {len(real_datasets['hg19_CHS'])}")
print(f"JPT Dataset size: {len(real_datasets['hg19_JPT'])}")
print(f"FIN Dataset size: {len(real_datasets['hg19_FIN'])}")

train_loader = DataLoader(train_dataset, batch_size=16, shuffle=True)
test_loader_larger200 = DataLoader(test_dataset_larger200, batch_size=16, shuffle=False)
test_loader_larger400 = DataLoader(test_dataset_larger400, batch_size=16, shuffle=False)
test_loader_larger600 = DataLoader(test_dataset_larger600, batch_size=16, shuffle=False)
test_loader_larger800 = DataLoader(test_dataset_larger800, batch_size=16, shuffle=False)
CEU_loader = DataLoader(CEU_dataset_meanstdReal, batch_size=16, shuffle=False)
real_loaders = {}
for i in real_datasets.keys():
    data_loader = DataLoader(real_datasets[i], batch_size=16, shuffle=False)
    real_loaders[i] = data_loader
# torch.isnan(inputs).sum() # torch.isnan(labels).sum()


# In[40]:


class CustomCNN(nn.Module):
    def __init__(self, num_classes=2):
        super(CustomCNN, self).__init__()
        # since the image is 1D (1*11)
        self.conv1 = nn.Conv1d(3, 32, kernel_size=3, stride=1, padding=1) 
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, stride=1, padding=1)
        self.conv3 = nn.Conv1d(64, 32, kernel_size=3, stride=1, padding=1)
        self.pool = nn.MaxPool1d(kernel_size=2, stride=2, padding=0)

        # Since we are using 1D convolutions, update the input shape accordingly
        self.num_features = self._get_conv_output(22)  # 20 features
        
        self.fc1 = nn.Linear(self.num_features, 512)
        self.fc2 = nn.Linear(512, num_classes) 
    
    def _get_conv_output(self, input_size):
        with torch.no_grad():
            # We assume that the input size is the size of the sequence.
            x = torch.zeros(1, 3, input_size)  # (batch_size, channels, sequence_length) #x = torch.zeros(1, *shape)
            x = self.pool(F.relu(self.conv1(x)))
            x = self.pool(F.relu(self.conv2(x)))
            x = self.pool(F.relu(self.conv3(x)))
            return int(np.prod(x.size()[1:]))  

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = self.pool(F.relu(self.conv3(x)))
        x = x.view(x.size(0), -1)  
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return torch.sigmoid(x) # Use sigmoid for binary classification; # F.log_softmax(x, dim=1) # Using log_softmax for classification; # x

# Initialize weights
def weights_init(m):
    if isinstance(m, nn.Conv1d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight)

model = CustomCNN(num_classes=2) # Initialize the model
model.apply(weights_init)
# device = torch.device('cpu')
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') # Select the device: checks if a GPU (CUDA) is available. If it is, it selects the CUDA device (i.e., GPU), otherwise, it selects the CPU.
model = model.to(device) # Move the model to the device(GPU (CUDA))

criterion = nn.CrossEntropyLoss()
optimizer = Adam(model.parameters(), lr=1e-6) # Reduced learning rate to avoid gradient exploding # optimizer = Adam(model.parameters(), lr=0.0001)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.1)

model_save_path = '/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features'
os.makedirs(model_save_path, exist_ok=True)

def train_model(model, train_loader,test_loader_larger200,test_loader_larger400,test_loader_larger600,test_loader_larger800,real_loaders,criterion, optimizer, scheduler, num_epochs, device,model_save_path):
    best_accuracy = 0.0
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for inputs, labels in tqdm(train_loader, desc=f'Epoch {epoch+1}/{num_epochs}', unit='batch'):
            inputs, labels = inputs.to(device), labels.to(device)

            # if torch.isnan(inputs).sum() > 0 or torch.isnan(labels).sum() > 0: # Check for NaN in inputs
            #     print("NaN detected in inputs or labels")

            optimizer.zero_grad()
            outputs = model(inputs)

            # if torch.isnan(outputs).sum() > 0:
            #     print("NaN detected in outputs")
                # print(outputs)
            
            loss = criterion(outputs, labels)
            loss.backward()

            # Gradient clipping: to avoid gradient exploding
            with contextlib.redirect_stdout(io.StringIO()): # Suppress the output of clip_grad_norm_ using contextlib.redirect_stdout
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)

            optimizer.step()
            running_loss += loss.item() * inputs.size(0)
        epoch_loss = running_loss / len(train_loader.dataset)
        
        model.eval()
        correct_train = 0
        total_train = 0
        with torch.no_grad():
            for inputs, labels in train_loader:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_train += labels.size(0)
                correct_train += (predicted == labels).sum().item()
        train_accuracy = correct_train / total_train

        correct_test = 0
        total_test = 0
        all_labels_200 = [] # true labels
        all_outputs_200 = []
        predictions_200 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger200:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_200.extend(labels.cpu().numpy())
                all_outputs_200.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_200.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_200 = correct_test / total_test

        correct_test = 0
        total_test = 0
        all_labels_400 = [] # true labels
        all_outputs_400 = []
        predictions_400 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger400:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_400.extend(labels.cpu().numpy())
                all_outputs_400.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_400.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_400 = correct_test / total_test

        correct_test = 0
        total_test = 0
        all_labels_600 = [] # true labels
        all_outputs_600 = []
        predictions_600 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger600:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_600.extend(labels.cpu().numpy())
                all_outputs_600.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_600.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_600 = correct_test / total_test

        correct_test = 0
        total_test = 0
        all_labels_800 = [] # true labels
        all_outputs_800 = []
        predictions_800 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger800:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_800.extend(labels.cpu().numpy())
                all_outputs_800.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_800.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_800 = correct_test / total_test
        print('test data is done.')
        
        # CEU_outputs_meanstdReal = []
        # CEU_predictions_meanstdReal = []
        # with torch.no_grad():
        #     for inputs in CEU_loader:
        #         inputs = inputs.to(device)
        #         outputs = model(inputs)
        #         _, predicted = torch.max(outputs.data, 1)
        #         CEU_outputs_meanstdReal.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
        #         CEU_predictions_meanstdReal.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        real_outputs_meanstdReal = {}
        real_predictions_meanstdReal = {}
        with torch.no_grad():
            for dataset_name in real_loaders.keys():
                loader = real_loaders[dataset_name]
                dataset_outputs = []
                dataset_predictions = []
                for inputs in loader:
                    inputs = inputs.to(device)
                    outputs = model(inputs)
                    _, predicted = torch.max(outputs.data, 1)
                    dataset_outputs.extend(outputs.cpu().numpy()[:, 1])  # 假设是二元分类
                    dataset_predictions.extend(predicted.cpu().numpy())
                real_outputs_meanstdReal[dataset_name] = dataset_outputs
                real_predictions_meanstdReal[dataset_name] = dataset_predictions


        if test_accuracy_200 > best_accuracy:
            best_test_accuracy_200 = test_accuracy_200 # Update best accuracy
            best_test_predictions_200 = predictions_200  # Update best predictions
            test_outputs_best_200 = all_outputs_200
            best_test_accuracy_400 = test_accuracy_400 # Update best accuracy
            best_test_predictions_400 = predictions_400  # Update best predictions
            test_outputs_best_400 = all_outputs_400
            best_test_accuracy_600 = test_accuracy_600 # Update best accuracy
            best_test_predictions_600 = predictions_600  # Update best predictions
            test_outputs_best_600 = all_outputs_600
            best_test_accuracy_800 = test_accuracy_800 # Update best accuracy
            best_test_predictions_800 = predictions_800  # Update best predictions
            test_outputs_best_800 = all_outputs_800
            # CEU_outputs_meanstdReal_best = CEU_outputs_meanstdReal
            # CEU_predictions_meanstdReal_best = CEU_predictions_meanstdReal
            real_outputs_meanstdReal_best = real_outputs_meanstdReal
            real_predictions_meanstdReal_best = real_predictions_meanstdReal
            torch.save(model.state_dict(), os.path.join(model_save_path, 'best_model.pth'))
        
        scheduler.step()

        print(f'Epoch {epoch+1}/{num_epochs}, Train Accuracy: {train_accuracy:.4f}, Test Accuracy (>=200): {test_accuracy_200:.4f}')
        print(f'Test Accuracy (>=400): {test_accuracy_400:.4f}, Test Accuracy (>=600): {test_accuracy_600:.4f}, Test Accuracy (>=800): {test_accuracy_800:.4f}')

    return best_test_accuracy_200,all_labels_200,test_outputs_best_200,best_test_predictions_200,best_test_accuracy_400,all_labels_400,test_outputs_best_400,best_test_predictions_400,best_test_accuracy_600,all_labels_600,test_outputs_best_600,best_test_predictions_600,best_test_accuracy_800,all_labels_800,test_outputs_best_800,best_test_predictions_800,real_outputs_meanstdReal_best,real_predictions_meanstdReal_best

num_epochs = 20
Results = train_model(model, train_loader,test_loader_larger200,test_loader_larger400,test_loader_larger600,test_loader_larger800,real_loaders,criterion, optimizer, scheduler, num_epochs, device, model_save_path)

# Save predictions
best_test_accuracy_200,all_labels_200,test_outputs_best_200,best_test_predictions_200,best_test_accuracy_400,all_labels_400,test_outputs_best_400,best_test_predictions_400,best_test_accuracy_600,all_labels_600,test_outputs_best_600,best_test_predictions_600,best_test_accuracy_800,all_labels_800,test_outputs_best_800,best_test_predictions_800,real_outputs_meanstdReal_best,real_predictions_meanstdReal_best = Results
print(f'Best accuracy (test>=200) = {best_test_accuracy_200}')
print(f'Best accuracy (test>=400) = {best_test_accuracy_400}')
print(f'Best accuracy (test>=600) = {best_test_accuracy_600}')
print(f'Best accuracy (test>=800) = {best_test_accuracy_800}')

test_df_200 = pd.DataFrame({'True_Labels':all_labels_200,'Predicted_Labels':best_test_predictions_200,'Prediction_Prob': test_outputs_best_200})
test_df_200.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger200.csv'), index=False)
test_df_400 = pd.DataFrame({'True_Labels':all_labels_400,'Predicted_Labels':best_test_predictions_400,'Prediction_Prob': test_outputs_best_400})
test_df_400.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger400.csv'), index=False)
test_df_600 = pd.DataFrame({'True_Labels':all_labels_600,'Predicted_Labels':best_test_predictions_600,'Prediction_Prob': test_outputs_best_600})
test_df_600.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger600.csv'), index=False)
test_df_800 = pd.DataFrame({'True_Labels':all_labels_800,'Predicted_Labels':best_test_predictions_800,'Prediction_Prob': test_outputs_best_800})
test_df_800.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger800.csv'), index=False)
# real_df = pd.DataFrame({'Predicted_Labels_meanstdReal': CEU_predictions_meanstdReal_best,'Prediction_Prob_meanstdReal': CEU_outputs_meanstdReal_best})
# real_df.to_csv(os.path.join(model_save_path, 'CEU_prediction_result.csv'), index=False)
for dataset_name in real_outputs_meanstdReal_best.keys():
    real_df = pd.DataFrame({
        'Predicted_Labels_meanstdReal': real_predictions_meanstdReal_best[dataset_name],
        'Prediction_Prob_meanstdReal': real_outputs_meanstdReal_best[dataset_name]
    })
    save_path = os.path.join(model_save_path, f'{dataset_name}_prediction_result.csv')
    real_df.to_csv(save_path, index=False)


# In[41]:


roc_data = []
for threshold in [200, 400, 600, 800]:
        file_path = f'/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/test_prediction_result_larger{threshold}.csv'
        sample_df = pd.read_csv(file_path)
        test_labels_best = sample_df['True_Labels']
        test_outputs_best = sample_df['Prediction_Prob']
        auc_roc = roc_auc_score(test_labels_best, test_outputs_best)
        fpr, tpr, _ = roc_curve(test_labels_best, test_outputs_best)
        print(f'auc (test>{threshold}) =',auc_roc)
        roc_data.append((fpr, tpr, auc_roc, f'test>{threshold}'))

plt.figure(figsize=(10, 8))
for fpr, tpr, auc_roc, label in roc_data:
    plt.plot(fpr, tpr, label=f'{label} (AUC = {auc_roc:.2f})')

# plt.plot([0, 1], [0, 1], 'k--', label='No Skill')  # add a line for no skill classifier
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curves for Different Thresholds')
plt.legend(loc='best')
plt.show()
plt.savefig(os.path.join(model_save_path,'roc_curve_cnn.png'))




# ===============================# SECTION 1.2: 1.2 A 2-layer CNN (Sigmoid) with gradient clipping# ===============================
# 1DConvolutionalArchitecture：
# * Conv1d(32 3x1 kernels) + ReLU + MaxPool1d(2x1)
# * Conv1d(64 3x1 kernels) + ReLU + MaxPool1d(2x1)
# * Conv1d(32 3x1 kernels) + ReLU + MaxPool1d(2x1)
# * Dense Layer (512 units) + ReLU
# * Dense Layer (number of classes units) + Sigmoid (assuming classification problem)

# In[ ]:


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

set_seed(1)

class WindowDataset(Dataset):
    def __init__(self, csv_file, mean_std_dir, transform=None):
        self.csv_file = csv_file
        self.transform = transform
        self.column_names = [
            'num_seg_p2','num_variant_window','exon_density', 'mean_introg_anc', 'exon_window', 'introg_anc_window', 'divergence_p3_p1', 'divergence_p3_p2', 'watterson_theta_p3', 'windowed_tajima_d_p3', 'D', 
            'Het', 'Q95', 'U0', 'U20', 'U50', 'U80', 'num_seg_p1', 'num_seg_p3', 'num_private_seg_p1', 'num_private_seg_p2', 'num_private_seg_p3'
        ]
        
        # Load entire dataset to fetch column indices
        all_data = np.genfromtxt(csv_file, delimiter=',', names=True, dtype=None)
        self.data = np.column_stack([all_data[name] for name in self.column_names])
        # Check if 'dominance' column exists
        if 'dominance' in all_data.dtype.names:
            self.labels = all_data['dominance'].astype(int)
            self.has_labels = True
        else:
            self.has_labels = False

        # Load means and stds for the selected columns only
        mean_std_data = pd.read_csv(mean_std_dir)
        filtered_mean_std_data = mean_std_data.set_index('Unnamed: 0').reindex(self.column_names).reset_index()
        self.means = filtered_mean_std_data['Mean'].to_numpy()
        self.stds = filtered_mean_std_data['Std'].to_numpy()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        x = self.data[idx].astype(float)
        x = (x - self.means) / self.stds  # Normalize x
        x = np.stack([x, x, x], axis=0)  # Duplicate the single channel to create 3 channels
        x = torch.tensor(x, dtype=torch.float32)
        if self.transform:
            x = self.transform(x)
        
        if self.has_labels:
            y = self.labels[idx]
            return x, torch.tensor(y, dtype=torch.long)
        else:
            return x

transform = transforms.Compose([
    transforms.Lambda(lambda x: x.unsqueeze(0))
])

# Emprical data
data_dir = '/home/zhulx/lab/summary_statistics/data/final/empirical data'
def create_datasets(data_dir):
    datasets = {}
    for file in os.listdir(data_dir): # Loop through all files in the directory
        if file.endswith('_all.csv'):
            base_name = file.replace('_empirical-stats_all.csv', '')
            meanstd_file = f'{base_name}_empirical-stats_all_meanstd.csv'
            if meanstd_file in os.listdir(data_dir):# Check if the meanstd file exists
                file_path = os.path.join(data_dir, file) # Full paths for files
                meanstd_path = os.path.join(data_dir, meanstd_file)
                datasets[base_name] = WindowDataset(file_path, meanstd_path) # Create WindowDataset object (assuming WindowDataset is a predefined class or function)
    return datasets

start_time = time.time()
train_dataset = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_train.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger200 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_200.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger400 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_400.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger600 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_600.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
test_dataset_larger800 = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/trainAll_test_larger_800.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/summary_statistics.csv')#,transform=transform)
data_dir = '/home/zhulx/lab/summary_statistics/data/final/empirical data'
CEU_dataset_meanstdReal = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_CEU_empirical-stats_all.csv',
                        '/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_CEU_empirical-stats_all_meanstd.csv')#,transform=transform)
real_datasets = create_datasets(data_dir) # print(datasets.keys())
end_time = time.time()
print(f"Loading data took {end_time - start_time} seconds.") 
print(f"Training Dataset size: {len(train_dataset)}; additive size: {len([i for i, label in enumerate(train_dataset.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(train_dataset.labels) if label == 1])}")
print(f"Testing>200 Dataset size: {len(test_dataset_larger200)}; additive size: {len([i for i, label in enumerate(test_dataset_larger200.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger200.labels) if label == 1])}") 
print(f"Testing>400 Dataset size: {len(test_dataset_larger400)}; additive size: {len([i for i, label in enumerate(test_dataset_larger400.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger400.labels) if label == 1])}") 
print(f"Testing>600 Dataset size: {len(test_dataset_larger600)}; additive size: {len([i for i, label in enumerate(test_dataset_larger600.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger600.labels) if label == 1])}") 
print(f"Testing>800 Dataset size: {len(test_dataset_larger800)}; additive size: {len([i for i, label in enumerate(test_dataset_larger800.labels) if label == 0])} ; ressisive size: {len([i for i, label in enumerate(test_dataset_larger800.labels) if label == 1])}") 
print(f"CEU Dataset size: {len(CEU_dataset_meanstdReal)}")
print(f"TSI Dataset size: {len(real_datasets['hg19_TSI'])}")
print(f"GBR Dataset size: {len(real_datasets['hg19_GBR'])}")
print(f"IBS Dataset size: {len(real_datasets['hg19_IBS'])}")
print(f"CHB Dataset size: {len(real_datasets['hg19_CHB'])}")
print(f"CHS Dataset size: {len(datasets['hg19_CHS'])}")
print(f"JPT Dataset size: {len(real_datasets['hg19_JPT'])}")
print(f"FIN Dataset size: {len(real_datasets['hg19_FIN'])}")

train_loader = DataLoader(train_dataset, batch_size=16, shuffle=True)
test_loader_larger200 = DataLoader(test_dataset_larger200, batch_size=16, shuffle=False)
test_loader_larger400 = DataLoader(test_dataset_larger400, batch_size=16, shuffle=False)
test_loader_larger600 = DataLoader(test_dataset_larger600, batch_size=16, shuffle=False)
test_loader_larger800 = DataLoader(test_dataset_larger800, batch_size=16, shuffle=False)
CEU_loader = DataLoader(CEU_dataset_meanstdReal, batch_size=16, shuffle=False)
real_loaders = {}
for i in datasets.keys():
    data_loader = DataLoader(real_datasets[i], batch_size=16, shuffle=False)
    real_loaders[i] = data_loader
# torch.isnan(inputs).sum() # torch.isnan(labels).sum()


# In[44]:


class CustomCNN(nn.Module):
    def __init__(self, num_classes=2):
        super(CustomCNN, self).__init__()
        # since the image is 1D (1*11)
        self.conv1 = nn.Conv1d(3, 32, kernel_size=3, stride=1, padding=1) 
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, stride=1, padding=1)
        # self.conv3 = nn.Conv1d(64, 32, kernel_size=3, stride=1, padding=1)
        self.pool = nn.MaxPool1d(kernel_size=2, stride=2, padding=0)

        # Since we are using 1D convolutions, update the input shape accordingly
        self.num_features = self._get_conv_output(20)  # 20 features
        
        self.fc1 = nn.Linear(self.num_features, 512)
        self.fc2 = nn.Linear(512, num_classes) 
    
    def _get_conv_output(self, input_size):
        with torch.no_grad():
            # We assume that the input size is the size of the sequence.
            x = torch.zeros(1, 3, input_size)  # (batch_size, channels, sequence_length) #x = torch.zeros(1, *shape)
            x = self.pool(F.relu(self.conv1(x)))
            x = self.pool(F.relu(self.conv2(x)))
            # x = self.pool(F.relu(self.conv3(x)))
            return int(np.prod(x.size()[1:]))  

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        # x = self.pool(F.relu(self.conv3(x)))
        x = x.view(x.size(0), -1)  
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return torch.sigmoid(x) # Use sigmoid for binary classification; # F.log_softmax(x, dim=1) # Using log_softmax for classification; # x

# Initialize weights
def weights_init(m):
    if isinstance(m, nn.Conv1d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight)

model = CustomCNN(num_classes=2) # Initialize the model
model.apply(weights_init)
# device = torch.device('cpu')
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') # Select the device: checks if a GPU (CUDA) is available. If it is, it selects the CUDA device (i.e., GPU), otherwise, it selects the CPU.
model = model.to(device) # Move the model to the device(GPU (CUDA))

criterion = nn.CrossEntropyLoss()
optimizer = Adam(model.parameters(), lr=1e-6) # Reduced learning rate to avoid gradient exploding # optimizer = Adam(model.parameters(), lr=0.0001)
scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.1)

model_save_path = '/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/2layerCNN_sigmoid_22features'
os.makedirs(model_save_path, exist_ok=True)

def train_model(model, train_loader,test_loader_larger200,test_loader_larger400,test_loader_larger600,test_loader_larger800,real_loaders,criterion, optimizer, scheduler, num_epochs, device,model_save_path):
    best_accuracy = 0.0
    for epoch in range(num_epochs):
        model.train()
        running_loss = 0.0
        for inputs, labels in tqdm(train_loader, desc=f'Epoch {epoch+1}/{num_epochs}', unit='batch'):
            inputs, labels = inputs.to(device), labels.to(device)

            # if torch.isnan(inputs).sum() > 0 or torch.isnan(labels).sum() > 0: # Check for NaN in inputs
            #     print("NaN detected in inputs or labels")

            optimizer.zero_grad()
            outputs = model(inputs)

            # if torch.isnan(outputs).sum() > 0:
            #     print("NaN detected in outputs")
                # print(outputs)
            
            loss = criterion(outputs, labels)
            loss.backward()

            # Gradient clipping: to avoid gradient exploding
            with contextlib.redirect_stdout(io.StringIO()): # Suppress the output of clip_grad_norm_ using contextlib.redirect_stdout
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)

            optimizer.step()
            running_loss += loss.item() * inputs.size(0)
        epoch_loss = running_loss / len(train_loader.dataset)
        
        model.eval()
        correct_train = 0
        total_train = 0
        with torch.no_grad():
            for inputs, labels in train_loader:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_train += labels.size(0)
                correct_train += (predicted == labels).sum().item()
        train_accuracy = correct_train / total_train

        correct_test = 0
        total_test = 0
        all_labels_200 = [] # true labels
        all_outputs_200 = []
        predictions_200 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger200:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_200.extend(labels.cpu().numpy())
                all_outputs_200.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_200.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_200 = correct_test / total_test

        correct_test = 0
        total_test = 0
        all_labels_400 = [] # true labels
        all_outputs_400 = []
        predictions_400 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger400:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_400.extend(labels.cpu().numpy())
                all_outputs_400.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_400.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_400 = correct_test / total_test

        correct_test = 0
        total_test = 0
        all_labels_600 = [] # true labels
        all_outputs_600 = []
        predictions_600 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger600:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_600.extend(labels.cpu().numpy())
                all_outputs_600.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_600.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_600 = correct_test / total_test

        correct_test = 0
        total_test = 0
        all_labels_800 = [] # true labels
        all_outputs_800 = []
        predictions_800 = []
        with torch.no_grad():
            for inputs, labels in test_loader_larger800:
                inputs, labels = inputs.to(device), labels.to(device)
                outputs = model(inputs)
                _, predicted = torch.max(outputs.data, 1)
                total_test += labels.size(0)
                correct_test += (predicted == labels).sum().item()
                all_labels_800.extend(labels.cpu().numpy())
                all_outputs_800.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
                predictions_800.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        test_accuracy_800 = correct_test / total_test
        print('test data is done.')
        
        # CEU_outputs_meanstdReal = []
        # CEU_predictions_meanstdReal = []
        # with torch.no_grad():
        #     for inputs in CEU_loader:
        #         inputs = inputs.to(device)
        #         outputs = model(inputs)
        #         _, predicted = torch.max(outputs.data, 1)
        #         CEU_outputs_meanstdReal.extend(outputs.cpu().numpy()[:, 1])  # Assuming binary classification and outputs are logits
        #         CEU_predictions_meanstdReal.extend(predicted.cpu().numpy())  # Storing predictions for each batch
        real_outputs_meanstdReal = {}
        real_predictions_meanstdReal = {}
        with torch.no_grad():
            for dataset_name in real_loaders.keys():
                loader = real_loaders[dataset_name]
                dataset_outputs = []
                dataset_predictions = []
                for inputs in loader:
                    inputs = inputs.to(device)
                    outputs = model(inputs)
                    _, predicted = torch.max(outputs.data, 1)
                    dataset_outputs.extend(outputs.cpu().numpy()[:, 1])  # 假设是二元分类
                    dataset_predictions.extend(predicted.cpu().numpy())
                real_outputs_meanstdReal[dataset_name] = dataset_outputs
                real_predictions_meanstdReal[dataset_name] = dataset_predictions


        if test_accuracy_200 > best_accuracy:
            best_test_accuracy_200 = test_accuracy_200 # Update best accuracy
            best_test_predictions_200 = predictions_200  # Update best predictions
            test_outputs_best_200 = all_outputs_200
            best_test_accuracy_400 = test_accuracy_400 # Update best accuracy
            best_test_predictions_400 = predictions_400  # Update best predictions
            test_outputs_best_400 = all_outputs_400
            best_test_accuracy_600 = test_accuracy_600 # Update best accuracy
            best_test_predictions_600 = predictions_600  # Update best predictions
            test_outputs_best_600 = all_outputs_600
            best_test_accuracy_800 = test_accuracy_800 # Update best accuracy
            best_test_predictions_800 = predictions_800  # Update best predictions
            test_outputs_best_800 = all_outputs_800
            # CEU_outputs_meanstdReal_best = CEU_outputs_meanstdReal
            # CEU_predictions_meanstdReal_best = CEU_predictions_meanstdReal
            real_outputs_meanstdReal_best = real_outputs_meanstdReal
            real_predictions_meanstdReal_best = real_predictions_meanstdReal
            torch.save(model.state_dict(), os.path.join(model_save_path, 'best_model.pth'))
        
        scheduler.step()

        print(f'Epoch {epoch+1}/{num_epochs}, Train Accuracy: {train_accuracy:.4f}, Test Accuracy (>=200): {test_accuracy_200:.4f}')
        print(f'Test Accuracy (>=400): {test_accuracy_400:.4f}, Test Accuracy (>=600): {test_accuracy_600:.4f}, Test Accuracy (>=800): {test_accuracy_800:.4f}')

    return best_test_accuracy_200,all_labels_200,test_outputs_best_200,best_test_predictions_200,best_test_accuracy_400,all_labels_400,test_outputs_best_400,best_test_predictions_400,best_test_accuracy_600,all_labels_600,test_outputs_best_600,best_test_predictions_600,best_test_accuracy_800,all_labels_800,test_outputs_best_800,best_test_predictions_800,real_outputs_meanstdReal_best,real_predictions_meanstdReal_best

num_epochs = 20
Results = train_model(model, train_loader,test_loader_larger200,test_loader_larger400,test_loader_larger600,test_loader_larger800,real_loaders,criterion, optimizer, scheduler, num_epochs, device, model_save_path)

# Save predictions
best_test_accuracy_200,all_labels_200,test_outputs_best_200,best_test_predictions_200,best_test_accuracy_400,all_labels_400,test_outputs_best_400,best_test_predictions_400,best_test_accuracy_600,all_labels_600,test_outputs_best_600,best_test_predictions_600,best_test_accuracy_800,all_labels_800,test_outputs_best_800,best_test_predictions_800,real_outputs_meanstdReal_best,real_predictions_meanstdReal_best = Results
print(f'Best accuracy (test>=200) = {best_test_accuracy_200}')
print(f'Best accuracy (test>=400) = {best_test_accuracy_400}')
print(f'Best accuracy (test>=600) = {best_test_accuracy_600}')
print(f'Best accuracy (test>=800) = {best_test_accuracy_800}')

test_df_200 = pd.DataFrame({'True_Labels':all_labels_200,'Predicted_Labels':best_test_predictions_200,'Prediction_Prob': test_outputs_best_200})
test_df_200.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger200.csv'), index=False)
test_df_400 = pd.DataFrame({'True_Labels':all_labels_400,'Predicted_Labels':best_test_predictions_400,'Prediction_Prob': test_outputs_best_400})
test_df_400.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger400.csv'), index=False)
test_df_600 = pd.DataFrame({'True_Labels':all_labels_600,'Predicted_Labels':best_test_predictions_600,'Prediction_Prob': test_outputs_best_600})
test_df_600.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger600.csv'), index=False)
test_df_800 = pd.DataFrame({'True_Labels':all_labels_800,'Predicted_Labels':best_test_predictions_800,'Prediction_Prob': test_outputs_best_800})
test_df_800.to_csv(os.path.join(model_save_path, 'test_prediction_result_larger800.csv'), index=False)
# real_df = pd.DataFrame({'Predicted_Labels_meanstdReal': CEU_predictions_meanstdReal_best,'Prediction_Prob_meanstdReal': CEU_outputs_meanstdReal_best})
# real_df.to_csv(os.path.join(model_save_path, 'CEU_prediction_result.csv'), index=False)
for dataset_name in real_outputs_meanstdReal_best.keys():
    real_df = pd.DataFrame({
        'Predicted_Labels_meanstdReal': real_predictions_meanstdReal_best[dataset_name],
        'Prediction_Prob_meanstdReal': real_outputs_meanstdReal_best[dataset_name]
    })
    save_path = os.path.join(model_save_path, f'{dataset_name}_prediction_result.csv')
    real_df.to_csv(save_path, index=False)


# In[45]:


roc_data = []
for threshold in [200, 400, 600, 800]:
        file_path = f'/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/2layerCNN_sigmoid_22features/test_prediction_result_larger{threshold}.csv'
        sample_df = pd.read_csv(file_path)
        test_labels_best = sample_df['True_Labels']
        test_outputs_best = sample_df['Prediction_Prob']
        auc_roc = roc_auc_score(test_labels_best, test_outputs_best)
        fpr, tpr, _ = roc_curve(test_labels_best, test_outputs_best)
        print(f'auc (test>{threshold}) =',auc_roc)
        roc_data.append((fpr, tpr, auc_roc, f'test>{threshold}'))

plt.figure(figsize=(10, 8))
for fpr, tpr, auc_roc, label in roc_data:
    plt.plot(fpr, tpr, label=f'{label} (AUC = {auc_roc:.2f})')

# plt.plot([0, 1], [0, 1], 'k--', label='No Skill')  # add a line for no skill classifier
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curves for Different Thresholds')
plt.legend(loc='best')
plt.show()
plt.savefig(os.path.join(model_save_path,'roc_curve_cnn.png'))


# 1. 假阳性率 (fpr)
# FPR (False Positive Rate) is the proportion of negative samples incorrectly labeled as positive: $FPR = \frac{FP}{FP + TN}$
# 
# 2. True Positive Rate (TPR)
# TPR, also called sensitivity, is the proportion of actual positives correctly identified: $TPR = \frac{TP}{TP + FN}$
# 
# 3. Length difference from True_Labels
# In roc_curve, fpr and tpr arrays are usually not the same length as True_Labels. This is because they are computed at various thresholds. Each unique prediction probability in data[i][label]['Prediction_Prob'] is treated as a threshold. FPR and TPR are calculated at each of these, plus usually an additional threshold of 0. So their length equals the number of unique prediction scores plus one.
# 
# 4. Use of ROC Curve
# ROC curve is a visual representation with FPR on the x-axis and TPR on the y-axis. It shows the model's ability to distinguish between classes across thresholds. The closer the curve is to the top-left corner, the better the model, indicating high sensitivity with low false positives.


# ===============================# SECTION 3: 3 Analysis# ===============================


# ===============================# SECTION 3.1: 3.1 Model Comparision# ===============================


# ===============================# SECTION 3.1: 3.1.1 Test data# ===============================


# ===============================# SECTION 3.1: 3.1.1.1 pairwise plots# ===============================

# In[56]:


# Updated file paths including both CSV and Excel files
file_paths = {
    "catboost":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_catboost.xlsx",
    "elastic_net":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_elasticnet.xlsx",
    "extra_trees":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_extratrees.xlsx",
    "gradient_boosting":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_gb.xlsx",
    "lasso":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_lasso.xlsx",
    "lightgbm":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_lightgbm.xlsx",
    "mlp":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_mlp.xlsx",
    "random_forest":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_rf.xlsx",
    "xgboost":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_xgboost.xlsx",
    "xgboost_allfea":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_xgboost_allfea.xlsx",
    "logistic": "/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_logistic.xlsx",
    "2layerCNN": "/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/2layerCNN_sigmoid_22features",
    "3layerCNN": "/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features",
}

# Load data into a dictionary, handling both CSV and Excel files
data = {}
i=1
exon_density_labels = ['test_exon_density_-1.3', 'test_exon_density_-0.74', 'test_exon_density_-0.18', 'test_exon_density_0.37']

for i in range(4):
    data[i] = {}
    for label, path in file_paths.items():
        if path.endswith('features'):
            data[i][label] = pd.read_csv(os.path.join(path,f'test_prediction_result_larger{(i+1)*200}.csv'))
            data[i][label]['PredictedLabels_thr0.9'] = (data[i][label]['Prediction_Prob'] >= 0.9).astype(int)
        elif path.endswith('.xlsx'):
            # Read specific columns from the Excel file
            data[i][label] = pd.read_excel(path, sheet_name=exon_density_labels[i], engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
            data[i][label].rename(columns={'Probability': 'Prediction_Prob','Predicted_Label':'Predicted_Labels'}, inplace=True)
            data[i][label]['PredictedLabels_thr0.9'] = (data[i][label]['Prediction_Prob'] >= 0.9).astype(int)


# In[8]:


for i in range(4):   
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30,65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {200*(i+1)})', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(data[i].keys(), 2)):
        label_x, label_y = pair
        ax.scatter(data[i][label_x]['Prediction_Prob'], data[i][label_y]['Prediction_Prob'], alpha=0.5)
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.1.2 Correlation Heatmap# ===============================

# In[17]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
final_prob = pd.DataFrame()
for i in range(4):
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = data[i][label]['Prediction_Prob']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density >= {200*(i+1)})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.1.3 Confusion matrix# ===============================

# In[18]:


for i in range(4):
    print(f'For exon density >= {200*(i+1)}')
    truelabels = data[i]['2layerCNN']['True_Labels']
    for label in file_paths:
        print(f"Confusion matrix for {label}:")
        print(confusion_matrix(truelabels, data[i][label]['Predicted_Labels']))
    print('--'*50)


# confusion matrix:  
# |                      | Predicted Negative | Predicted Positive |
# |----------------------|--------------------|--------------------|
# | **Actual Negative**  | TN                 | FP                 |
# | **Actual Positive**  | FN                 | TP                 |


# ===============================# SECTION 3.1: 3.1.1.4 ROC AUC# ===============================

# In[72]:


# using threshold = 0.5
for i in range(4):
    plt.figure(figsize=(10, 8))
    True_Labels = data[i]['2layerCNN']['True_Labels']

    auc_list = []
    for label in file_paths:
        auc_roc = roc_auc_score(True_Labels, data[i][label]['Prediction_Prob'])
        auc_list.append((label, auc_roc))
    auc_list.sort(key=lambda x: x[1], reverse=True) # sorting labels using AUC!!!
    
    max_label_length = max(len(label) for label, _ in auc_list)
    max_auc_length = max(len(f"{auc:.5f}") for _, auc in auc_list)
    for label, auc_roc in auc_list:
        auc_roc = roc_auc_score(True_Labels, data[i][label]['Prediction_Prob'])
        fpr, tpr, _ = roc_curve(True_Labels, data[i][label]['Prediction_Prob'])
        # plt.plot(fpr, tpr, label=f'{label} (AUC = {auc_roc:.5f})')
                
        intersection_index = np.where(fpr >= 0.01)[0][0]
        x_intersection = fpr[intersection_index]
        y_intersection = tpr[intersection_index]
        # plt.scatter(x_intersection, y_intersection, color='red')
        # plt.text(x_intersection, y_intersection, f'({x_intersection:.2f}, {y_intersection:.2f})', 
        #         fontsize=10, verticalalignment='bottom', horizontalalignment='right')

        # plt.plot(fpr, tpr, label=f'{label} (AUC = {auc_roc:.5f}, TPR = {y_intersection:.2f} at FPR=0.01)')
        label_fmt = "{:<{label_width}}  AUC = {:<{auc_width}.5f},  TPR = {:<.3f} at FPR=0.01"
        formatted_label = label_fmt.format(label, auc_roc, y_intersection, label_width=max_label_length, auc_width=max_auc_length)
        plt.plot(fpr, tpr, label=formatted_label)

        print(label,'when fpr=0.01')
        print('prob =',data[i][label]['Prediction_Prob'][intersection_index])
        print('tpr=',y_intersection)

    plt.axvline(x=0.01, color='gray', linestyle='--', linewidth=1)
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.title(f'ROC Curve for exon density >= {200*(i+1)}') #Receiver Operating Characteristic (ROC) Curve
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right', labelspacing=0.5) 
    plt.grid(True) 
    plt.show()


# In[73]:


# using threshold = 0.5
label_list = ['catboost','xgboost','xgboost_allfea']
for i in range(4):
    plt.figure(figsize=(10, 8))
    True_Labels = data[i]['2layerCNN']['True_Labels']

    auc_list = []
    for label in label_list:
        auc_roc = roc_auc_score(True_Labels, data[i][label]['Prediction_Prob'])
        auc_list.append((label, auc_roc))
    auc_list.sort(key=lambda x: x[1], reverse=True) # sorting labels using AUC!!!
    
    max_label_length = max(len(label) for label, _ in auc_list)
    max_auc_length = max(len(f"{auc:.5f}") for _, auc in auc_list)
    for label, auc_roc in auc_list:
        auc_roc = roc_auc_score(True_Labels, data[i][label]['Prediction_Prob'])
        fpr, tpr, _ = roc_curve(True_Labels, data[i][label]['Prediction_Prob'])
        # plt.plot(fpr, tpr, label=f'{label} (AUC = {auc_roc:.5f})')
                
        intersection_index = np.where(fpr >= 0.01)[0][0]
        x_intersection = fpr[intersection_index]
        y_intersection = tpr[intersection_index]
        # plt.scatter(x_intersection, y_intersection, color='red')
        # plt.text(x_intersection, y_intersection, f'({x_intersection:.2f}, {y_intersection:.2f})', 
        #         fontsize=10, verticalalignment='bottom', horizontalalignment='right')

        # plt.plot(fpr, tpr, label=f'{label} (AUC = {auc_roc:.5f}, TPR = {y_intersection:.2f} at FPR=0.01)')
        label_fmt = "{:<{label_width}}  AUC = {:<{auc_width}.5f},  TPR = {:<.3f} at FPR=0.01"
        formatted_label = label_fmt.format(label, auc_roc, y_intersection, label_width=max_label_length, auc_width=max_auc_length)
        plt.plot(fpr, tpr, label=formatted_label)

        print(label,'when fpr=0.01')
        print('prob =',data[i][label]['Prediction_Prob'][intersection_index])
        print('tpr=',y_intersection)

    plt.axvline(x=0.01, color='gray', linestyle='--', linewidth=1)
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.title(f'ROC Curve for exon density >= {200*(i+1)}') #Receiver Operating Characteristic (ROC) Curve
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right', labelspacing=0.5) 
    plt.grid(True) 
    plt.show()




# ===============================# SECTION 3.1: 3.1.2 Real data# ===============================


# ===============================# SECTION 3.1: 3.1.2.1 pairwise plots# ===============================

# In[21]:


file_paths = {
    "catboost":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_catboost.xlsx",
    "elastic_net":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_elasticnet.xlsx",
    "extra_trees":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_extratrees.xlsx",
    "gradient_boosting":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_gb.xlsx",
    "lasso":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_lasso.xlsx",
    "lightgbm":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_lightgbm.xlsx",
    "mlp":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_mlp.xlsx",
    "random_forest":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_rf.xlsx",
    "xgboost":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_xgboost.xlsx",
    "xgboost_allfea":"/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_xgboost_allfea.xlsx",
    "logistic": "/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features/comparison/result_logistic.xlsx",
    "2layerCNN": "/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/2layerCNN_sigmoid_22features",
    "3layerCNN": "/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features",
}

# Load data into a dictionary, handling both CSV and Excel files
CEU_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_CEU_empirical-stats_all.csv')
CHB_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_CHB_empirical-stats_all.csv')
CHS_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_CHS_empirical-stats_all.csv')
FIN_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_FIN_empirical-stats_all.csv')
GBR_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_GBR_empirical-stats_all.csv')
IBS_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_IBS_empirical-stats_all.csv')
JPT_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_JPT_empirical-stats_all.csv')
TSI_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/hg19_TSI_empirical-stats_all.csv')

# CEU
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_CEU_prediction_result.csv'))
        CEU_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        CEU_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='CEU_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        CEU_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        CEU_data[f'PredictedProb_{label}'] = a['Probability']
outliers_given = pd.read_csv('/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/comparison/Extra tree classifier outlier.txt',delimiter='\t')
outliers_given['ExtraTreeClassifier_outlier_identified'] = 1
CEU_data = pd.merge(CEU_data, outliers_given[['chrom', 'start', 'end', 'ExtraTreeClassifier_outlier_identified']], 
                       on=['chrom', 'start', 'end'], how='left')
CEU_data['ExtraTreeClassifier_outlier_identified'].fillna(0, inplace=True)
CEU_data['ExtraTreeClassifier_outlier_identified'] = CEU_data['ExtraTreeClassifier_outlier_identified'].astype(int)

# CHB
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_CHB_prediction_result.csv'))
        CHB_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        CHB_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='CHB_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        CHB_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        CHB_data[f'PredictedProb_{label}'] = a['Probability']

# CHS
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_CHS_prediction_result.csv'))
        CHS_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        CHS_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='CHS_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        CHS_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        CHS_data[f'PredictedProb_{label}'] = a['Probability']

# FIN
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_FIN_prediction_result.csv'))
        FIN_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        FIN_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='FIN_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        FIN_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        FIN_data[f'PredictedProb_{label}'] = a['Probability']

# GBR
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_GBR_prediction_result.csv'))
        GBR_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        GBR_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='GBR_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        GBR_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        GBR_data[f'PredictedProb_{label}'] = a['Probability']

# IBS
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_IBS_prediction_result.csv'))
        IBS_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        IBS_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='IBS_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        IBS_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        IBS_data[f'PredictedProb_{label}'] = a['Probability']

# JPT
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_JPT_prediction_result.csv'))
        JPT_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        JPT_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='JPT_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        JPT_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        JPT_data[f'PredictedProb_{label}'] = a['Probability']

# TSI
for label, path in file_paths.items():
    if path.endswith('features'):
        a = pd.read_csv(os.path.join(path,'hg19_TSI_prediction_result.csv'))
        TSI_data[f'PredictedLabels_{label}'] = a['Predicted_Labels_meanstdReal']
        TSI_data[f'PredictedProb_{label}'] = a['Prediction_Prob_meanstdReal']
    elif path.endswith('.xlsx'):
        # Read specific columns from the Excel file
        a = pd.read_excel(path, sheet_name='TSI_empirical', engine='openpyxl', usecols=['Predicted_Label', 'Probability'])
        TSI_data[f'PredictedLabels_{label}'] = a['Predicted_Label']
        TSI_data[f'PredictedProb_{label}'] = a['Probability']


# In[28]:


CEU_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/CEU_data_final_result.csv')
CHB_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/CHB_data_final_result.csv')
CHS_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/CHS_data_final_result.csv')
FIN_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/FIN_data_final_result.csv')
GBR_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/GBR_data_final_result.csv')
IBS_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/IBS_data_final_result.csv')
JPT_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/JPT_data_final_result.csv')
TSI_data.to_csv('/home/zhulx/lab/summary_statistics/data/final/TSI_data_final_result.csv')




# ===============================# SECTION 3.1: 3.1.2.1.1 CEU# ===============================

# In[29]:


# CEU
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for CEU', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = CEU_data[CEU_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.1.2 CHB# ===============================

# In[30]:


# CHB
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for CHB', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = CHB_data[CHB_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.1.3 CHS# ===============================

# In[31]:


# CHS
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for CHS', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = CHS_data[CHS_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.1.4 FIN# ===============================

# In[32]:


# FIN
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for FIN', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = FIN_data[FIN_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.1.5 GBR# ===============================

# In[33]:


# GBR
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for GBR', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = GBR_data[GBR_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.1.6 IBS# ===============================

# In[34]:


# IBS
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for IBS', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = IBS_data[IBS_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.1.7 JPT# ===============================

# In[35]:


# JPT
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for JPT', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = JPT_data[JPT_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.1.8 TSI# ===============================

# In[36]:


# TSI
for threshold in [200,400,600,800]:
    # Create scatter plots for each pair of datasets
    fig, axs = plt.subplots(13,6, figsize=(30, 65))  # Adjust subplot grid size based on the number of combinations
    axs = axs.flatten()
    fig.suptitle(f'Pairwise Comparison of Prediction_Prob (exon density >= {threshold}) for TSI', fontsize=40,y=0.95)

    # Generate scatter plots for each combination of datasets
    for ax, (pair) in zip(axs, combinations(file_paths.keys(), 2)):
        label_x, label_y = pair
        filtered_real_data = TSI_data[TSI_data['exon_density']>=threshold]
        scatter_normal = ax.scatter(filtered_real_data[f'PredictedProb_{label_x}'], filtered_real_data[f'PredictedProb_{label_y}'], alpha=0.5, label='Normal')
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)
        ax.set_title(f'{label_x} vs {label_y}')

        # outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1]
        # scatter_outlier = ax.scatter(outlier_data[f'PredictedProb_{label_x}'], outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='red', edgecolor='black', marker='o', label='ExtraTreeClassifier_outlier_identified')

        # new_outlier_data = filtered_real_data[filtered_real_data['outlier_identified'] == 1][filtered_real_data['chrom']==6][filtered_real_data['start']>=26400000][filtered_real_data['end']<=34200000]
        # new_scatter_outlier = ax.scatter(new_outlier_data[f'PredictedProb_{label_x}'], new_outlier_data[f'PredictedProb_{label_y}'], alpha=0.7, color='yellow', edgecolor='black', marker='o', label='New Outlier Identified')

        # Set aspect ratio to be equal, making each subplot a square
        ax.set_aspect('equal')

        # Add y=x line for reference
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])
        ]
        ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.legend(handles=[scatter_normal], loc='upper right')# ax.legend(handles=[scatter_normal, scatter_outlier,new_scatter_outlier], loc='upper right')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    
# Remove any unused subplots
for j in range(78, len(axs)): # 78 = num of plots
    fig.delaxes(axs[j])




# ===============================# SECTION 3.1: 3.1.2.2 Correlation Heatmaps# ===============================


# ===============================# SECTION 3.1: 3.1.2.2.1 CEU# ===============================

# In[22]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = CEU_data[CEU_data['exon_density']>threshold]
    for label in file_paths: 
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.2.2.2 CHB# ===============================

# In[23]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = CHB_data[CHB_data['exon_density']>threshold]
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.2.2.3 CHS# ===============================

# In[24]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = CHS_data[CHS_data['exon_density']>threshold]
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.2.2.4 FIN# ===============================

# In[25]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = FIN_data[FIN_data['exon_density']>threshold]
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.2.2.5 GBR# ===============================

# In[26]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = GBR_data[GBR_data['exon_density']>threshold]
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.2.2.6 IBS# ===============================

# In[27]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = IBS_data[IBS_data['exon_density']>threshold]
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.2.2.7 JPT# ===============================

# In[28]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = JPT_data[JPT_data['exon_density']>threshold]
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.1: 3.1.2.2.8 TSI# ===============================

# In[29]:


# colors = ["#6e7f80", "#708090", "#8f9779",  "#536872"]# colors = ["#F7FF3C", "#DC143C"]# colors = ["yellow","red"]
cmap = sns.light_palette("#9f162c", as_cmap=True) # cmap = sns.diverging_palette(0, 25, as_cmap=True,s=80,l=50)# cmap = sns.diverging_palette(230, 25, as_cmap=True,s=80,l=50)# cmap = LinearSegmentedColormap.from_list("mycmap", colors) #cmap = sns.diverging_palette(230, 20, as_cmap=True)
for threshold in [200,400,600,800]:
    final_prob = pd.DataFrame()
    filtered_real_data = TSI_data[TSI_data['exon_density']>threshold]
    for label in file_paths:
        final_prob[f'PredictedProb_{label}'] = filtered_real_data[f'PredictedProb_{label}']
    correlation_matrix = final_prob.corr()
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap=cmap, fmt='.2f',mask=mask,linewidths=.5, linecolor='white', annot_kws={"size": 10},square=True)
    plt.title(f'Correlation Matrix Heatmap for prediction prob (exon density > {threshold})')
    plt.show()




# ===============================# SECTION 3.2: 3.2 verification result# ===============================

# In[31]:


##### load data #####
class WindowDataset(Dataset):
    def __init__(self, csv_file, mean_std_dir, transform=None):
        self.csv_file = csv_file
        self.transform = transform
        self.column_names = [
            'num_seg_p2','num_variant_window','exon_density', 'mean_introg_anc', 'exon_window', 'introg_anc_window', 'divergence_p3_p1', 'divergence_p3_p2', 'watterson_theta_p3', 'windowed_tajima_d_p3', 'D', 
            'Het', 'Q95', 'U0', 'U20', 'U50', 'U80', 'num_seg_p1', 'num_seg_p3', 'num_private_seg_p1', 'num_private_seg_p2', 'num_private_seg_p3'
        ]
        
        # Load entire dataset to fetch column indices
        all_data = np.genfromtxt(csv_file, delimiter=',', names=True, dtype=None)
        self.data = np.column_stack([all_data[name] for name in self.column_names])
        # Check if 'dominance' column exists
        if 'dominance' in all_data.dtype.names:
            self.labels = all_data['dominance'].astype(int)
            self.has_labels = True
        else:
            self.has_labels = False

        # Load means and stds for the selected columns only
        mean_std_data = pd.read_csv(mean_std_dir)
        filtered_mean_std_data = mean_std_data.set_index('Unnamed: 0').reindex(self.column_names).reset_index()
        self.means = filtered_mean_std_data['Mean'].to_numpy()
        self.stds = filtered_mean_std_data['Std'].to_numpy()

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        x = self.data[idx].astype(float)
        x = (x - self.means) / self.stds  # Normalize x
        x = np.stack([x, x, x], axis=0)  # Duplicate the single channel to create 3 channels
        x = torch.tensor(x, dtype=torch.float32)
        if self.transform:
            x = self.transform(x)
        
        if self.has_labels:
            y = self.labels[idx]
            return x, torch.tensor(y, dtype=torch.long)
        else:
            return x

transform = transforms.Compose([
    transforms.Lambda(lambda x: x.unsqueeze(0))
])


start_time = time.time()
verification_dataset_meanstdReal = WindowDataset('/home/zhulx/lab/summary_statistics/data/final/empirical data/24520combined_varyh_test_stats.csv',
                    '/home/zhulx/lab/summary_statistics/data/final/empirical data/24520combined_varyh_test_stats_meanstd.csv')#,transform=transform)
end_time = time.time()
print(f"Loading data took {end_time - start_time} seconds.") 
print(f"Verification Dataset size: {len(verification_dataset_meanstdReal)}")

Verification_loader = DataLoader(verification_dataset_meanstdReal, batch_size=16, shuffle=False)




# ===============================# SECTION 3.2: 3.2.1 two-layer CNN# ===============================

# In[32]:


##### load model #####
class CustomCNN(nn.Module):
    def __init__(self, num_classes=2):
        super(CustomCNN, self).__init__()
        # since the image is 1D (1*11)
        self.conv1 = nn.Conv1d(3, 32, kernel_size=3, stride=1, padding=1) 
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, stride=1, padding=1)
        # self.conv3 = nn.Conv1d(64, 32, kernel_size=3, stride=1, padding=1)
        self.pool = nn.MaxPool1d(kernel_size=2, stride=2, padding=0)

        # Since we are using 1D convolutions, update the input shape accordingly
        self.num_features = self._get_conv_output(20)  # 20 features
        
        self.fc1 = nn.Linear(self.num_features, 512)
        self.fc2 = nn.Linear(512, num_classes) 
    
    def _get_conv_output(self, input_size):
        with torch.no_grad():
            # We assume that the input size is the size of the sequence.
            x = torch.zeros(1, 3, input_size)  # (batch_size, channels, sequence_length) #x = torch.zeros(1, *shape)
            x = self.pool(F.relu(self.conv1(x)))
            x = self.pool(F.relu(self.conv2(x)))
            # x = self.pool(F.relu(self.conv3(x)))
            return int(np.prod(x.size()[1:]))  

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        # x = self.pool(F.relu(self.conv3(x)))
        x = x.view(x.size(0), -1)  
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return torch.sigmoid(x) # Use sigmoid for binary classification; # F.log_softmax(x, dim=1) # Using log_softmax for classification; # x

# Initialize weights
def weights_init(m):
    if isinstance(m, nn.Conv1d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight)

model = CustomCNN(num_classes=2) # Initialize the model

model_save_path = '/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/2layerCNN_sigmoid_22features'
model.load_state_dict(torch.load(os.path.join(model_save_path, 'best_model.pth')))

##### fit verification data set #####
def evaluate_model_on_verification(model, Verification_loader, device):
    all_labels = []
    all_outputs = []
    predictions = []
    
    with torch.no_grad():  # No gradient computation needed during evaluation
        for inputs, labels in Verification_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            _, predicted = torch.max(outputs.data, 1)  # Get predicted class for each input
            
            all_labels.extend(labels.cpu().numpy())
            all_outputs.extend(outputs.cpu().numpy()[:, 1])  # Get predicted probabilities (binary classification assumed)
            predictions.extend(predicted.cpu().numpy())  # Store predicted class
    
    return all_labels, all_outputs, predictions

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') # Select the device: checks if a GPU (CUDA) is available. If it is, it selects the CUDA device (i.e., GPU), otherwise, it selects the CPU.
model = model.to(device) # Move the model to the device(GPU (CUDA))
verification_labels, verification_outputs, verification_predictions = evaluate_model_on_verification(model, Verification_loader, device)
verification_df = pd.DataFrame({
    'Predicted_Labels': verification_predictions,
    'Prediction_Prob': verification_outputs
})
verification_save_path = os.path.join(model_save_path, 'verification_prediction_results.csv')
verification_df.to_csv(verification_save_path, index=False)




# ===============================# SECTION 3.2: 3.2.2 three-layer CNN# ===============================

# In[39]:


##### load model #####
class CustomCNN(nn.Module):
    def __init__(self, num_classes=2):
        super(CustomCNN, self).__init__()
        # since the image is 1D (1*11)
        self.conv1 = nn.Conv1d(3, 32, kernel_size=3, stride=1, padding=1) 
        self.conv2 = nn.Conv1d(32, 64, kernel_size=3, stride=1, padding=1)
        self.conv3 = nn.Conv1d(64, 32, kernel_size=3, stride=1, padding=1)
        self.pool = nn.MaxPool1d(kernel_size=2, stride=2, padding=0)

        # Since we are using 1D convolutions, update the input shape accordingly
        self.num_features = self._get_conv_output(22)  # 20 features
        
        self.fc1 = nn.Linear(self.num_features, 512)
        self.fc2 = nn.Linear(512, num_classes) 
    
    def _get_conv_output(self, input_size):
        with torch.no_grad():
            # We assume that the input size is the size of the sequence.
            x = torch.zeros(1, 3, input_size)  # (batch_size, channels, sequence_length) #x = torch.zeros(1, *shape)
            x = self.pool(F.relu(self.conv1(x)))
            x = self.pool(F.relu(self.conv2(x)))
            x = self.pool(F.relu(self.conv3(x)))
            return int(np.prod(x.size()[1:]))  

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = self.pool(F.relu(self.conv3(x)))
        x = x.view(x.size(0), -1)  
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return torch.sigmoid(x) # Use sigmoid for binary classification; # F.log_softmax(x, dim=1) # Using log_softmax for classification; # x

# Initialize weights
def weights_init(m):
    if isinstance(m, nn.Conv1d) or isinstance(m, nn.Linear):
        nn.init.kaiming_normal_(m.weight)

model = CustomCNN(num_classes=2) # Initialize the model

model_save_path = '/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features'
model.load_state_dict(torch.load(os.path.join(model_save_path, 'best_model.pth')))


##### fit verification data set #####
def evaluate_model_on_verification(model, Verification_loader, device):
    all_labels = []
    all_outputs = []
    predictions = []
    
    with torch.no_grad():  # No gradient computation needed during evaluation
        for inputs, labels in Verification_loader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            _, predicted = torch.max(outputs.data, 1)  # Get predicted class for each input
            
            all_labels.extend(labels.cpu().numpy())
            all_outputs.extend(outputs.cpu().numpy()[:, 1])  # Get predicted probabilities (binary classification assumed)
            predictions.extend(predicted.cpu().numpy())  # Store predicted class
    
    return all_labels, all_outputs, predictions

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu') # Select the device: checks if a GPU (CUDA) is available. If it is, it selects the CUDA device (i.e., GPU), otherwise, it selects the CPU.
model = model.to(device) # Move the model to the device(GPU (CUDA))
verification_labels, verification_outputs, verification_predictions = evaluate_model_on_verification(model, Verification_loader, device)

verification_df = pd.DataFrame({
    'Predicted_Labels': verification_predictions,
    'Prediction_Prob': verification_outputs
})
verification_save_path = os.path.join(model_save_path, 'verification_prediction_results.csv')
verification_df.to_csv(verification_save_path, index=False)


# In[42]:


label = "2layerCNN"
Verify_data = pd.read_csv('/home/zhulx/lab/summary_statistics/data/final/empirical data/24520combined_varyh_test_stats.csv')
a = pd.read_csv(os.path.join('/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/2layerCNN_sigmoid_22features','verification_prediction_results.csv'))
Verify_data[f'PredictedProb_{label}'] = a['Prediction_Prob']
Verify_data[f'PredictedLabels_{label}'] = a['Predicted_Labels']
Verify_data[f'PredictedLabels_thr0.9_{label}'] = (Verify_data[f'PredictedProb_{label}'] >= 0.9).astype(int)


label = "3layerCNN"
a = pd.read_csv(os.path.join('/home/zhulx/lab/summary_statistics/model_cnn_1row/feature_selection_stratified/trainAll_testOthers/3layerCNN_sigmoid_22features','verification_prediction_results.csv'))
Verify_data[f'PredictedProb_{label}'] = a['Prediction_Prob']
Verify_data[f'PredictedLabels_{label}'] = a['Predicted_Labels']
Verify_data[f'PredictedLabels_thr0.9_{label}'] = (Verify_data[f'PredictedProb_{label}'] >= 0.9).astype(int)

Verify_data


# In[51]:


import pandas as pd
import matplotlib.pyplot as plt

# Step 1: Get unique values in 'dominance' column and sort
sorted_dominance_values = sorted(Verify_data['dominance'].unique())

# Step 2: Calculate proportion where PredictedLabels_2layerCNN == 1
proportions = []
for value in sorted_dominance_values:
    group = Verify_data[Verify_data['dominance'] == value]  # Get group for each dominance value
    proportion = (group['PredictedLabels_2layerCNN'] == 1).mean()  # 计算 PredictedLabels_2layerCNN == 1 的比例
    proportions.append(proportion)

# Step 3: Create scatter plot
plt.figure(figsize=(10, 6))
plt.scatter(sorted_dominance_values, proportions, color='b', marker='o')
plt.plot(sorted_dominance_values, proportions, color='b', linestyle='-', label='Line Connection')  # Connect points with line


# Step 4: Set plot labels and title
plt.xlabel('Dominance (Sorted Unique Values)', fontsize=12)
plt.ylabel('Proportion of PredictedLabels_2layerCNN == 1', fontsize=12)
plt.title('Proportion of PredictedLabels_2layerCNN == 1 for Each Dominance Value', fontsize=14)

# Show grid lines
plt.grid(True)

# Step 5: Display the plot
plt.show()


# In[82]:


import pandas as pd
import matplotlib.pyplot as plt

# Get unique values in 'dominance' column and sort
sorted_dominance_values = sorted(Verify_data['dominance'].unique())

# Initialize dictionary to store proportions for each condition
proportions = {
    'PredictedLabels_2layerCNN': [],
    'PredictedLabels_thr0.9_2layerCNN': [],
    'PredictedLabels_3layerCNN': [],
    'PredictedLabels_thr0.9_3layerCNN': []
}

# Calculate proportion where PredictedLabels == 0 for each dominance value
for value in sorted_dominance_values:
    group = Verify_data[Verify_data['dominance'] == value]  # Get group for each dominance value
    proportions['PredictedLabels_2layerCNN'].append((group['PredictedLabels_2layerCNN'] == 1).mean())
    proportions['PredictedLabels_thr0.9_2layerCNN'].append((group['PredictedLabels_thr0.9_2layerCNN'] == 1).mean())
    proportions['PredictedLabels_3layerCNN'].append((group['PredictedLabels_3layerCNN'] == 1).mean())
    proportions['PredictedLabels_thr0.9_3layerCNN'].append((group['PredictedLabels_thr0.9_3layerCNN'] == 1).mean())

# Create scatter plot
plt.figure(figsize=(10, 6))
colors = ['b', 'g', 'r', 'c']  # Different colors represent different columns
labels = [
    'PredictedLabels_2layerCNN',
    'PredictedLabels_thr0.9_2layerCNN',
    'PredictedLabels_3layerCNN',
    'PredictedLabels_thr0.9_3layerCNN'
]

for i, label in enumerate(labels):
    plt.scatter(sorted_dominance_values, proportions[label], color=colors[i], marker='o', label=f'Proportion of {label} == 0')
    plt.plot(sorted_dominance_values, proportions[label], color=colors[i], linestyle='-')

# Set plot labels and title
plt.xlabel('Dominance (Sorted Unique Values)', fontsize=12)
plt.ylabel('Proportion of Predicted Labels == 1', fontsize=12)
plt.title('Proportion of Zero Predicted Labels for Each Dominance Value by Model', fontsize=14)

# Show grid lines
plt.grid(True)
plt.legend()

# Display the plot
plt.show()



plt.figure(figsize=(10, 6))
colors = ['b', 'g', 'r', 'c']  # Different colors represent different columns
labels = [
    'PredictedLabels_2layerCNN',
    'PredictedLabels_thr0.9_2layerCNN',
    'PredictedLabels_3layerCNN',
    'PredictedLabels_thr0.9_3layerCNN'
]

for i, label in enumerate(labels):
    plt.scatter(np.log(sorted_dominance_values), proportions[label], color=colors[i], marker='o', label=f'Proportion of {label} == 0')
    plt.plot(np.log(sorted_dominance_values), proportions[label], color=colors[i], linestyle='-')

# Set plot labels and title
plt.xlabel('log(Dominance) (Sorted Unique Values)', fontsize=12)
plt.ylabel('Proportion of Predicted Labels == 1', fontsize=12)
plt.title('Proportion of Zero Predicted Labels for Each log(Dominance) Value by Model', fontsize=14)

# Show grid lines
plt.grid(True)
plt.legend()

# Display the plot
plt.show()


# In[83]:


colors = ['b', 'g']  # Limit colors to two per plot
labels_2layer = [
    'PredictedLabels_2layerCNN',
    'PredictedLabels_thr0.9_2layerCNN'
]
labels_3layer = [
    'PredictedLabels_3layerCNN',
    'PredictedLabels_thr0.9_3layerCNN'
]



# Plot figures related to 2-layer CNN
plt.figure(figsize=(10, 6))
for i, label in enumerate(labels_2layer):
    plt.scatter(sorted_dominance_values, proportions[label], color=colors[i], marker='o', label=f'Proportion of {label} == 0')
    plt.plot(sorted_dominance_values, proportions[label], color=colors[i], linestyle='-')
plt.xlabel('Dominance (Sorted Unique Values)', fontsize=12)
plt.ylabel('Proportion of Predicted Labels == 1', fontsize=12)
plt.title('Proportion of Zero Predicted Labels for 2layer Models by Dominance Value', fontsize=14)
plt.grid(True)
plt.legend()
plt.show()
plt.figure(figsize=(10, 6))
for i, label in enumerate(labels_2layer):
    plt.scatter(np.log(sorted_dominance_values), proportions[label], color=colors[i], marker='o', label=f'Proportion of {label} == 0')
    plt.plot(np.log(sorted_dominance_values), proportions[label], color=colors[i], linestyle='-')
plt.xlabel('log(Dominance) (Sorted Unique Values)', fontsize=12)
plt.ylabel('Proportion of Predicted Labels == 1', fontsize=12)
plt.title('Proportion of Zero Predicted Labels for 2layer Models by log(Dominance) Value', fontsize=14)
plt.grid(True)
plt.legend()
plt.show()

# Plot figures related to 3-layer CNN
plt.figure(figsize=(10, 6))
for i, label in enumerate(labels_3layer):
    plt.scatter(sorted_dominance_values, proportions[label], color=colors[i], marker='o', label=f'Proportion of {label} == 0')
    plt.plot(sorted_dominance_values, proportions[label], color=colors[i], linestyle='-')
plt.xlabel('Dominance (Sorted Unique Values)', fontsize=12)
plt.ylabel('Proportion of Predicted Labels == 1', fontsize=12)
plt.title('Proportion of Zero Predicted Labels for 3layer Models by Dominance Value', fontsize=14)
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
for i, label in enumerate(labels_3layer):
    plt.scatter(np.log(sorted_dominance_values), proportions[label], color=colors[i], marker='o', label=f'Proportion of {label} == 0')
    plt.plot(np.log(sorted_dominance_values), proportions[label], color=colors[i], linestyle='-')
plt.xlabel('log(Dominance) (Sorted Unique Values)', fontsize=12)
plt.ylabel('Proportion of Predicted Labels == 1', fontsize=12)
plt.title('Proportion of Zero Predicted Labels for 3layer Models by log(Dominance) Value', fontsize=14)
plt.grid(True)
plt.legend()
plt.show()


# 

# 
