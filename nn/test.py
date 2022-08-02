import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset
from torchvision import datasets
from torchvision.transforms import ToTensor

plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')


data_folder = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-04-19'
invivo_folder = r'C:/Users/BRoehrich/Desktop/EIS-EAB data/2022-05-11/rat 1 -330mV'
train_csv = os.path.join(data_folder, 'labels.csv')
test_csv  = os.path.join(data_folder, 'test.csv')
invivo_csv  = os.path.join(invivo_folder, 'test.csv')

def extract_titration(d_folder):
    concs = []
    data = []
    
    files = []
    test_files = []
    test_concs = []
    
    for folder in os.listdir(d_folder):
        if folder.endswith('.csv'):
            continue
        
        _, conc = folder.split('_')
        
        conc = float(conc)
        
        if conc != 0:
            conc = np.log10(conc)/10
            for file in os.listdir(os.path.join(d_folder, folder)):
                if file.endswith('s.txt'):
                                        
                    fi = os.path.join(d_folder, folder, file)
                    f, re, im = np.loadtxt(fi, skiprows=1, unpack=True)
                    
                    # # Sometimes 1Hz re is < 0 due to noise, throw out that data file
                    if not any(([i<0 for i in re] + [i>0 for i in im])):
                        
                        data.append([np.hstack((np.log10(re)/10, np.log10(-im)/10))])
                        
                        if file.endswith('0002s.txt'):
                            test_files.append(os.path.join(folder, file))
                            test_concs.append(conc)
                            
                        else:
                            files.append(os.path.join(folder, file))
                            concs.append(conc)
                            
    
                
    data  = np.array(data)
    concs = np.array(concs)
    
    return data, concs, files, test_files, test_concs



def extract_invivo(d_folder):
    times = []
    data  = []
    files = []
    
    time_file = os.path.join(d_folder, '0000_time_list.txt')
    ts = np.loadtxt(time_file)
    
    j = 0
    for file in os.listdir(d_folder):
        if file.endswith('s.txt'):
            file = os.path.join(d_folder, file)
            
            f, re, im = np.loadtxt(file, skiprows=1, unpack=True)
            
            if not any(([i<0 for i in re] + [i>0 for i in im])):
                times.append(ts[j])
                files.append(file)
                data.append([np.hstack((np.log10(re)/10, np.log10(-im)/10))])
        
            j += 1
            
    
    return data, times, files


data, concs, files, test_files, test_concs = extract_titration(data_folder)

df = pd.DataFrame(np.array([files, concs]).T)
df.to_csv(train_csv, header=False, index=False)

df = pd.DataFrame(np.array([test_files, test_concs]).T)
df.to_csv(test_csv, header=False, index=False)


data, times, files = extract_invivo(invivo_folder)
df = pd.DataFrame(np.array([files, times]).T)
df.to_csv(invivo_csv, header=False, index=False)



class CustomDataset(Dataset):
    def __init__(self, annotations_file, data_dir, transform=None, 
                 target_transform=None):
        self.data_labels = pd.read_csv(annotations_file)
        self.data_dir = data_dir
        self.transform = transform
        self.target_transform = target_transform


    def __len__(self):
        return len(self.data_labels)


    def __getitem__(self, idx):
        data_path = os.path.join(self.data_dir, self.data_labels.iloc[idx, 0])
        f, re, im = np.loadtxt(data_path, unpack=True, skiprows=1)
        data = np.vstack((np.log10(re)/10, np.log10(-im)/10))        
        # data = np.vstack((re, im))
        data = torch.from_numpy(data).double()
        label = self.data_labels.iloc[idx, 1]
        return data, label
    


class NeuralNetwork(nn.Module):
    def __init__(self):
        super(NeuralNetwork, self).__init__()
        self.flatten = nn.Flatten()
        self.linear_relu_stack = nn.Sequential(
            nn.Linear(44, 90),
            nn.ReLU(),
            nn.Linear(90,50),
            nn.ReLU(),
            nn.Linear(50, 1),
            nn.Flatten(0,1)
        )
        self.double()

    def forward(self, x):
        x = self.flatten(x)
        logits = self.linear_relu_stack(x)
        return logits
    

    
def train_loop(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    for batch, (X, y) in enumerate(dataloader):
        # Compute prediction and loss
        pred = model(X)
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch % 100 == 0:
            loss, current = loss.item(), batch * len(X)
            # print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")



def test_loop(dataloader, model, loss_fn, losses=[], p=False, epoch=0):
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    test_loss, correct = 0, 0

    with torch.no_grad():
        for X, y in dataloader:
            pred = model(X)
            test_loss += loss_fn(pred, y).item()
            # correct += (pred.argmax(1) == y).type(torch.float).sum().item()

    test_loss /= num_batches
    correct /= size
    losses.append(test_loss)
    
    if p:
        print(f'Epoch {epoch}. loss: {test_loss}')
    
    return losses





training_dataset = CustomDataset(train_csv, data_folder)
test_dataset = CustomDataset(test_csv, data_folder)
invivo_dataset = CustomDataset(invivo_csv, invivo_folder)

train_dataloader = DataLoader(training_dataset, 
                              batch_size=50, 
                              shuffle=True)
test_dataloader = DataLoader(test_dataset, 
                             batch_size=len(test_dataset), 
                             shuffle=False)
invivo_dataloader = DataLoader(invivo_dataset, 
                               batch_size=len(invivo_dataset), 
                               shuffle=False)


# model = torch.load('model.pth')




model = NeuralNetwork().to('cpu')

learning_rate = 1e-3

# # loss_fn = nn.CrossEntropyLoss()
loss_fn = nn.MSELoss(reduction='sum')
optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)

epochs = 5000
losses = []
for t in range(epochs):
    train_loop(train_dataloader, model, loss_fn, optimizer)
    if t%100 == 0:
        losses = test_loop(test_dataloader, model, loss_fn, losses,
                            p=True, epoch=t)
    else:
        losses = test_loop(test_dataloader, model, loss_fn, losses)
print("Done!")


def test(dataloader, model):
    with torch.no_grad():
        for X, y in dataloader:
            pred = model(X)
            plt.plot(10**(10*y), color='k')
            plt.plot(10**(10*pred), color='red')
            plt.yscale('log')

test(test_dataloader, model)


fig, ax = plt.subplots()
ax.plot(losses)
ax.set_yscale('log')
ax.set_ylabel('log loss')

# torch.save(model.state_dict(), 'model_weights.pth')
# torch.save(model, 'model.pth')

# with torch.no_grad():
#     for X, y in invivo_dataloader:
#         print(X)
#         pred = model(X)
#         plt.plot(10**(10*pred))




