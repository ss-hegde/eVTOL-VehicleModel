"""
Element Spatial Convolutional Neural Network model
Number of convolutional layers - 4
Number of fully connected layers - 2
"""

import torch    # type:ignore
import torch.nn as nn   # type:ignore

# Model for lift coefficient
class ESCNN_Cl(nn.Module):
    def __init__(self):
        super(ESCNN_Cl, self).__init__()
        
        # Conv1: Assume 1D Convolution
        self.conv1 = nn.Conv1d(in_channels=1, out_channels=200, kernel_size=5, stride=5)
        self.relu1 = nn.ReLU()

        #conv2
        self.conv2 = nn.Conv1d(in_channels=200, out_channels=200, kernel_size=1)
        self.relu2 = nn.ReLU()

        # Conv3
        self.conv3 = nn.Conv1d(in_channels=200, out_channels=100, kernel_size=1)
        self.relu3 = nn.ReLU()
        
        # Conv4
        self.conv4 = nn.Conv1d(in_channels=100, out_channels=1, kernel_size=5, padding=2)
        self.relu4 = nn.ReLU()
        
        # Final fully connected layer to output scalar
        self.fc1 = nn.Linear(in_features=34, out_features=34)
        self.relu5 = nn.ReLU()

        self.fc2 = nn.Linear(in_features=34, out_features=1)
    
    def forward(self, x):
        # Reshape input if necessary, ensure it's in the shape (batch_size, channels, elements)
        x = x.view(-1, 1, 170)  # Reshape to (batch_size, channel=1, elements=170)
        
        x = self.conv1(x)  
        x = self.relu1(x)
        
        x = self.conv2(x)  
        x = self.relu2(x)

        x = self.conv3(x)  
        x = self.relu3(x)
        
        x = self.conv4(x)  
        x = self.relu4(x)
        
        x = torch.flatten(x, 1)
        
        x = self.fc1(x)  
        x = self.relu5(x)

        x = self.fc2(x)
        
        return x