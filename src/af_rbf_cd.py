import torch # type: ignore
import torch.nn as nn   # type: ignore
import torch.nn.functional as F #type: ignore

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Cd Model
class RBFLayer_cd(nn.Module):
    """
    Radial Basis Function (RBF) Layer with centers determined by KMeans or explicitly provided.
    Args:
        in_features (int): Number of input features.
        out_features (int): Number of output features (number of RBF centers).
        centers (Tensor, optional): Predefined centers for the RBFs. If None, centers are initialized randomly.
    Attributes:
        in_features (int): Number of input features.
        out_features (int): Number of output features (number of RBF centers).
        centers (Parameter): Tensor containing the centers of the RBFs.
        log_beta (Parameter): Logarithm of the beta (width) parameter for the RBFs.
    Methods:
        beta:
            Property that returns the positive beta (width) parameter using the softplus function.
        forward(input):
            Computes the output of the RBF layer for the given input.
            Args:
                input (Tensor): Input tensor of shape (batch_size, in_features).
            Returns:
                Tensor: Output tensor of shape (batch_size, out_features) containing the RBF activations.
    """
    def __init__(self, in_features, out_features, centers=None):
        super(RBFLayer_cd, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        
        # Centers initialized using KMeans or passed explicitly
        if centers is None:
            self.centers = nn.Parameter(torch.Tensor(out_features, in_features))
            nn.init.uniform_(self.centers, -1, 1)  # Random if not initialized with KMeans
        else:
            self.centers = nn.Parameter(torch.Tensor(centers))  # Set centers from KMeans
        
        # Initialize the beta (width) parameter, positive constraint with softplus
        self.log_beta = nn.Parameter(torch.ones(out_features) * torch.log(torch.tensor(18.2)))

    @property
    def beta(self):
        # Ensures beta is positive using softplus
        return F.softplus(self.log_beta)
    
    def forward(self, input):
        # Compute distances between inputs and centers
        x = input.unsqueeze(1).expand(-1, self.out_features, self.in_features)
        c = self.centers.unsqueeze(0).expand(input.size(0), -1, -1)
        
        # Squared Euclidean distance
        distances = torch.sum((x - c) ** 2, dim=-1).to(device)
        
        # Apply Gaussian RBF
        return torch.exp(-self.beta.unsqueeze(0) * distances)

class RBFNet_cd(nn.Module):
    def __init__(self, input_size, rbf_units, output_size, centers=None):
        super(RBFNet_cd, self).__init__()
        self.rbf = RBFLayer_cd(input_size, rbf_units, centers)
        self.fc = nn.Linear(rbf_units, output_size)
        
        # Initialize weights for the linear layer
        nn.init.xavier_uniform_(self.fc.weight)
        nn.init.zeros_(self.fc.bias)
        
    def forward(self, x):
        rbf_out = self.rbf(x)
        return self.fc(rbf_out)
     