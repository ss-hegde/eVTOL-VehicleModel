import torch # type: ignore
import torch.nn as nn   # type: ignore

# LSTM Model Definition - static wing data 
class LSTMNet_static(nn.Module):
    """
    A static LSTM network model for sequence prediction.

    Args:
        input_size (int): The number of expected features in the input `x`.
        hidden_size (int): The number of features in the hidden state `h`.
        output_size (int): The number of features in the output.
        num_layers (int): Number of recurrent layers.

    Attributes:
        lstm (nn.LSTM): LSTM layer for processing input sequences.
        fc (nn.Linear): Fully connected layer for transforming LSTM output.
        fc2 (nn.Linear): Additional fully connected layer for further transformation.

    Methods:
        forward(x):
            Defines the forward pass of the network.
            Args:
                x (torch.Tensor): The input tensor of shape (batch_size, sequence_length, input_size).
            Returns:
                torch.Tensor: The output tensor of shape (batch_size, sequence_length, output_size).
    """
    def __init__(self, input_size, hidden_size, output_size, num_layers):
        super(LSTMNet_static, self).__init__()
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_size, output_size)
        self.fc2 = nn.Linear(output_size, output_size)

    def forward(self, x):
        lstm_out, (hn, cn) = self.lstm(x)
        out = self.fc(lstm_out)
        out2 = self.fc2(out)
        return out2