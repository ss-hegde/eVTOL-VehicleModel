import torch.nn as nn

# LSTM Model Definition
class LSTMNet_rotor(nn.Module):
    """
    A Long Short-Term Memory (LSTM) network for rotor modeling.

    Args:
        input_size (int): The number of expected features in the input.
        hidden_size (int): The number of features in the hidden state.
        output_size (int): The number of features in the output.
        num_layers (int): The number of recurrent layers in the LSTM.

    Attributes:
        lstm (nn.LSTM): The LSTM layer.
        fc (nn.Linear): The fully connected layer that maps the LSTM output to the desired output size.

    Methods:
        forward(x):
            Defines the forward pass of the network.
            Args:
                x (torch.Tensor): The input tensor of shape (batch_size, sequence_length, input_size).
            Returns:
                torch.Tensor: The output tensor of shape (batch_size, sequence_length, output_size).
    """
    def __init__(self, input_size, hidden_size, output_size, num_layers):
        super(LSTMNet_rotor, self).__init__()
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=True)
        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        lstm_out, (hn, cn) = self.lstm(x)
        out = self.fc(lstm_out)
        return out