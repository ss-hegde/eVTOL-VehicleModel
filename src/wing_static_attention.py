import torch
import torch.nn as nn

# LSTM Model Definition - static wing data 
class LSTMNetWithAttention(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, num_layers, num_heads=2, bias=0.0):
        super(LSTMNetWithAttention, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.bias = nn.Parameter(torch.tensor(bias))  # Learnable bias term

        # LSTM layer
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, batch_first=True)

        # Multi-head attention layer
        self.multihead_attention = nn.MultiheadAttention(hidden_size, num_heads, batch_first=True)

        # Fully connected output layer
        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        """
        x: Input tensor of shape (batch_size, seq_length, input_size)
        """
        # Pass through LSTM
        lstm_out, _ = self.lstm(x)  # lstm_out: (batch_size, seq_length, hidden_size)

        # Multi-head attention: Query, Key, and Value are all lstm_out
        attention_output, attention_weights = self.multihead_attention(lstm_out, lstm_out, lstm_out)

        # Pass attention output through the fully connected layer
        output = self.fc(attention_output)  # (batch_size, seq_length, output_size)

        # Add bias to the output
        output = output + self.bias

        return output