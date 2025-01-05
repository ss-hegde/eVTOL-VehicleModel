import torch # type: ignore
import torch.nn as nn   # type: ignore

# LSTM Model Definition - static wing data 
class LSTMNet_static(nn.Module):
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