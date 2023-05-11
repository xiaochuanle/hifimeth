import torch 
import torch.nn as nn
import attention
from attention import Attention

# Using BS-seq modfreqs as gold standard to train a AggrAttRNN regression model, no softmax
class AggrAttRNN(nn.Module):
    def __init__(self, 
                kmer_size,
                num_features,
                num_layers = 3,
                dropout_rate = 0,
                hidden_size = 32,
                device = "cpu"):
        super(AggrAttRNN, self).__init__()
        self.model_type = "attbigru"
        self.device = device

        self.seq_len = kmer_size
        self.num_layers = num_layers
        self.num_classes = 1
        self.hidden_size = hidden_size

        self.num_features = num_features
        self.feas_ccs = num_features
        if self.model_type == "attbilstm":
            self.rnn_cell = "lstm"
            self.rnn = nn.LSTM(self.feas_ccs, self.hidden_size, self.num_layers,
                               dropout=0, batch_first=True, bidirectional=True)
        elif self.model_type == "attbigru":
            self.rnn_cell = "gru"
            self.rnn = nn.GRU(self.feas_ccs, self.hidden_size, self.num_layers,
                              dropout=0, batch_first=True, bidirectional=True)
        else:
            raise ValueError("--model_type not set right!")

        self._att3 = Attention(self.hidden_size * 2, self.hidden_size * 2, self.hidden_size)

        self.dropout1 = nn.Dropout(p=dropout_rate)
        self.fc1 = nn.Linear(self.hidden_size * 2, self.num_classes)  # 2 for bidirection

    def get_model_type(self):
        return self.model_type

    def init_hidden(self, batch_size, num_layers, hidden_size):
        # Set initial states
        h0 = torch.randn(num_layers * 2, batch_size, hidden_size, requires_grad=True)
        if self.device != torch.device("cpu"):
            h0 = h0.cuda(self.device)
        if self.rnn_cell == "lstm":
            c0 = torch.randn(num_layers * 2, batch_size, hidden_size, requires_grad=True)
            if self.device != torch.device("cpu"):
                c0 = c0.cuda(self.device)
            return h0, c0
        return h0

    def forward(self, fwd_x):

        out, n_states = self.rnn(fwd_x, self.init_hidden(fwd_x.size(0),
                                                       self.num_layers,
                                                       self.hidden_size))  # (N, L, nhid*2)
        # attention_net3 ======
        # h_n: (num_layer * 2, N, nhid), h_0, c_0 -> h_n, c_n not affected by batch_first
        # h_n (last layer) = out[:, -1, :self.hidden_size] concats out1[:, 0, self.hidden_size:]
        h_n = n_states[0] if self.rnn_cell == "lstm" else n_states
        h_n = h_n.reshape(self.num_layers, 2, -1, self.hidden_size)[-1]  # last layer (2, N, nhid)
        h_n = h_n.transpose(0, 1).reshape(-1, 1, 2 * self.hidden_size)
        out, att_weights = self._att3(h_n, out)

        out = self.dropout1(out)
        out = self.fc1(out)
        out = torch.sigmoid(out)

        return out