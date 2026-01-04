
from pathlib import Path

import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM


def load_nttx():
    model_path = Path("/data/nvme/models/nucleotide-transformer-500m-human-ref").resolve()
    
    tokenizer = AutoTokenizer.from_pretrained(model_path, trust_remote_code = True)
    model = AutoModelForMaskedLM.from_pretrained(model_path, trust_remote_code = True)
    
    return model, tokenizer

def run_nttx(model, tkn, sequence, max_length = 128):
    
    
    tokens_ids = tkn.encode(sequence, return_tensors="pt", padding="max_length", max_length = max_length)
    
    attention_mask = tokens_ids != tkn.pad_token_id
    torch_outs = model(
        tokens_ids,
        attention_mask=attention_mask,
        encoder_attention_mask=attention_mask,
        output_hidden_states=True
    )
    logits = torch_outs['logits'][0]
    probs = (torch.tanh(logits/2) + 1)/2
    choices = torch.argmax(probs, dim = -1).to(dtype=torch.int)
    # max_probs = [probs[i, choices[i]] for i in range(probs.shape[0])]
    
    inf = tkn.decode(choices)
    
    entrs = torch.sum(-probs*torch.log(probs), -1)
    entrs[torch.isnan(entrs)] = 0.0
    
    # seqs_inf = ["".join(s).removeprefix("<cls> ").replace(" ","") for s in inf]
    seqs_inf = "".join(inf)
    return seqs_inf, entrs

def run_nttx_batch(model, tkn, sequences, max_length = 128):
    
    
    tokens_ids = tkn.batch_encode_plus(sequences, return_tensors="pt", padding="max_length", max_length = max_length)["input_ids"]
    
    attention_mask = tokens_ids != tkn.pad_token_id
    torch_outs = model(
        tokens_ids,
        attention_mask=attention_mask,
        encoder_attention_mask=attention_mask,
        output_hidden_states=True
    )
    logits = torch_outs['logits']
    choices = torch.argmax(logits, dim = -1).to(dtype=torch.int)

    inf = tkn.batch_decode(choices)
    
    entrs = torch.sum(logits[...,choices]*torch.exp(logits[...,choices]), -1)
    seqs_inf = ["".join(s) for s in inf]
    return seqs_inf, entrs


