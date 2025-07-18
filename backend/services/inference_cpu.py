import os
import sys
import time
import json
from datetime import datetime

import torch
from tqdm import tqdm
from peft import PeftModel
from transformers import AutoModelForCausalLM, AutoTokenizer, pipeline, BitsAndBytesConfig
from datasets import Dataset, load_dataset

from utils.prompter import Prompter
from utils.config import DEFAULT_MAX_NEW_TOKENS, TARGET_TASKS

os.environ["TOKNIZERS_PARALLELISM"] = "false"

cuda_available = torch.cuda.is_available()
device = "cuda" if cuda_available else "cpu"
torch_type = torch.bfloat16 if cuda_available else torch.float32


def generate(
    load_in_8bit: bool = False,
    use_lora: bool = True,
    base_model: str = "",
    lora_weights: str = "",
    opt_type: str = "thresh",
    task: str = "",
    setting: str = "seen",
    data_path: str = "",
    output_dir: str = "output",
    batch_size: int = 16,
    num_return_sequences: int = 1,
):
    print (torch_type)
    response_split = "%%% Response:"
    prompter = Prompter(opt_type, response_split)

    # Task mapping (if not registered, use as-is)
    task = TARGET_TASKS.get(task, task)

    # Handle local JSON or HuggingFace dataset
    if data_path.endswith(".json"):
        print(f"Loading local JSON file: {data_path}")
        with open(data_path, "r") as f:
            data = json.load(f)
        test_data = Dataset.from_list(data["test"])
    else:
        print(f"Loading HuggingFace Dataset: {data_path}")
        test_data = load_dataset(data_path)["test"]

    # Filtering
    test_data = test_data.filter(
        lambda example: example['property_comb'] == task and example['instr_setting'] == setting
    )

    print(f"Loaded {len(test_data)} samples for task '{task}' with setting '{setting}'")

    # Handle model name
    if base_model == 'Llama-3.1-8B-Instruct' or base_model == 'Llama-3.1-70B-Instruct':
        base_model = f'meta-llama/{base_model}'
    elif base_model == 'Mistral-7B-Instruct-v0.3':
        base_model = f'mistralai/{base_model}'

    if "Llama-3.1-70B" in base_model:
        print("⚡ Using 70B model, forcing 8-bit loading")
        load_in_8bit = True
        batch_size = 4

    # Get Hugging Face token
    hf_token = os.getenv("HF_TOKEN")
    if not hf_token:
        print("Hugging Face token not set.")
        print("In PowerShell, enter the following:")
        print('$env:HF_TOKEN="hf_xxxYourTokenxxx"')
        sys.exit(1)

    # Generate prompts
    instructions = [
        prompter.generate_prompt(data_point, prompt_type='ins', add_response=False) for data_point in test_data
    ]
    print(f"Example prompt:\n{instructions[0]}")

    # Load model
    tokenizer = AutoTokenizer.from_pretrained(base_model, trust_remote_code=True, token=hf_token)
    print ("Tokenizer complete")
    quant_config = BitsAndBytesConfig(load_in_8bit=load_in_8bit) if cuda_available else None
    model = AutoModelForCausalLM.from_pretrained(
        base_model,
        torch_dtype=torch_type,
        quantization_config=quant_config,
        device_map = "auto" if cuda_available else {"": "cpu"},
        trust_remote_code=True,
        token=hf_token
    )
    if not cuda_available: 
        model.to("cpu")
    print ("Model complete")

    if use_lora and device == "cuda":
        print(f"Backend: {device}, Dtype: {torch_type}, 8bit: {load_in_8bit}")
        model = PeftModel.from_pretrained(
            model, 
            lora_weights, 
            torch_dtype=torch_type, 
            device_map=device
        )
        print ("Lora complete")

    if cuda_available and not load_in_8bit:
        model.bfloat16()

    model.eval()
    if cuda_available and torch.__version__ >= "2" and sys.platform != "win32":
        model = torch.compile(model)

    tokenizer.pad_token_id = model.config.eos_token_id
    tokenizer.padding_side = 'left'

    pipe = pipeline(
        "text-generation",
        model=model,
        tokenizer=tokenizer,
        torch_dtype=torch_type,
        device_map=0 if cuda_available else -1
    )

    print ("Pipeline complete")

    results = []
    for i in range(0, len(instructions), batch_size):
        batch = instructions[i:i + batch_size]
        print(f"Processing batch {i // batch_size + 1}/{(len(instructions) - 1) // batch_size + 1}")
        start_time = time.time()

        outputs = evaluate(prompter, batch, tokenizer, pipe, batch_size, num_return_sequences)
        results.extend(outputs)

        elapsed = time.time() - start_time
        print(f"Batch {i // batch_size + 1} completed in {elapsed:.2f} seconds")

    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"{task}_response_{timestamp}.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    return output_file


def evaluate(prompter, prompts, tokenizer, pipe, batch_size, num_return_sequences=1,
             max_new_tokens=DEFAULT_MAX_NEW_TOKENS, do_sample=False, temperature=1):
    print ("Pipe started")
    if cuda_available : 
        outputs = pipe(
            prompts,
            do_sample=do_sample,
            max_new_tokens=max_new_tokens,
            temperature=temperature,
            num_return_sequences=num_return_sequences,
            num_beams=num_return_sequences,
            pad_token_id=tokenizer.eos_token_id,
            batch_size=batch_size
        )
    else : 
        outputs = pipe(
            prompts,
            do_sample=do_sample,
            max_new_tokens=32,
            temperature=temperature,
            num_return_sequences=num_return_sequences,
            num_beams=num_return_sequences,
            pad_token_id=tokenizer.eos_token_id,
        )

    print ("Pipe complete")
    batch_results = []
    for i, result in enumerate(outputs):
        responses = [prompter.get_response(r['generated_text']) for r in result]
        batch_results.append({'prompt': prompts[i], 'response': responses})
    return batch_results