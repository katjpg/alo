import os
import sys
import time
import json
from datetime import datetime

# import fire
import torch
from peft import PeftModel
from transformers import AutoModelForCausalLM, AutoTokenizer, pipeline, BitsAndBytesConfig
from datasets import Dataset, load_dataset

from utils.prompter import Prompter
from utils.config import DEFAULT_MAX_NEW_TOKENS, TARGET_TASKS

os.environ["TOKNIZERS_PARALLELISM"] = "false"

device = "cuda" if torch.cuda.is_available() else "cpu"

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
    response_split = "%%% Response:"
    prompter = Prompter(opt_type, response_split)

    # ✅ Task mapping (if not registered, use as-is)
    task = TARGET_TASKS.get(task, task)

    # ✅ Handle local JSON or HuggingFace dataset
    if data_path.endswith(".json"):
        print(f"🔹 Loading local JSON file: {data_path}")
        with open(data_path, "r") as f:
            data = json.load(f)
        test_data = Dataset.from_list(data["test"])
    else:
        print(f"🔹 Loading HuggingFace Dataset: {data_path}")
        test_data = load_dataset(data_path)["test"]

    # ✅ Filtering
    test_data = test_data.filter(
        lambda example: example['property_comb'] == task and example['instr_setting'] == setting
    )

    print(f"✅ Loaded {len(test_data)} samples for task '{task}' with setting '{setting}'")

    # ✅ Handle model name
    if base_model == 'Llama-3.1-8B-Instruct' or base_model == 'Llama-3.1-70B-Instruct':
        base_model = f'meta-llama/{base_model}'
    elif base_model == 'Mistral-7B-Instruct-v0.3':
        base_model = f'mistralai/{base_model}'

    if "Llama-3.1-70B" in base_model:
        print("⚡ Using 70B model, forcing 8-bit loading")
        load_in_8bit = True
        batch_size = 4

    # ✅ Get Hugging Face token
    hf_token = os.getenv("HF_TOKEN")
    if not hf_token:
        print("❌ Hugging Face token not set.")
        print("✅ In PowerShell, enter the following:")
        print('$env:HF_TOKEN="hf_xxxYourTokenxxx"')
        sys.exit(1)

    # ✅ Generate prompts
    instructions = [
        prompter.generate_prompt(data_point, prompt_type='ins', add_response=False) for data_point in test_data
    ]
    print(f"🔹 Example prompt:\n{instructions[0]}")

    # ✅ Load model
    tokenizer = AutoTokenizer.from_pretrained(base_model, trust_remote_code=True, token=hf_token)
    quant_config = BitsAndBytesConfig(load_in_8bit=load_in_8bit)
    model = AutoModelForCausalLM.from_pretrained(
        base_model,
        torch_dtype=torch.bfloat16,
        quantization_config=quant_config,
        device_map="auto",
        trust_remote_code=True,
        token=hf_token
    )

    if use_lora:
        model = PeftModel.from_pretrained(model, lora_weights, torch_dtype=torch.bfloat16, device_map="auto")

    if not load_in_8bit:
        model.bfloat16()

    model.eval()
    if torch.__version__ >= "2" and sys.platform != "win32":
        model = torch.compile(model)

    tokenizer.pad_token_id = model.config.eos_token_id
    tokenizer.padding_side = 'left'

    pipe = pipeline(
        "text-generation",
        model=model,
        tokenizer=tokenizer,
        torch_dtype=torch.float16,
        device_map="auto"
    )

    results = []
    for i in range(0, len(instructions), batch_size):
        batch = instructions[i:i + batch_size]
        print(f"🚀 Processing batch {i // batch_size + 1}/{(len(instructions) - 1) // batch_size + 1}")
        start_time = time.time()

        outputs = evaluate(prompter, batch, tokenizer, pipe, batch_size, num_return_sequences)
        results.extend(outputs)

        elapsed = time.time() - start_time
        print(f"✅ Batch {i // batch_size + 1} completed in {elapsed:.2f} seconds")

    # ✅ 결과 저장
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(output_dir, f"{task}_response_{timestamp}.json")
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)

    return output_file


def evaluate(prompter, prompts, tokenizer, pipe, batch_size, num_return_sequences=1,
             max_new_tokens=DEFAULT_MAX_NEW_TOKENS, do_sample=False, temperature=1):
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
    batch_results = []
    for i, result in enumerate(outputs):
        responses = [prompter.get_response(r['generated_text']) for r in result]
        batch_results.append({'prompt': prompts[i], 'response': responses})
    return batch_results


# if __name__ == "__main__":
#     torch.cuda.empty_cache()
#     fire.Fire(main)



    # // "instructions": [
    # //     "Adjust the molecular structure to ensure that each specified property reaches the corresponding threshold listed in <THRESHOLD> </THRESHOLD>. Minimize structural changes and try to maintain the core scaffold. Return the resulting molecule using <SMILES> </SMILES> tags. You must keep the input <SMILES>.",
    # //     "Alter the molecule to satisfy the provided property thresholds in <THRESHOLD> </THRESHOLD>. Preserve the core scaffold and make as few structural changes as possible. Output the SMILES of the new molecule, enclosed in <SMILES> </SMILES>. Keep the input <SMILES> as much as possible.",
    # //     "Update the given molecule so that the specified properties fall within acceptable ranges defined by the values in <THRESHOLD> </THRESHOLD>. Maintain as much of the original structure as possible. Output only the modified molecule enclosed in <SMILES> </SMILES> tags. Keep the input <SMILES> as much as possible.",
    # //     "Modify the molecule to bring its properties to at least the levels defined in <THRESHOLD> </THRESHOLD>. Avoid excessive modifications and preserve the core scaffold. Output only the resulting molecule's SMILES wrapped in <SMILES> </SMILES>. Keep the input <SMILES> as much as possible.",
    # //     "Edit the molecular structure so that all required properties match or exceed the threshold values defined in <THRESHOLD> </THRESHOLD>. Try to retain the core scaffold. Output only the SMILES representation of the optimized molecule enclosed in <SMILES> </SMILES>. Keep the input <SMILES> as much as possible."
    # // ],