#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/sample_reads.py

from pathlib import Path

import torch
from smf_net.model import SMFNet
from smf_net.evaluation import (
    SamplingConfig, batch_regions, regions_from_bed, regions_from_fasta,
    sample_regions, write_samples_h5,
)

#Set paths and input type before running; reference FASTA is needed only for BED.
CHECKPOINT = None
SOURCE = None
INPUT_FILES = []
REFERENCE_FASTA = None
OUTPUT_DIR = Path("output/read_sampling/samples")
DEVICE = "cuda"
COMPUTE_DTYPE = torch.bfloat16
READS_PER_REGION = 512
READS_PER_CHUNK = 512
REGIONS_PER_BATCH = 32
SAMPLE_STEPS = 32
SEED = 0


def main():
    if CHECKPOINT is None:
        raise ValueError("Set CHECKPOINT to Path(\"trained_models/your_model.ckpt\")")
    if SOURCE not in {"fasta", "bed"}:
        raise ValueError("SOURCE must be fasta or bed")
    if not INPUT_FILES or len({p.stem for p in INPUT_FILES}) != len(INPUT_FILES):
        raise ValueError("Provide input files with distinct stems")
    required = [CHECKPOINT, *INPUT_FILES]
    if SOURCE == "bed":
        if REFERENCE_FASTA is None:
            raise ValueError("Set REFERENCE_FASTA to the genome FASTA for BED input")
        required.append(REFERENCE_FASTA)
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    for path in INPUT_FILES:
        destination = OUTPUT_DIR / path.stem / "samples.h5"
        if destination.exists():
            raise FileExistsError(destination)

    checkpoint = torch.load(CHECKPOINT, map_location="cpu", weights_only=False)
    model = SMFNet(**checkpoint["hyper_parameters"])
    #Compiled checkpoints carry an extra wrapper in their parameter names.
    state = {key.replace("_orig_mod.", ""): value
             for key, value in checkpoint["state_dict"].items()}
    model.load_state_dict(state, strict=True)
    model = model.to(device=DEVICE, dtype=torch.float32).eval()
    config = SamplingConfig(
        reads_per_region=READS_PER_REGION, reads_per_chunk=READS_PER_CHUNK,
        regions_per_batch=REGIONS_PER_BATCH, sample_steps=SAMPLE_STEPS, seed=SEED,
    )
    for path in INPUT_FILES:
        if SOURCE == "fasta":
            regions = regions_from_fasta(path, sequence_length=model.cfg.input_length)
        else:
            regions = regions_from_bed(
                path, REFERENCE_FASTA, sequence_length=model.cfg.input_length,
            )
        destination = OUTPUT_DIR / path.stem / "samples.h5"
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(".h5.tmp")
        metadata = {
            "checkpoint": str(CHECKPOINT), "input": str(path), "source": SOURCE,
            "seed": SEED, "reads_per_region": READS_PER_REGION,
            "reads_per_chunk": READS_PER_CHUNK, "regions_per_batch": REGIONS_PER_BATCH,
            "sample_steps": SAMPLE_STEPS, "device": DEVICE,
            "dtype": str(COMPUTE_DTYPE), "parameter_dtype": "float32",
            "model_input_length": model.cfg.input_length,
        }
        if SOURCE == "bed":
            metadata["reference_fasta"] = str(REFERENCE_FASTA)
        with torch.inference_mode():
            batches = sample_regions(
                model, batch_regions(regions, REGIONS_PER_BATCH), config,
                total_regions=regions.total_regions, compute_dtype=COMPUTE_DTYPE,
            )
            count = write_samples_h5(temporary, batches, metadata=metadata)
        temporary.replace(destination)
        print(f"Wrote {count} regions to {destination}")


if __name__ == "__main__":
    main()
