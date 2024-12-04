import requests
import os
from pathlib import Path
import tarfile
import sys
from artic.utils import clair3_manifest
from clint.textui import colored


def download_file(url: str, local_path: Path):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


def get_model(model_dir: Path, model_fname: str, model_url: str):

    model_path = Path(model_dir, model_fname)

    os.makedirs(model_dir, exist_ok=True)

    download_file(model_url, model_path)

    with tarfile.open(model_path, "r") as tar:
        tar.extractall(model_dir)

    os.remove(model_path)

    return model_path


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=f"{os.getenv('CONDA_PREFIX')}/bin/models/",
        help="Directory to download the model to, default is: %(default)s",
    )
    args = parser.parse_args()

    if not os.getenv("CONDA_PREFIX"):
        print(
            f"CONDA_PREFIX is not set, this probably means you are not running this inside a conda environment, if you have not provided a model path argument '--model-dir' the models might be downloaded somewhere you don't want them to be.",
            file=sys.stderr,
        )

    model_manifest = clair3_manifest()
    models = model_manifest.models

    for model in models:
        if not model["rerio"]:
            continue

        if not os.path.exists(Path(args.model_dir, model["name"])):
            get_model(
                model_dir=args.model_dir,
                model_fname=model["model_fname"],
                model_url=model["model_url"],
            )
            print(f"Downloaded model: {model['name']}", file=sys.stderr)

        else:
            print(
                f"Model {model['name']} already exists, skipping download",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
