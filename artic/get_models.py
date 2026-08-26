import requests
import os
from pathlib import Path
import tarfile
import sys
import shutil
from clint.textui import colored
from artic.utils import CLAIR3_MANIFEST, CLAIR3_MODEL_BACKUP_URL


def download_file(url: str, local_path: Path):
    with requests.get(url, stream=True, timeout=(10, 60)) as r:
        r.raise_for_status()
        with open(local_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)


def get_model(
    model_dir: Path, model_fname: str, model_url: str, model_name: str = None
):

    model_path = Path(model_dir, model_fname)

    os.makedirs(model_dir, exist_ok=True)

    download_file(model_url, model_path)

    with tarfile.open(model_path, "r") as tar:
        paths = [Path(x) for x in tar.getnames()]
        root_paths = [str(x.parent) for x in paths if str(x.parent) != "."]
        if len(set(root_paths)) != 1:
            raise ValueError(
                f"The Clair3 model tarfile {model_fname} contains multiple root directories (there can only be one), please check the tar file."
            )

        tar.extractall(model_dir)

        if model_name:
            if root_paths[0] != model_name:
                shutil.move(Path(model_dir, root_paths[0]), Path(model_dir, model_name))

    os.remove(model_path)

    return model_path


def get_pytorch_model(model_dir: Path, model_name: str, model_url: str):
    model_path = Path(model_dir, model_name)
    os.makedirs(model_path, exist_ok=True)

    base_url = model_url.rstrip("/")
    for fname in ["pileup.pt", "full_alignment.pt"]:
        download_file(f"{base_url}/{fname}", Path(model_path, fname))

    return model_path


def _model_backup_url(model: dict, backup_base_url: str) -> str:
    return f"{backup_base_url.rstrip('/')}/{model['name']}/"


def _cleanup_partial_download(model_dir: Path, model: dict):
    if model.get("pytorch"):
        shutil.rmtree(Path(model_dir, model["name"]), ignore_errors=True)
    else:
        Path(model_dir, model["model_fname"]).unlink(missing_ok=True)


def fetch_model(model: dict, model_dir: Path, backup_base_url: str = None):
    """Download a clair3 model, trying model["model_url"] first and falling
    back to backup_base_url (if given) if the primary URL fails."""

    candidate_urls = [model["model_url"]]
    if backup_base_url:
        candidate_urls.append(_model_backup_url(model, backup_base_url))

    last_error = None
    for i, model_url in enumerate(candidate_urls):
        try:
            if model.get("pytorch"):
                get_pytorch_model(
                    model_dir=model_dir,
                    model_name=model["name"],
                    model_url=model_url,
                )
            else:
                get_model(
                    model_dir=model_dir,
                    model_fname=model["model_fname"],
                    model_url=model_url,
                    model_name=model["name"],
                )
            return
        except (requests.RequestException, tarfile.TarError, ValueError) as e:
            last_error = e
            _cleanup_partial_download(model_dir, model)
            is_last_candidate = i == len(candidate_urls) - 1
            print(
                colored.yellow(
                    f"Failed to download model '{model['name']}' from {model_url}: {e}"
                ),
                file=sys.stderr,
            )
            if not is_last_candidate:
                print(
                    colored.yellow(f"Retrying '{model['name']}' from backup source"),
                    file=sys.stderr,
                )

    raise RuntimeError(
        f"Failed to download model '{model['name']}' from all sources: {candidate_urls}"
    ) from last_error


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=f"{os.getenv('CONDA_PREFIX')}/bin/models/",
        help="Directory to download the model to, default is: %(default)s",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        metavar="MODEL",
        help="Only download the specified model(s) by name. Downloads all models if omitted.",
    )
    parser.add_argument(
        "--backup-url",
        default=CLAIR3_MODEL_BACKUP_URL,
        help="Base URL of a backup model source, tried if the primary model URL fails. Default is: %(default)s",
    )
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Do not fall back to the backup model source if the primary download fails.",
    )
    args = parser.parse_args()

    backup_base_url = None if args.no_backup else args.backup_url

    if not os.getenv("CONDA_PREFIX"):
        print(
            "CONDA_PREFIX is not set, this probably means you are not running this inside a conda environment, if you have not provided a model path argument '--model-dir' the models might be downloaded somewhere you don't want them to be.",
            file=sys.stderr,
        )

    models = CLAIR3_MANIFEST
    if args.models:
        models = [m for m in models if m["name"] in args.models]
        unknown = set(args.models) - {m["name"] for m in models}
        for name in unknown:
            print(f"Warning: unknown model '{name}' — skipping", file=sys.stderr)

    for model in models:

        model_path = Path(args.model_dir, model["name"])

        if model.get("pytorch"):
            needs_download = not all(
                Path(model_path, f).exists()
                for f in ["pileup.pt", "full_alignment.pt"]
            )
        else:
            needs_download = (
                not model_path.exists() or len(os.listdir(model_path)) == 0
            )

        if needs_download:
            try:
                fetch_model(model, args.model_dir, backup_base_url)
            except RuntimeError as e:
                print(colored.red(str(e)), file=sys.stderr)
                raise SystemExit(1)
            print(f"Downloaded model: {model['name']}", file=sys.stderr)
        else:
            print(
                f"Model {model['name']} already exists, skipping download",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
