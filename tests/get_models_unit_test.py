"""Unit tests for artic/get_models.py — focuses on needs_download logic."""

from unittest.mock import MagicMock, patch

import pytest

from artic.get_models import main


_PYTORCH_MODEL = {
    "name": "test_pytorch_model",
    "model_url": "https://example.com/clair3_models_pytorch/test_pytorch_model/",
    "pytorch": True,
}

_TF_MODEL = {
    "name": "test_tf_model",
    "model_fname": "test_tf_model.tar.gz",
    "model_url": "https://example.com/clair3/test_tf_model.tar.gz",
}

# TF checkpoint file names that Clair3 writes when a model is extracted
_TF_FILES = [
    "pileup.data-00000-of-00001",
    "pileup.index",
    "full_alignment.data-00000-of-00001",
    "full_alignment.index",
]


def _run_main(tmp_path, manifest, mock_get_pytorch, mock_get_model):
    with (
        patch("artic.get_models.CLAIR3_MANIFEST", manifest),
        patch("artic.get_models.get_pytorch_model", mock_get_pytorch),
        patch("artic.get_models.get_model", mock_get_model),
        patch("sys.argv", ["artic_get_models", "--model-dir", str(tmp_path)]),
    ):
        main()


class TestPytorchNeedsDownload:
    def test_missing_directory_triggers_download(self, tmp_path):
        mock_dl = MagicMock()
        _run_main(tmp_path, [_PYTORCH_MODEL], mock_dl, MagicMock())
        mock_dl.assert_called_once()

    def test_empty_directory_triggers_download(self, tmp_path):
        (tmp_path / _PYTORCH_MODEL["name"]).mkdir()
        mock_dl = MagicMock()
        _run_main(tmp_path, [_PYTORCH_MODEL], mock_dl, MagicMock())
        mock_dl.assert_called_once()

    def test_tf_only_directory_triggers_download(self, tmp_path):
        """A directory containing only TF checkpoint files must still be downloaded."""
        model_dir = tmp_path / _PYTORCH_MODEL["name"]
        model_dir.mkdir()
        for fname in _TF_FILES:
            (model_dir / fname).write_text("")

        mock_dl = MagicMock()
        _run_main(tmp_path, [_PYTORCH_MODEL], mock_dl, MagicMock())
        mock_dl.assert_called_once()

    def test_partial_pt_files_triggers_download(self, tmp_path):
        """Only pileup.pt present (full_alignment.pt missing) → still needs download."""
        model_dir = tmp_path / _PYTORCH_MODEL["name"]
        model_dir.mkdir()
        (model_dir / "pileup.pt").write_text("")

        mock_dl = MagicMock()
        _run_main(tmp_path, [_PYTORCH_MODEL], mock_dl, MagicMock())
        mock_dl.assert_called_once()

    def test_tf_files_plus_partial_pt_triggers_download(self, tmp_path):
        """TF files + only one .pt file → still needs download."""
        model_dir = tmp_path / _PYTORCH_MODEL["name"]
        model_dir.mkdir()
        for fname in _TF_FILES:
            (model_dir / fname).write_text("")
        (model_dir / "pileup.pt").write_text("")

        mock_dl = MagicMock()
        _run_main(tmp_path, [_PYTORCH_MODEL], mock_dl, MagicMock())
        mock_dl.assert_called_once()

    def test_complete_pt_files_skips_download(self, tmp_path):
        """Both pileup.pt and full_alignment.pt present → skip download."""
        model_dir = tmp_path / _PYTORCH_MODEL["name"]
        model_dir.mkdir()
        (model_dir / "pileup.pt").write_text("")
        (model_dir / "full_alignment.pt").write_text("")

        mock_dl = MagicMock()
        _run_main(tmp_path, [_PYTORCH_MODEL], mock_dl, MagicMock())
        mock_dl.assert_not_called()

    def test_mixed_tf_and_pt_files_skips_download(self, tmp_path):
        """TF checkpoint files alongside both .pt files → already complete, skip download."""
        model_dir = tmp_path / _PYTORCH_MODEL["name"]
        model_dir.mkdir()
        for fname in _TF_FILES + ["pileup.pt", "full_alignment.pt"]:
            (model_dir / fname).write_text("")

        mock_dl = MagicMock()
        _run_main(tmp_path, [_PYTORCH_MODEL], mock_dl, MagicMock())
        mock_dl.assert_not_called()


class TestTfNeedsDownload:
    def test_missing_directory_triggers_download(self, tmp_path):
        mock_dl = MagicMock()
        _run_main(tmp_path, [_TF_MODEL], MagicMock(), mock_dl)
        mock_dl.assert_called_once()

    def test_empty_directory_triggers_download(self, tmp_path):
        (tmp_path / _TF_MODEL["name"]).mkdir()
        mock_dl = MagicMock()
        _run_main(tmp_path, [_TF_MODEL], MagicMock(), mock_dl)
        mock_dl.assert_called_once()

    def test_non_empty_directory_skips_download(self, tmp_path):
        model_dir = tmp_path / _TF_MODEL["name"]
        model_dir.mkdir()
        for fname in _TF_FILES:
            (model_dir / fname).write_text("")

        mock_dl = MagicMock()
        _run_main(tmp_path, [_TF_MODEL], MagicMock(), mock_dl)
        mock_dl.assert_not_called()
