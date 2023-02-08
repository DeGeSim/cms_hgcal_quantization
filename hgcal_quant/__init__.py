from pathlib import Path

from omegaconf import OmegaConf

conf = OmegaConf.load(Path("~/hgcal_quant/conf.yaml").expanduser())
