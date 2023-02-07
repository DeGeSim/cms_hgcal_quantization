from omegaconf import OmegaConf
from pathlib import Path

conf = OmegaConf.load(Path("../conf.yaml").expanduser())
