from fairchem.core.units import mlip_unit
from fairchem.core import FAIRChemCalculator

MODEL = {
    'small': 'model/UMA/uma-s-1p1',
    'medium': 'model/UMA/uma-m-1p1',
    'large': ''
}

class PredictorLoader:
    def __init__(self, model_size='small', device='cuda'):
        if not MODEL[model_size]:
            raise ValueError(f"Model path for {model_size} is empty or invalid.")
        self.predictor = mlip_unit.load_predict_unit(path=MODEL[model_size], device=device)
        self.calc = FAIRChemCalculator(self.predictor, task_name='omol')
