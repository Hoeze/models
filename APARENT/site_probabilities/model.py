from kipoi.model import BaseModel
from keras.models import load_model
import numpy as np


class APARENTModel(BaseModel):

    def __init__(self, weights):
        self.weights = weights
        self.model = load_model(weights)

    def predict_on_batch(self, inputs):
        batch_size = inputs.shape[0]

        input_1 = np.expand_dims(inputs, -1)
        input_2 = np.zeros((batch_size, 13))
        input_3 = np.ones((batch_size, 1))

        _, out = self.model.predict_on_batch([input_1, input_2, input_3])

        return out
