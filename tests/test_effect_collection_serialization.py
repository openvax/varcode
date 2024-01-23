from nose.tools import eq_
import pickle
from varcode import EffectCollection
from .data import tcga_ov_variants, ov_wustle_variants

tcga_ov_effects = tcga_ov_variants.effects()
ov_wustle_effects = ov_wustle_variants.effects()

def test_tcga_effect_collection_to_dict():
    eq_(
        tcga_ov_effects,
        EffectCollection.from_dict(tcga_ov_effects.to_dict()))

def test_wustle_effect_collection_to_dict():
    eq_(
        ov_wustle_effects,
        EffectCollection.from_dict(ov_wustle_effects.to_dict()))

def test_tcga_effect_collection_to_json():
    eq_(tcga_ov_effects, EffectCollection.from_json(tcga_ov_effects.to_json()))

def test_wustle_effect_collection_to_json():
    eq_(
        ov_wustle_effects,
        EffectCollection.from_json(ov_wustle_effects.to_json()))

def test_tcga_effect_collection_pickling():
    reconstructed = pickle.loads(pickle.dumps(tcga_ov_effects))
    eq_(tcga_ov_effects, reconstructed)

def test_wustle_effect_collection_pickling():
    reconstructed = pickle.loads(pickle.dumps(ov_wustle_effects))
    eq_(ov_wustle_effects, reconstructed)