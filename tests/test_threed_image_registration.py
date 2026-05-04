from __future__ import annotations

import numpy as np
import pytest
from scipy.ndimage import shift as ndi_shift
from skimage.draw import ellipse
from skimage.transform import rotate

from histoseg.threed import (
    apply_similarity_to_points,
    estimate_radon_rotation,
    estimate_translation,
)


def test_estimate_radon_rotation_recovers_rotated_tissue_mask():
    fixed = _ellipse_mask()
    moving = rotate(
        fixed,
        17.0,
        resize=False,
        order=1,
        mode="constant",
        cval=0.0,
        preserve_range=True,
    )

    result = estimate_radon_rotation(fixed, moving, angle_step=0.5)

    assert result.rotation_degrees == pytest.approx(-17.0, abs=1.0)


def test_estimate_translation_recovers_shifted_image():
    fixed = _ellipse_mask()
    moving = ndi_shift(fixed, shift=(7, -5), order=0, mode="constant", cval=0.0)

    result = estimate_translation(fixed, moving)

    assert result.shift_y == pytest.approx(-7.0, abs=1.0)
    assert result.shift_x == pytest.approx(5.0, abs=1.0)


def test_apply_similarity_to_points_uses_xy_translation_and_origin():
    coords = np.array([[1.0, 0.0]])

    transformed = apply_similarity_to_points(
        coords,
        rotation_degrees=90.0,
        translate_xy=(10.0, -2.0),
        origin=(0.0, 0.0),
    )

    np.testing.assert_allclose(transformed, [[10.0, -1.0]], atol=1e-6)


def _ellipse_mask() -> np.ndarray:
    image = np.zeros((128, 128), dtype=float)
    rr, cc = ellipse(64, 64, 18, 42, rotation=np.deg2rad(28), shape=image.shape)
    image[rr, cc] = 1.0
    return image
