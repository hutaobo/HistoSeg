from __future__ import annotations

import numpy as np
import pytest
from scipy.ndimage import shift as ndi_shift
from skimage.draw import ellipse, polygon
from skimage.transform import rotate

from histoseg.threed import (
    apply_similarity_to_points,
    estimate_coda_image_registration,
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


def test_coda_image_registration_resolves_radon_half_turn_ambiguity():
    fixed = _asymmetric_mask()
    moved = rotate(
        fixed,
        120.0,
        resize=False,
        order=1,
        mode="constant",
        cval=0.0,
        preserve_range=True,
    )
    moving = ndi_shift(moved, shift=(8, -11), order=0, mode="constant", cval=0.0)

    radon_only = estimate_radon_rotation(fixed, moving, angle_step=1.0)
    result = estimate_coda_image_registration(fixed, moving)

    assert radon_only.rotation_degrees == pytest.approx(60.0, abs=1.0)
    assert result.rotation.rotation_degrees == pytest.approx(-120.0, abs=1.0)
    assert len(result.orientation_candidates) == 2
    assert max(candidate["translation_iou"] for candidate in result.orientation_candidates) > 0.9


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


def _asymmetric_mask() -> np.ndarray:
    image = np.zeros((192, 192), dtype=float)
    rr, cc = polygon(
        [60, 35, 70, 120, 150, 118, 82],
        [44, 100, 154, 142, 76, 28, 16],
        shape=image.shape,
    )
    image[rr, cc] = 1.0
    rr, cc = polygon([82, 96, 116, 102], [78, 112, 100, 68], shape=image.shape)
    image[rr, cc] = 0.0
    rr, cc = polygon([28, 48, 52, 32], [132, 152, 176, 164], shape=image.shape)
    image[rr, cc] = 1.0
    return image
