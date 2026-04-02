from __future__ import annotations

from pathlib import Path
import json

import pandas as pd

from .artifact_loader import build_case_bundle
from .config import SpatialPathologistConfig
from .heuristics import (
    build_case_summary,
    build_cluster_review,
    build_structure_review,
)
from .openai_client import OpenAIJsonClient, OpenAISettings
from .report import write_html_report


def _cluster_review_schema() -> dict[str, object]:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "standardized_label": {"type": "string"},
            "cell_family": {"type": "string"},
            "confidence": {"type": "number", "minimum": 0, "maximum": 1},
            "review_priority": {"type": "string", "enum": ["low", "medium", "high"]},
            "summary": {"type": "string"},
            "evidence": {"type": "array", "items": {"type": "string"}, "maxItems": 8},
            "concerns": {"type": "array", "items": {"type": "string"}, "maxItems": 8},
        },
        "required": [
            "standardized_label",
            "cell_family",
            "confidence",
            "review_priority",
            "summary",
            "evidence",
            "concerns",
        ],
    }


def _structure_review_schema() -> dict[str, object]:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "title": {"type": "string"},
            "assigned_label": {"type": "string"},
            "behavior": {"type": "string"},
            "confidence": {"type": "number", "minimum": 0, "maximum": 1},
            "review_priority": {"type": "string", "enum": ["low", "medium", "high"]},
            "discovery_flag": {"type": "boolean"},
            "summary": {"type": "string"},
            "top_celltype_summary": {"type": "string"},
            "key_evidence": {"type": "array", "items": {"type": "string"}, "maxItems": 10},
            "differential": {"type": "array", "items": {"type": "string"}, "maxItems": 8},
            "recommended_checks": {"type": "array", "items": {"type": "string"}, "maxItems": 8},
        },
        "required": [
            "title",
            "assigned_label",
            "behavior",
            "confidence",
            "review_priority",
            "discovery_flag",
            "summary",
            "top_celltype_summary",
            "key_evidence",
            "differential",
            "recommended_checks",
        ],
    }


def _case_summary_schema() -> dict[str, object]:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "headline": {"type": "string"},
            "overall_impression": {"type": "string"},
            "key_findings": {"type": "array", "items": {"type": "string"}, "maxItems": 8},
            "review_priorities": {"type": "array", "items": {"type": "string"}, "maxItems": 8},
            "discovery_candidates": {"type": "array", "items": {"type": "string"}, "maxItems": 8},
        },
        "required": [
            "headline",
            "overall_impression",
            "key_findings",
            "review_priorities",
            "discovery_candidates",
        ],
    }


def _cluster_prompt(cluster: dict, study_context: str) -> tuple[str, str]:
    system = (
        "You are an expert spatial transcriptomics reviewer. "
        "Return JSON only with keys: standardized_label, cell_family, confidence, "
        "review_priority, summary, evidence, concerns."
    )
    user = json.dumps(
        {"study_context": study_context, "cluster": cluster},
        ensure_ascii=False,
        indent=2,
    )
    return system, user


def _structure_prompt(structure: dict, study_context: str) -> tuple[str, str]:
    system = (
        "You are an expert pathologist reviewing structure-level spatial transcriptomics evidence. "
        "Return JSON only with keys: title, assigned_label, behavior, confidence, review_priority, "
        "discovery_flag, summary, top_celltype_summary, key_evidence, differential, recommended_checks."
    )
    user = json.dumps(
        {"study_context": study_context, "structure": structure},
        ensure_ascii=False,
        indent=2,
    )
    return system, user


def _case_prompt(case_bundle: dict, structure_reviews: list[dict]) -> tuple[str, str]:
    system = (
        "You are an expert spatial pathology summarizer. "
        "Return JSON only with keys: headline, overall_impression, key_findings, "
        "review_priorities, discovery_candidates."
    )
    user = json.dumps(
        {
            "case_name": case_bundle["case_name"],
            "study_context": case_bundle["study_context"],
            "structures": structure_reviews,
        },
        ensure_ascii=False,
        indent=2,
    )
    return system, user


def _merge_review(default: dict, llm_result: dict, *, engine_name: str) -> dict:
    merged = dict(default)
    for key, value in llm_result.items():
        if value is None:
            continue
        merged[key] = value
    merged["engine"] = engine_name
    return merged


def run_spatial_pathologist(cfg: SpatialPathologistConfig) -> dict[str, str]:
    case_bundle = build_case_bundle(cfg)
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    llm_client = None
    settings = OpenAISettings(
        enabled=cfg.openai_enabled,
        api_key_env=cfg.openai_api_key_env,
        model=cfg.openai_model,
        reasoning_effort=cfg.openai_reasoning_effort,
        store=cfg.openai_store,
    )
    if settings.available:
        llm_client = OpenAIJsonClient(settings)

    cluster_reviews = []
    for cluster in case_bundle["clusters"]:
        heuristic_review = build_cluster_review(cluster)
        if llm_client is not None:
            try:
                system, user = _cluster_prompt(cluster, case_bundle["study_context"])
                llm_review = llm_client.generate_json(
                    system_prompt=system,
                    user_prompt=user,
                    schema_name="cluster_review",
                    schema=_cluster_review_schema(),
                )
                cluster_reviews.append(
                    _merge_review(heuristic_review, llm_review, engine_name=f"openai:{cfg.openai_model}")
                )
                continue
            except Exception:
                pass
        cluster_reviews.append(heuristic_review)

    structure_reviews = []
    for structure in case_bundle["structures"]:
        heuristic_review = build_structure_review(
            structure,
            low_confidence_threshold=cfg.low_confidence_threshold,
            ambiguity_margin_threshold=cfg.ambiguity_margin_threshold,
        )
        if llm_client is not None:
            try:
                system, user = _structure_prompt(structure, case_bundle["study_context"])
                llm_review = llm_client.generate_json(
                    system_prompt=system,
                    user_prompt=user,
                    schema_name="structure_review",
                    schema=_structure_review_schema(),
                )
                structure_reviews.append(
                    _merge_review(heuristic_review, llm_review, engine_name=f"openai:{cfg.openai_model}")
                )
                continue
            except Exception:
                pass
        structure_reviews.append(heuristic_review)

    case_summary = build_case_summary(case_bundle, structure_reviews)
    if llm_client is not None:
        try:
            system, user = _case_prompt(case_bundle, structure_reviews)
            llm_summary = llm_client.generate_json(
                system_prompt=system,
                user_prompt=user,
                schema_name="case_summary",
                schema=_case_summary_schema(),
            )
            case_summary = _merge_review(case_summary, llm_summary, engine_name=f"openai:{cfg.openai_model}")
        except Exception:
            pass

    (output_dir / "cluster_reviews.json").write_text(
        json.dumps(cluster_reviews, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    pd.DataFrame(cluster_reviews).to_csv(
        output_dir / "cluster_reviews.csv",
        index=False,
    )
    (output_dir / "structure_reviews.json").write_text(
        json.dumps(structure_reviews, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    pd.DataFrame(structure_reviews).to_csv(
        output_dir / "structure_reviews.csv",
        index=False,
    )
    (output_dir / "case_summary.json").write_text(
        json.dumps(case_summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (output_dir / "run_metadata.json").write_text(
        json.dumps(
            {
                "openai_requested": cfg.openai_enabled,
                "openai_api_key_env": cfg.openai_api_key_env,
                "openai_api_key_present": settings.api_key is not None,
                "openai_dependency_available": settings.available,
                "openai_model": cfg.openai_model,
                "openai_reasoning_effort": cfg.openai_reasoning_effort,
                "openai_store": cfg.openai_store,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    report_path = write_html_report(
        output_dir=output_dir,
        case_bundle=case_bundle,
        cluster_reviews=cluster_reviews,
        structure_reviews=structure_reviews,
        case_summary=case_summary,
    )

    return {
        "output_dir": str(output_dir),
        "report_html": str(report_path),
        "cluster_reviews_json": str(output_dir / "cluster_reviews.json"),
        "structure_reviews_json": str(output_dir / "structure_reviews.json"),
        "case_summary_json": str(output_dir / "case_summary.json"),
    }
