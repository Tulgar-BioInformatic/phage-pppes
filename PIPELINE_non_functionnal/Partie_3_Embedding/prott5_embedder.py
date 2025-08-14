#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prott5_embedder.py
------------------
But : calculer les embeddings ProtT5 pour une sélection d'ID issues d'un FASTA,
en gérant automatiquement les erreurs de mémoire GPU (découpage de batch, puis
fallback CPU si une séquence est trop grosse).

Entrées principales :
  - un FASTA d'acides aminés (argument -i/--input)
  - une liste d'IDs à conserver (fichier texte, 1 ID/ligne)
      -> ici : 'Partie_4_UMAP/input/IDs_sup15_all_or_ref.txt'

Sortie :
  - écrit *progressivement* un fichier TSV de vecteurs (ID \t v0 v1 v2 ...)
    dans le chemin donné par -o/--output (NB : malgré le nommage interne
    'out_npz', la fonction utilisée ci-dessous sauvegarde un **TSV**).

Exemple :
  python prott5_embedder.py -i mon.faa -o out/embeddings.tsv \\
         --per_protein --max_residues 12000 --max_batch 256 \\
         --save_every 20 --auto_shrink
"""

from __future__ import annotations
import argparse
import time
from pathlib import Path
from typing import List, Tuple, Dict

import torch
from transformers import T5EncoderModel, T5Tokenizer
import numpy as np

# --- Sélection du device (GPU si dispo) --------------------------------------
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# --- Utilitaires de sauvegarde ------------------------------------------------
def save_ids_and_embeddings(emb_dict: Dict[str, np.ndarray], npz_path: Path) -> None:
    """Sauvegarde en NPZ (utile si on souhaite un format compressé NumPy)."""
    if not emb_dict:
        return
    ids = sorted(emb_dict.keys())
    embs = [emb_dict[_id] for _id in ids]
    try:
        arr = np.stack(embs)
    except ValueError:
        # tableaux de tailles variables → ragged array
        arr = np.array(embs, dtype=object)
    np.savez_compressed(npz_path, ids=np.array(ids, dtype=object),
                        embeddings=arr, allow_pickle=True)

def save_ids_and_embeddings_tsv(emb_dict: Dict[str, np.ndarray], tsv_path: Path) -> None:
    """Sauvegarde en TSV : <seqID> \t <v0 v1 v2 ...>"""
    if not emb_dict:
        return
    with tsv_path.open("w") as fh:
        for uid in sorted(emb_dict.keys()):
            vec = " ".join(map(str, emb_dict[uid].reshape(-1)))
            fh.write(f"{uid}\t{vec}\n")

# --- IDs à garder (fichier texte, 1 par ligne) --------------------------------
IDs: set[str] = set()
with open("Partie_4_UMAP/input/IDs_sup15_all_or_ref.txt") as f:
    IDs = {l.strip() for l in f}
    print(f"{len(IDs)} IDs chargés")

# --- Lecture FASTA (ne charge que les IDs utiles) -----------------------------
def read_fasta(path: Path, ids_keep: set[str]) -> Dict[str, str]:
    """Ne lit que les séquences dont l’ID (mot juste après '>') appartient à ids_keep."""
    seqs: Dict[str, str] = {}
    keep = False
    uid = ""
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                uid = line[1:].split()[0]
                keep = uid in ids_keep
                if keep:
                    # éviter caractères spéciaux qui cassent les noms de fichiers/colonnes
                    uid = uid.replace("/", "_").replace(".", "_")
                    seqs[uid] = ""
            elif keep:
                seqs[uid] += line.strip().upper().replace("-", "")
    return seqs

# --- Chargement du modèle -----------------------------------------------------
def load_model(cache: Path | None, ckpt: str = "Rostlab/prot_t5_xl_half_uniref50-enc"):
    """Charge l'encodeur ProtT5 + tokenizer."""
    print(f"Loading model: {ckpt}")
    model = T5EncoderModel.from_pretrained(ckpt, cache_dir=cache)
    tok = T5Tokenizer.from_pretrained(ckpt, do_lower_case=False, cache_dir=cache)
    if device.type == "cpu":
        model = model.to(torch.float32)  # en CPU, rester en FP32
    model = model.to(device).eval()
    return model, tok

# --- Embedding sécurisé (split auto en cas d'OOM (out of memory)) ----------------------------
def safe_embed(
    model: T5EncoderModel,
    tok: T5Tokenizer,
    ids: List[str],
    seqs_spaced: List[str],
    lens: List[int],
    per_protein: bool,
    emb_dict: Dict[str, np.ndarray],
    depth: int = 0,
) -> None:
    """
    Tente d'embedder un sous-batch. Si OOM GPU :
      - on scinde récursivement en deux si batch>1 (auto_shrink),
      - si batch==1, on bascule sur CPU pour cette séquence.
    """
    try:
        toks = tok(seqs_spaced, add_special_tokens=True, padding=True, return_tensors="pt")
        inp  = toks["input_ids"].to(device)
        attn = toks["attention_mask"].to(device)
        with torch.no_grad():
            reps = model(inp, attention_mask=attn).last_hidden_state  # [B, L, H]
        for b, uid in enumerate(ids):
            emb = reps[b, : lens[b]]          # tronquer au vrai length
            if per_protein:
                emb = emb.mean(dim=0)         # pooling moyen par protéine
            emb_dict[uid] = emb.cpu().numpy().squeeze()
    except RuntimeError as e:
        if "out of memory" in str(e).lower():
            torch.cuda.empty_cache()
            if len(ids) > 1:
                mid = len(ids) // 2
                safe_embed(model, tok, ids[:mid], seqs_spaced[:mid], lens[:mid], per_protein, emb_dict, depth + 1)
                safe_embed(model, tok, ids[mid:], seqs_spaced[mid:], lens[mid:], per_protein, emb_dict, depth + 1)
            else:
                uid = ids[0]
                print(f"⚠️  {uid} trop gros pour le GPU → fallback CPU")
                toks_cpu = tok(seqs_spaced, add_special_tokens=True, padding=True, return_tensors="pt")
                with torch.no_grad():
                    reps_cpu = model.to("cpu", dtype=torch.float32)(
                        toks_cpu["input_ids"], attention_mask=toks_cpu["attention_mask"]
                    ).last_hidden_state
                emb = reps_cpu[0, : lens[0]]
                if per_protein:
                    emb = emb.mean(dim=0)
                emb_dict[uid] = emb.numpy().squeeze()
                model.to(device)  # remettre le modèle sur GPU pour la suite
        else:
            raise

# --- Pipeline principal -------------------------------------------------------
def embed_fasta(
    fasta_path: Path,
    out_npz: Path,                 #  voir note : on sauvegarde en TSV plus bas
    cache_dir: Path | None,
    per_protein: bool,
    *,
    max_residues: int,
    max_seq_len: int,
    max_batch: int,
    save_every: int,
    auto_shrink: bool,
) -> None:
    # Lire uniquement les séquences utiles
    seqs = read_fasta(fasta_path, IDs)
    if not seqs:
        raise ValueError("Empty FASTA or no IDs matched")

    model, tok = load_model(cache_dir)

    # Trier par longueur décroissante (pratique pour équilibrer les batchs)
    entries = sorted(seqs.items(), key=lambda kv: len(kv[1]), reverse=True)
    emb_dict: Dict[str, np.ndarray] = {}

    batch: List[Tuple[str, str, int]] = []
    batch_idx = 0
    start = time.time()

    for i, (pid, seq) in enumerate(entries, 1):
        prev_len = len(emb_dict)
        # Remplacements classiques des AA non standards
        seq = seq.replace("U", "X").replace("Z", "X").replace("O", "X")
        batch.append((pid, " ".join(seq), len(seq)))

        # Déclenche le flush si on dépasse les limites
        res_ct = sum(l for _, _, l in batch)
        flush = (len(batch) >= max_batch) or (res_ct >= max_residues) or (i == len(entries)) or (len(seq) > max_seq_len)
        if not flush:
            continue

        ids, spaced, lens = zip(*batch)
        batch.clear()
        batch_idx += 1

        # Passage par safe_embed (gère auto_shrink si demandé)
        safe_embed(model, tok, list(ids), list(spaced), list(lens), per_protein, emb_dict)

        print(f"Batch {batch_idx}: +{len(emb_dict) - prev_len}  total={len(emb_dict)}/{len(entries)}")
        # Sauvegarde périodique (ici en TSV — voir NOTE)
        if batch_idx % save_every == 0 or i == len(entries):
            save_ids_and_embeddings_tsv(emb_dict, out_npz)  # ← TSV volontaire

    dur = time.time() - start
    print(f"Finished in {dur/60:.1f} min")

# --- Arguments CLI ------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="ProtT5 incremental embedder (auto-shrink)")
    p.add_argument("-i", "--input", required=True, help="FASTA file")
    p.add_argument("-o", "--output", required=True, help="Output file (TSV conseillé)")
    p.add_argument("--model", default=None, help="HF cache directory")
    p.add_argument("--per_protein", action="store_true", help="Mean-pooled per protein")

    p.add_argument("--max_residues", type=int, default=4000, help="Seuil total d'acides aminés par batch")
    p.add_argument("--max_seq_len", type=int, default=1000, help="Longueur max déclenchant un flush")
    p.add_argument("--max_batch", type=int, default=256, help="Taille max du batch (nb de séquences)")
    p.add_argument("--save_every", type=int, default=10000, help="Flush écritures toutes n batches")
    p.add_argument("--auto_shrink", action="store_true", help="Split récursif si OOM GPU")
    return p

# --- Main ---------------------------------------------------------------------
def main() -> None:
    args = build_parser().parse_args()

    fasta = Path(args.input).expanduser().resolve()
    out_dir = Path(args.output).expanduser().resolve()
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    out_path = out_dir  # on traite 'output' comme un chemin de fichier TSV

    cache = Path(args.model).expanduser().resolve() if args.model else None

    embed_fasta(
        fasta_path=fasta,
        out_npz=out_path,                 # NOTE : chemin de sortie TSV ici
        cache_dir=cache,
        per_protein=args.per_protein,
        max_residues=args.max_residues,
        max_seq_len=args.max_seq_len,
        max_batch=args.max_batch,
        save_every=args.save_every,
        auto_shrink=args.auto_shrink,
    )

if __name__ == "__main__":
    main()
