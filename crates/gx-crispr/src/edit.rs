//! Simple edit simulation: NHEJ indels and HDR replacement at/around a cut site.

use gx_plasmid::translate;

#[derive(Debug, Clone)]
pub struct EditPreview {
    pub description: String,
    /// A window of the original sequence around the edit.
    pub before: String,
    /// The same window after the edit.
    pub after: String,
    /// `Some(true)` if the edit shifts a reading frame (indel not ÷3).
    pub frameshift: Option<bool>,
    /// Protein around the edit before/after (frame +0 from window start).
    pub protein_before: String,
    pub protein_after: String,
}

const FLANK: usize = 15;

fn window(seq: &[u8], lo: usize, hi: usize) -> (usize, usize) {
    (lo.saturating_sub(FLANK), (hi + FLANK).min(seq.len()))
}

/// Model an NHEJ deletion of `del_len` bp immediately 3' of the cut site.
pub fn nhej_deletion(seq: &[u8], cut: usize, del_len: usize) -> EditPreview {
    let n = seq.len();
    let start = cut.min(n);
    let end = (start + del_len).min(n);
    let mut edited = Vec::with_capacity(n - (end - start));
    edited.extend_from_slice(&seq[..start]);
    edited.extend_from_slice(&seq[end..]);

    let (w0, w1) = window(seq, start, end);
    let before = &seq[w0..w1];
    let ew1 = w1.saturating_sub(end - start).min(edited.len());
    let after = &edited[w0.min(edited.len())..ew1];

    EditPreview {
        description: format!("NHEJ — {del_len} bp deletion at cut site {}", cut + 1),
        before: String::from_utf8_lossy(before).into_owned(),
        after: String::from_utf8_lossy(after).into_owned(),
        frameshift: Some(del_len % 3 != 0),
        protein_before: translate(before),
        protein_after: translate(after),
    }
}

/// Model an HDR event that replaces `[start, end)` with a donor sequence.
pub fn hdr_replace(seq: &[u8], start: usize, end: usize, donor: &[u8]) -> EditPreview {
    let n = seq.len();
    let start = start.min(n);
    let end = end.clamp(start, n);
    let mut edited = Vec::new();
    edited.extend_from_slice(&seq[..start]);
    edited.extend_from_slice(donor);
    edited.extend_from_slice(&seq[end..]);

    let (w0, w1) = window(seq, start, end);
    let before = &seq[w0..w1];
    let ew1 = (w1 + donor.len()).saturating_sub(end - start).min(edited.len());
    let after = &edited[w0.min(edited.len())..ew1];
    let net = donor.len() as isize - (end - start) as isize;

    EditPreview {
        description: format!(
            "HDR — replace {}..{} ({} bp) with donor ({} bp)",
            start + 1,
            end,
            end - start,
            donor.len()
        ),
        before: String::from_utf8_lossy(before).into_owned(),
        after: String::from_utf8_lossy(after).into_owned(),
        frameshift: Some(net % 3 != 0),
        protein_before: translate(before),
        protein_after: translate(after),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deletion_changes_frame_when_not_multiple_of_three() {
        let seq = b"AAAAAAAAAAAAAAAACCCCGGGGTTTTAAAAAAAAAAAAAAAA";
        let e1 = nhej_deletion(seq, 20, 1);
        assert_eq!(e1.frameshift, Some(true));
        let e3 = nhej_deletion(seq, 20, 3);
        assert_eq!(e3.frameshift, Some(false));
        assert!(e1.after.len() < e1.before.len());
    }

    #[test]
    fn hdr_inserts_donor() {
        let seq = b"AAAAAAAAAAAAAAAATTTTAAAAAAAAAAAAAAAA";
        let e = hdr_replace(seq, 16, 20, b"GGGGGG");
        assert!(e.after.contains("GGGGGG"));
    }
}
