# ViralSafeTarget — מדריך קצר בעברית

ViralSafeTarget הוא כלי מחקר חישובי שמתחיל מגנומים של וירוס ומפיק מועמדים
מתועדים לבחינה נוספת: שימור בין זנים, מיפוי לגנים, סריקת דמיון למארח,
השערות שיבוש חלבון וראיות מהספרות עם ביקורת אנושית.

[README באנגלית](README.md) · [מחברות](notebooks/README.md) ·
[תוצאות HSV-2](reports/README.md) · [כל התיעוד](docs/README.md)

> זהו כלי ליצירת השערות מחקריות. הוא אינו מוכיח חיתוך, בטיחות, השבתת וירוס,
> ניקוי לטנטיות, יעילות טיפולית או ריפוי, ואינו כולל פרוטוקול מעבדה.

## איפה מתחילים?

| המטרה | נקודת הכניסה |
|---|---|
| ללמוד בעברית | [`notebooks/learning/he/00_START_HERE.ipynb`](notebooks/learning/he/00_START_HERE.ipynb) |
| לראות מה כבר מצאנו | [`reports/README.md`](reports/README.md) |
| להריץ דמו קטן | `make demo` |
| לשחזר את מחקר HSV-2 | `vst reproduce hsv2` |
| לבדוק וירוס חדש | `vst project init ...` |
| לחפש ראיות בספרות | `vst evidence discover ...` |

פקודת `vst` היא הממשק הראשי לחוקר. הסקריפטים משמשים לשחזור case studies
ולתחזוקה ומתועדים ב־[`scripts/README.md`](scripts/README.md).

## התקנה

```bash
conda env create -f environment.yml
conda activate viral-safe-target
python -m pip install -e .
pytest -q
ruff check .
```

לאחר מכן:

```bash
make demo
make notebook
```

## וירוס חדש

```bash
vst project init \
  --id my-virus \
  --display-name "My virus" \
  --reference-accession REF_ACCESSION \
  --out-dir projects/my-virus

vst project validate --project projects/my-virus/project.yaml
vst project run --project projects/my-virus/project.yaml
vst project status --project projects/my-virus/project.yaml
```

הפרויקט מכיל profiles נפרדים עבור הווירוס, המארח והעורך. מוסיפים FASTA,
GFF ויישור זנים, בלי לקודד שמות גנים ספציפיים בתוך Python. מדריך מלא:
[`docs/getting-started/NEW_VIRUS_WORKFLOW.md`](docs/getting-started/NEW_VIRUS_WORKFLOW.md).

## שחזור HSV-2

```bash
vst reproduce hsv2             # מציג תוכנית בלבד
vst reproduce hsv2 --execute   # מריץ או ממשיך מ-cache
```

התוצאות הציבוריות נמצאות ב־[`reports/`](reports/README.md). יש שם שתי ריצות
בעלות עומק דגימה שונה: showcase מאוזן וריצה exhaustive מאוחרת. אין להשוות את
הדירוגים בלי לקרוא את ההסבר המצורף.

## איך מפרשים את הפלט?

המערכת מפרידה בין:

1. איכות הרצף כמטרה חישובית;
2. ראיות ביולוגיות מאושרות ומצוטטות;
3. השערת שיבוש ברמת הרצף והחלבון;
4. כיסוי הראיות והמידע החסר.

ראיה חסרה נשארת `unknown`. תוצאה של אפס פגיעות מארח חזויות תקפה רק למודל
החיפוש שהוגדר ואינה הוכחת בטיחות. ראיית HSV-1 נשמרת כראיית ortholog ואינה
מוצגת כהוכחה ישירה ב־HSV-2.

## Evidence Agent

```bash
vst evidence discover --project project.yaml
```

הפקודה יוצרת `review_queue.tsv`. כל שורה מתחילה כ־`pending`; החוקר בודק את
המאמר, סוג הניסוי, מודל הניסוי, מין הווירוס והמשפט הרלוונטי. רק שורות שאושרו
עם שם בודק, תאריך ומקור נכנסות לטבלת הראיות:

```bash
vst evidence apply --project project.yaml \
  --review-queue results/evidence/review_queue.tsv
```

מדריך: [`docs/workflows/EVIDENCE_AGENT.md`](docs/workflows/EVIDENCE_AGENT.md).

## מפת המאגר

- [`notebooks/README.md`](notebooks/README.md) — סדר המחברות ומה כל אחת עושה.
- [`docs/README.md`](docs/README.md) — תיעוד לפי נושא.
- [`scripts/README.md`](scripts/README.md) — entry points מול helpers פנימיים.
- [`configs/README.md`](configs/README.md) — profiles וקונפיגורציות.
- [`reports/README.md`](reports/README.md) — תוצאות שכבר הורצו.

למגבלות המדעיות קראו את
[`docs/research/KNOWN_LIMITATIONS.md`](docs/research/KNOWN_LIMITATIONS.md), את
[`DISCLAIMER.md`](DISCLAIMER.md) ואת [`SECURITY.md`](SECURITY.md).
