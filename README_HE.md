# ViralSafeTarget — סביבת מחקר חישובית לווירוסים

המחברת המלאה החדשה כתובה באנגלית:

```bash
jupyter lab notebooks/08_RUN_FULL_PIPELINE_EN.ipynb
```

הפרויקט מקבל אוסף גנומים של וירוס, מחפש אתרים שמורים בין הזנים, ממפה אותם למידע ביולוגי, מכין סריקת סיכון מול גנום המארח ומחשב **מה היה משתנה ברצף** אילו חיתוך אידאלי התרחש.

> זהו כלי ליצירת השערות מחקריות. הוא לא מוכיח שהווירוס נשבר, לא מוכיח בטיחות, לא מסלק HSV מהגוף ואינו פרוטוקול מעבדה.

## התמונה הפשוטה

```text
FASTA מיושר של הרבה זנים
        ↓
איפה הרצף נשמר?
        ↓
איפה יש guide + PAM שמתאימים לעורך?
        ↓
מה אומר ה-GFF על האזור הזה?
        ↓
האם קיימים אתרים דומים באדם?
        ↓
רשימת מועמדים מוסברת
        ↓
סימולציית רצף: מה יימחק אם שני חיתוכים אכן יתרחשו?
```

## התחלה מיידית

```bash
conda env create -f environment.yml
conda activate viral-safe-target
python -m pip install -e .
pytest -q
ruff check .
bash scripts/run_demo.sh
```

לאחר מכן פתח את:

```text
reports/demo/report.html
reports/demo/candidates.csv
reports/demo/simulated_pairs.csv
```

## ממשק העלאת קבצים

```bash
streamlit run app.py
```

בממשק מעלים:

1. `aligned_virus.fasta` — כמה גנומי וירוס שכבר עברו alignment.
2. `reference.gff3` — אופציונלי, מפת הגנים של גנום הייחוס.
3. FASTA קטן של מארח — רק להדגמה. עבור כל GRCh38 משתמשים בכלי חיצוני כמו Cas-OFFinder או CRISPRitz.

## הרצה על HSV-2 אמיתי

```bash
bash scripts/run_real_hsv2.sh --sample-size 25
```

הסקריפט:

- מוריד גנומי HSV-2 ציבוריים מ-NCBI;
- שומר את `NC_001798.2` כ-reference;
- מסנן כפילויות ורצפים חסרים;
- מיישר את הדגימות באמצעות MAFFT;
- מחפש אתרי SpCas9 שמורים;
- ממפה אותם ל-GFF;
- מפיק טבלת מועמדים וסימולציית מחיקות אידאלית;
- שומר `run_manifest.json` עם checksums ופרמטרים.

להכנת סריקה מול האדם:

```bash
bash scripts/run_real_hsv2.sh --with-human --sample-size 25
```

הפקודה מורידה גם GRCh38 ומכינה קלט ל-Cas-OFFinder. היא עדיין לא מריצה הוכחת בטיחות.

לפיילוט הממוקד והניתן לשחזור של UL19 ו-UL30:

```bash
bash scripts/run_hsv2_pilot.sh
```

הבחירה נעשית אחרי דירוג, סינון כפילויות ו-stratification לפי גן, עם עד 100
מועמדים לכל גן. פרטי ההרצה הדו-שלבית נמצאים ב-[`docs/HSV2_PILOT.md`](docs/HSV2_PILOT.md).

## פירוש הציונים

- `pre_human_score` משלב רכיבי שימור, ייחודיות בווירוס, GC, מורכבות רצף,
  annotation וראיות אוצרות אם הן קיימות.
- ראיה ביולוגית חסרה נשארת `null`; היא אינה הופכת אוטומטית ל-1.
- `post_human_score` מחושב בנפרד רק אחרי סיכום תוצאות Cas-OFFinder.
- היעדר פגיעה חזויה עד סף ה-mismatches שהוגדר אינו הוכחת בטיחות.
- לזוג מטרות בגנים שונים לא מיוחסת מחיקה פיזית אחת; זו השערת multi-target בלבד.

## מה הסימולציה כן עושה

לכל guide היא מחשבת את נקודת החיתוך הקנונית המשוערת של SpCas9. לזוג guides היא מחשבת:

- מיקום שני החיתוכים;
- אורך הקטע שהיה נמחק אילו שניהם התרחשו;
- באילו features של ה-GFF המחיקה חופפת;
- איזה חלק מה-feature נמחק ברמת הקואורדינטות;
- בכמה מהזנים שני האתרים קיימים בדיוק.

## מה הסימולציה לא יודעת

היא אינה יודעת:

- האם כלי העריכה הגיע לתא העצב;
- האם ה-DNA הלטנטי נגיש;
- האם באמת התרחש חיתוך;
- מה תהיה התפלגות תיקון ה-DNA;
- האם הפגיעה משביתה את הווירוס;
- האם יש נזק לתא או לאדם;
- האם הווירוס יכול להתעורר שוב.

בנוסף, שימור בין רצפים אינו הוכחת יעילות ניסויית, ושיבוש מתמטי של רצף אינו
הוכחה להשבתת הווירוס. נגישות בלטנטיות, delivery, יעילות עריכה, תיקון DNA,
רעילות ו-reactivation נמצאים מחוץ למודל.

כדי לדעת מה התרחש בפועל משתמשים בנתוני sequencing ובכלים כמו CRISPResso2. כדי להוכיח פגיעה בווירוס צריך ניסויי וירולוגיה מתאימים.

## איך מעלים ל-GitHub

הוראות מדויקות נמצאות ב-[`docs/GITHUB_UPLOAD_HE.md`](docs/GITHUB_UPLOAD_HE.md). הכנתי גם CI, קובץ citation, תבניות Issues, רישיון ומסמכי תרומה.

## איפה מתחילים ללמוד

- [`docs/CONCEPTS_HE.md`](docs/CONCEPTS_HE.md) — FASTA, GFF, guide, PAM ו-off-target.
- [`docs/DATA_FORMATS.md`](docs/DATA_FORMATS.md) — חוזה הקלט והפלט.
- [`docs/SIMULATION_LIMITS.md`](docs/SIMULATION_LIMITS.md) — מה אפשר להסיק ומה לא.
- [`docs/EXISTING_TOOLS.md`](docs/EXISTING_TOOLS.md) — כלים קיימים ומה אנחנו מוסיפים.
- [`docs/RESEARCH_PLAN.md`](docs/RESEARCH_PLAN.md) — מסלול להפיכת הפרויקט למאמר.
