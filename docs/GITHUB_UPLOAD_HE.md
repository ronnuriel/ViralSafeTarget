# העלאה ל-GitHub

## דרך האתר

1. צור repository חדש בשם `ViralSafeTarget`.
2. אל תסמן יצירת README או LICENSE, כי הם כבר קיימים.
3. פתח Terminal בתוך תיקיית הפרויקט והריץ את הפקודות שמופיעות בעמוד GitHub.

## דרך Terminal

```bash
cd ViralSafeTarget
git init
git add .
git commit -m "Initial public research release"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/ViralSafeTarget.git
git push -u origin main
```

אם מותקן GitHub CLI:

```bash
cd ViralSafeTarget
git init
git add .
git commit -m "Initial public research release"
gh repo create ViralSafeTarget --public --source=. --remote=origin --push
```

## לפני פרסום ציבורי

- החלף `OWNER` ב-`pyproject.toml` וב-`CITATION.cff`.
- הכנס את שמך ל-`CITATION.cff` ולרישיון אם תרצה.
- ודא שאין ב-`data/raw/` קבצים גדולים או נתונים פרטיים.
- הרץ `pytest` ו-`ruff check .`.
- פתח Issue ראשון שמגדיר את benchmark של HSV-2.
