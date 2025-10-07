# Test Data Input/Output Examples

**⚠️ Auto-generated from `test_generate_readme.py` - Do not edit manually!**

## Input → Output Transformations

*Generated from actual test runs showing how photon-mosaic transforms input
folder structures:*

### Basic Dataset

**Input:**
```
├── 001
│   ├── type_1_01.tif
│   ├── type_1_02.tif
│   └── type_2_02.tif
├── 002
│   ├── type_1_01.tif
│   ├── type_1_02.tif
│   ├── type_2_01.tif
│   └── type_2_02.tif
└── 003
    └── imaging
        ├── type_1_01.tif
        └── type_1_02.tif
```

**Output:**
```
derivatives/
├── sub-001_id-001/
│   └── ses-001/
│       └── funcimg/
│           └── {output_pattern}type_1_01.tif
│           └── {output_pattern}type_1_02.tif
├── sub-001_id-001/
│   └── ses-002/
│       └── funcimg/
│           └── {output_pattern}type_2_02.tif
├── sub-002_id-002/
│   └── ses-001/
│       └── funcimg/
│           └── {output_pattern}type_1_01.tif
│           └── {output_pattern}type_1_02.tif
├── sub-002_id-002/
│   └── ses-002/
│       └── funcimg/
│           └── {output_pattern}type_2_01.tif
│           └── {output_pattern}type_2_02.tif
├── sub-003_id-003/
│   └── ses-001/
│       └── funcimg/
│           └── {output_pattern}type_1_01.tif
│           └── {output_pattern}type_1_02.tif
```

### Custom Metadata Dataset

**Input:**
```
└── mouse-001_genotype-WT_age-P60_treatment-saline
    └── session-001_condition-baseline_paradigm-open-field
        └── recording.tif
```

**Output:**
```
derivatives/
├── sub-001_id-mouse-001_genotype-WT_age-P60_treatment-saline/
│   └── ses-001_condition-baseline_paradigm-open-field/
│       └── funcimg/
│           └── {output_pattern}recording.tif
```

### NeuroBlueprint Format Dataset

**Input:**
```
├── sub-001_strain-C57BL6_sex-M
│   └── ses-001_date-20250225_protocol-training
│       └── recording.tif
└── sub-002_strain-C57BL6_sex-F
    └── ses-001_date-20250226_protocol-testing
        └── recording.tif
```

**Output:**
```
derivatives/
├── sub-001_strain-C57BL6_sex-M/
│   └── ses-001_date-20250225_protocol-training/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-002_strain-C57BL6_sex-F/
│   └── ses-001_date-20250226_protocol-testing/
│       └── funcimg/
│           └── {output_pattern}recording.tif
```

### Non-continuous NeuroBlueprint Dataset

**Input:**
```
├── sub-005_strain-BALBC_sex-M
│   ├── ses-001_date-20250221_protocol-test
│   │   └── recording.tif
│   ├── ses-003_date-20250223_protocol-test
│   │   └── recording.tif
│   └── ses-007_date-20250227_protocol-test
│       └── recording.tif
├── sub-010_strain-BALBC_sex-F
│   ├── ses-002_date-20250222_protocol-test
│   │   └── recording.tif
│   └── ses-005_date-20250225_protocol-test
│       └── recording.tif
└── sub-025_strain-BALBC_sex-M
    ├── ses-001_date-20250221_protocol-test
    │   └── recording.tif
    ├── ses-004_date-20250224_protocol-test
    │   └── recording.tif
    ├── ses-008_date-20250228_protocol-test
    │   └── recording.tif
    └── ses-009_date-20250220_protocol-test
        └── recording.tif
```

**Output:**
```
derivatives/
├── sub-005_strain-BALBC_sex-M/
│   └── ses-001_date-20250221_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-005_strain-BALBC_sex-M/
│   └── ses-003_date-20250223_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-005_strain-BALBC_sex-M/
│   └── ses-007_date-20250227_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-010_strain-BALBC_sex-F/
│   └── ses-002_date-20250222_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-010_strain-BALBC_sex-F/
│   └── ses-005_date-20250225_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-025_strain-BALBC_sex-M/
│   └── ses-001_date-20250221_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-025_strain-BALBC_sex-M/
│   └── ses-004_date-20250224_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-025_strain-BALBC_sex-M/
│   └── ses-008_date-20250228_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
├── sub-025_strain-BALBC_sex-M/
│   └── ses-009_date-20250220_protocol-test/
│       └── funcimg/
│           └── {output_pattern}recording.tif
```

## Key Rules

- **Custom formats** → `sub-{INDEX:03d}_id-{original-name}` + `ses-000`, `ses-001`...
- **NeuroBlueprint formats** → Preserve original names + session IDs
- **All outputs** → `derivatives/{subject}/{session}/funcimg/`

**Regenerate:** `pytest tests/test_unit/test_generate_readme.py::test_generate_readme -v`
