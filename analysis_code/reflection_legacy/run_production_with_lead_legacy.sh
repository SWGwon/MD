#!/bin/bash

# 1. 사용할 스크립트 경로 설정
PRODUCTION_SCRIPT="/home/sgwon/MD/analysis_scripts/calculate_charge.C"

# 2. 1부터 9까지 루프 시작
for i in {2..9}
do
    # 파일 경로 생성 (test_1.root, test_2.root ...)
    FILE="/mnt/H54T/MD/17cm/caen/with_lead/test_with_lead_v2_${i}.root"


    echo "------------------------------------------------------"
    echo "Running analysis for: test_${i}.root"
    echo "------------------------------------------------------"

    # 3. ROOT 실행
    # 큰따옴표 내부에 변수를 넣을 때 역슬래시(\)를 써서 ROOT가 경로를 문자열로 인식하게 합니다.
    root -b -q "${PRODUCTION_SCRIPT}(\"${FILE}\", 2)"

done

echo "All tasks are completed!"
