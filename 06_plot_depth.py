import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
import re

# HTML 파일 읽기
with open('version-2_integrated/table_summary_stats/sample-frequency-detail.html', 'r') as f:
    html_content = f.read()

# JSON 데이터 추출
json_match = re.search(r'<script id="table-data" type="application/json">(.*?)</script>', html_content, re.DOTALL)
if json_match:
    json_data = json.loads(json_match.group(1))
    # Frequency 딕셔너리를 DataFrame으로 변환
    df = pd.DataFrame(list(json_data['Frequency'].items()), columns=['SampleID', 'Frequency'])
    df = df.sort_values('Frequency')
else:
    print("Error: Could not find JSON data in HTML")
    exit(1)

# 시각화 설정
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# 1. 히스토그램 (서열 개수 분포)
ax1.hist(df['Frequency'], bins=50, color='skyblue', edgecolor='black', alpha=0.7)
ax1.set_title('Distribution of Sequence Counts per Sample', fontsize=15, fontweight='bold')
ax1.set_xlabel('Number of Sequences', fontsize=12)
ax1.set_ylabel('Number of Samples', fontsize=12)
ax1.grid(axis='y', linestyle='--', alpha=0.5)

# 2. 샘플 유지 곡선 (Sample Retention Curve)
depths = np.sort(df['Frequency'].unique())
retention = [sum(df['Frequency'] >= d) for d in depths]

ax2.plot(depths, retention, linestyle='-', color='salmon', linewidth=2)
ax2.set_title('Sample Retention by Sequencing Depth', fontsize=15, fontweight='bold')
ax2.set_xlabel('Sequencing Depth (Threshold)', fontsize=12)
ax2.set_ylabel('Number of Remaining Samples', fontsize=12)

# 주요 통계량
median_val = df['Frequency'].median()
min_val = df['Frequency'].min()
total_samples = len(df)

# 샘플의 90%를 유지하는 최대 Depth 찾기
target_retention = int(total_samples * 0.9)
recommended_depth = df.iloc[total_samples - target_retention]['Frequency']

ax2.axvline(x=recommended_depth, color='blue', linestyle='-.', label=f'90% Retention Depth: {int(recommended_depth):,}')
ax2.axvline(x=median_val, color='green', linestyle='--', label=f'Median: {int(median_val):,}')
ax2.axvline(x=min_val, color='red', linestyle=':', label=f'Min: {int(min_val):,}')

ax2.legend(fontsize=11)
ax2.grid(linestyle='--', alpha=0.5)

# 하단에 텍스트 정보 추가
info_text = (f"Total Samples: {total_samples}\n"
             f"Min: {int(min_val):,} | Max: {int(df['Frequency'].max()):,}\n"
             f"Median: {int(median_val):,}\n"
             f"Recommended Depth (90% samples): {int(recommended_depth):,}")
plt.figtext(0.5, 0.01, info_text, ha="center", fontsize=12, bbox={"facecolor":"orange", "alpha":0.2, "pad":5})

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.savefig('version-2_integrated/sequencing_depth_distribution.png')

# 콘솔 출력용
print(info_text)
print("\nSamples with very low frequency (< 1000):")
low_freq = df[df['Frequency'] < 1000]
print(low_freq)
