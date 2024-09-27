import os
import argparse
import sys
import pandas as pd
import subprocess
import time
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../monitor'))
#sys.path.append('/opt/software/tagFinder-inhouse/bin/monitor')
#sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'D:\Project_Pharmablock\/01.tagFinder-inhouse\/bin\monitor'))
import brc

def tag_dwar(count_file, ratio_file):
    #计算counts的比例，enrichment，以及生成dwar文件
    count_f = pd.read_csv(count_file, sep='\t', dtype=object, header=0)
    count_f['Raw_count'] = count_f['Raw_count'].astype(int)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    count_f_mean = count_f['Raw_count'].mean()
    count_f_sum = count_f['Raw_count'].sum()
    count_f['RAW_ratio'] = count_f['Raw_count']/count_f_sum
    count_f['RAW_ratio2'] = (count_f['Raw_count'] - count_f_mean)/count_f_sum
    count_f['RAW_z_score'] = brc.z_score(count_f['Raw_count'], 0)
    count_f_dedup_mean = count_f['Dedup_count'].mean()
    count_f_dedup_sum = count_f['Dedup_count'].sum()
    count_f['Dedup_ratio'] = count_f['Dedup_count'] / count_f_dedup_sum
    count_f['Dedup_ratio2'] = (count_f['Dedup_count'] - count_f_dedup_mean) / count_f_dedup_sum
    count_f['Dedup_z_score'] = brc.z_score(count_f['Dedup_count'], 0)
    print('the raw median counts is {}'.format(count_f['Raw_count'].median()))
    print('the deduped median counts is {}'.format(count_f['Dedup_count'].median()))
    print('the Raw total counts is {}'.format(count_f_sum))
    print('the deduped total counts is {}'.format(count_f_dedup_sum))
    #display index as compound ID
    #count_f.sort_values(by='Dedup_count', ascending=False).round(4).to_csv(ratio_file, sep='\t', index=True)
    data = count_f.sort_values(by='Dedup_count', ascending=False).round(4)
    data.reset_index(drop=True, inplace=True)
    data.index = range(1, len(data)+1)
    data.to_csv(ratio_file, sep='\t', index=True)
    return count_f

def dedup_filter(ratio_file):
    count_f = pd.read_csv(ratio_file, sep='\t', dtype=object, header=0)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    dedup_max = count_f['Dedup_count'].max()
    BBmin = 0
    BB1, BB2, BB3 = [], [], []
    BB_lst = []
    dedup_filter_threshold = 1
    test = count_f[count_f['Dedup_count'] >= 1]
    BB1 = set(test['TAG2'].tolist())
    BB2 = set(test['TAG3'].tolist())
    BB3 = set(test['TAG4'].tolist())
    BBmin_pre = min(len(BB1), len(BB2), len(BB3))
    dedup_max_pre = count_f['Dedup_count'].max()
    # filtering for displaying dot, line, plane
    while (dedup_max >= 1):
        new_dataframe = count_f[count_f['Dedup_count'] >= (dedup_max)]
        BB1 = set(new_dataframe['TAG2'].tolist())
        BB2 = set(new_dataframe['TAG3'].tolist())
        BB3 = set(new_dataframe['TAG4'].tolist())
        BBmin = min(len(BB1), len(BB2), len(BB3))
        if BBmin < 50:
            dedup_max = dedup_max / 2 + 0.1
        elif BBmin > 100:
            dedup_max = (dedup_max + dedup_max_pre) / 2
            if (dedup_max_pre - dedup_max) <= 1:
                dedup_filter_threshold = dedup_max
                dedup_max = 0.5
            dedup_max_pre = dedup_max
        elif BBmin >= 50:
            dedup_filter_threshold = dedup_max
            dedup_max = 0.5
    return dedup_filter_threshold


def txt2dwar(ratio_file, dwar_file):
    count_f = pd.read_csv(ratio_file, sep='\t', dtype=object, header=0)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    dedup_max = count_f['Dedup_count'].max()
    dedup_filter_threshold = dedup_filter(ratio_file)

#定义dwar文件的头尾
    tag_header = '''<datawarrior-fileinfo>
<version="3.2">
<created="1559631713662">
<rowcount="{0}">
</datawarrior-fileinfo>
'''.format(count_f.shape[0])

#dwar格式要求特别严，另外文本处理时对tab格式必须从外部指定，内部输入tab时在输出时变为空格。
#这里设置dedup的count数过滤参数，这样大概得到点线面的关系图
    tag_tail = '''<datawarrior properties>
<axisColumn_2D View_0="Column 1">
<axisColumn_2D View_1="Raw_count">
<axisColumn_3D Dedup_counts_0="TAG2">
<axisColumn_3D Dedup_counts_1="TAG3">
<axisColumn_3D Dedup_counts_2="TAG4">
<chartType_2D View="scatter">
<chartType_3D Dedup_counts="scatter">
<chartType_3D Raw_counts="scatter">
<columnWidth_Table_Dedup_count="80">
<columnWidth_Table_Dedup_ratio="80">
<columnWidth_Table_Dedup_ratio2="80">
<columnWidth_Table_Dedup_z_score="80">
<columnWidth_Table_Library_ID="80">
<columnWidth_Table_Protein_CC_ID="80">
<columnWidth_Table_RAW_ratio="80">
<columnWidth_Table_RAW_ratio2="80">
<columnWidth_Table_RAW_z_score="80">
<columnWidth_Table_Raw_count="80">
<columnWidth_Table_TAG2="80">
<columnWidth_Table_TAG2_seq="80">
<columnWidth_Table_TAG2_structure="80">
<columnWidth_Table_TAG3="80">
<columnWidth_Table_TAG3_seq="80">
<columnWidth_Table_TAG3_structure="80">
<columnWidth_Table_TAG4="80">
<columnWidth_Table_TAG4_seq="80">
<columnWidth_Table_TAG4_structure="80">
<detailView="height[Data]=1">
<faceColor3D_3D Dedup_counts="-1250054">
<faceColor3D_3D Raw_counts="-1250054">
<fastRendering_2D View="true">
<fastRendering_3D Dedup_counts="true">
<filter0="#browser#	#disabled#	TAG2	1.001">
<filter1="#double#	TAG2">
<filter10="#double#	Dedup_count{0}{1}{2}{3}">
<filter11="#double#	Dedup_ratio">
<filter12="#double#	Dedup_z_score">
<filter13="#double#	Raw_count">
<filter14="#double#	RAW_ratio">
<filter15="#double#	RAW_z_score">
<filter2="#string#	TAG2_seq">
<filter3="#string#	TAG2_structure">
<filter4="#double#	TAG3">
<filter5="#string#	TAG3_seq">
<filter6="#string#	TAG3_structure">
<filter7="#double#	TAG4">
<filter8="#string#	TAG4_seq">
<filter9="#string#	TAG4_structure">
<filter16="#double#	TAG4">
<filter17="#string#	TAG4_seq">
<filter18="#string#	TAG4_structure">
<filter19="#string#	TAG2_DLP">
<filter20="#string#	TAG3_DLP">
<filter21="#string#	TAG4_DLP">
<filter22="#string#	TAG_DOT">
<filterAnimation0="state=stopped delay=500">
<mainSplitting="0.7227">
<mainView="3D Raw_counts">
<mainViewCount="4">
<mainViewDockInfo0="root">
<mainViewDockInfo1="Table	bottom	0.5">
<mainViewDockInfo2="2D View	center">
<mainViewDockInfo3="3D Raw_counts	right	0.497">
<mainViewName0="Table">
<mainViewName1="2D View">
<mainViewName2="3D Raw_counts">
<mainViewName3="3D Dedup_counts">
<mainViewType0="tableView">
<mainViewType1="2Dview">
<mainViewType2="3Dview">
<mainViewType3="3Dview">
<masterView_3D Raw_counts="3D Dedup_counts">
<rightSplitting="0.62829">
<rotationMatrix_3D Dedup_counts00="0.99021">
<rotationMatrix_3D Dedup_counts01="-0.082155">
<rotationMatrix_3D Dedup_counts02="-0.11283">
<rotationMatrix_3D Dedup_counts10="0.0938">
<rotationMatrix_3D Dedup_counts11="0.99034">
<rotationMatrix_3D Dedup_counts12="0.10212">
<rotationMatrix_3D Dedup_counts20="0.10335">
<rotationMatrix_3D Dedup_counts21="-0.11171">
<rotationMatrix_3D Dedup_counts22="0.98836">
<rowHeight_Table="16">
<showNaNValues_2D View="true">
<showNaNValues_3D Dedup_counts="true">
<showNaNValues_3D Raw_counts="true">
<sizeColumn_3D Dedup_counts="Dedup_count">
<sizeColumn_3D Raw_counts="Raw_count">
</datawarrior properties>

'''.format('\t', dedup_filter_threshold, '\t', count_f['Dedup_count'].max())
    #生成dwar文件
    with open(ratio_file, "r") as f:
        old = f.read()
        with open(dwar_file, 'w') as dwar:
            dwar.seek(0)
            dwar.write(tag_header)
            dwar.write(old)
            dwar.write(tag_tail)
    #temp = pd.read_csv(dwar_file, sep='\t', dtype=object, header=None)
    #temp.to_csv('D:\Project_Pharmablock\/01.tagFinder-inhouse\/bin\pipeline', sep='\t', index=False)
    # return dedup_filter

def DotLinePlane(ratio_file, ratio_file_filter):
    count_f = pd.read_csv(ratio_file, sep='\t', dtype=object, header=0)
    count_f['Dedup_count'] = count_f['Dedup_count'].astype(int)
    dedup_filter_threshold = dedup_filter(ratio_file)
    dedup_mean = count_f['Dedup_count'].mean()
    count_f['TAG2_DLP'] = str(0)
    count_f['TAG3_DLP'] = str(0)
    count_f['TAG4_DLP'] = str(0)
    count_f['TAG_DOT'] = str(0)
    count_f_filter = count_f[count_f['Dedup_count'] >= dedup_filter_threshold]
    BB1 = set(count_f_filter['TAG2'].tolist())
    BB2 = set(count_f_filter['TAG3'].tolist())
    BB3 = set(count_f_filter['TAG4'].tolist())
    BB1_len = len(BB1)
    BB2_len = len(BB2)
    BB3_len = len(BB3)
    count_f_filter.loc[(count_f_filter['Dedup_count'] >= 1000 * dedup_mean), 'TAG_DOT'] = str(2)  # 判断清晰的点
    count_f_filter.loc[(count_f_filter['Dedup_count'] < 1000 * dedup_mean) & (count_f_filter['Dedup_count'] >= 100 * dedup_mean), 'TAG_DOT'] = str(
        1)  # 判断点，1.10表示tag1，点，清晰程度是0

    # print(BB1)
    # print(BB2)
    # print(BB3)
    for Tag1 in BB1:
        count_f_filter_tag1 = count_f_filter[count_f_filter['TAG2'] == Tag1]
        BB1_f1 = count_f_filter_tag1.shape[0] # 同一个BB1的行数，即
        BB2_f1 = set(count_f_filter_tag1['TAG3'].tolist())
        BB3_f1 = set(count_f_filter_tag1['TAG4'].tolist())
        if (BB1_f1 == 1):
            count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['Dedup_count'] >= 100*dedup_mean), 'TAG2_DLP'] = str(2.11) #判断清晰的点
            count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['Dedup_count'] < 100*dedup_mean),'TAG2_DLP'] = str(2.10) # 判断点，1.10表示tag1，点，清晰程度是0
        elif (len(BB2_f1) - BB1_f1 == 0) and (len(BB3_f1) - BB1_f1 == 0):
            # 任意两点都不在一条线上就认为不成面（即分散的点），如果2个点以上才算
            count_f_filter.loc[count_f_filter['TAG2'] == Tag1, 'TAG2_DLP'] = str(0)
        else:
            if len(BB3_f1) >= (0.9 * BB3_len) and len(BB2_f1) >= (0.9 * BB2_len):
                count_f_filter.loc[count_f_filter['TAG2'] == Tag1,'TAG2_DLP'] = str(2.34)
                # 判断面
            elif (len(BB3_f1) >= (0.9 * BB3_len) and len(BB2_f1) >= (0.8 * BB2_len)) or (len(BB3_f1) >= (0.8 * BB3_len) and len(BB2_f1) >= (0.9 * BB2_len)):
                count_f_filter.loc[count_f_filter['TAG2'] == Tag1,'TAG2_DLP'] = str(2.33)
            elif len(BB3_f1) >= (0.8 * BB3_len) and len(BB2_f1) >= (0.8 * BB2_len):
                count_f_filter.loc[count_f_filter['TAG2'] == Tag1,'TAG2_DLP'] = str(2.32)
            elif len(BB2_f1) <= (0.2 * BB2_len) and len(BB3_f1) >= (0.8 * BB3_len):
                for Tag2 in BB2_f1:
                    BB3_f2 = count_f_filter_tag1[count_f_filter_tag1['TAG3']==Tag2].shape[0]
                    # 当BB1和BB2固定时，判断BB3的个数
                    if BB3_f2 >= (0.8 * BB3_len):
                        count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG3'] == Tag2),'TAG2_DLP'] = str(2.24)
                        # 判断线，很明显的线.1.21：1表示tag1,2表示线，1表示线明显程度。明显程度逐渐递增（0,1,2）
                    elif BB3_f2 >= (0.8 * len(BB3_f1)):
                        count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG3'] == Tag2),'TAG2_DLP'] = str(2.23)
                        # 判断线，相对没那么明显的线
            elif len(BB3_f1) <= (0.2 * BB3_len) and len(BB2_f1) >= (0.8 * BB2_len):
                for Tag3 in BB3_f1:
                    BB2_f2 = count_f_filter_tag1[count_f_filter_tag1['TAG4']==Tag3].shape[0]
                    # 当BB1和BB3固定时，判断BB2的个数
                    if BB2_f2 >= (0.8 * BB2_len):
                        count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG4'] == Tag3),'TAG2_DLP'] = str(2.24)
                        # 判断线，很明显的线.1.21：1表示tag1,2表示线，1表示线明显程度。明显程度逐渐递增（0,1,2）
                    elif BB2_f2 >= (0.8 * len(BB2_f1)):
                        count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG4'] == Tag3),'TAG2_DLP'] = str(2.23)
                        # 判断线，相对没那么明显的线
            elif (len(BB2_f1) >= max(3, (0.5 * BB2_len))) or (len(BB3_f1) >= max(3, (0.5 * BB3_len))):
                count_f_filter.loc[(count_f_filter['TAG2'] == Tag1), 'TAG2_DLP'] = str(2.22)

            # 如果面和线都不明显，一条线上只有3个点或以上时，DLP对应的值都为零时，缩小条件重新判定点线面
            elif BB1_f1 >= 3:
                # 判断不清晰的面
                if len(BB2_f1) >= max(3, 0.2*BB1_f1) and len(BB3_f1) >= max(3, 0.2*BB1_f1):
                    count_f_filter.loc[count_f_filter['TAG2'] == Tag1, 'TAG2_DLP'] = str(2.31)
                elif len(BB2_f1) >= max(3, 0.2*BB1_f1) or len(BB3_f1) >= max(3, 0.2*BB1_f1):
                    count_f_filter.loc[count_f_filter['TAG2'] == Tag1, 'TAG2_DLP'] = str(2.30)
                elif len(BB2_f1) < max(3, 0.2*BB1_f1):
                    for Tag2 in BB2_f1:
                        BB3_f3 = count_f_filter_tag1[count_f_filter_tag1['TAG3']==Tag2].shape[0]
                        if BB3_f3 >= 3:
                            count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG3'] == Tag2), 'TAG2_DLP'] = str(2.21)
                        else:
                            count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG3'] == Tag2), 'TAG2_DLP'] = str(2.20)
                elif len(BB3_f1) < max(3, 0.2*BB1_f1):
                    for Tag3 in BB2_f1:
                        BB2_f3 = count_f_filter_tag1[count_f_filter_tag1['TAG4']==Tag3].shape[0]
                        if BB2_f3 >= 3:
                            count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG4'] == Tag3), 'TAG2_DLP'] = str(2.21)
                        else:
                            count_f_filter.loc[(count_f_filter['TAG2'] == Tag1) & (count_f_filter['TAG4'] == Tag3), 'TAG2_DLP'] = str(2.20)

######################################################################################################################


    for Tag2 in BB2:
        count_f_filter_tag2 = count_f_filter[count_f_filter['TAG3'] == Tag2]
        BB2_f1 = count_f_filter_tag2.shape[0] # 同一个BB2的行数，即
        BB1_f1 = set(count_f_filter_tag2['TAG2'].tolist())
        BB3_f1 = set(count_f_filter_tag2['TAG4'].tolist())
        if BB2_f1 == 1:
            count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (count_f_filter['Dedup_count'] >= 100 * dedup_mean), 'TAG3_DLP'] = str(3.11)
            count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (count_f_filter['Dedup_count'] < 100 * dedup_mean),'TAG3_DLP'] = str(3.10) # 判断点，2.10表示tag2，点，清晰程度是0
        elif len(BB1_f1) == BB2_f1 and len(BB3_f1) == BB2_f1:
            # 任意两点都不在一条线上就认为不成面
            count_f_filter.loc[count_f_filter['TAG3'] == Tag2, 'TAG3_DLP'] = str(0)
        else:
            if len(BB3_f1) >= (0.9 * BB3_len) and len(BB1_f1) >= (0.9 * BB1_len):
                count_f_filter.loc[count_f_filter['TAG3'] == Tag2,'TAG3_DLP'] = str(3.34)
                # 判断面
            elif (len(BB3_f1) >= (0.9 * BB3_len) and len(BB1_f1) >= (0.8 * BB1_len)) or (len(BB3_f1) >= (0.8 * BB3_len) and len(BB1_f1) >= (0.9 * BB1_len)):
                count_f_filter.loc[count_f_filter['TAG3'] == Tag2,'TAG3_DLP'] = str(3.33)
            elif len(BB3_f1) >= (0.8 * BB3_len) and len(BB1_f1) >= (0.8 * BB1_len):
                count_f_filter.loc[count_f_filter['TAG3'] == Tag2,'TAG3_DLP'] = str(3.32)
            elif len(BB1_f1) <= (0.2 * BB1_len) and len(BB3_f1) >= (0.8 * BB3_len):
                for Tag1 in BB1_f1:
                    BB3_f2 = count_f_filter_tag2[count_f_filter_tag2['TAG2']==Tag1].shape[0]
                    # 当BB1和BB2固定时，判断BB3的个数
                    if BB3_f2 >= (0.8 * BB3_len):
                        count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (count_f_filter['TAG2'] == Tag1),'TAG3_DLP'] = str(3.24)
                        # 判断线，很明显的线.2.21：2表示tag2,2表示线，1表示线明显程度。明显程度逐渐递增（0,1,2）
                    elif BB3_f2 >= (0.8 * len(BB3_f1)):
                        count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (count_f_filter['TAG2'] == Tag1),'TAG3_DLP'] = str(3.23)
                        # 判断线，相对没那么明显的线
            elif len(BB3_f1) <= (0.2 * BB3_len) and len(BB1_f1) >= (0.8 * BB1_len):
                for Tag3 in BB3_f1:
                    BB1_f2 = count_f_filter_tag2[count_f_filter_tag2['TAG4']==Tag3].shape[0]
                    # 当BB2和BB3固定时，判断BB2的个数
                    if BB1_f2 >= (0.8 * BB1_len):
                        count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (count_f_filter['TAG4'] == Tag3),'TAG3_DLP'] = str(3.24)
                        # 判断线，很明显的线.1.21：1表示tag1,2表示线，1表示线明显程度。明显程度逐渐递增（0,1,2）
                    elif BB1_f2 >= (0.8 * len(BB1_f1)):
                        count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (count_f_filter['TAG4'] == Tag3),'TAG3_DLP'] = str(3.23)
                        # 判断线，相对没那么明显的线
            elif (len(BB1_f1) >= max(3, (0.5 * BB1_len))) or (len(BB3_f1) >= max(3, (0.5 * BB3_len))):
                count_f_filter.loc[(count_f_filter['TAG3'] == Tag2), 'TAG3_DLP'] = str(3.22)

            # 如果面和线都不明显，一条线上只有3个点或以上时，DLP对应的值都为零时，缩小条件重新判定点线面
            elif BB2_f1 >= 3:
                # 判断不清晰的面
                if len(BB1_f1) >= max(3, 0.2 * BB2_f1) and len(BB3_f1) >= max(3, 0.2 * BB2_f1):
                    count_f_filter.loc[count_f_filter['TAG3'] == Tag2, 'TAG3_DLP'] = str(3.31)
                elif len(BB1_f1) >= max(3, 0.2 * BB2_f1) or len(BB3_f1) >= max(3, 0.2 * BB2_f1):
                    count_f_filter.loc[count_f_filter['TAG3'] == Tag2, 'TAG3_DLP'] = str(3.30)
                elif len(BB1_f1) < max(3, 0.2 * BB2_f1):
                    for Tag1 in BB1_f1:
                        BB3_f3 = count_f_filter_tag2[count_f_filter_tag2['TAG2'] == Tag1].shape[0]
                        if BB3_f3 >= 3:
                            count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (
                                        count_f_filter['TAG2'] == Tag1), 'TAG3_DLP'] = str(3.21)
                        else:
                            count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (
                                        count_f_filter['TAG2'] == Tag1), 'TAG3_DLP'] = str(3.20)
                elif len(BB3_f1) < max(3, 0.2 * BB2_f1):
                    for Tag3 in BB1_f1:
                        BB1_f3 = count_f_filter_tag2[count_f_filter_tag2['TAG4'] == Tag3].shape[0]
                        if BB1_f3 >= 3:
                            count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (
                                        count_f_filter['TAG4'] == Tag3), 'TAG3_DLP'] = str(3.21)
                        else:
                            count_f_filter.loc[(count_f_filter['TAG3'] == Tag2) & (
                                        count_f_filter['TAG4'] == Tag3), 'TAG3_DLP'] = str(3.20)

    ######################################################################################################################


    for Tag3 in BB3:
        count_f_filter_tag3 = count_f_filter[count_f_filter['TAG4'] == Tag3]
        BB3_f1 = count_f_filter_tag3.shape[0] # 同一个BB1的行数，即
        BB2_f1 = set(count_f_filter_tag3['TAG3'].tolist())
        BB1_f1 = set(count_f_filter_tag3['TAG2'].tolist())
        if BB3_f1 == 1:
            count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (count_f_filter['Dedup_count'] >= 100 * dedup_mean), 'TAG4_DLP'] = str(4.11)
            count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (count_f_filter['Dedup_count'] < 100 * dedup_mean),'TAG4_DLP'] = str(4.10) # 判断点，1.10表示tag1，点，清晰程度是0
        elif len(BB2_f1) == BB3_f1 and len(BB1_f1) == BB3_f1:
            # 任意两点都不在一条线上就认为不成面
            count_f_filter.loc[count_f_filter['TAG4'] == Tag3, 'TAG4_DLP'] = str(0)
        else:
            if len(BB1_f1) >= (0.9 * BB1_len) and len(BB2_f1) >= (0.9 * BB2_len):
                count_f_filter.loc[count_f_filter['TAG4'] == Tag3,'TAG4_DLP'] = str(4.32)
                # 判断面
            elif (len(BB1_f1) >= (0.9 * BB1_len) and len(BB2_f1) >= (0.8 * BB2_len)) or (len(BB1_f1) >= (0.8 * BB1_len) and len(BB2_f1) >= (0.9 * BB2_len)):
                count_f_filter.loc[count_f_filter['TAG4'] == Tag3,'TAG4_DLP'] = str(4.31)
            elif len(BB1_f1) >= (0.8 * BB1_len) and len(BB2_f1) >= (0.8 * BB2_len):
                count_f_filter.loc[count_f_filter['TAG4'] == Tag3,'TAG4_DLP'] = str(4.30)
            elif len(BB2_f1) <= (0.2 * BB2_len) and len(BB1_f1) >= (0.8 * BB1_len):
                for Tag2 in BB2_f1:
                    BB1_f2 = count_f_filter_tag3[count_f_filter_tag3['TAG3']==Tag2].shape[0]
                    # 当BB1和BB2固定时，判断BB3的个数
                    if BB1_f2 >= (0.8 * BB1_len):
                        count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (count_f_filter['TAG3'] == Tag2),'TAG4_DLP'] = str(4.21)
                        # 判断线，很明显的线.1.21：1表示tag1,2表示线，1表示线明显程度。明显程度逐渐递增（0,1,2）
                    elif BB1_f2 >= (0.8 * len(BB1_f1)):
                        count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (count_f_filter['TAG3'] == Tag2),'TAG4_DLP'] = str(4.20)
                        # 判断线，相对没那么明显的线
            elif len(BB1_f1) <= (0.2 * BB1_len) and len(BB2_f1) >= (0.8 * BB2_len):
                for Tag1 in BB1_f1:
                    BB2_f2 = count_f_filter_tag3[count_f_filter_tag3['TAG2']==Tag1].shape[0]
                    # 当BB1和BB3固定时，判断BB2的个数
                    if BB2_f2 >= (0.8 * BB2_len):
                        count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (count_f_filter['TAG2'] == Tag1),'TAG4_DLP'] = str(4.21)
                        # 判断线，很明显的线.1.21：1表示tag1,2表示线，1表示线明显程度。明显程度逐渐递增（0,1,2）
                    elif BB2_f2 >= (0.8 * len(BB2_f1)):
                        count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (count_f_filter['TAG2'] == Tag1),'TAG4_DLP'] = str(4.20)
                        # 判断线，相对没那么明显的线

            elif (len(BB2_f1) >= max(3, (0.5 * BB2_len))) or (len(BB1_f1) >= max(3, (0.5 * BB1_len))):
                count_f_filter.loc[(count_f_filter['TAG4'] == Tag3), 'TAG4_DLP'] = str(4.22)

            # 如果面和线都不明显，一条线上只有3个点或以上时，DLP对应的值都为零时，缩小条件重新判定点线面
            elif BB3_f1 >= 3:
                # 判断不清晰的面
                if len(BB1_f1) >= max(3, 0.2 * BB3_f1) and len(BB2_f1) >= max(3, 0.2 * BB3_f1):
                    count_f_filter.loc[count_f_filter['TAG4'] == Tag3, 'TAG4_DLP'] = str(4.31)
                elif len(BB1_f1) >= max(3, 0.2 * BB3_f1) or len(BB2_f1) >= max(3, 0.2 * BB3_f1):
                    count_f_filter.loc[count_f_filter['TAG4'] == Tag3, 'TAG4_DLP'] = str(4.30)
                elif len(BB1_f1) < max(3, 0.2 * BB3_f1):
                    for Tag1 in BB1_f1:
                        BB2_f3 = count_f_filter_tag3[count_f_filter_tag3['TAG2'] == Tag1].shape[0]
                        if BB2_f3 >= 3:
                            count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (
                                        count_f_filter['TAG2'] == Tag1), 'TAG4_DLP'] = str(4.21)
                        else:
                            count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (
                                        count_f_filter['TAG2'] == Tag1), 'TAG4_DLP'] = str(4.20)
                elif len(BB2_f1) < max(3, 0.2 * BB3_f1):
                    for Tag2 in BB1_f1:
                        BB1_f3 = count_f_filter_tag3[count_f_filter_tag3['TAG3'] == Tag2].shape[0]
                        if BB1_f3 >= 3:
                            count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (
                                        count_f_filter['TAG3'] == Tag2), 'TAG4_DLP'] = str(4.21)
                        else:
                            count_f_filter.loc[(count_f_filter['TAG4'] == Tag3) & (
                                        count_f_filter['TAG3'] == Tag2), 'TAG4_DLP'] = str(4.20)

    ######################################################################################################################


    count_f_filter.to_csv(ratio_file_filter, sep='\t', index=False)

def datawarrior_structure(EnrichmentFile3, DwarFile4):
    FilterFileName = os.path.basename(EnrichmentFile3)
    libraryName = FilterFileName.split('.')[0]
    dwam_file = os.path.dirname(EnrichmentFile3) + '/' + libraryName + '.dwam'
    template = """<macro name="{0}">
<task name="openFile">
fileName={1}
</task>
<task name="selectWindow">
viewName={2}
</task>
<task name="assignOrZoomAxes">
column1=TAG3
column0=TAG2
viewName=3D View
high0=0.0
high1=0.0
column2=TAG4
high2=0.0
low2=1.0
millis=1000
low0=1.0
low1=1.0
</task>
<task name="setMarkerSize">
inverse=false
adaptive=true
viewName=3D View
size=1.0
column=Dedup_count
proportional=false
</task>
<task name="saveFileAs">
fileName={3}
</task>
<task name="exitProgram">
saveChanges=ask
</task>
</macro>
""".format(libraryName, EnrichmentFile3, FilterFileName, DwarFile4)
    with open(dwam_file, 'w') as dwam:
        dwam.write(template)
    cmd = "/opt/datawarrior/datawarrior {}".format(dwam_file)
    subprocess.check_call(cmd, shell=True)

def main(EnrichmentFile, EnrichmentFile3, DwarFile4):
    # tag_dwar(CountFile, EnrichmentFile)
    DotLinePlane(EnrichmentFile, EnrichmentFile3)
    # txt2dwar(EnrichmentFile, DwarFile)
    txt2dwar(EnrichmentFile3, DwarFile4) # 生成dwar文件
    # datawarrior_structure_windows(EnrichmentFile3) # 仅生成dwam文件，用的时候在H盘中生成新的带结构信息的dwar文件
    # datawarrior_structure(EnrichmentFile3, DwarFile4) # 生成带有结构的dwar文件


if __name__ == "__main__":
    data_start = time.time()
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    parser = argparse.ArgumentParser(description='Analyze data and Conver count file to dwar file')
    # parser.add_argument('--input1', '-i1', dest="CountFile", help="Count file containing raw and deduped counts")
    parser.add_argument('--input', '-i', dest="EnrichmentFile", help="Enrichment analysis file output")
    # parser.add_argument('--output_dwar', '-o2', dest="DwarFile", help="Dwar file output")
    parser.add_argument('--output3', '-o3', dest="EnrichmentFile3", help="Enrichment analysis file output by filtering")
    parser.add_argument('--output_dwar4', '-o4', dest="DwarFile4", help="Dwar file output by filtering")
    #parser.add_argument('--config', '-c', help='config file')
    args = parser.parse_args()
    main(args.EnrichmentFile, args.EnrichmentFile3, args.DwarFile4)

    #tag_dwar(args.CountFile, args.EnrichmentFile, args.DwarFile)
    #txt2dwar(args.EnrichmentFile, args.DwarFile)
    #txt2dwar('0001865-196.count.freq.txt', '0001865-196.count.freq.dwar')
    #DotLinePlane('PBDEL0014-2-CP0785.count.freq.txt', 'PBDEL0014-2-CP0785.count.freq_filter.txt')

