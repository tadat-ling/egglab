################################################################################
#      EGGLab: Processing Electroglottography Signal with Praat (ver 1.0)      #
#                                                                              #
#                    Written in April 2025 by Dat-Danny Ta                     #
#                           Last updated: 2025.06.18                           #
#                This script is released under the GNU GPL v3.0                #
#                           © 2025.04.30 Ta Van Dat                            #
#                tavandat@stu.pku.edu.cn | datbg.0702@gmail.com                #
#                                                                              #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                                                              #
#                                 INTRODUCTION                                 #
#                                                                              #
# This program was inspired by EGG tools such as praatdet, peakdet and PrEgg.  #
# I began developing this script after completing a fieldwork trip in Northern #
# Vietnam. The script extracts several glottal parameters from an EGG signal.  #
# Durations are calculated in seconds. Multiply by 1000 manually if millisecond#
# units are needed.                                                            #
#                                                                              #
# Methodology: Rothenberg method (Rothenberg & Mahshie, 1988), dEGG (Henrich   #
# et al., 2004), First central difference (Herbst et al., 2014) and hybrid     #
# approaches proposed by Howard (1995) and Henry Tehrani (2009).               #
# Processing speed: 45-50 (upto 80) files per minute.                          #
# Users can manually adjust methods and thresholds via the Settings window.    #
#                                                                              #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                                                              #
#                                   CITATION                                   #
#                                                                              #
# Ta Van Dat. 2025. EGGLab: Processing Electroglottography Signal with Praat.  #
#	Retrieved on (date) from https://tadat-ling.github.io/resource-EGGLab.html #
#                                                                              #
#                                 BIBLIOGRAPHY                                 #
#                                                                              #
# [1] C. Herbst, J. Lohscheller, J. Švec, N. Henrich, G. Weissengruber & W.    #
#		Tecumseh Fitch (2014): “Glottal opening and closing events investigated#
#		by electroglottography and super-high-speed video recordings”, The     #
#		Journal of Experimental Biology 217: 955–963.                          #
# [2] Chan, May Pik Yu & Jianjing Kuang. 2024. PrEgg: A free and open source   #
#		Praat script for electroglottography measurements, JASA, 155, A336.    #
#		DOI: doi.org/10.1121/10.0027728.                                       #
# [3] Chan, May Pik Yu & Jianjing Kuang. 2024. PrEgg:A free and open source    #
#		Praat script for electroglottography measurements (v0.1). Retrieved    #
#		on 10 March, 2025 from https://github.com/maypychan/praat-egg.         #
# [4] Christian T. Herbst. Electroglottography – An Update, Journal of Voice,  #
#		34, (4), 503-526.                                                      #
# [5] David M. Howard, Geoffrey A. Lindsey, Bridget Allen. 1990. Toward the    #
# 		quantification  of vocal efficiency, Journal of Voice, 4, (3), 205-212.#
# [6] Howard, D. M. (1995). Variation of electrolaryngographically derived     #
#		closed quotient for trained and untrained adult female singers,        #
#		J. Voice 9, 163–72.                                                    #
# [7] Kirby, James. 2020. Praatdet: Praat-based tools for EGG analysis (v0.3). #
#		Retrieved on 20 March, 2025 from https://github.com/kirbyj/praatdet.   #
# [8] Michaud, A. 2004. Final consonants and glottalization: new perspectives  #
# 		from Hanoi Vietnamese, Phonetica 61, 119–146.                          #
# [9] Nathalie Henrich et al. 2004. On the use of the derivative of electro-   #
#		glottographic signals for characterization of nonpathological          #
#		phonation, JASA, 115 (3): 1321–1332. https://doi.org/10.1121/1.1646401 #
# [10] Rothenberg, M., & Mahshie, J. J. (1988). Monitoring vocal fold abduction#
#       through vocal fold contact area. Journal of Speech, Language, and Hear-#
#       ing Research, 31(3), 338-351.                                          #
# [11] Tehrani., H. EggWorks. from                                             #
#		http://phonetics.linguistics.ucla.edu/facilities/physiology/EGG.htm.   #
#                                                                              #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                                                                              #
#                              OUTPUT PARAMETERS                               #
#                                                                              #
# method: Method of analysis (dEGG, hybrid, threshold, Henry Tehrani)          #
# filename: Name of the analyzed wav file                                      #
# speaker: Speaker code                                                        #
# tier: Reference tier                                                         #
# label: Reference label                                                       #
# duration: Duration of the analyzed label                                     #
# period: Number of period within the analyzed label                           #
# period_dur: Duration of the period                                           #
# cp: Duration of closed phase                                                 #
# op: Duration of open phase                                                   #
# cq: Contact Quotient                                                         #
# oq: Open Quotient                                                            #
# cp_op_ratio: Ratio of CP and OP                                              #
# p_pascal: Amplitude of the EGG peak within the period                        #
# v_pascal: Amplitude of the EGG valley within the period                      #
# pic: Peak of increased contact                                               #
# pic_time: Time of PIC (measured from the beginning of sound file)            #
# pdc: Peak of decreased contact                                               #
# pdc_time: Time of PDC (measured from the beginning of sound file)            #
# f0_egg: F0 calculated by 1/cycle's duration                                  #
# f0_auto: F0 tracked automatically by Praat's algorithm                       #
# pv_slope: Slope of the EGG signal's peak-to-valley                           #
# closing: duration of closing phase                                           #
#          smaller closing & bigger sq -> more leftward-skewed EGG peak        #
# opening: duration of opening phase                                           #
# sq: Speed Quotient                                                           #
# op_slope: Slope of opening phase                                             #
# sq_m: Skew Quotient (Speed Quotient) - The ratio between closing duration and#
# opening duration measured at a 10% threshold (Marasek, 1997)                 #
# closing_marasek: Closing duration measured at a 10% threshold                #
# opening_marasek: Opening duration measured at a 10% threshold                #
# contact_dur: Glottis contacted phase at 10% threshold                        #
# decontact_dur: Glottis decontacted phase at 10% threshold                    #
# con_dec_r: Ratio of contact_dur and decontact_dur (contact_dur/decontact_dur)#
#                                                                              #
################################################################################

form EGG Lab v1.0
	#comment ***********************  Welcome to EGG Lab ver 1.0!  ***********************
  	comment Pitch and Filter Range:
 	positive pitch_floor 40
	natural pitch_ceiling 450
	natural low_pass 5000
	natural high_pass 20

	comment Advanced Pitch Settings:
	optionmenu pitch_method 1
		option Filter autocorrelation
		option Raw cross-correlation
		option Raw autocorrelation
		option Filter cross-correlation
	real time_step 0
	natural max_candidates 15
	boolean very_accurate 0
  	positive silence_threshold 0.03
  	positive voicing_threshold 0.01
  	positive octave_cost 0.01
  	positive octave_jump_cost 0.15
  	positive voiced_cost 0.14

	comment EGG Settings:
	optionmenu egg_channel 2
		option Channel 1
		option Channel 2
	natural	smoothing 100
	positive closing_threshold 0.30
	#positive peak_threshold 0.05
	positive hybrid_threshold 0.25
	optionmenu derivative_method 1
		option Derivative
		option First central difference
	optionmenu smooth_method 4
        option No smoothing
        option Narrow (5-point)
        option Medium (7-point)
        option Wide (21-point)
        option Loess-like smoothing
	comment Interpolation of Amplitude:
	optionmenu interpolation 4
		option nearest
		option linear
		option cubic
		option sinc70
		option sinc700
endform

at_least_one_selected = 0
while !at_least_one_selected
	beginPause: "Directory and Analysis method"
		folder: "Directory", "/Browse/your/WAV/and/TextGrid/input/folder/(please save wav and TextGrid files in the same folder)"
		comment: "Enter speaker code (The data will be saved as 'speaker-egg.csv'.):"
		sentence: "Speaker", "M23"
	    comment: "Enter reference tier for analysing (specify of enter [0] if there is no ref tier.):"
		real: "Reference tier", 2
		comment: "Specify or enter [all] to analyze all of the labeled intervals within the ref tier."
		sentence: "Interval label", "all"
		real: "Origin label tier", 3
		
	    comment: "Please select at least one method to process your EGG signals:"
		boolean: "Threshold method", 0
		boolean: "dEGG method", 0
		boolean: "Hybrid method", 1
		boolean: "Henry Tehrani method", 0
		positive: "dEGG crossing threshold", 0.25
		positive: "dEGG crossing time step",  0.002
		real: "Slope step", 15
		positive: "Skew quotient threshold", 0.10
		boolean: "Analyse spectrum", 0
		
		comment: "Other settings:"
		#comment: "Would you like to enable a notification alert when processing finishes?"
		boolean: "Completion alert", 0
		#comment: "Would you like to manually check the closing-opening intervals?"
		boolean: "Manual check", 0
		#comment: "Would you like to remove aperiodic intervals? (On by default, turn off for creaky voice analysis)"
		boolean: "Remove aperiodic intervals", 1
		#comment: "Automatically remove temporary objects from the Object list after processing?"
		boolean: "Remove temporary objects", 1
		boolean: "Save EGG TextGrids", 1
	endPause: "Continue", 1
    if threshold_method or dEGG_method or hybrid_method or henry_Tehrani_method
        at_least_one_selected = 1
    else
        beginPause: "Error - No Method Selected"
            comment: "Please select at least one method!"
            comment: "Click 'Continue' to reselect."
        endPause: "Continue", 1
    endif
endwhile


# Get processing duration
stopwatch

# Praat version checking
if (praatVersion < 6401)
	printline Requires Praat version 6.4 or higher. Please upgrade your Praat version 
	exit
endif

if very_accurate = 0
	very_accurate_status$ = "off"
else
	very_accurate_status$ = "on"
endif


# Directory format
if !endsWith(directory$, "/") && !endsWith(directory$, "\")
    directory$ = directory$ + if macintosh or unix then "/" else "\" fi
endif

textGrid_directory$ = directory$
if !endsWith(textGrid_directory$, "/") && !endsWith(textGrid_directory$, "\")
    textGrid_directory$ = textGrid_directory$ + if macintosh or unix then "/" else "\" fi
endif

# Create output folder
output_directory$ = directory$ + "output/"
createDirectory: output_directory$

if save_EGG_TextGrids == 1
	output_wav_textgrid$ = directory$ + "output_wav_textgrid/"
	createDirectory: output_wav_textgrid$
endif


if pitch_method = 1
	f0_method$ = "ac"
elif pitch_method = 2
	f0_method$ = "raw_cc"
elif pitch_method = 3
	f0_method$ = "raw_ac"
else
	f0_method$ = "cc"
endif

if analyse_spectrum == 1
	output_file$ = output_directory$ + speaker$ + "-egg.csv"
	writeFile: output_file$
	appendFile: output_file$, "method", ",", "filename", ",", "origin_label", ",", "speaker", ",", "tier", ",", "label", ",", "duration", ",", "period", ",", "period_dur", ",", "cp", ",", "op", ",", "cq", ",", "oq", ",", "cp_op_ratio", ",", "p_pascal", ",", "v_pascal", ",", "pic", ",", "pic_time", ",", "pdc", ",", "pdc_time", ",", "f0_egg", ",", "f0_auto", ",", "pv_slope", ",", "closing", ",", "opening", ",", "sq", ",", "op_slope", ",", "sq_m", ",", "closing_marasek", ",", "opening_marasek", ",", "contact_dur", ",", "decontact_dur", ",", "con_dec_r", ",", "h1h2c", newline$
else
	output_file$ = output_directory$ + speaker$ + "-egg.csv"
	writeFile: output_file$
	appendFile: output_file$, "method", ",", "filename", ",", "origin_label", ",", "speaker", ",", "tier", ",", "label", ",", "duration", ",", "period", ",", "period_dur", ",", "cp", ",", "op", ",", "cq", ",", "oq", ",", "cp_op_ratio", ",", "p_pascal", ",", "v_pascal", ",", "pic", ",", "pic_time", ",", "pdc", ",", "pdc_time", ",", "f0_egg", ",", "f0_auto", ",", "pv_slope", ",", "closing", ",", "opening", ",", "sq", ",", "op_slope", ",", "sq_m", ",", "closing_marasek", ",", "opening_marasek", ",", "contact_dur", ",", "decontact_dur", ",", "con_dec_r", newline$
endif


filestrings = Create Strings as file list: "list", directory$ + "/*.wav"
filenumbers = Get number of strings
clearinfo
for numWAV from 1 to filenumbers
	selectObject: "Strings list"
	current_token$ = Get string: numWAV
	soundName$ = directory$ + current_token$
	wavfile_raw = Read from file: soundName$
    filename$ = selected$ ("Sound", 1)
	#appendInfoLine: newline$, numWAV, "/", filenumbers, ": ", filename$ + ".wav"
	appendInfo: newline$, "Now processing ", numWAV, "/", filenumbers, ": ", filename$ + ".wav"

	#pitch object
	selectObject: wavfile_raw
	if pitch_method = 1
		pitch = To Pitch (ac): time_step, pitch_floor, max_candidates, very_accurate_status$, silence_threshold, voicing_threshold, octave_cost, octave_jump_cost, voiced_cost, pitch_ceiling
	elif pitch_method = 2
		pitch = To Pitch (raw cc): time_step, pitch_floor, pitch_ceiling, max_candidates, very_accurate_status$, silence_threshold, voicing_threshold, octave_cost, octave_jump_cost, voiced_cost
	elif pitch_method = 3
		pitch = To Pitch (raw ac): time_step, pitch_floor, pitch_ceiling, max_candidates, very_accurate_status$, silence_threshold, voicing_threshold, octave_cost, octave_jump_cost, voiced_cost
	else
		pitch = To Pitch (cc): time_step, pitch_floor, max_candidates, very_accurate_status$, silence_threshold, silence_threshold, octave_cost, octave_jump_cost, voiced_cost, pitch_ceiling
	endif

	# If textgrid exists, open it
	if reference_tier > 0
		textgrid_raw = Read from file: textGrid_directory$ + filename$ + ".TextGrid"
	endif

	# 开始分析
	selectObject: wavfile_raw
	egg = Extract Electroglottogram: egg_channel, "no"
	eggfilter = High-pass filter: high_pass, smoothing
	selectObject: egg
	egg_raw = To Sound
	selectObject: egg_raw
	dc_component = Filter (pass Hann band): 0, 20, smoothing

	selectObject: eggfilter
	#Derivative Signal
	if derivative_method = 1
		deriv = Derivative: low_pass, smoothing, 0
	else
		deriv = First central difference: 0
	endif

	# smooth
	selectObject: deriv
	procedure validateSound: .minSamples
	    .nFrames = Get number of samples
	    if .nFrames < .minSamples
	        pause: "Sound too short (needs ≥" + string$(.minSamples) + " samples)."
	        exit
	    endif
	endproc
	duration = Get total duration
	nFrames = Get number of samples
	if smooth_method = 1
	    deriv_smoothed = Copy: "original_copy_" + filename$
	    kernel$ = "No smoothing"
	    smooth_notice$ = "Original sound preserved (no smoothing applied)"
	else
	    deriv_smoothed = Copy: "smoothed_" + filename$
	    selectObject: deriv_smoothed
	
	    if smooth_method = 2
	        @validateSound: 5
	        Formula: "if col < 3 || col > ncol-2 then self else (1*(self[col-2]+self[col+2]) + 2*(self[col-1]+self[col+1]) + 3*self[col])/9 fi"
	        kernel$ = "5-point narrow [9]"
	        smooth_notice$ = "Applied 5-point narrow smoothing kernel"
	
	    elif smooth_method = 3
	        @validateSound: 7
	        Formula: "if col < 4 || col > ncol-3 then self else (1*(self[col-3]+self[col+3]) + 2*(self[col-2]+self[col+2]) + 3*(self[col-1]+self[col+1]) + 3*self[col])/15 fi"
	        kernel$ = "7-point medium [15]"
	        smooth_notice$ = "Applied 7-point medium smoothing kernel"
	
	    elif smooth_method = 4
	        @validateSound: 21
	        Formula: "if col < 11 || col > ncol-10 then self else (1*(self[col-10]+self[col+10]) + 2*(self[col-9]+self[col+9]) + 3*(self[col-8]+self[col+8]) + 4*(self[col-7]+self[col+7]) + 5*(self[col-6]+self[col+6]) + 6*(self[col-5]+self[col+5]) + 7*(self[col-4]+self[col+4]) + 8*(self[col-3]+self[col+3]) + 9*(self[col-2]+self[col+2]) + 10*(self[col-1]+self[col+1]) + 11*self[col])/121 fi"
	        kernel$ = "21-point wide [121]"
	        smooth_notice$ = "Applied 21-point wide smoothing kernel"
	
	    elif smooth_method = 5
	        @validateSound: 21
	        Formula: "if col < 11 || col > ncol-10 then self else (self[col-10]*0.005 + self[col-9]*0.01 + self[col-8]*0.02 + self[col-7]*0.03 + self[col-6]*0.05 + self[col-5]*0.07 + self[col-4]*0.1 + self[col-3]*0.12 + self[col-2]*0.15 + self[col-1]*0.17 + self[col]*0.18 + self[col+1]*0.17 + self[col+2]*0.15 + self[col+3]*0.12 + self[col+4]*0.1 + self[col+5]*0.07 + self[col+6]*0.05 + self[col+7]*0.03 + self[col+8]*0.02 + self[col+9]*0.01 + self[col+10]*0.005) fi"
	        kernel$ = "Loess-like 21-point"
	        smooth_notice$ = "Applied loess-like super smooth (21-point weighted average)"
	    endif
	endif

	#Save as WAV file: directory$ + filename$ + "_Deriv.wav"

	selectObject: deriv_smoothed
	# Create point tier
	tg_degg500 = To TextGrid: "lv500", "lv500"
	selectObject: deriv_smoothed
	tg_pic = To TextGrid: "pic", "pic"
	selectObject: deriv_smoothed
	tg_pdc = To TextGrid: "pdc", "pdc"

	# 选择平滑后的音频
	# Để tránh tính giá trị sai khi cắt nhầm: tính giá trị cao nhất trong khoảng 0.0005 (0.5ms đầu)
	# đến dur_total - 0.0005 (0.5ms trước điểm kết thúc)
	selectObject: deriv_smoothed
	start_time_degg = Get start time
	end_time_degg = Get end time
	degg_peak_detect_start = start_time_degg + 0.001
	degg_peak_detect_end = end_time_degg - 0.001
	degg_max_peak = Get maximum: degg_peak_detect_start, degg_peak_detect_end, interpolation$
	### degg_max_peak = Get maximum: 0, 0, interpolation$
	cross_lvl = degg_max_peak*dEGG_crossing_threshold
	# Get the first level crossing
	t_current = Get nearest level crossing: 1, 0.001, cross_lvl, "right"
	selectObject: deriv_smoothed
	total_duration = Get total duration
	# 选择EGG标注层并插入第一个点
	selectObject: tg_degg500
	Insert point: 1, t_current, "l"
	# 初始化t_next
	t_next = t_current + dEGG_crossing_time_step
	due_duration = total_duration - dEGG_crossing_time_step
	while t_current < due_duration and t_next != undefined
	    search_start = t_current + dEGG_crossing_time_step
	    selectObject: deriv_smoothed
	    t_next = Get nearest level crossing: 1, search_start, cross_lvl, "right"
	    # 检查是否找到有效点
	    if t_next != undefined and t_next > t_current
	        # 插入点
	        selectObject: tg_degg500
	        Insert point: 1, t_next, "l"
	        # 更新当前时间&准备下一次搜索
	        t_current = t_next
	        t_next = t_current + dEGG_crossing_time_step
	    endif
	endwhile

	# 获取所有点数
	selectObject: tg_degg500
	n_points_cross = Get number of points: 1
	# 从后往前循环，避免索引错乱
	i = n_points_cross - 1
	while i >= 1
	    t1 = Get time of point: 1, i
	    t2 = Get time of point: 1, i + 1
	    interval = t2 - t1 - 0.001
	    frequency = 1 / interval
	    if frequency > pitch_ceiling
	        # 查找最大值时间
	        selectObject: deriv_smoothed
	        time_max = Get time of maximum: t1, t2, interpolation$
	        # 获取新的穿越点
	        t_cross_new = Get nearest level crossing: 1, time_max, cross_lvl, "left"
	        # 删除原始两个点
	        selectObject: tg_degg500
	        Remove point: 1, i + 1
	        Remove point: 1, i
	        # 插入新的点
			if t_cross_new != undefined
				selectObject: tg_degg500
				Insert point: 1, t_cross_new, "l"
			endif
	        # 点数变化，更新总点数
	        n_points_cross = Get number of points: 1
	        # 因为删除了两个，插入了一个，往前退一个点继续判断
	        i = i - 1
	    else
	        i = i - 1
	    endif
	endwhile

	# 获取 tg_degg500 中所有 "l" 点的数量
	selectObject: tg_degg500
	num_l_points = Get number of points: 1
	# 确保至少有两个 "l" 点才能计算区间
	if num_l_points < 2
	    exitScript: "需要至少两个'l'点才能进行计算"
	endif
	# 遍历所有相邻的 "l" 点对
	for i from 1 to num_l_points - 1
	    # 获取当前和下一个 "l" 点的时间
		selectObject: tg_degg500
	    t1 = Get time of point: 1, i
	    t2 = Get time of point: 1, i+1
	    # 在平滑后的音频中查找这个时间范围内的极值
	    selectObject: deriv_smoothed
	    d1 = Get time of maximum: t1, t2, interpolation$
	    d2 = Get time of minimum: d1, t2, interpolation$
		if d1 < d2
			# 在标注层中插入极值点
			selectObject: tg_pic
			Insert point: 1, d1, "d1"
			selectObject: tg_pdc
			Insert point: 1, d2, "d2"
		endif
	endfor


	###### Filtered and Smoothes EGG signal ######
	selectObject: eggfilter
	wavfile = To Sound
	# smooth
	selectObject: wavfile
	procedure validateSound: .minSamples
	    .nFrames = Get number of samples
	    if .nFrames < .minSamples
	        pause: "Sound too short (needs ≥" + string$(.minSamples) + " samples)."
	        exit
	    endif
	endproc
	duration = Get total duration
	nFrames = Get number of samples
	if smooth_method = 1
	    wavfile_smoothed = Copy: "original_copy_" + filename$
	    kernel$ = "No smoothing"
	    smooth_notice$ = "Original sound preserved (no smoothing applied)"
	else
	    wavfile_smoothed = Copy: "smoothed_" + filename$
	    selectObject: wavfile_smoothed
	    if smooth_method = 2
	        @validateSound: 5
	        Formula: "if col < 3 || col > ncol-2 then self else (1*(self[col-2]+self[col+2]) + 2*(self[col-1]+self[col+1]) + 3*self[col])/9 fi"
	        kernel$ = "5-point narrow [9]"
	        smooth_notice$ = "Applied 5-point narrow smoothing kernel"
	
	    elif smooth_method = 3
	        @validateSound: 7
	        Formula: "if col < 4 || col > ncol-3 then self else (1*(self[col-3]+self[col+3]) + 2*(self[col-2]+self[col+2]) + 3*(self[col-1]+self[col+1]) + 3*self[col])/15 fi"
	        kernel$ = "7-point medium [15]"
	        smooth_notice$ = "Applied 7-point medium smoothing kernel"
	
	    elif smooth_method = 4
	        @validateSound: 21
	        Formula: "if col < 11 || col > ncol-10 then self else (1*(self[col-10]+self[col+10]) + 2*(self[col-9]+self[col+9]) + 3*(self[col-8]+self[col+8]) + 4*(self[col-7]+self[col+7]) + 5*(self[col-6]+self[col+6]) + 6*(self[col-5]+self[col+5]) + 7*(self[col-4]+self[col+4]) + 8*(self[col-3]+self[col+3]) + 9*(self[col-2]+self[col+2]) + 10*(self[col-1]+self[col+1]) + 11*self[col])/121 fi"
	        kernel$ = "21-point wide [121]"
	        smooth_notice$ = "Applied 21-point wide smoothing kernel"
	
	    elif smooth_method = 5
	        @validateSound: 21
	        Formula: "if col < 11 || col > ncol-10 then self else (self[col-10]*0.005 + self[col-9]*0.01 + self[col-8]*0.02 + self[col-7]*0.03 + self[col-6]*0.05 + self[col-5]*0.07 + self[col-4]*0.1 + self[col-3]*0.12 + self[col-2]*0.15 + self[col-1]*0.17 + self[col]*0.18 + self[col+1]*0.17 + self[col+2]*0.15 + self[col+3]*0.12 + self[col+4]*0.1 + self[col+5]*0.07 + self[col+6]*0.05 + self[col+7]*0.03 + self[col+8]*0.02 + self[col+9]*0.01 + self[col+10]*0.005) fi"
	        kernel$ = "Loess-like 21-point"
	        smooth_notice$ = "Applied loess-like super smooth (21-point weighted average)"
	    endif
	endif


######### Preparing for TextGrid checking
	# Resample
	if egg_channel = 1
		selectObject: wavfile_raw
		audio = Extract one channel: 2
	elif egg_channel = 2
		selectObject: wavfile_raw
		audio = Extract one channel: 1
	endif
	selectObject: audio
	audio_resampled = Resample: 22050, 50
	selectObject: wavfile_smoothed
	wavfile_resampled = Resample: 22050, 50
	# Resample deriv
	selectObject: deriv_smoothed
	deriv_resampled = Resample: 22050, 50
	# Merge them together!
	selectObject: audio_resampled, wavfile_resampled, deriv_resampled
	stereo_combined = Combine to stereo
	if save_EGG_TextGrids == 1
		selectObject: stereo_combined
		nowarn Save as binary file: output_wav_textgrid$ + speaker$ + "_" + filename$ + "_EGG" + ".Sound"
	endif

### Merge and compelet TextGrids manufacturing
	if manual_check == 1
		cq_min = 1
		cq_max = 0
		f0_min = pitch_ceiling
		f0_max = pitch_floor
		h1h2c_min = 1000
		h1h2c_max = -1000
		index_max = 0
	
		selectObject: tg_pic
		n_point_pic = Get number of points: 1
		n_row = n_point_pic * 2
		# Initialize matrix with n_clmn rows and 7 columns (index, cq_d, cq_h, cq_t, cq_ht, f0_degg, f0_thres, h1h2c, f0_auto)
		analysis_matrix = Create simple Matrix: "analysis_matrix", n_row, 9, "row*col"
	endif

	### dEGG Method
	if dEGG_method == 1
		selectObject: tg_pic
		plusObject: tg_pdc
		tg_degg = Merge
		# Rename tiers
		selectObject: tg_degg
		Set tier name: 1, "dEGGpeak"
		Set tier name: 2, "dEGGvalley"
		# Add tier 3 (dEGGCI) with intervals bounded by tier1 and tier2 points
		selectObject: tg_degg
		Insert interval tier: 3, "dEGGCI"
	
		# Get all points from tier 1 (dEGGpeak) and tier 2 (dEGGvalley)
		selectObject: tg_degg
		numPointsTier1 = Get number of points: 1
		numPointsTier2 = Get number of points: 2
		# Create lists of all points
		d1_times# = zero#(numPointsTier1)
		d2_times# = zero#(numPointsTier2)
		for i from 1 to numPointsTier1
		    d1_times#[i] = Get time of point: 1, i
		endfor
		for i from 1 to numPointsTier2
		    d2_times#[i] = Get time of point: 2, i
		endfor
		# Check if first d2 comes before first d1
		if numPointsTier2 > 0 and numPointsTier1 > 0
		    if d2_times#[1] < d1_times#[1]
		        # Create new array without the first d2 point
		        temp_d2_times# = zero#(numPointsTier2 - 1)
		        for i from 2 to numPointsTier2
		            temp_d2_times#[i-1] = d2_times#[i]
		        endfor
		        d2_times# = temp_d2_times#
		        numPointsTier2 = numPointsTier2 - 1
		    endif
		endif
		# Manually combine the two arrays
		boundaryTimes# = zero#(numPointsTier1 + numPointsTier2)
		for i from 1 to numPointsTier1
		    boundaryTimes#[i] = d1_times#[i]
		endfor
		for i from 1 to numPointsTier2
		    boundaryTimes#[numPointsTier1 + i] = d2_times#[i]
		endfor
		# Sort the boundary times
		boundaryTimes# = sort#(boundaryTimes#)
		# Add boundaries to tier 3
		for i from 1 to size(boundaryTimes#)
		    time = boundaryTimes#[i]
		    Insert boundary: 3, time
		endfor
	
		# Tiers now:
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. dEGGCI (closed intvl)
	
		### Add tiers dEGGTime
		selectObject: tg_degg
		Insert point tier: 4, "dEGGTime"
		# Get all points from tier 1 (dEGGpeak), tier 2 (dEGGvalley)
		selectObject: tg_degg
		numPointsTier1 = Get number of points: 1
		numPointsTier2 = Get number of points: 2
		# Create lists of all points
		d1_times# = zero#(numPointsTier1)
		d2_times# = zero#(numPointsTier2)
		for i from 1 to numPointsTier1
		    d1_times#[i] = Get time of point: 1, i
		endfor
		for i from 1 to numPointsTier2
		    d2_times#[i] = Get time of point: 2, i
		endfor
		# Tiers now:
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. dEGGCI (closed intvl)
		# 4. dEGGTime (points)
		### Adjust TIME tiers (dEGGTime-Tier4)
		### Taking numbers from the TextGrid (tier 3) and adjust dEGG (tier 4)	
		selectObject: tg_degg
		nInt = do ("Get number of intervals...", 3)
	    # Assuming that the second interval annotates the first glottal period's close phase
		for i from 2 to nInt
			time_dEGGCI = Get start time of interval: 3, i
			endtime_dEGGCI = Get end time of interval: 3, i
			# Checking whether the point is dEGG peak or valley
			selectObject: deriv_smoothed
		    check1 = Get value at time: 1, time_dEGGCI, interpolation$
		    check2 = Get value at time: 1, endtime_dEGGCI, interpolation$
			selectObject: tg_degg
			if check1 > check2
				Insert point: 4, time_dEGGCI, "t1"
			elif check1 < check2
				Insert point: 4, time_dEGGCI, "t2"
			else
				Insert point: 4, time_dEGGCI, "error"
			endif
		endfor
	
		### Label closing intervals and opening intervals of dEGGCI (3), time tier: 4
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. dEGGCI (closed intvl)
		# 4. dEGGTime (points)
	
		selectObject: tg_degg
		nInt_dEGG = do ("Get number of intervals...", 3)
	
		for i from 2 to nInt_dEGG
			starttime = Get start time of interval: 3, i
			endtime = Get end time of interval: 3, i
	
			ptlabel_start$ = Get label of point: 4, i-1
	
			if i != nInt_dEGG
				ptlabel_end$ = Get label of point: 4, i
			endif
	
			# if the time of the peak is between t1 and t2, then it's a closed interval.
			if ptlabel_start$ == "t1" and ptlabel_end$ == "t2"
				Set interval text: 3, i, "cp"
			# if the start of the interval corresponds to a t2 label, and the end of the interval has a t1 label,
			# it's probably an open interval
			elif ptlabel_start$ == "t2" and ptlabel_end$ == "t1"
				Set interval text: 3, i, "op"
			endif
		endfor
	
		selectObject: tg_degg
		Insert point tier: 5, "check"
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. dEGGCI (closed intvl)
		# 4. dEGGTime (points)
		# 5. dEGGcheck
		Remove tier: 1
		Remove tier: 1
		# Tiers now: 1. dEGGCI, 2. dEGGTime, 3. dEGGCheck
		Remove tier: 2
		# Tiers now: 1. dEGGCI, 2. dEGGCheck

		# Pitch checking
		selectObject: tg_degg
		nInt_tocheck = do ("Get number of intervals...", 1)
		for i from 1 to nInt_tocheck
			selectObject: tg_degg
			intval$ = Get label of interval: 1, i
		    nextval$ = "NA"
		    if i < nInt_tocheck
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
				selectObject: tg_degg
				start_time_pitch_check = Get start time of interval: 1, i
				end_time_pitch_check = Get end time of interval: 1, i+1
				duration = end_time_pitch_check - start_time_pitch_check
        			freq = 1 / duration
				if remove_aperiodic_intervals == 1
					if freq > pitch_ceiling || freq < pitch_floor
						selectObject: tg_degg
						Set interval text: 1, i, ""
						Set interval text: 1, i+1, ""
					endif
				endif
			endif
		endfor

		############ Manual Check ############
		if manual_check = 1
			selectObject: stereo_combined
			plusObject: tg_degg
			Edit
			pauseScript: "Edit the TextGrid as needed, then click Continue to proceed."
			#appendInfo: " Checked"
		endif

############ CALCULATION ############
		method$ = "dEGG"
		selectObject: tg_degg	
		# dEGG OQ, CQ, F0 measurements
		nInt_t1 = do ("Get number of intervals...", 1)
		line_count = 1
		for i from 1 to nInt_t1
			selectObject: tg_degg
			intval$ = Get label of interval: 1, i
			#nextval$ = Get label of interval: 1, i+1
		    # Safely get next interval label with boundary check
		    nextval$ = "NA"
		    if i < nInt_t1
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
		        # Get timing information
		        cp_starttime = Get start time of interval: 1, i
		        cp_endtime = Get end time of interval: 1, i
		        op_endtime = Get end time of interval: 1, i+1
		        # Calculate durations
		        cp_dur = cp_endtime - cp_starttime
		        op_dur = op_endtime - cp_endtime
		        cycle_dur = cp_dur + op_dur
		        # Calculate parameters with safety checks
		        cq = if cycle_dur > 0 then cp_dur / cycle_dur else undefined fi
		        oq = if cycle_dur > 0 then op_dur / cycle_dur else undefined fi
		        cqoq = if cycle_dur > 0 then cp_dur / op_dur else undefined fi
		        f0_egg = if cycle_dur > 0 then 1 / cycle_dur else undefined fi
		        midpoint1 = if cp_dur > 0 then cp_dur/2 + cp_starttime else undefined fi
		        midpoint2 = if op_dur > 0 then op_dur/2 + cp_endtime else undefined fi
		        
				selectObject: pitch
				f0_auto = Get value at time: cp_starttime, "Hertz", "linear"
				if f0_auto == undefined
		            f0_auto = f0_egg
		        endif

				selectObject: wavfile_smoothed
				p_pascal = Get maximum: cp_starttime, cp_endtime, interpolation$
				peak_time = Get time of maximum: cp_starttime, cp_endtime, interpolation$
				v_pascal = Get minimum: peak_time, op_endtime, interpolation$
				valley_time = Get time of minimum: peak_time, op_endtime, interpolation$
				# or: cp_starttime

				# SQ, contacting duration, decontacting duration, decontacting slope
				c_starttime = cp_starttime
		        c_endtime = peak_time
		        c_dur = c_endtime - c_starttime
		        o_starttime = peak_time
		        o_endtime = cp_endtime
		        o_dur = o_endtime - o_starttime
		        sq = o_dur / c_dur

				# 10% threshold (Marasek, 1997) Speed Quotient: sq_m
				thres_for_cend_ostart = p_pascal - (p_pascal - v_pascal) * skew_quotient_threshold
				c_end = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "left"
				o_start = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "right"
				thres_for_oend = v_pascal + (p_pascal - v_pascal) * skew_quotient_threshold
				o_end = Get nearest level crossing: 1, valley_time, thres_for_oend, "left"
				c_start_valley = Get nearest level crossing: 1, valley_time, thres_for_oend, "right"
				# Get c_start
				if i > 1
					selectObject: tg_degg
					prev_op_starttime = Get start time of interval: 1, i-1
				else
					selectObject: tg_degg
					prev_op_starttime = cp_starttime - 0.002
				endif
				curr_cp_starttime = cp_starttime
				selectObject: wavfile_smoothed
				prev_v_pascal = Get minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				prev_valley_time = Get time of minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				thres_for_cstart = prev_v_pascal + (p_pascal - prev_v_pascal) * skew_quotient_threshold
				c_start = Get nearest level crossing: 1, prev_valley_time, thres_for_cstart, "right"
				c_dur_marasek = c_end - c_start
				o_dur_marasek = o_end - o_start
				sq_m = c_dur_marasek / o_dur_marasek
				#sq_m = o_dur_marasek / c_dur_marasek
				# contact vs decontact ratio (10% threshold)
				contact_dur = o_start - c_end
				decontact_dur = c_start_valley - o_end
				con_dec_r = contact_dur / decontact_dur
		
		        # For checking purposes, a point on tier four refers to where the sq was taken
				selectObject: tg_degg
		        Insert point: 2, o_endtime, fixed$(sq,3)
				
				## Calculate Decontact Skewness: op_slope
				sound_part_dur = o_dur
				sum_diff = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = o_starttime + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = o_starttime + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_diff += sk2 - sk1
				endfor
				dec_skew = sum_diff / (slope_step - 1)

				## Calculate EGG Skewness: pv_slope
				selectObject: wavfile_smoothed
				sound_part_dur = valley_time - peak_time
				sum_differences = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = peak_time + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = peak_time + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_differences += sk2 - sk1
				endfor
				skewness = sum_differences / (slope_step - 1)
	
				#selectObject: deriv_smoothed
				#pic = Get value at time: 0, cp_starttime, interpolation$
				##pic_time = cp_starttime - c_start
				#pic_time = cp_starttime
				#pdc = Get value at time: 0, cp_endtime, interpolation$
				##pdc_time = cp_endtime - c_start
				#pdc_time = cp_endtime

				# Get value from raw degg signal
				selectObject: deriv
				pic_starttime = cp_starttime - 0.0005
				pdc_endtime = op_endtime - 0.0005
				pic = Get maximum: pic_starttime, cp_endtime, interpolation$
				pic_time = Get time of maximum: pic_starttime, cp_endtime, interpolation$
				pdc = Get minimum: pic_time, pdc_endtime, interpolation$
				pdc_time = Get time of minimum: pic_time, pdc_endtime, interpolation$

				# Spectral
				if analyse_spectrum == 1
					resample = 11025
					selectObject: audio
					audio_11kHz = Resample: resample, 50
					selectObject: audio_11kHz
					lpc_audio_11kHz = To LPC (burg): 20, 0.025, 0.005, 34.45
					selectObject: audio_11kHz, lpc_audio_11kHz
					wavfile_invrsfltr = Filter (inverse)
					selectObject: wavfile_invrsfltr
					specgram = To Spectrogram: 0.03, resample/2, 0.01, 1, "Gaussian"
					h1h2c_sum = 0
					n_step = 10
					for x from 1 to n_step
						step1 = x * cycle_dur/n_step
						selectObject: specgram
						spectrum_slice1 = To Spectrum (slice): cp_starttime+step1
						selectObject: spectrum_slice1
						f0_spectrum_1 = Get frequency of nearest maximum: f0_egg
						h1c_1 = Get sound pressure level of nearest maximum: f0_spectrum_1
						h2c_freq_1 = f0_spectrum_1 * 2
						h2c_1 = Get sound pressure level of nearest maximum: h2c_freq_1
						h1h2c_1 = h1c_1 - h2c_1
						h1h2c_sum += h1h2c_1
						removeObject: spectrum_slice1
					endfor
					h1h2c = h1h2c_sum / n_step
					removeObject: specgram
					removeObject: audio_11kHz, lpc_audio_11kHz, wavfile_invrsfltr
				endif
		
		        # Insert points with error handling
				selectObject: tg_degg
		        if cq != undefined
		            nocheck Insert point: 2, cp_starttime, fixed$(cq, 3)
		        else
		            nocheck Insert point: 2, cp_starttime, "NA"
		        endif
		        selectObject: tg_degg
		        if midpoint1 != undefined and f0_auto != undefined
		            nocheck Insert point: 2, midpoint1, fixed$(f0_auto, 0)
		        else
		            nocheck Insert point: 2, cp_starttime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif
				selectObject: tg_degg
		        if midpoint2 != undefined and f0_egg != undefined
		            nocheck Insert point: 2, midpoint2, fixed$(f0_egg, 0)
		        else
		            nocheck Insert point: 2, cp_endtime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif

				if origin_label_tier > 0
					selectObject: textgrid_raw
					current_intvl = Get interval at time: origin_label_tier, cp_starttime
					if current_intvl > 0
				        origin_label$ = Get label of interval: origin_label_tier, current_intvl
				        if origin_label$ = ""
				            origin_label$ = "NA"
				        endif
				    else
				        origin_label$ = "NA"
				    endif
				else
				    origin_label$ = "NA"
				endif
				
				if reference_tier > 0
				    selectObject: textgrid_raw
				    current_interval = Get interval at time: reference_tier, cp_starttime
				    if current_interval > 0
				        intvl_label$ = Get label of interval: reference_tier, current_interval
				        # 修改判断条件：空字符串或与interval_label$不同时设为NA
				        if intvl_label$ = "" or (interval_label$ != "all" and intvl_label$ != interval_label$)
				            intvl_label$ = "NA"
				        endif
				        # 获取当前interval的开始时间、结束时间和duration
				        intvl_start = Get start time of interval: reference_tier, current_interval
				        intvl_end = Get end time of interval: reference_tier, current_interval
				        intvl_duration = intvl_end - intvl_start
				    else
				        intvl_label$ = "NA"
				        intvl_start = undefined
				        intvl_end = undefined
				        selectObject: wavfile_smoothed
				    		intvl_duration = Get total duration
				    endif
				else
				    current_interval = 0
				    intvl_label$ = "NA"
				    intvl_start = undefined
				    intvl_end = undefined
				    selectObject: wavfile_smoothed
				    intvl_duration = Get total duration
				endif
	
				if reference_tier > 0
					if intvl_label$ != "NA"
						if analyse_spectrum == 1
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
							line_count = line_count + 1
						else
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
							line_count = line_count + 1
						endif
						# Save for manual check
						if manual_check = 1
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index
					        Set value: line_count, 1, line_count
					        Set value: line_count, 2, cq
					        Set value: line_count, 6, f0_egg
							Set value: line_count, 9, f0_auto
							if cq < cq_min
								cq_min = cq
							endif
							if cq > cq_max
								cq_max = cq
							endif
							if f0_egg < f0_min
								f0_min = f0_egg
							endif
							if f0_egg > f0_max
								f0_max = f0_egg
							endif
							if f0_auto < f0_min
								f0_min = f0_auto
							endif
							if f0_auto > f0_max
								f0_max = f0_auto
							endif
							if line_count > index_max
								index_max = line_count
							endif
							if analyse_spectrum == 1
								if h1h2c > h1h2c_max
									h1h2c_max = h1h2c
								endif
								if h1h2c < h1h2c_min
									h1h2c_min = h1h2c
								endif
								#appendInfoLine: cq
								# Store in matrix
						        selectObject: analysis_matrix
								# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
								Set value: line_count, 1, line_count
						        Set value: line_count, 8, h1h2c
							endif
						endif
					endif
				else
					if analyse_spectrum == 1
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
						line_count = line_count + 1
					else
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
						line_count = line_count + 1
					endif
					# Save for manual check
					if manual_check = 1
						#appendInfoLine: cq
						# Store in matrix
				        selectObject: analysis_matrix
						# index
				        Set value: line_count, 1, line_count
				        Set value: line_count, 2, cq
				        Set value: line_count, 6, f0_egg
						Set value: line_count, 9, f0_auto
						if cq < cq_min
							cq_min = cq
						endif
						if cq > cq_max
							cq_max = cq
						endif
						if f0_egg < f0_min
							f0_min = f0_egg
						endif
						if f0_egg > f0_max
							f0_max = f0_egg
						endif
						if f0_auto < f0_min
								f0_min = f0_auto
						endif
						if f0_auto > f0_max
							f0_max = f0_auto
						endif
						if line_count > index_max
							index_max = line_count
						endif
						if analyse_spectrum == 1
							if h1h2c > h1h2c_max
								h1h2c_max = h1h2c
							endif
							if h1h2c < h1h2c_min
								h1h2c_min = h1h2c
							endif
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
							Set value: line_count, 1, line_count
					        Set value: line_count, 8, h1h2c
						endif
					endif
				endif
			endif
		endfor
		appendInfo: " dEGG;"
	endif
	
	
### Hybrid Method
	if hybrid_method == 1
		# Get hybrid EGG closure
		selectObject: wavfile_smoothed
		tg_25threshold = To TextGrid: "25cl", "25cl"
		selectObject: tg_pic
		num_pic = Get number of points: 1
		for i from 1 to num_pic - 1
		    # 获取当前peak点和下一个peak点的时间
			selectObject: tg_pic
		    peak1_time = Get time of point: 1, i
		    peak2_time = Get time of point: 1, i+1
			peakmid_time = peak1_time + ((peak2_time - peak1_time)/2)
		    # 获取peak点的振幅值
		    selectObject: wavfile_smoothed
			egg_peak_time = Get time of maximum: peak1_time, peakmid_time, interpolation$
			p = Get maximum: peak1_time, peakmid_time, interpolation$
			egg_valley_time = Get time of minimum: egg_peak_time, peak2_time, interpolation$
		    v = Get minimum: egg_peak_time, peak2_time, interpolation$
		    # 计算阈值
		    cr_25threshold = v + ((p - v) * hybrid_threshold)
			selectObject: wavfile_smoothed
		    c2_25 = Get nearest level crossing: 1, egg_valley_time, cr_25threshold, "left"
		    ##c1_25 = Get nearest level crossing: 1, egg_valley_time, cr_25threshold, "right"
			if c2_25 < peak2_time
				selectObject: tg_25threshold
				Insert point: 1, c2_25, "t2"
				##Insert point: 1, c1_25, "t1"
			endif
		endfor
	
		# Add: hybrid: tg_25threshold
		selectObject: tg_pic
		plusObject: tg_pdc
		plusObject: tg_25threshold
		tg_hybrid = Merge
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. tg_25threshold
		selectObject: tg_hybrid
		Insert interval tier: 4, "hybridCI"
		# Get all points from tier 1 (dEGGpeak) and tier 3 (25% Threshold)
		selectObject: tg_hybrid
		numPointsTier1 = Get number of points: 1
		numPointsTier2 = Get number of points: 2
		numPointsTier3 = Get number of points: 3
		
		# Create lists of all points
		d1_times# = zero#(numPointsTier1)
		d2_times# = zero#(numPointsTier2)
		d25_times# = zero#(numPointsTier3)
	
		for i from 1 to numPointsTier1
		    d1_times#[i] = Get time of point: 1, i
		endfor
		for i from 1 to numPointsTier2
		    d2_times#[i] = Get time of point: 2, i
		endfor
		for i from 1 to numPointsTier3
		    d25_times#[i] = Get time of point: 3, i
		endfor
	
		# Combine d1 and d25 times manually
		boundaryTimes# = zero#(numPointsTier1 + numPointsTier3)
		for i from 1 to numPointsTier1
		    boundaryTimes#[i] = d1_times#[i]
		endfor
		for i from 1 to numPointsTier3
		    boundaryTimes#[numPointsTier1 + i] = d25_times#[i]
		endfor
		
		# Sort the boundary times
		boundaryTimes# = sort#(boundaryTimes#)
		
		# Initialize final boundaries array
		finalBoundaryCount = 1
		finalBoundaries# = zero#(finalBoundaryCount)
		finalBoundaries#[1] = boundaryTimes#[1]
		
		# Check each interval between sorted boundary times
		for i from 1 to size(boundaryTimes#)-1
		    startTime = boundaryTimes#[i]
		    endTime = boundaryTimes#[i+1]
		    
		    # Add the end boundary to final array
		    finalBoundaryCount = finalBoundaryCount + 1
		    tempArray# = zero#(finalBoundaryCount)
		    for j from 1 to finalBoundaryCount-1
		        tempArray#[j] = finalBoundaries#[j]
		    endfor
		    tempArray#[finalBoundaryCount] = endTime
		    finalBoundaries# = tempArray#
		    
		    # Check if both boundaries are d1 points (no d25 in between)
		    is_d1_start = 0
		    is_d1_end = 0
		    for j from 1 to numPointsTier1
		        if abs(d1_times#[j] - startTime) < 0.0001
		            is_d1_start = 1
		        endif
		        if abs(d1_times#[j] - endTime) < 0.0001
		            is_d1_end = 1
		        endif
		    endfor
		    
		    # If both are d1 points, find the closest d2 point in between
		    if is_d1_start and is_d1_end
		        found_d2 = 0
		        closest_d2_time = 0
		        min_distance = 1000  
				# Large initial value
		        
		        for j from 1 to numPointsTier2
		            d2_time = d2_times#[j]
		            if d2_time > startTime and d2_time < endTime
		                distance = abs(d2_time - (startTime + endTime)/2)
		                if distance < min_distance
		                    min_distance = distance
		                    closest_d2_time = d2_time
		                    found_d2 = 1
		                endif
		            endif
		        endfor
		        
		        # If found a d2 point in between, add it as boundary
		        if found_d2
		            # Create new array with inserted value
		            finalBoundaryCount = finalBoundaryCount + 1
		            tempArray# = zero#(finalBoundaryCount)
		            inserted = 0
		            for k from 1 to finalBoundaryCount-1
		                if not inserted and closest_d2_time < finalBoundaries#[k]
		                    tempArray#[k] = closest_d2_time
		                    inserted = 1
		                    tempArray#[k+1] = finalBoundaries#[k]
		                else
		                    tempArray#[k+inserted] = finalBoundaries#[k]
		                endif
		            endfor
		            if not inserted
		                tempArray#[finalBoundaryCount] = closest_d2_time
		            endif
		            finalBoundaries# = sort#(tempArray#)
		        endif
		    endif
		endfor
	
		# Add boundaries to tier 4 (hybridCI)
		selectObject: tg_hybrid
		for i from 1 to size(finalBoundaries#)
		    time = finalBoundaries#[i]
		    Insert boundary: 4, time
		endfor
		### Add tiers hybridTime
		selectObject: tg_hybrid
		Insert point tier: 5, "hybridTime"
	
		# Get all points from tier 1 (dEGGpeak), tier 3 (25%)
		selectObject: tg_hybrid
		numPointsTier1 = Get number of points: 1
		numPointsTier3 = Get number of points: 3
		
		# Create lists of all points
		d1_times# = zero#(numPointsTier1)
		d25_times# = zero#(numPointsTier3)
		
		for i from 1 to numPointsTier1
		    d1_times#[i] = Get time of point: 1, i
		endfor
		for i from 1 to numPointsTier3
		    d25_times#[i] = Get time of point: 3, i
		endfor
	
		# Tiers now:
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. tg_25threshold
		# 4. hybridCI
		# 5. hybridTime
		### Adjust TIME tiers (hybridTime)
		### Taking numbers from the TextGrid (tier 4) and adjust hybridTime (tier 5)	
		selectObject: tg_hybrid
		nInt = do ("Get number of intervals...", 4)
	    # Assuming that the second interval annotates the first glottal period's close phase
		for i from 2 to nInt
			time_hybridCI = Get start time of interval: 4, i
			endtime_hybridCI = Get end time of interval: 4, i
			# Checking whether the wave is going up or down
			# Using 10% of the duration of the cycle to as following / previous check points
			# if the following interval is too far away, it takes the dur as 0.01
			hybrid_dur = (endtime_hybridCI - time_hybridCI)
			if hybrid_dur < 0.01
				hybrid_dur = hybrid_dur
			else
				hybrid_dur = 0.01
			endif
			fol_10perc_hybrid = time_hybridCI + hybrid_dur*0.1
			prev_10perc_hybrid = time_hybridCI - hybrid_dur*0.1
			selectObject: wavfile_smoothed
			point_y = Get value at time: 0, time_hybridCI, interpolation$
			fol_y = Get value at time: 0, fol_10perc_hybrid, interpolation$
			prev_y = Get value at time: 0, prev_10perc_hybrid, interpolation$
			selectObject: tg_hybrid
			if fol_y > point_y
				Insert point: 5, time_hybridCI, "t1"
			elif prev_y > point_y
				Insert point: 5, time_hybridCI, "t2"
			else
				Insert point: 5, time_hybridCI, "error"
			endif
		endfor
	
		### Label closing intervals and opening intervals of hybridCI (4), time tier: 5
		# Tiers now:
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. tg_25threshold
		# 4. hybridCI
		# 5. hybridTime
		selectObject: tg_hybrid
		nInt_hybrid = do ("Get number of intervals...", 4)
		for i from 2 to nInt_hybrid
			starttime = Get start time of interval: 4, i
			endtime = Get end time of interval: 4, i
			ptlabel_start$ = Get label of point: 5, i-1
			if i != nInt_hybrid
				ptlabel_end$ = Get label of point: 5, i
			endif
			if ptlabel_start$ == "t1" and ptlabel_end$ == "t2"
				Set interval text: 4, i, "cp"
			elif ptlabel_start$ == "t2" and ptlabel_end$ == "t1"
				Set interval text: 4, i, "op"
			endif
		endfor
		
		selectObject: tg_hybrid
		Insert point tier: 6, "check"
		# Tiers now:
		# 1. dEGGpeak
		# 2. dEGGvalley
		# 3. tg_25threshold
		# 4. hybridCI
		# 5. hybridTime
		# 6. hybridCheck
		# nowarn Save as text file: output_wav_textgrid$ + "testing_" + filename$ + "_EGG" + ".TextGrid"
		Remove tier: 1
		Remove tier: 1
		Remove tier: 1
		# Tiers now: 1. hybridCI, 2. hybridTime, 3. hybridCheck
		Remove tier: 2
		# Tiers now: 1. hybridCI, 2. hybridCheck

		# Pitch checking
		selectObject: tg_hybrid
		nInt_tocheck = do ("Get number of intervals...", 1)
		for i from 1 to nInt_tocheck
			selectObject: tg_hybrid
			intval$ = Get label of interval: 1, i
		    nextval$ = "NA"
		    if i < nInt_tocheck
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
				selectObject: tg_hybrid
				start_time_pitch_check = Get start time of interval: 1, i
				end_time_pitch_check = Get end time of interval: 1, i+1
				duration = end_time_pitch_check - start_time_pitch_check
        			freq = 1 / duration
				if remove_aperiodic_intervals == 1
					if freq > pitch_ceiling || freq < pitch_floor
						selectObject: tg_hybrid
						Set interval text: 1, i, ""
						Set interval text: 1, i+1, ""
					endif
				endif
			endif
		endfor
		############ Manual Check ############
		if manual_check = 1
			selectObject: stereo_combined
			plusObject: tg_hybrid
			Edit
			pauseScript: "Edit the TextGrid as needed, then click Continue to proceed."
			#appendInfo: " Checked"
		endif

		#### CALCULATION
		method$ = "hybrid"
		selectObject: tg_hybrid	
		nInt_t1 = do ("Get number of intervals...", 1)
		line_count = 1
		for i from 1 to nInt_t1
			selectObject: tg_hybrid
			intval$ = Get label of interval: 1, i
			#nextval$ = Get label of interval: 1, i+1
		    # Safely get next interval label with boundary check
		    nextval$ = "NA"
		    if i < nInt_t1
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
		        # Get timing information
		        cp_starttime = Get start time of interval: 1, i
		        cp_endtime = Get end time of interval: 1, i
		        op_endtime = Get end time of interval: 1, i+1
		        # Calculate durations
		        cp_dur = cp_endtime - cp_starttime
		        op_dur = op_endtime - cp_endtime
		        cycle_dur = cp_dur + op_dur
		        # Calculate parameters with safety checks
		        cq = if cycle_dur > 0 then cp_dur / cycle_dur else undefined fi
		        oq = if cycle_dur > 0 then op_dur / cycle_dur else undefined fi
		        cqoq = if cycle_dur > 0 then cp_dur / op_dur else undefined fi
		        f0_egg = if cycle_dur > 0 then 1 / cycle_dur else undefined fi
		        midpoint1 = if cp_dur > 0 then cp_dur/2 + cp_starttime else undefined fi
		        midpoint2 = if op_dur > 0 then op_dur/2 + cp_endtime else undefined fi
		        
				selectObject: pitch
				f0_auto = Get value at time: cp_starttime, "Hertz", "linear"
				#f0_auto = Get mean: cp_starttime, op_endtime, "Hertz"
				# or: midpoint1
				if f0_auto == undefined
		            f0_auto = f0_egg
		        endif

				selectObject: wavfile_smoothed
				p_pascal = Get maximum: cp_starttime, cp_endtime, interpolation$
				peak_time = Get time of maximum: cp_starttime, cp_endtime, interpolation$
				v_pascal = Get minimum: peak_time, op_endtime, interpolation$
				valley_time = Get time of minimum: peak_time, op_endtime, interpolation$
				# or: cp_starttime

				# SQ, contacting duration, decontacting duration, decontacting slope
				c_starttime = cp_starttime
		        c_endtime = peak_time
		        c_dur = c_endtime - c_starttime
		        o_starttime = peak_time
		        o_endtime = cp_endtime
		        o_dur = o_endtime - o_starttime
		        sq = o_dur / c_dur

				# 10% threshold (Marasek, 1997) Speed Quotient: sq_m
				thres_for_cend_ostart = p_pascal - (p_pascal - v_pascal) * skew_quotient_threshold
				c_end = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "left"
				o_start = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "right"
				thres_for_oend = v_pascal + (p_pascal - v_pascal) * skew_quotient_threshold
				o_end = Get nearest level crossing: 1, valley_time, thres_for_oend, "left"
				c_start_valley = Get nearest level crossing: 1, valley_time, thres_for_oend, "right"
				# Get c_start
				if i > 1
					selectObject: tg_hybrid
					prev_op_starttime = Get start time of interval: 1, i-1
				else
					selectObject: tg_hybrid
					prev_op_starttime = cp_starttime - 0.002
				endif
				curr_cp_starttime = cp_starttime
				selectObject: wavfile_smoothed
				prev_v_pascal = Get minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				prev_valley_time = Get time of minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				thres_for_cstart = prev_v_pascal + (p_pascal - prev_v_pascal) * skew_quotient_threshold
				c_start = Get nearest level crossing: 1, prev_valley_time, thres_for_cstart, "right"
				c_dur_marasek = c_end - c_start
				o_dur_marasek = o_end - o_start
				sq_m = c_dur_marasek / o_dur_marasek
				#sq_m = o_dur_marasek / c_dur_marasek
				# contact vs decontact ratio (10% threshold)
				contact_dur = o_start - c_end
				decontact_dur = c_start_valley - o_end
				con_dec_r = contact_dur / decontact_dur

		        # For checking purposes, a point on tier four refers to where the sq was taken
				selectObject: tg_hybrid
		        Insert point: 2, o_endtime, fixed$(sq,3)
				
				## Calculate Decontact Skewness: op_slope
				sound_part_dur = o_dur
				sum_diff = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = o_starttime + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = o_starttime + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_diff += sk2 - sk1
				endfor
				dec_skew = sum_diff / (slope_step - 1)

				## Calculate EGG Skewness: pv_slope
				selectObject: wavfile_smoothed
				sound_part_dur = valley_time - peak_time
				sum_differences = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = peak_time + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = peak_time + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_differences += sk2 - sk1
				endfor
				skewness = sum_differences / (slope_step - 1)
	
				#selectObject: deriv_smoothed
				#pic = Get value at time: 0, cp_starttime, interpolation$
				##pic_time = cp_starttime - c_start
				#pic_time = cp_starttime
				#pdc = Get minimum: cp_starttime, valley_time, interpolation$
				#pdc_instant = Get time of minimum: cp_starttime, valley_time, interpolation$
				##pdc_time = pdc_instant - c_start
				#pdc_time = pdc_instant

				# Get value from raw degg signal
				selectObject: deriv
				pic_starttime = cp_starttime - 0.0005
				pdc_endtime = op_endtime - 0.0005
				pic = Get maximum: pic_starttime, cp_endtime, interpolation$
				pic_time = Get time of maximum: pic_starttime, cp_endtime, interpolation$
				pdc = Get minimum: pic_time, pdc_endtime, interpolation$
				pdc_time = Get time of minimum: pic_time, pdc_endtime, interpolation$

				# Spectral
				if analyse_spectrum == 1
					resample = 11025
					selectObject: audio
					audio_11kHz = Resample: resample, 50
					selectObject: audio_11kHz
					lpc_audio_11kHz = To LPC (burg): 20, 0.025, 0.005, 34.45
					selectObject: audio_11kHz, lpc_audio_11kHz
					wavfile_invrsfltr = Filter (inverse)
					selectObject: wavfile_invrsfltr
					specgram = To Spectrogram: 0.03, resample/2, 0.01, 1, "Gaussian"
					h1h2c_sum = 0
					n_step = 10
					for x from 1 to n_step
						step1 = x * cycle_dur/n_step
						selectObject: specgram
						spectrum_slice1 = To Spectrum (slice): cp_starttime+step1
						selectObject: spectrum_slice1
						f0_spectrum_1 = Get frequency of nearest maximum: f0_egg
						h1c_1 = Get sound pressure level of nearest maximum: f0_spectrum_1
						h2c_freq_1 = f0_spectrum_1 * 2
						h2c_1 = Get sound pressure level of nearest maximum: h2c_freq_1
						h1h2c_1 = h1c_1 - h2c_1
						h1h2c_sum += h1h2c_1
						removeObject: spectrum_slice1
					endfor
					h1h2c = h1h2c_sum / n_step
					removeObject: specgram
					removeObject: audio_11kHz, lpc_audio_11kHz, wavfile_invrsfltr
				endif
		
		        # Insert points with error handling
				selectObject: tg_hybrid
		        if cq != undefined
		            nocheck Insert point: 2, cp_starttime, fixed$(cq, 3)
		        else
		            nocheck Insert point: 2, cp_starttime, "NA"
		        endif
		        selectObject: tg_hybrid
		        if midpoint1 != undefined and f0_auto != undefined
		            nocheck Insert point: 2, midpoint1, fixed$(f0_auto, 0)
		        else
		            nocheck Insert point: 2, cp_starttime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif
				selectObject: tg_hybrid
		        if midpoint2 != undefined and f0_egg != undefined
		            nocheck Insert point: 2, midpoint2, fixed$(f0_egg, 0)
		        else
		            nocheck Insert point: 2, cp_endtime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif

				if origin_label_tier > 0
					selectObject: textgrid_raw
					current_intvl = Get interval at time: origin_label_tier, cp_starttime
					if current_intvl > 0
				        origin_label$ = Get label of interval: origin_label_tier, current_intvl
				        if origin_label$ = ""
				            origin_label$ = "NA"
				        endif
				    else
				        origin_label$ = "NA"
				    endif
				else
				    origin_label$ = "NA"
				endif

				if reference_tier > 0
				    selectObject: textgrid_raw
				    current_interval = Get interval at time: reference_tier, cp_starttime
				    if current_interval > 0
				        intvl_label$ = Get label of interval: reference_tier, current_interval
				        # 修改判断条件：空字符串或与interval_label$不同时设为NA
				        if intvl_label$ = "" or (interval_label$ != "all" and intvl_label$ != interval_label$)
				            intvl_label$ = "NA"
				        endif
				        # 获取当前interval的开始时间、结束时间和duration
				        intvl_start = Get start time of interval: reference_tier, current_interval
				        intvl_end = Get end time of interval: reference_tier, current_interval
				        intvl_duration = intvl_end - intvl_start
				    else
				        intvl_label$ = "NA"
				        intvl_start = undefined
				        intvl_end = undefined
				        selectObject: wavfile_smoothed
				    		intvl_duration = Get total duration
				    endif
				else
				    current_interval = 0
				    intvl_label$ = "NA"
				    intvl_start = undefined
				    intvl_end = undefined
				    selectObject: wavfile_smoothed
				    intvl_duration = Get total duration
				endif
	
				if reference_tier > 0
					if intvl_label$ != "NA"
						if analyse_spectrum == 1
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
							line_count = line_count + 1
						else
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
							line_count = line_count + 1
						endif
						# Save datafor manual check
						if manual_check = 1
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index
					        Set value: line_count, 1, line_count
					        Set value: line_count, 3, cq
							Set value: line_count, 6, f0_egg
							Set value: line_count, 9, f0_auto
							if cq < cq_min
								cq_min = cq
							endif
							if cq > cq_max
								cq_max = cq
							endif
							if f0_egg < f0_min
								f0_min = f0_egg
							endif
							if f0_egg > f0_max
								f0_max = f0_egg
							endif
							if f0_auto < f0_min
								f0_min = f0_auto
							endif
							if f0_auto > f0_max
								f0_max = f0_auto
							endif
							if line_count > index_max
								index_max = line_count
							endif
							if analyse_spectrum == 1
								if h1h2c > h1h2c_max
									h1h2c_max = h1h2c
								endif
								if h1h2c < h1h2c_min
									h1h2c_min = h1h2c
								endif
								#appendInfoLine: cq
								# Store in matrix
						        selectObject: analysis_matrix
								# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
								Set value: line_count, 1, line_count
						        Set value: line_count, 8, h1h2c
							endif
						endif
					endif
				else
					if analyse_spectrum == 1
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
						line_count = line_count + 1
					else
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
						line_count = line_count + 1
					endif
					# Save datafor manual check
					if manual_check = 1
						#appendInfoLine: cq
						# Store in matrix
				        selectObject: analysis_matrix
						# index
				        Set value: line_count, 1, line_count
				        Set value: line_count, 3, cq
						Set value: line_count, 6, f0_egg
						Set value: line_count, 9, f0_auto
						if cq < cq_min
							cq_min = cq
						endif
						if cq > cq_max
							cq_max = cq
						endif
						if f0_egg < f0_min
							f0_min = f0_egg
						endif
						if f0_egg > f0_max
							f0_max = f0_egg
						endif
						if f0_auto < f0_min
							f0_min = f0_auto
						endif
						if f0_auto > f0_max
							f0_max = f0_auto
						endif
						if line_count > index_max
							index_max = line_count
						endif
						if analyse_spectrum == 1
							if h1h2c > h1h2c_max
								h1h2c_max = h1h2c
							endif
							if h1h2c < h1h2c_min
								h1h2c_min = h1h2c
							endif
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
							Set value: line_count, 1, line_count
					        Set value: line_count, 8, h1h2c
						endif
					endif
				endif
			endif
		endfor
		appendInfo: " Hybrid;"
	endif
	
	
### Threshold method
	if threshold_method == 1
		selectObject: wavfile_smoothed
		tg_37threshold = To TextGrid: "37cl", "37cl"
		# 开头补点，利用dEGG
		selectObject: tg_pic
		c1_add = Get time of point: 1, 1
		selectObject: tg_37threshold
		Insert point: 1, c1_add, "t1"

		selectObject: tg_pic
		num_pic = Get number of points: 1
		# 处理每个dEGG peak区间
		for i from 1 to num_pic - 1
		    # 获取当前peak点和下一个peak点的时间
			selectObject: tg_pic
		    peak1_time = Get time of point: 1, i
		    peak2_time = Get time of point: 1, i+1
			peakmid_time = peak1_time + ((peak2_time - peak1_time)/2)
		    # 获取peak点的振幅值
		    selectObject: wavfile_smoothed
			egg_peak_time = Get time of maximum: peak1_time, peakmid_time, interpolation$
			p = Get maximum: peak1_time, peakmid_time, interpolation$
			egg_valley_time = Get time of minimum: egg_peak_time, peak2_time, interpolation$
		    v = Get minimum: egg_peak_time, peak2_time, interpolation$
		    # 计算阈值
		    cr_threshold = v + ((p - v) * closing_threshold)
			selectObject: wavfile_smoothed
		    c2 = Get nearest level crossing: 1, egg_valley_time, cr_threshold, "left"
		    c1 = Get nearest level crossing: 1, egg_valley_time, cr_threshold, "right"
		    # 在标注层中插入极值点
		    selectObject: tg_37threshold
		    Insert point: 1, c2, "t2"
		    Insert point: 1, c1, "t1"
		endfor
	
		selectObject: tg_37threshold
		Insert interval tier: 2, "37CI"
		# Tiers now:
		# 1. tg_37threshold
		# 2. 37CI
		# Get all points from tier 1 (3/7 Threshold)
		selectObject: tg_37threshold
		numPointsTier37 = Get number of points: 1
		for i from 1 to numPointsTier37
			selectObject: tg_37threshold
		    time37 = Get time of point: 1, i
			selectObject: tg_37threshold
			Insert boundary: 2, time37
		endfor
	
		### Label closing intervals and opening intervals of 37CI (2), time tier: 1
		selectObject: tg_37threshold
		nInt_37 = do ("Get number of intervals...", 2)
		for i from 2 to nInt_37
			starttime = Get start time of interval: 2, i
			endtime = Get end time of interval: 2, i
			ptlabel_start$ = Get label of point: 1, i-1
			if i != nInt_37
				ptlabel_end$ = Get label of point: 1, i
			endif
			if ptlabel_start$ == "t1" and ptlabel_end$ == "t2"
				Set interval text: 2, i, "cp"
			elif ptlabel_start$ == "t2" and ptlabel_end$ == "t1"
				Set interval text: 2, i, "op"
			endif
		endfor
		
		selectObject: tg_37threshold
		Insert point tier: 3, "check"
		Remove tier: 1
			# tg_37threshold -> 37CI becomes tier 1, 37check becomes tier 2

		# Pitch checking
		selectObject: tg_37threshold
		nInt_tocheck = do ("Get number of intervals...", 1)
		for i from 1 to nInt_tocheck
			selectObject: tg_37threshold
			intval$ = Get label of interval: 1, i
		    nextval$ = "NA"
		    if i < nInt_tocheck
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
				selectObject: tg_37threshold
				start_time_pitch_check = Get start time of interval: 1, i
				end_time_pitch_check = Get end time of interval: 1, i+1
				duration = end_time_pitch_check - start_time_pitch_check
        			freq = 1 / duration
				if remove_aperiodic_intervals == 1
					if freq > pitch_ceiling || freq < pitch_floor
						selectObject: tg_37threshold
						Set interval text: 1, i, ""
						Set interval text: 1, i+1, ""
					endif
				endif
			endif
		endfor
		############ Manual Check ############
		if manual_check = 1
			selectObject: stereo_combined
			plusObject: tg_37threshold
			Edit
			pauseScript: "Edit the TextGrid as needed, then click Continue to proceed."
			#appendInfo: " Checked"
		endif

		#### CALCULATION
		method$ = "threshold"
		selectObject: tg_37threshold	
		nInt_t1 = do ("Get number of intervals...", 1)
		line_count = 1
		for i from 1 to nInt_t1
			selectObject: tg_37threshold
			intval$ = Get label of interval: 1, i
			#nextval$ = Get label of interval: 1, i+1
		    # Safely get next interval label with boundary check
		    nextval$ = "NA"
		    if i < nInt_t1
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
		        # Get timing information
		        cp_starttime = Get start time of interval: 1, i
		        cp_endtime = Get end time of interval: 1, i
		        op_endtime = Get end time of interval: 1, i+1
		        # Calculate durations
		        cp_dur = cp_endtime - cp_starttime
		        op_dur = op_endtime - cp_endtime
		        cycle_dur = cp_dur + op_dur
		        # Calculate parameters with safety checks
		        cq = if cycle_dur > 0 then cp_dur / cycle_dur else undefined fi
		        oq = if cycle_dur > 0 then op_dur / cycle_dur else undefined fi
		        cqoq = if cycle_dur > 0 then cp_dur / op_dur else undefined fi
		        f0_egg = if cycle_dur > 0 then 1 / cycle_dur else undefined fi
		        midpoint1 = if cp_dur > 0 then cp_dur/2 + cp_starttime else undefined fi
		        midpoint2 = if op_dur > 0 then op_dur/2 + cp_endtime else undefined fi
		        
				selectObject: pitch
				f0_auto = Get value at time: cp_starttime, "Hertz", "linear"
				#f0_auto = Get mean: cp_starttime, op_endtime, "Hertz"
				# or: midpoint1
				if f0_auto == undefined
		            f0_auto = f0_egg
		        endif

				selectObject: wavfile_smoothed
				p_pascal = Get maximum: cp_starttime, cp_endtime, interpolation$
				peak_time = Get time of maximum: cp_starttime, cp_endtime, interpolation$
				v_pascal = Get minimum: peak_time, op_endtime, interpolation$
				valley_time = Get time of minimum: peak_time, op_endtime, interpolation$
				# or: cp_starttime

				# SQ, contacting duration, decontacting duration, decontacting slope
				c_starttime = cp_starttime
		        c_endtime = peak_time
		        c_dur = c_endtime - c_starttime
		        o_starttime = peak_time
		        o_endtime = cp_endtime
		        o_dur = o_endtime - o_starttime
		        sq = o_dur / c_dur

				# 10% threshold (Marasek, 1997) Speed Quotient: sq_m
				thres_for_cend_ostart = p_pascal - (p_pascal - v_pascal) * skew_quotient_threshold
				c_end = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "left"
				o_start = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "right"
				thres_for_oend = v_pascal + (p_pascal - v_pascal) * skew_quotient_threshold
				o_end = Get nearest level crossing: 1, valley_time, thres_for_oend, "left"
				c_start_valley = Get nearest level crossing: 1, valley_time, thres_for_oend, "right"
				# Get c_start
				if i > 1
					selectObject: tg_37threshold
					prev_op_starttime = Get start time of interval: 1, i-1
				else
					selectObject: tg_37threshold
					prev_op_starttime = cp_starttime - 0.002
				endif
				curr_cp_starttime = cp_starttime
				selectObject: wavfile_smoothed
				prev_v_pascal = Get minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				prev_valley_time = Get time of minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				thres_for_cstart = prev_v_pascal + (p_pascal - prev_v_pascal) * skew_quotient_threshold
				c_start = Get nearest level crossing: 1, prev_valley_time, thres_for_cstart, "right"
				c_dur_marasek = c_end - c_start
				o_dur_marasek = o_end - o_start
				sq_m = c_dur_marasek / o_dur_marasek
				#sq_m = o_dur_marasek / c_dur_marasek
				# contact vs decontact ratio (10% threshold)
				contact_dur = o_start - c_end
				decontact_dur = c_start_valley - o_end
				con_dec_r = contact_dur / decontact_dur
		
		        # For checking purposes, a point on tier four refers to where the sq was taken
				selectObject: tg_37threshold
		        Insert point: 2, o_endtime, fixed$(sq,3)
				
				## Calculate Decontact Skewness: op_slope
				sound_part_dur = o_dur
				sum_diff = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = o_starttime + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = o_starttime + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_diff += sk2 - sk1
				endfor
				dec_skew = sum_diff / (slope_step - 1)

				## Calculate EGG Skewness: pv_slope
				selectObject: wavfile_smoothed
				sound_part_dur = valley_time - peak_time
				sum_differences = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = peak_time + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = peak_time + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_differences += sk2 - sk1
				endfor
				skewness = sum_differences / (slope_step - 1)
	
				#selectObject: deriv_smoothed
				#pic_start = peak_time - 0.002
				#pic = Get maximum: pic_start, peak_time, interpolation$
				#pic_instant = Get time of maximum: pic_start, peak_time, interpolation$
				##pic_time = pic_instant - c_start
				#pic_time = pic_instant
				#pdc = Get minimum: peak_time, valley_time, interpolation$
				#pdc_instant = Get time of minimum: peak_time, valley_time, interpolation$
				##pdc_time = pdc_instant - c_start
				#pdc_time = pdc_instant

				# Get value from raw degg signal
				selectObject: deriv
				pic_starttime = cp_starttime - 0.0005
				pdc_endtime = op_endtime - 0.0005
				pic = Get maximum: pic_starttime, cp_endtime, interpolation$
				pic_time = Get time of maximum: pic_starttime, cp_endtime, interpolation$
				pdc = Get minimum: pic_time, pdc_endtime, interpolation$
				pdc_time = Get time of minimum: pic_time, pdc_endtime, interpolation$
		
				# Spectral
				if analyse_spectrum == 1
					resample = 11025
					selectObject: audio
					audio_11kHz = Resample: resample, 50
					selectObject: audio_11kHz
					lpc_audio_11kHz = To LPC (burg): 20, 0.025, 0.005, 34.45
					selectObject: audio_11kHz, lpc_audio_11kHz
					wavfile_invrsfltr = Filter (inverse)
					selectObject: wavfile_invrsfltr
					specgram = To Spectrogram: 0.03, resample/2, 0.01, 1, "Gaussian"
					h1h2c_sum = 0
					n_step = 10
					for x from 1 to n_step
						step1 = x * cycle_dur/n_step
						selectObject: specgram
						spectrum_slice1 = To Spectrum (slice): cp_starttime+step1
						selectObject: spectrum_slice1
						f0_spectrum_1 = Get frequency of nearest maximum: f0_egg
						h1c_1 = Get sound pressure level of nearest maximum: f0_spectrum_1
						h2c_freq_1 = f0_spectrum_1 * 2
						h2c_1 = Get sound pressure level of nearest maximum: h2c_freq_1
						h1h2c_1 = h1c_1 - h2c_1
						h1h2c_sum += h1h2c_1
						removeObject: spectrum_slice1
					endfor
					h1h2c = h1h2c_sum / n_step
					removeObject: specgram
					removeObject: audio_11kHz, lpc_audio_11kHz, wavfile_invrsfltr
				endif

		        # Insert points with error handling
				selectObject: tg_37threshold
		        if cq != undefined
		            nocheck Insert point: 2, cp_starttime, fixed$(cq, 3)
		        else
		            nocheck Insert point: 2, cp_starttime, "NA"
		        endif
		        selectObject: tg_37threshold
		        if midpoint1 != undefined and f0_auto != undefined
		            nocheck Insert point: 2, midpoint1, fixed$(f0_auto, 0)
		        else
		            nocheck Insert point: 2, cp_starttime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif
				selectObject: tg_37threshold
		        if midpoint2 != undefined and f0_egg != undefined
		            nocheck Insert point: 2, midpoint2, fixed$(f0_egg, 0)
		        else
		            nocheck Insert point: 2, cp_endtime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif

				if origin_label_tier > 0
					selectObject: textgrid_raw
					current_intvl = Get interval at time: origin_label_tier, cp_starttime
					if current_intvl > 0
				        origin_label$ = Get label of interval: origin_label_tier, current_intvl
				        if origin_label$ = ""
				            origin_label$ = "NA"
				        endif
				    else
				        origin_label$ = "NA"
				    endif
				else
				    origin_label$ = "NA"
				endif

				if reference_tier > 0
				    selectObject: textgrid_raw
				    current_interval = Get interval at time: reference_tier, cp_starttime
				    if current_interval > 0
				        intvl_label$ = Get label of interval: reference_tier, current_interval
				        # 修改判断条件：空字符串或与interval_label$不同时设为NA
				        if intvl_label$ = "" or (interval_label$ != "all" and intvl_label$ != interval_label$)
				            intvl_label$ = "NA"
				        endif
				        # 获取当前interval的开始时间、结束时间和duration
				        intvl_start = Get start time of interval: reference_tier, current_interval
				        intvl_end = Get end time of interval: reference_tier, current_interval
				        intvl_duration = intvl_end - intvl_start
				    else
				        intvl_label$ = "NA"
				        intvl_start = undefined
				        intvl_end = undefined
				        selectObject: wavfile_smoothed
				    		intvl_duration = Get total duration
				    endif
				else
				    current_interval = 0
				    intvl_label$ = "NA"
				    intvl_start = undefined
				    intvl_end = undefined
				    selectObject: wavfile_smoothed
				    intvl_duration = Get total duration
				endif
	
				if reference_tier > 0
					if intvl_label$ != "NA"
						if analyse_spectrum == 1
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
							line_count = line_count + 1
						else
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
							line_count = line_count + 1
						endif
						# Save data for manual check
						if manual_check = 1
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index
					        Set value: line_count, 1, line_count
					        Set value: line_count, 4, cq
					        Set value: line_count, 7, f0_egg
							Set value: line_count, 9, f0_auto
							if cq < cq_min
								cq_min = cq
							endif
							if cq > cq_max
								cq_max = cq
							endif
							if f0_egg < f0_min
								f0_min = f0_egg
							endif
							if f0_egg > f0_max
								f0_max = f0_egg
							endif
							if f0_auto < f0_min
								f0_min = f0_auto
							endif
							if f0_auto > f0_max
								f0_max = f0_auto
							endif
							if line_count > index_max
								index_max = line_count
							endif
							if analyse_spectrum == 1
								if h1h2c > h1h2c_max
									h1h2c_max = h1h2c
								endif
								if h1h2c < h1h2c_min
									h1h2c_min = h1h2c
								endif
								#appendInfoLine: cq
								# Store in matrix
						        selectObject: analysis_matrix
								# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
								Set value: line_count, 1, line_count
						        Set value: line_count, 8, h1h2c
							endif
						endif
					endif
				else
					if analyse_spectrum == 1
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
						line_count = line_count + 1
					else
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
						line_count = line_count + 1
					endif
					# Save data for manual check
					if manual_check = 1
						#appendInfoLine: cq
						# Store in matrix
				        selectObject: analysis_matrix
						# index
				        Set value: line_count, 1, line_count
				        Set value: line_count, 4, cq
				        Set value: line_count, 7, f0_egg
						Set value: line_count, 9, f0_auto
						if cq < cq_min
							cq_min = cq
						endif
						if cq > cq_max
							cq_max = cq
						endif
						if f0_egg < f0_min
							f0_min = f0_egg
						endif
						if f0_egg > f0_max
							f0_max = f0_egg
						endif
						if f0_auto < f0_min
							f0_min = f0_auto
						endif
						if f0_auto > f0_max
							f0_max = f0_auto
						endif
						if line_count > index_max
							index_max = line_count
						endif
						if analyse_spectrum == 1
							if h1h2c > h1h2c_max
								h1h2c_max = h1h2c
							endif
							if h1h2c < h1h2c_min
								h1h2c_min = h1h2c
							endif
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
							Set value: line_count, 1, line_count
					        Set value: line_count, 8, h1h2c
						endif
					endif
				endif
			endif
		endfor
		appendInfo: " Threshold;"
	endif



##### Henry Tehrani method (dEGG + EGG & DC Component's cross)
	if henry_Tehrani_method == 1
		###DC Component
		#selectObject: wavfile_smoothed
		# Invert DC Component (for EGG + -DC)
		selectObject: dc_component
		#Formula: -self
		Multiply: -1
		selectObject: wavfile_smoothed, dc_component
		egg_dc1 = Combine to stereo
		selectObject: egg_dc1
		egg_dc = Convert to mono

		selectObject: egg_dc
		# Create point tier
		tg_HTtime = To TextGrid: "HTtime", "HTtime"
		selectObject: egg_dc
		tg_HTCI = To TextGrid: "HTCI", ""
		selectObject: tg_pic
		plusObject: tg_HTtime
		plusObject: tg_HTCI
		tg_HT = Merge
		Insert point tier: 4, "HTcheck"
		num_pic_tier = Get number of points: 1
		for i from 1 to num_pic_tier - 1
			selectObject: tg_HT
		    t1 = Get time of point: 1, i
		    t2 = Get time of point: 1, i+1
		    selectObject: egg_dc
		    valley_to_seek = Get time of minimum: t1, t2, interpolation$
			peak_to_seek = Get time of maximum: t1, valley_to_seek, interpolation$
			#d2 = Get nearest level crossing: 1, valley_to_seek, 0, "left"
			d2 = Get nearest level crossing: 1, peak_to_seek, 0, "right"
			if d2 >= t2
				selectObject: wavfile_smoothed
		   		valley_to_seek_d2 = Get time of minimum: t1, t2, interpolation$
				selectObject: deriv_smoothed
			    d2 = Get time of minimum: t1, valley_to_seek_d2, interpolation$
			else
				d2 = d2
			endif
			selectObject: tg_HT
			Insert point: 2, t1, "d1"
			Insert point: 2, d2, "d2"
		endfor
		# Get all points from tier 1
		selectObject: tg_HT
		Remove tier: 1
		selectObject: tg_HT
		numPointsTierHT = Get number of points: 1
		for i from 1 to numPointsTierHT
			selectObject: tg_HT
		    timeHT = Get time of point: 1, i
			selectObject: tg_HT
			Insert boundary: 2, timeHT
		endfor
		### Label closing intervals and opening intervals
		selectObject: tg_HT
		nInt_HT = do ("Get number of intervals...", 2)
		for i from 2 to nInt_HT
			starttime = Get start time of interval: 2, i
			endtime = Get end time of interval: 2, i
			ptlabel_start$ = Get label of point: 1, i-1
			if i != nInt_HT
				ptlabel_end$ = Get label of point: 1, i
			endif
			if ptlabel_start$ == "d1" and ptlabel_end$ == "d2"
				Set interval text: 2, i, "cp"
			elif ptlabel_start$ == "d2" and ptlabel_end$ == "d1"
				Set interval text: 2, i, "op"
			endif
		endfor
		selectObject: tg_HT
		Remove tier: 1
		# Tiers now: 1.HTCI; 2.HTcheck

		# Pitch checking
		selectObject: tg_HT
		nInt_tocheck = do ("Get number of intervals...", 1)
		for i from 1 to nInt_tocheck
			selectObject: tg_HT
			intval$ = Get label of interval: 1, i
		    nextval$ = "NA"
		    if i < nInt_tocheck
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
				selectObject: tg_HT
				start_time_pitch_check = Get start time of interval: 1, i
				end_time_pitch_check = Get end time of interval: 1, i+1
				duration = end_time_pitch_check - start_time_pitch_check
        			freq = 1 / duration
				if remove_aperiodic_intervals == 1
					if freq > pitch_ceiling || freq < pitch_floor
						selectObject: tg_HT
						Set interval text: 1, i, ""
						Set interval text: 1, i+1, ""
					endif
				endif
			endif
		endfor

		############ Manual Check ############
		if manual_check = 1
			selectObject: egg_dc
			egg_dc_resampled = Resample: 22050, 50
			selectObject: stereo_combined, egg_dc_resampled
			check_HT_sound = Combine to stereo
			selectObject: check_HT_sound
			plusObject: tg_HT
			Edit
			pauseScript: "Edit the TextGrid as needed, then click Continue to proceed."
			#appendInfo: " Checked"
			removeObject: egg_dc_resampled
			removeObject: check_HT_sound
		endif
		
		############ CALCULATION ############
		#### CALCULATION
		method$ = "HenryTehrani"
		selectObject: tg_HT
		nInt_t1 = do ("Get number of intervals...", 1)
		line_count = 1
		for i from 1 to nInt_t1
			selectObject: tg_HT
			intval$ = Get label of interval: 1, i
			#nextval$ = Get label of interval: 1, i+1
		    # Safely get next interval label with boundary check
		    nextval$ = "NA"
		    if i < nInt_t1
		        nextval$ = Get label of interval: 1, i+1
		    endif
		    if intval$ == "cp" and nextval$ == "op"
		        # Get timing information
		        cp_starttime = Get start time of interval: 1, i
		        cp_endtime = Get end time of interval: 1, i
		        op_endtime = Get end time of interval: 1, i+1
		        # Calculate durations
		        cp_dur = cp_endtime - cp_starttime
		        op_dur = op_endtime - cp_endtime
		        cycle_dur = cp_dur + op_dur
		        # Calculate parameters with safety checks
		        cq = if cycle_dur > 0 then cp_dur / cycle_dur else undefined fi
		        oq = if cycle_dur > 0 then op_dur / cycle_dur else undefined fi
		        cqoq = if cycle_dur > 0 then cp_dur / op_dur else undefined fi
		        f0_egg = if cycle_dur > 0 then 1 / cycle_dur else undefined fi
		        midpoint1 = if cp_dur > 0 then cp_dur/2 + cp_starttime else undefined fi
		        midpoint2 = if op_dur > 0 then op_dur/2 + cp_endtime else undefined fi
		        
				selectObject: pitch
				f0_auto = Get value at time: cp_starttime, "Hertz", "linear"
				#f0_auto = Get mean: cp_starttime, op_endtime, "Hertz"
				# or: midpoint1
				if f0_auto == undefined
		            f0_auto = f0_egg
		        endif

				selectObject: wavfile_smoothed
				p_pascal = Get maximum: cp_starttime, cp_endtime, interpolation$
				peak_time = Get time of maximum: cp_starttime, cp_endtime, interpolation$
				v_pascal = Get minimum: peak_time, op_endtime, interpolation$
				valley_time = Get time of minimum: peak_time, op_endtime, interpolation$
				# or: cp_starttime

				# SQ, contacting duration, decontacting duration, decontacting slope
				c_starttime = cp_starttime
		        c_endtime = peak_time
		        c_dur = c_endtime - c_starttime
		        o_starttime = peak_time
		        o_endtime = cp_endtime
		        o_dur = o_endtime - o_starttime
		        sq = o_dur / c_dur

				# 10% threshold (Marasek, 1997) Speed Quotient: sq_m
				thres_for_cend_ostart = p_pascal - (p_pascal - v_pascal) * skew_quotient_threshold
				c_end = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "left"
				o_start = Get nearest level crossing: 1, peak_time, thres_for_cend_ostart, "right"
				thres_for_oend = v_pascal + (p_pascal - v_pascal) * skew_quotient_threshold
				o_end = Get nearest level crossing: 1, valley_time, thres_for_oend, "left"
				c_start_valley = Get nearest level crossing: 1, valley_time, thres_for_oend, "right"
				# Get c_start
				if i > 1
					selectObject: tg_HT
					prev_op_starttime = Get start time of interval: 1, i-1
				else
					selectObject: tg_HT
					prev_op_starttime = cp_starttime - 0.002
				endif
				curr_cp_starttime = cp_starttime
				selectObject: wavfile_smoothed
				prev_v_pascal = Get minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				prev_valley_time = Get time of minimum: prev_op_starttime, curr_cp_starttime, interpolation$
				thres_for_cstart = prev_v_pascal + (p_pascal - prev_v_pascal) * skew_quotient_threshold
				c_start = Get nearest level crossing: 1, prev_valley_time, thres_for_cstart, "right"
				c_dur_marasek = c_end - c_start
				o_dur_marasek = o_end - o_start
				sq_m = c_dur_marasek / o_dur_marasek
				#sq_m = o_dur_marasek / c_dur_marasek
				# contact vs decontact ratio (10% threshold)
				contact_dur = o_start - c_end
				decontact_dur = c_start_valley - o_end
				con_dec_r = contact_dur / decontact_dur

		        # For checking purposes, a point on tier four refers to where the sq was taken
				selectObject: tg_HT
		        Insert point: 2, o_endtime, fixed$(sq,3)
				
				## Calculate Decontact Skewness: op_slope
				sound_part_dur = o_dur
				sum_diff = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = o_starttime + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = o_starttime + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_diff += sk2 - sk1
				endfor
				dec_skew = sum_diff / (slope_step - 1)

				## Calculate EGG Skewness: pv_slope
				selectObject: wavfile_smoothed
				sound_part_dur = valley_time - peak_time
				sum_differences = 0
				for s from 0 to slope_step - 1
					skew_timestep1 = s * sound_part_dur/slope_step
					skew_time1 = peak_time + skew_timestep1
					selectObject: wavfile_smoothed
					sk1 = Get value at time: 1, skew_time1, interpolation$
					skew_timestep2 = (s+1) * sound_part_dur/slope_step
					skew_time2 = peak_time + skew_timestep2
					selectObject: wavfile_smoothed
					sk2 = Get value at time: 1, skew_time2, interpolation$
					sum_differences += sk2 - sk1
				endfor
				skewness = sum_differences / (slope_step - 1)
	
				#selectObject: deriv_smoothed
				#pic = Get value at time: 0, cp_starttime, interpolation$
				##pic_time = cp_starttime - c_start
				#pic_time = cp_starttime
				#pdc = Get minimum: cp_starttime, valley_time, interpolation$
				#pdc_instant = Get time of minimum: cp_starttime, valley_time, interpolation$
				##pdc_time = pdc_instant - c_start
				#pdc_time = pdc_instant

				# Get value from raw degg signal
				selectObject: deriv
				pic_starttime = cp_starttime - 0.0005
				pdc_endtime = op_endtime - 0.0005
				pic = Get maximum: pic_starttime, cp_endtime, interpolation$
				pic_time = Get time of maximum: pic_starttime, cp_endtime, interpolation$
				pdc = Get minimum: pic_time, pdc_endtime, interpolation$
				pdc_time = Get time of minimum: pic_time, pdc_endtime, interpolation$

				# Spectral
				if analyse_spectrum == 1
					resample = 11025
					selectObject: audio
					audio_11kHz = Resample: resample, 50
					selectObject: audio_11kHz
					lpc_audio_11kHz = To LPC (burg): 20, 0.025, 0.005, 34.45
					selectObject: audio_11kHz, lpc_audio_11kHz
					wavfile_invrsfltr = Filter (inverse)
					selectObject: wavfile_invrsfltr
					specgram = To Spectrogram: 0.03, resample/2, 0.01, 1, "Gaussian"
					h1h2c_sum = 0
					n_step = 10
					for x from 1 to n_step
						step1 = x * cycle_dur/n_step
						selectObject: specgram
						spectrum_slice1 = To Spectrum (slice): cp_starttime+step1
						selectObject: spectrum_slice1
						f0_spectrum_1 = Get frequency of nearest maximum: f0_egg
						h1c_1 = Get sound pressure level of nearest maximum: f0_spectrum_1
						h2c_freq_1 = f0_spectrum_1 * 2
						h2c_1 = Get sound pressure level of nearest maximum: h2c_freq_1
						h1h2c_1 = h1c_1 - h2c_1
						h1h2c_sum += h1h2c_1
						removeObject: spectrum_slice1
					endfor
					h1h2c = h1h2c_sum / n_step
					removeObject: specgram
					removeObject: audio_11kHz, lpc_audio_11kHz, wavfile_invrsfltr
				endif

		        # Insert points with error handling
				selectObject: tg_HT
		        if cq != undefined
		            nocheck Insert point: 2, cp_starttime, fixed$(cq, 3)
		        else
		            nocheck Insert point: 2, cp_starttime, "NA"
		        endif
		        selectObject: tg_HT
		        if midpoint1 != undefined and f0_auto != undefined
		            nocheck Insert point: 2, midpoint1, fixed$(f0_auto, 0)
		        else
		            nocheck Insert point: 2, cp_starttime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif
				selectObject: tg_HT
		        if midpoint2 != undefined and f0_egg != undefined
		            nocheck Insert point: 2, midpoint2, fixed$(f0_egg, 0)
		        else
		            nocheck Insert point: 2, cp_endtime + 0.001, "NA"
					# Small offset if midpoint undefined
		        endif

				if origin_label_tier > 0
					selectObject: textgrid_raw
					current_intvl = Get interval at time: origin_label_tier, cp_starttime
					if current_intvl > 0
				        origin_label$ = Get label of interval: origin_label_tier, current_intvl
				        if origin_label$ = ""
				            origin_label$ = "NA"
				        endif
				    else
				        origin_label$ = "NA"
				    endif
				else
				    origin_label$ = "NA"
				endif

				if reference_tier > 0
				    selectObject: textgrid_raw
				    current_interval = Get interval at time: reference_tier, cp_starttime
				    if current_interval > 0
				        intvl_label$ = Get label of interval: reference_tier, current_interval
				        # 修改判断条件：空字符串或与interval_label$不同时设为NA
				        if intvl_label$ = "" or (interval_label$ != "all" and intvl_label$ != interval_label$)
				            intvl_label$ = "NA"
				        endif
				        # 获取当前interval的开始时间、结束时间和duration
				        intvl_start = Get start time of interval: reference_tier, current_interval
				        intvl_end = Get end time of interval: reference_tier, current_interval
				        intvl_duration = intvl_end - intvl_start
				    else
				        intvl_label$ = "NA"
				        intvl_start = undefined
				        intvl_end = undefined
				        selectObject: wavfile_smoothed
				    		intvl_duration = Get total duration
				    endif
				else
				    current_interval = 0
				    intvl_label$ = "NA"
				    intvl_start = undefined
				    intvl_end = undefined
				    selectObject: wavfile_smoothed
				    intvl_duration = Get total duration
				endif
	
				if reference_tier > 0
					if intvl_label$ != "NA"
						if analyse_spectrum == 1
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
							line_count = line_count + 1
						else
							appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
							line_count = line_count + 1
						endif
						# Save data for manual check
						if manual_check = 1
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index
					        Set value: line_count, 1, line_count
					        Set value: line_count, 5, cq
							Set value: line_count, 6, f0_egg
							Set value: line_count, 9, f0_auto
							if cq < cq_min
								cq_min = cq
							endif
							if cq > cq_max
								cq_max = cq
							endif
							if f0_egg < f0_min
								f0_min = f0_egg
							endif
							if f0_egg > f0_max
								f0_max = f0_egg
							endif
							if f0_auto < f0_min
								f0_min = f0_auto
							endif
							if f0_auto > f0_max
								f0_max = f0_auto
							endif
							if line_count > index_max
								index_max = line_count
							endif
							if analyse_spectrum == 1
								if h1h2c > h1h2c_max
									h1h2c_max = h1h2c
								endif
								if h1h2c < h1h2c_min
									h1h2c_min = h1h2c
								endif
								#appendInfoLine: cq
								# Store in matrix
						        selectObject: analysis_matrix
								# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
								Set value: line_count, 1, line_count
						        Set value: line_count, 8, h1h2c
							endif
						endif
					endif
				else
					if analyse_spectrum == 1
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, ",", h1h2c, newline$
						line_count = line_count + 1
					else
						appendFile: output_file$, method$, ",", filename$, ",", origin_label$, ",", speaker$, ",", current_interval, ",", intvl_label$, ",", intvl_duration, ",", line_count, ",", cycle_dur, ",", cp_dur, ",", op_dur, ",", cq, ",", oq, ",", cqoq, ",", p_pascal, ",", v_pascal, ",", pic, ",", pic_time, ",", pdc, ",", pdc_time, ",", f0_egg, ",", f0_auto, ",", skewness, ",", c_dur, ",", o_dur, ",", sq, ",", dec_skew, ",", sq_m, ",", c_dur_marasek, ",", o_dur_marasek, ",", contact_dur, ",", decontact_dur, ",", con_dec_r, newline$
						line_count = line_count + 1
					endif
					# Save data for manual check
					if manual_check = 1
						#appendInfoLine: cq
						# Store in matrix
				        selectObject: analysis_matrix
						# index
				        Set value: line_count, 1, line_count
				        Set value: line_count, 5, cq
						Set value: line_count, 6, f0_egg
						Set value: line_count, 9, f0_auto
						if cq < cq_min
							cq_min = cq
						endif
						if cq > cq_max
							cq_max = cq
						endif
						if f0_egg < f0_min
							f0_min = f0_egg
						endif
						if f0_egg > f0_max
							f0_max = f0_egg
						endif
						if f0_auto < f0_min
							f0_min = f0_auto
						endif
						if f0_auto > f0_max
							f0_max = f0_auto
						endif
						if line_count > index_max
							index_max = line_count
						endif
						if analyse_spectrum == 1
							if h1h2c > h1h2c_max
								h1h2c_max = h1h2c
							endif
							if h1h2c < h1h2c_min
								h1h2c_min = h1h2c
							endif
							#appendInfoLine: cq
							# Store in matrix
					        selectObject: analysis_matrix
							# index (vẫn thêm vào, phòng trường hợp số period nhiều hơn các methods khác)
							Set value: line_count, 1, line_count
					        Set value: line_count, 8, h1h2c
						endif
					endif
				endif
			endif
		endfor

		appendInfo: " Henry Tehrani;"
		
		if save_EGG_TextGrids == 1
			selectObject: egg_dc1
			egg_dc1_resampled = Resample: 22050, 50
			selectObject: egg_dc
			egg_dc_resampled = Resample: 22050, 50
			# Merge them together!
			selectObject: egg_dc1_resampled, egg_dc_resampled
			egg_dc_combined = Combine to stereo
			selectObject: egg_dc_combined
			nowarn Save as binary file: output_wav_textgrid$ + speaker$ + "_" + filename$ + "_EGG-DC" + ".Sound"
			removeObject: egg_dc1_resampled
			removeObject: egg_dc_resampled
			removeObject: egg_dc_combined
		endif
		removeObject: egg_dc1
		removeObject: egg_dc
		removeObject: tg_HTtime
		removeObject: tg_HTCI
	endif


##### CALCULATION COMPLETE #####
	# Save TextGrid
	if save_EGG_TextGrids == 1
		if dEGG_method == 1
			selectObject: tg_degg
		elif hybrid_method == 1
			selectObject: tg_hybrid
		elif threshold_method == 1
			selectObject: tg_37threshold
		elif henry_Tehrani_method == 1
			selectObject: tg_HT
		endif
		if dEGG_method == 1
			plusObject: tg_degg
		endif
		if hybrid_method == 1
			plusObject: tg_hybrid
		endif
		if threshold_method == 1
			plusObject: tg_37threshold
		endif
		if henry_Tehrani_method == 1
			plusObject: tg_HT
		endif
		tg_tosave = Merge
		selectObject: tg_tosave
		nowarn Save as text file: output_wav_textgrid$ + speaker$ + "_" + filename$ + "_EGG" + ".TextGrid"
		removeObject: tg_tosave
	endif
	# Remove temporary objects
	if remove_temporary_objects == 1
		# Clean up textgrid objects
		removeObject: tg_pic
		removeObject: tg_pdc
		removeObject: tg_degg500
		if dEGG_method == 1
			removeObject: tg_degg
		endif
		if hybrid_method == 1
			removeObject: tg_25threshold, tg_hybrid
		endif
		if threshold_method == 1
			removeObject: tg_37threshold
		endif
		if henry_Tehrani_method == 1
			removeObject: tg_HT
		endif
		# Clean up sound objects
		removeObject: wavfile_raw
		removeObject: deriv
		removeObject: wavfile
		removeObject: stereo_combined
		removeObject: audio
		removeObject: audio_resampled
		removeObject: egg
		removeObject: eggfilter
		removeObject: deriv_resampled
		removeObject: wavfile_resampled
		removeObject: wavfile_smoothed
		removeObject: deriv_smoothed
		removeObject: dc_component
		removeObject: egg_raw
		removeObject: pitch
		if reference_tier > 0
			removeObject: textgrid_raw
		endif
	endif

	if manual_check == 1
	    ########################## [1] Compute y ranges ##########################
	    selectObject: analysis_matrix
		
	    # Add some padding
		#y_min1 = cq_min * 0.9
	    #y_max1 = cq_max * 1.1
		y_max1 = 1
		y_min1 = 0
	    y_min2 = f0_min * 0.9
	    y_max2 = f0_max * 1.1
	    
	    ########################## [2] Plot ##########################
	    Erase all
	    ###### Top plot: CQ and RET/PV ######
	    Select inner viewport: 1, 5, 1, 3.5
		Black
	    #Text top: "yes", "Plot of CQ and RET/PV"
		#Text left: "yes", "CQ and RET/PV"
		Text top: "yes", "Plot of CQ detected by EGGLab"
		Text left: "yes", "CQ"
	    Text bottom: "yes", "Period"
		line_max = index_max - 1
		pos1 = line_max*1/5
		pos2 = line_max*2/5
		pos3 = line_max*3/5
		pos4 = line_max*4/5
		ver_pos = 0.9

		if dEGG_method == 1
			Red
			selectObject: analysis_matrix
		    Scatter plot: 1, 2, 1, line_max, y_min1, y_max1, 1.2, "+", "no"
			#One mark top: pos1, "no", "no", "no", "+ CQ_d_E_G_G"
			Text: pos1, "left", ver_pos, "half", "+ CQ_d_E_G_G"
		endif
		if hybrid_method == 1
			Blue
			selectObject: analysis_matrix
		    Scatter plot: 1, 3, 1, line_max, y_min1, y_max1, 1.2, "o", "no"
			#One mark top: pos2, "no", "no", "no", "o CQ_H_o_w_a_r_d"
			Text: pos2, "left", ver_pos, "half", "o CQ_H_o_w_a_r_d"
		endif
		if threshold_method == 1
			Magenta
			selectObject: analysis_matrix
		    Scatter plot: 1, 4, 1, line_max, y_min1, y_max1, 1.2, "x", "no"
			#One mark top: pos3, "no", "no", "no", "x CQ_T_h_r_e_s_h_o_l_d"
			Text: pos3, "left", ver_pos, "half", "x CQ_T_h_r_e_s_h_o_l_d"
		endif
		if henry_Tehrani_method == 1
			Black
			selectObject: analysis_matrix
		    Scatter plot: 1, 4, 1, line_max, y_min1, y_max1, 1.2, ".", "no"
			#One mark top: pos4, "no", "no", "no", ". CQ_H_T"
			Text: pos4, "left", ver_pos, "half", ". CQ_H_T"
		endif
		Black
		Draw inner box
	    Marks left: 5, "yes", "yes", "no"
	    Marks bottom every: 1, 5, "yes", "yes", "no"
		One mark bottom: 1, "yes", "yes", "no", ""
		One mark bottom: line_max, "yes", "yes", "no", ""
	
	    ###### Bottom plot: F0 ######
	    Select inner viewport: 1, 5, 4.5, 7
		Black
	    Text top: "yes", "Plot of F0 detected by EGGLab"
		Text left: "yes", "F0 (Hz)"
	    Text bottom: "yes", "Period"
		pos1 = line_max*1/4
		pos2 = line_max*2/4
		pos3 = line_max*3/4
		ver_pos = (y_max2 - y_min2) * 0.9 + y_min2

		if dEGG_method == 1 || hybrid_method == 1 || henry_Tehrani_method == 1
		    Blue
			selectObject: analysis_matrix
		    Scatter plot: 1, 6, 1, line_max, y_min2, y_max2, 1.2, "x", "no"
			#One mark top: pos1, "no", "no", "no", "x F0_d_E_G_G"
			Text: pos1, "left", ver_pos, "half", "x F0_d_E_G_G"
		endif
		if threshold_method == 1
			Red
			selectObject: analysis_matrix
		    Scatter plot: 1, 7, 1, line_max, y_min2, y_max2, 1.2, "o", "no"
			#One mark top: pos2, "no", "no", "no", "o F0_T_h_r_e_s_h_o_l_d"
			Text: pos2, "left", ver_pos, "half", "o F0_T_h_r_e_s_h_o_l_d"
		endif
		Magenta
		selectObject: analysis_matrix
		Scatter plot: 1, 9, 1, line_max, y_min2, y_max2, 1.2, "+", "no"
		#One mark top: pos3, "no", "no", "no", "+ F0_T_h_r_e_s_h_o_l_d"
		Text: pos3, "left", ver_pos, "half", "o F0_P_r_a_a_t"

	    Black
		Draw inner box
	    Marks left: 5, "yes", "yes", "no"
	    Marks bottom every: 1, 5, "yes", "yes", "no"
		One mark bottom: 1, "yes", "yes", "no", ""
		One mark bottom: line_max, "yes", "yes", "no", ""

		if analyse_spectrum == 1
		    ###### H1*-H2* ######
		    Select inner viewport: 1, 5, 8, 10.5
			Black
		    Text top: "yes", "Plot of H1*-H2* detected by EGGLab"
			Text left: "yes", "H1*-H2* (dB)"
		    Text bottom: "yes", "Period"
			pos = line_max/2
			y_min3 = h1h2c_min * 0.9
			y_max3 = h1h2c_max * 1.1
	
		    Blue
			selectObject: analysis_matrix
		    Scatter plot: 1, 8, 1, line_max, y_min3, y_max3, 1.2, "o", "yes"
		    Black
		    Marks left: 5, "yes", "yes", "no"
		    Marks bottom every: 1, 5, "yes", "yes", "no"
		endif

		removeObject: analysis_matrix

		beginPause: "Please check CQ, F0, H1*-H2* (if any) detected by EGGLab."
			comment: "1. Check whether the parameters are correct."
		    comment: "2. Then click 'Continue'."
		endPause: "Continue", 1
	endif

endfor

removeObject: filestrings

execution_time = stopwatch
minutes = execution_time div 60
seconds = execution_time mod 60
if minutes > 1
	min$ = " mins "
else
	min$ = " min "
endif
if seconds > 1
	sec$ = " secs"
else
	sec$ = " sec"
endif
if save_EGG_TextGrids == 1
	# Completion message
	appendInfoLine: newline$, newline$, "=== PROCESSING COMPLETE ===", newline$,
	... "Complete at: ", date$(), newline$,
	... "Execution time: ", fixed$(minutes,0), min$, fixed$(seconds,0), sec$, newline$,
	... "EGG parameters extraction successful for all files!", newline$,
	... "Data saved to: ", directory$ + "output/", newline$,
	... "WAV and TextGrid for checking saved to: ", output_wav_textgrid$, newline$,
	... "Total files processed: ", filenumbers, newline$,
	... newline$,
	... "Note:", newline$,
	... "Parameters in the TextGrid files are shown in this order:", newline$,
	... "  CQ, f0 (Praat-tracked), SQ, and f0 (1/EGG-cycle).", newline$,
	... "==========================="
else
	# Completion message
	appendInfoLine: newline$, newline$, "=== PROCESSING COMPLETE ===", newline$,
	... "Complete at: ", date$(), newline$,
	... "Execution time: ", fixed$(minutes,0), min$, fixed$(seconds,0), sec$, newline$,
	... "EGG parameters extraction successful for all files!", newline$,
	... "Data saved to: ", directory$ + "output/", newline$,
	... "Total files processed: ", filenumbers, newline$,
	... "==========================="
endif
if completion_alert = 1
	# Optional sound notification (uncomment if desired)
	if windows
	    runSystem: "powershell [console]::beep(1000,300)"
	elsif macintosh
	    runSystem: "afplay /System/Library/Sounds/Glass.aiff"
	endif
endif
