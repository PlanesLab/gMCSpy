\ Model MCS
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  zp_g5_ag_on_redflux54 + zp_g8_ag_on_redflux56 + zp_g3_ag_on_redflux58
   + zp_g7_ag_on_redflux60 + zp_g2_ag_on_redflux62 + zp_g4_ag_on_redflux64
   + zp_g6_ag_on_redflux66 + 0 zw_67
Subject To
 Ndual_constraint0: u_0 - u_1 >= 0
 Ndual_constraint1: u_0 - u_2 >= 0
 Ndual_constraint2: - u_1 + u_3 >= 0
 Ndual_constraint3: u_3 - u_4 >= 0
 Ndual_constraint4: - u_0 - u_3 - u_5 + u_6 >= 0
 Ndual_constraint5: u_7 - u_8 >= 0
 Ndual_constraint6: u_8 - u_9 - u_10 >= 0
 Ndual_constraint7: u_7 - u_11 >= 0
 Ndual_constraint8: - u_9 + u_11 - u_12 >= 0
 Ndual_constraint9: u_7 - u_13 >= 0
 Ndual_constraint10: u_5 - u_14 >= 0
 Ndual_constraint11: u_5 - u_15 >= 0
 Ndual_constraint12: u_13 - u_16 - u_17 >= 0
 Ndual_constraint13: u_14 - u_18 >= 0
 Ndual_constraint14: u_16 - u_19 >= 0
 Ndual_constraint15: u_2 - u_20 >= 0
 Ndual_constraint16: u_2 - u_21 >= 0
 Ndual_constraint17: u_10 - u_22 - u_23 >= 0
 Ndual_constraint18: u_22 - u_24 >= 0
 Ndual_constraint19: u_20 - u_25 >= 0
 Ndual_constraint20: u_1 - u_26 >= 0
 Ndual_constraint21: - u_5 - u_22 + u_26 >= 0
 Ndual_constraint22: u_1 - u_27 >= 0
 Ndual_constraint23: - u_13 + u_28 >= 0
 Ndual_constraint24: - u_20 + u_28 >= 0
 Ndual_constraint25: u_9 - u_28 - u_29 >= 0
 Ndual_constraint26: u_4 - u_30 >= 0
 Ndual_constraint27: - u_22 + u_30 - u_31 >= 0
 Ndual_constraint28: u_4 - u_32 >= 0
 Ndual_constraint29: u_33 - u_34 >= 0
 Ndual_constraint30: - u_20 + u_33 >= 0
 Ndual_constraint31: u_12 - u_33 - u_35 >= 0
 Ndual_constraint32: u_31 - u_36 >= 0
 Ndual_constraint33: u_34 - u_37 >= 0
 Ndual_constraint34: u_15 + vp_38 >= 0
 Ndual_constraint35: u_17 + vp_39 >= 0
 Ndual_constraint36: u_18 + vp_40 >= 0
 Ndual_constraint37: u_18 >= 0
 Ndual_constraint38: u_19 + vp_41 >= 0
 Ndual_constraint39: u_21 + vp_42 >= 0
 Ndual_constraint40: u_23 + vp_43 >= 0
 Ndual_constraint41: u_24 + vp_44 >= 0
 Ndual_constraint42: u_24 >= 0
 Ndual_constraint43: u_25 + vp_45 >= 0
 Ndual_constraint44: u_27 + vp_46 >= 0
 Ndual_constraint45: u_29 + vp_47 >= 0
 Ndual_constraint46: u_32 + vp_48 >= 0
 Ndual_constraint47: u_35 + vp_49 >= 0
 Ndual_constraint48: u_36 + vp_50 >= 0
 Ndual_constraint49: u_36 >= 0
 Ndual_constraint50: u_37 + vp_51 >= 0
 Ndual_constraint51: - u_7 - w_52 >= 0
 Ndual_constraint52: - u_6 >= 0
 forceBioCons: - 0.001 w_52 <= -0.001
 forceAtLeastOne: zp_g5_ag_on_redflux54 + zp_g8_ag_on_redflux56
   + zp_g3_ag_on_redflux58 + zp_g7_ag_on_redflux60 + zp_g2_ag_on_redflux62
   + zp_g4_ag_on_redflux64 + zp_g6_ag_on_redflux66 >= 1
 relationship_g5: zp_g5_ag_off_redflux53 + zp_g5_ag_on_redflux54 = 1
 relationship_g3: zp_g3_ag_off_redflux57 + zp_g3_ag_on_redflux58 = 1
 relationship_g2: zp_g2_ag_off_redflux61 + zp_g2_ag_on_redflux62 = 1
 relationship_g7: zp_g7_ag_off_redflux59 + zp_g7_ag_on_redflux60 = 1
 relationship_g8: zp_g8_ag_off_redflux55 + zp_g8_ag_on_redflux56 = 1
 relationship_g6: zp_g6_ag_off_redflux65 + zp_g6_ag_on_redflux66 = 1
 relationship_g4: zp_g4_ag_off_redflux63 + zp_g4_ag_on_redflux64 = 1
 indicatorZOn_0: zp_g5_ag_off_redflux53 = 1 -> vp_38 >= 1
 indicatorZOn_1: zp_g5_ag_on_redflux54 = 1 -> vp_39 >= 1
 indicatorZOn_2: zp_g8_ag_off_redflux55 = 1 -> vp_40 >= 1
 indicatorZOn_3: zp_g8_ag_on_redflux56 = 1 -> vp_41 >= 1
 indicatorZOn_4: zp_g3_ag_off_redflux57 = 1 -> vp_42 >= 1
 indicatorZOn_5: zp_g3_ag_on_redflux58 = 1 -> vp_43 >= 1
 indicatorZOn_6: zp_g7_ag_off_redflux59 = 1 -> vp_44 >= 1
 indicatorZOn_7: zp_g7_ag_on_redflux60 = 1 -> vp_45 >= 1
 indicatorZOn_8: zp_g2_ag_off_redflux61 = 1 -> vp_46 >= 1
 indicatorZOn_9: zp_g2_ag_on_redflux62 = 1 -> vp_47 >= 1
 indicatorZOn_10: zp_g4_ag_off_redflux63 = 1 -> vp_48 >= 1
 indicatorZOn_11: zp_g4_ag_on_redflux64 = 1 -> vp_49 >= 1
 indicatorZOn_12: zp_g6_ag_off_redflux65 = 1 -> vp_50 >= 1
 indicatorZOn_13: zp_g6_ag_on_redflux66 = 1 -> vp_51 >= 1
 indicatorZOn_14: zw_67 = 1 -> w_52 >= 1
 indicatorZOff_0: zp_g5_ag_off_redflux53 = 0 -> vp_38 <= 0
 indicatorZOff_1: zp_g5_ag_on_redflux54 = 0 -> vp_39 <= 0
 indicatorZOff_2: zp_g8_ag_off_redflux55 = 0 -> vp_40 <= 0
 indicatorZOff_3: zp_g8_ag_on_redflux56 = 0 -> vp_41 <= 0
 indicatorZOff_4: zp_g3_ag_off_redflux57 = 0 -> vp_42 <= 0
 indicatorZOff_5: zp_g3_ag_on_redflux58 = 0 -> vp_43 <= 0
 indicatorZOff_6: zp_g7_ag_off_redflux59 = 0 -> vp_44 <= 0
 indicatorZOff_7: zp_g7_ag_on_redflux60 = 0 -> vp_45 <= 0
 indicatorZOff_8: zp_g2_ag_off_redflux61 = 0 -> vp_46 <= 0
 indicatorZOff_9: zp_g2_ag_on_redflux62 = 0 -> vp_47 <= 0
 indicatorZOff_10: zp_g4_ag_off_redflux63 = 0 -> vp_48 <= 0
 indicatorZOff_11: zp_g4_ag_on_redflux64 = 0 -> vp_49 <= 0
 indicatorZOff_12: zp_g6_ag_off_redflux65 = 0 -> vp_50 <= 0
 indicatorZOff_13: zp_g6_ag_on_redflux66 = 0 -> vp_51 <= 0
 indicatorZOff_14: zw_67 = 0 -> w_52 <= 0
Bounds
 u_0 free
 u_1 free
 u_2 free
 u_3 free
 u_4 free
 u_5 free
 u_6 free
 u_7 free
 u_8 free
 u_9 free
 u_10 free
 u_11 free
 u_12 free
 u_13 free
 u_14 free
 u_15 free
 u_16 free
 u_17 free
 u_18 free
 u_19 free
 u_20 free
 u_21 free
 u_22 free
 u_23 free
 u_24 free
 u_25 free
 u_26 free
 u_27 free
 u_28 free
 u_29 free
 u_30 free
 u_31 free
 u_32 free
 u_33 free
 u_34 free
 u_35 free
 u_36 free
 u_37 free
Binaries
 zp_g5_ag_off_redflux53 zp_g5_ag_on_redflux54 zp_g8_ag_off_redflux55
 zp_g8_ag_on_redflux56 zp_g3_ag_off_redflux57 zp_g3_ag_on_redflux58
 zp_g7_ag_off_redflux59 zp_g7_ag_on_redflux60 zp_g2_ag_off_redflux61
 zp_g2_ag_on_redflux62 zp_g4_ag_off_redflux63 zp_g4_ag_on_redflux64
 zp_g6_ag_off_redflux65 zp_g6_ag_on_redflux66 zw_67
End
