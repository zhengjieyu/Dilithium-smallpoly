.align 2
.global copy_64x8
.type copy_64x8, %function
copy_64x8:
    array1 .req x0
    array2  .req x1

    ld1 {V0.2D, V1.2D, V2.2D, V3.2D}, [array2]
    st1 {V0.2D, V1.2D, V2.2D, V3.2D}, [array1]

    .unreq array1
    .unreq array2

    ret

.align 2
.global arrayAcc_arm
.type arrayAcc_arm, %function
arrayAcc_arm:
    array1 .req x0
    array2  .req x1
    counter .req x2
    
    mov counter, #32 
    
    acc_loop:
        ld1 {V0.2D, V1.2D, V2.2D, V3.2D}, [array1]
        ld1 {V4.2D, V5.2D, V6.2D, V7.2D}, [array2], #64

        add V0.2D, V0.2D, V4.2D
        add V1.2D, V1.2D, V5.2D
        add V2.2D, V2.2D, V6.2D
        add V3.2D, V3.2D, V7.2D

        st1 {V0.2D, V1.2D, V2.2D, V3.2D}, [array1], #64
        sub counter, counter, #1
        cbnz counter, acc_loop

    .unreq array1
    .unreq array2
    .unreq counter

    ret

.align 2
.global cal_table
.type cal_table, %function
cal_table:
    coeffs .req x0
    add_num .req w1
    table .req x2
    offset .req x3

    ld1 {v1.2s, v2.2s, v3.2s, v4.2s}, [coeffs]
    dup v0.2s, add_num
    saddl v1.2d, v1.2s, v0.2s
    saddl v2.2d, v2.2s, v0.2s
    saddl v3.2d, v3.2s, v0.2s
    saddl v4.2d, v4.2s, v0.2s

    //se_table[k+N] = (se_table[k+N]<<8) | (temp);
    ld1 {v5.2d, v6.2d, v7.2d, v8.2d}, [table]
    dup v0.2d, offset
    ushl v5.2d, v5.2d, v0.2d
    ushl v6.2d, v6.2d, v0.2d
    ushl v7.2d, v7.2d, v0.2d
    ushl v8.2d, v8.2d, v0.2d

    orr v5.16b, v5.16b, v1.16b
    orr v6.16b, v6.16b, v2.16b
    orr v7.16b, v7.16b, v3.16b
    orr v8.16b, v8.16b, v4.16b

    st1 {v5.2d, v6.2d, v7.2d, v8.2d}, [table]

    .unreq coeffs
    .unreq add_num
    .unreq table
    .unreq offset

    ret

.align 2
.global st_table
.type st_table, %function
st_table:
    mask .req x0
    table .req x1

    dup v0.2d, mask

    ld1 {v5.2d, v6.2d, v7.2d, v8.2d}, [table]
    sub v5.2d, v0.2d, v5.2d
    sub v6.2d, v0.2d, v6.2d
    sub v7.2d, v0.2d, v7.2d
    sub v8.2d, v0.2d, v8.2d

    sub table, table, #256*8
    st1 {v5.2d, v6.2d, v7.2d, v8.2d}, [table]
    add table, table, #256*8*2
    st1 {v5.2d, v6.2d, v7.2d, v8.2d}, [table]

    .unreq mask
    .unreq table

    ret


.macro recover_evaluate temp, coeffs, and_num, sub_num, offset
    ld1 {v1.2d, v2.2d, v3.2d, v4.2d}, [\temp]
    dup v0.2d, \and_num
    and v5.16b, v1.16b, v0.16b
    and v6.16b, v2.16b, v0.16b
    and v7.16b, v3.16b, v0.16b
    and v8.16b, v4.16b, v0.16b

    dup v0.2d, \sub_num
    sub v5.2d, v5.2d, v0.2d
    sub v6.2d, v6.2d, v0.2d
    sub v7.2d, v7.2d, v0.2d
    sub v8.2d, v8.2d, v0.2d
    xtn v9.2s, v5.2d
    xtn v10.2s, v6.2d
    xtn v11.2s, v7.2d
    xtn v12.2s, v8.2d

    ushr v1.2d, v1.2d, \offset
    ushr v2.2d, v2.2d, \offset
    ushr v3.2d, v3.2d, \offset
    ushr v4.2d, v4.2d, \offset

    st1 {v1.2d, v2.2d, v3.2d, v4.2d}, [\temp]
    st1 {v9.2s, v10.2s, v11.2s, v12.2s}, [\coeffs]
.endm

.align 2
.global recover_cs_arm
.type recover_cs_arm, %function
recover_cs_arm:
    temp .req x0
    cs_coeffs .req x1
    and_num .req x2
    sub_num .req x3

    recover_evaluate temp, cs_coeffs, and_num, sub_num, #8

    .unreq temp
    .unreq cs_coeffs
    .unreq and_num
    .unreq sub_num

    ret

.align 2
.global evaluate_ct0_arm
.type evaluate_ct0_arm, %function
evaluate_ct0_arm:
    temp .req x0
    target_coeffs .req x1
    and_num .req x2
    sub_num .req x3

    recover_evaluate temp, target_coeffs, and_num, sub_num, #19

    .unreq temp
    .unreq target_coeffs
    .unreq and_num
    .unreq sub_num

    ret


.macro reduce32 a, t, add_num, mul_num, shr_num
    add \t, \a, \add_num
    sshr \t, \t, \shr_num
    mla \a, \t, \mul_num
.endm

.macro recover_evaluate_ct1 temp, coeffs, and_num, sub_num, offset
    ld1 {v1.2d, v2.2d, v3.2d, v4.2d}, [\temp]
    dup v0.2d, \and_num
    and v5.16b, v1.16b, v0.16b
    and v6.16b, v2.16b, v0.16b
    and v7.16b, v3.16b, v0.16b
    and v8.16b, v4.16b, v0.16b

    dup v0.2d, \sub_num
    sub v5.2d, v5.2d, v0.2d
    sub v6.2d, v6.2d, v0.2d
    sub v7.2d, v7.2d, v0.2d
    sub v8.2d, v8.2d, v0.2d
    xtn v9.2s, v5.2d
    xtn v10.2s, v6.2d
    xtn v11.2s, v7.2d
    xtn v12.2s, v8.2d


    ushr v1.2d, v1.2d, \offset
    ushr v2.2d, v2.2d, \offset
    ushr v3.2d, v3.2d, \offset
    ushr v4.2d, v4.2d, \offset


    shl v9.2s, v9.2s, #13
    shl v10.2s, v10.2s, #13
    shl v11.2s, v11.2s, #13
    shl v12.2s, v12.2s, #13

    reduce32 v9.2s, v13.2s, v14.2s, v15.2s, #23
    reduce32 v10.2s, v13.2s, v14.2s, v15.2s, #23
    reduce32 v11.2s, v13.2s, v14.2s, v15.2s, #23
    reduce32 v12.2s, v13.2s, v14.2s, v15.2s, #23

    st1 {v1.2d, v2.2d, v3.2d, v4.2d}, [\temp]
    st1 {v9.2s, v10.2s, v11.2s, v12.2s}, [\coeffs]
.endm

.align 2
.global evaluate_ct1_arm
.type evaluate_ct1_arm, %function
evaluate_ct1_arm:
    temp .req x0
    target_coeffs .req x1
    and_num .req x2
    sub_num .req x3
    negative_Q .req w4 

    two_22 .req w5
    
    mov two_22, #4194304 
    dup v14.2s, two_22
    dup v15.2s, negative_Q 

    recover_evaluate_ct1 temp, target_coeffs, and_num, sub_num, #17
    
    .unreq temp
    .unreq target_coeffs
    .unreq and_num
    .unreq sub_num
    .unreq negative_Q
    .unreq two_22

    ret