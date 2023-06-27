	.file	"test.cpp"
	.text
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"int func1(int)"
.LC1:
	.string	"test.cpp"
.LC2:
	.string	"b > 0"
	.text
	.globl	_Z5func1i
	.type	_Z5func1i, @function
_Z5func1i:
.LFB0:
	.cfi_startproc
	endbr64
	movl	%edi, %eax
	imull	%edi, %eax
	imull	%edi, %eax
	testl	%eax, %eax
	jle	.L6
	ret
.L6:
	subq	$8, %rsp
	.cfi_def_cfa_offset 16
	leaq	.LC0(%rip), %rcx
	movl	$9, %edx
	leaq	.LC1(%rip), %rsi
	leaq	.LC2(%rip), %rdi
	call	__assert_fail@PLT
	.cfi_endproc
.LFE0:
	.size	_Z5func1i, .-_Z5func1i
	.globl	_Z5func2i
	.type	_Z5func2i, @function
_Z5func2i:
.LFB1:
	.cfi_startproc
	endbr64
	movl	%edi, %eax
	imull	%edi, %eax
	imull	%edi, %eax
	ret
	.cfi_endproc
.LFE1:
	.size	_Z5func2i, .-_Z5func2i
	.ident	"GCC: (Ubuntu 11.3.0-1ubuntu1~22.04.1) 11.3.0"
	.section	.note.GNU-stack,"",@progbits
	.section	.note.gnu.property,"a"
	.align 8
	.long	1f - 0f
	.long	4f - 1f
	.long	5
0:
	.string	"GNU"
1:
	.align 8
	.long	0xc0000002
	.long	3f - 2f
2:
	.long	0x3
3:
	.align 8
4:
