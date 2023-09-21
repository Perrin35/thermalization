OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83672268) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(-0.81122938) q[0];
rz(-3.0468416) q[2];
sx q[2];
rz(-1.0031327) q[2];
sx q[2];
rz(1.809158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58798446) q[1];
sx q[1];
rz(-1.6349287) q[1];
sx q[1];
rz(1.9967805) q[1];
rz(-pi) q[2];
rz(0.26478404) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(0.36987723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-3.1100173) q[2];
rz(1.8850373) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(-2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8388222) q[0];
sx q[0];
rz(-1.4571723) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(-2.6020715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3300433) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(3.0898068) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9794481) q[2];
sx q[2];
rz(-0.89102972) q[2];
sx q[2];
rz(-2.0603927) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62944618) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(-2.9260103) q[1];
x q[2];
rz(2.5856421) q[3];
sx q[3];
rz(-2.2009938) q[3];
sx q[3];
rz(-2.6549908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-0.24307069) q[2];
rz(-2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3617525) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(2.6932122) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6305144) q[0];
sx q[0];
rz(-0.96195463) q[0];
sx q[0];
rz(-1.2468673) q[0];
x q[1];
rz(2.0492378) q[2];
sx q[2];
rz(-1.97176) q[2];
sx q[2];
rz(-1.5329597) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33830723) q[1];
sx q[1];
rz(-2.3803664) q[1];
sx q[1];
rz(-2.0558946) q[1];
rz(-1.6085298) q[3];
sx q[3];
rz(-1.747526) q[3];
sx q[3];
rz(2.7023774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-2.4352303) q[2];
sx q[2];
rz(-2.5668872) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213585) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(0.25948778) q[0];
rz(-1.150594) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(0.73192275) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61011945) q[0];
sx q[0];
rz(-1.1024794) q[0];
sx q[0];
rz(2.5584695) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.133693) q[2];
sx q[2];
rz(-2.4198654) q[2];
sx q[2];
rz(2.0267817) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5922028) q[1];
sx q[1];
rz(-1.376774) q[1];
sx q[1];
rz(0.26280304) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51283522) q[3];
sx q[3];
rz(-2.8898015) q[3];
sx q[3];
rz(0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-2.990492) q[2];
rz(-0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(-2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5916409) q[0];
sx q[0];
rz(-2.0786091) q[0];
sx q[0];
rz(1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(0.79777065) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3691468) q[0];
sx q[0];
rz(-1.4506842) q[0];
sx q[0];
rz(0.0025047501) q[0];
rz(-pi) q[1];
rz(-1.351379) q[2];
sx q[2];
rz(-1.2741538) q[2];
sx q[2];
rz(-2.9079633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7792015) q[1];
sx q[1];
rz(-2.5070842) q[1];
sx q[1];
rz(1.4125376) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5161414) q[3];
sx q[3];
rz(-2.6126385) q[3];
sx q[3];
rz(-1.4048502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(2.8395555) q[2];
rz(1.9942412) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(1.1557895) q[0];
rz(2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(0.070080431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72567155) q[0];
sx q[0];
rz(-1.7416746) q[0];
sx q[0];
rz(2.9861949) q[0];
x q[1];
rz(2.5885133) q[2];
sx q[2];
rz(-2.2773909) q[2];
sx q[2];
rz(-1.7318219) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97092123) q[1];
sx q[1];
rz(-1.0075924) q[1];
sx q[1];
rz(1.7621653) q[1];
rz(-2.0134986) q[3];
sx q[3];
rz(-3.0590995) q[3];
sx q[3];
rz(-3.1174297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(0.62136674) q[2];
rz(1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086534111) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(2.281718) q[0];
rz(-1.2043918) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(3.133657) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746349) q[0];
sx q[0];
rz(-1.8263706) q[0];
sx q[0];
rz(-2.5711683) q[0];
x q[1];
rz(-0.87848778) q[2];
sx q[2];
rz(-1.9108859) q[2];
sx q[2];
rz(1.0483339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.24212813) q[1];
sx q[1];
rz(-2.2367396) q[1];
sx q[1];
rz(0.66023402) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2357335) q[3];
sx q[3];
rz(-0.48478904) q[3];
sx q[3];
rz(-2.0743899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(2.725214) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(-0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96173441) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(0.36488786) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-2.6049728) q[1];
sx q[1];
rz(1.6392802) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84425981) q[0];
sx q[0];
rz(-1.5808006) q[0];
sx q[0];
rz(-2.8914333) q[0];
rz(-pi) q[1];
x q[1];
rz(2.084923) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(-2.4664997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8870526) q[1];
sx q[1];
rz(-2.2215543) q[1];
sx q[1];
rz(-0.73144967) q[1];
rz(-pi) q[2];
rz(1.3650465) q[3];
sx q[3];
rz(-1.660941) q[3];
sx q[3];
rz(0.84116018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(1.0650744) q[2];
rz(-2.8403357) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97312462) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(-2.705943) q[0];
rz(-1.3849974) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(-0.41697821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5676253) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(-2.2249939) q[0];
x q[1];
rz(-0.99879361) q[2];
sx q[2];
rz(-2.1527094) q[2];
sx q[2];
rz(-1.7172161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7855362) q[1];
sx q[1];
rz(-1.0724663) q[1];
sx q[1];
rz(-2.9836125) q[1];
rz(-pi) q[2];
rz(-1.4407518) q[3];
sx q[3];
rz(-1.4308883) q[3];
sx q[3];
rz(-3.1062982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(-2.9337692) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41480961) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(-1.6171932) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029862558) q[0];
sx q[0];
rz(-2.7435281) q[0];
sx q[0];
rz(1.8882621) q[0];
rz(-1.1216713) q[2];
sx q[2];
rz(-1.3438864) q[2];
sx q[2];
rz(-2.801148) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10510124) q[1];
sx q[1];
rz(-2.0836012) q[1];
sx q[1];
rz(2.3829616) q[1];
rz(-pi) q[2];
rz(0.15354746) q[3];
sx q[3];
rz(-2.2205177) q[3];
sx q[3];
rz(-0.49680199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(-2.1255169) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(-1.7779508) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(0.50921847) q[2];
sx q[2];
rz(-1.5794532) q[2];
sx q[2];
rz(3.0184359) q[2];
rz(-0.090311269) q[3];
sx q[3];
rz(-1.8009381) q[3];
sx q[3];
rz(1.2931852) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];