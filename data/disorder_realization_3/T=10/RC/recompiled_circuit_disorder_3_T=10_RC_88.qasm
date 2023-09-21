OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(-3.1103818) q[0];
sx q[0];
rz(-2.6565235) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038923351) q[0];
sx q[0];
rz(-1.9063213) q[0];
sx q[0];
rz(-1.194792) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4235052) q[2];
sx q[2];
rz(-0.5746595) q[2];
sx q[2];
rz(-1.5073843) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1297076) q[1];
sx q[1];
rz(-1.9958479) q[1];
sx q[1];
rz(0.070406291) q[1];
x q[2];
rz(0.36177735) q[3];
sx q[3];
rz(-0.28218111) q[3];
sx q[3];
rz(-0.85229814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9238613) q[2];
sx q[2];
rz(-1.2552746) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3027705) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(-2.6020715) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44885143) q[0];
sx q[0];
rz(-2.6896547) q[0];
sx q[0];
rz(-1.6777722) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4194854) q[2];
sx q[2];
rz(-1.8849843) q[2];
sx q[2];
rz(2.9177641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62944618) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(-2.9260103) q[1];
rz(-pi) q[2];
rz(0.94445618) q[3];
sx q[3];
rz(-0.8144905) q[3];
sx q[3];
rz(-0.32523793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(-2.4754751) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(0.4483805) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(2.8853436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8921593) q[0];
sx q[0];
rz(-1.306635) q[0];
sx q[0];
rz(0.63412068) q[0];
x q[1];
rz(-2.3149812) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(2.5342864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33830723) q[1];
sx q[1];
rz(-2.3803664) q[1];
sx q[1];
rz(-2.0558946) q[1];
x q[2];
rz(2.9647397) q[3];
sx q[3];
rz(-1.5336509) q[3];
sx q[3];
rz(1.1249441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8024575) q[2];
sx q[2];
rz(-0.70636237) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67128348) q[0];
sx q[0];
rz(-2.0844315) q[0];
sx q[0];
rz(1.0259823) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6875793) q[2];
sx q[2];
rz(-2.2846966) q[2];
sx q[2];
rz(-1.2920979) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5922028) q[1];
sx q[1];
rz(-1.376774) q[1];
sx q[1];
rz(-2.8787896) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9210806) q[3];
sx q[3];
rz(-1.6933428) q[3];
sx q[3];
rz(-1.8054655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(-2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.939643) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(1.4506838) q[0];
rz(-1.7902137) q[2];
sx q[2];
rz(-1.8674388) q[2];
sx q[2];
rz(0.23362939) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0610173) q[1];
sx q[1];
rz(-1.6643545) q[1];
sx q[1];
rz(2.1993125) q[1];
rz(-2.5161414) q[3];
sx q[3];
rz(-2.6126385) q[3];
sx q[3];
rz(-1.4048502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(0.30203715) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.4619504) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7323332) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(1.0844768) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(-0.070080431) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2698343) q[0];
sx q[0];
rz(-1.7239128) q[0];
sx q[0];
rz(-1.7437177) q[0];
x q[1];
rz(-2.5885133) q[2];
sx q[2];
rz(-0.86420176) q[2];
sx q[2];
rz(1.4097708) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8225704) q[1];
sx q[1];
rz(-0.59148568) q[1];
sx q[1];
rz(-2.849008) q[1];
rz(-pi) q[2];
rz(-3.1061884) q[3];
sx q[3];
rz(-1.4962713) q[3];
sx q[3];
rz(-2.7217334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.30248102) q[2];
sx q[2];
rz(-1.182686) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(1.7012043) q[3];
sx q[3];
rz(-0.50783235) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086534111) q[0];
sx q[0];
rz(-1.7250412) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(3.133657) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7746349) q[0];
sx q[0];
rz(-1.3152221) q[0];
sx q[0];
rz(0.5704244) q[0];
rz(0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(1.0483339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.36564) q[1];
sx q[1];
rz(-2.0740293) q[1];
sx q[1];
rz(2.3535437) q[1];
rz(-pi) q[2];
rz(-2.2357335) q[3];
sx q[3];
rz(-0.48478904) q[3];
sx q[3];
rz(2.0743899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69592151) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(-2.725214) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-2.8680153) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(-0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.6392802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72909268) q[0];
sx q[0];
rz(-1.3206498) q[0];
sx q[0];
rz(1.581122) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6697568) q[2];
sx q[2];
rz(-2.2484357) q[2];
sx q[2];
rz(-0.0080646947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9598436) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(0.77397857) q[1];
x q[2];
rz(-1.9873398) q[3];
sx q[3];
rz(-0.22437469) q[3];
sx q[3];
rz(0.32240543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6442948) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(1.0650744) q[2];
rz(2.8403357) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97312462) q[0];
sx q[0];
rz(-0.081806101) q[0];
sx q[0];
rz(-0.43564963) q[0];
rz(1.3849974) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(2.7246144) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57396736) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(-0.91659878) q[0];
rz(2.4531104) q[2];
sx q[2];
rz(-0.79198972) q[2];
sx q[2];
rz(-2.288523) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6777842) q[1];
sx q[1];
rz(-0.52075547) q[1];
sx q[1];
rz(-1.8522472) q[1];
rz(-0.74420332) q[3];
sx q[3];
rz(-0.19072285) q[3];
sx q[3];
rz(-0.71803367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(-2.9337692) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.726783) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(-1.5243994) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(-1.3226002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1117301) q[0];
sx q[0];
rz(-0.39806453) q[0];
sx q[0];
rz(1.8882621) q[0];
x q[1];
rz(1.1216713) q[2];
sx q[2];
rz(-1.3438864) q[2];
sx q[2];
rz(-0.34044468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0304071) q[1];
sx q[1];
rz(-2.2135995) q[1];
sx q[1];
rz(-0.91099693) q[1];
rz(-pi) q[2];
rz(-1.3721458) q[3];
sx q[3];
rz(-0.66505265) q[3];
sx q[3];
rz(-0.24634758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3502729) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-2.1255169) q[2];
rz(1.2223876) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(-0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9983457) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.7779508) q[1];
sx q[1];
rz(-1.2294055) q[1];
sx q[1];
rz(1.3235863) q[1];
rz(-1.5608816) q[2];
sx q[2];
rz(-1.0615988) q[2];
sx q[2];
rz(1.4524729) q[2];
rz(1.2033403) q[3];
sx q[3];
rz(-0.24693476) q[3];
sx q[3];
rz(0.91528391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
