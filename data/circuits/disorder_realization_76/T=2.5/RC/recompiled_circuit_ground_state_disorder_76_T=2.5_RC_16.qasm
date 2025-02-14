OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.051209) q[0];
sx q[0];
rz(-0.78855711) q[0];
sx q[0];
rz(2.8330044) q[0];
rz(-2.7170972) q[1];
sx q[1];
rz(-1.0841882) q[1];
sx q[1];
rz(-2.7660313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.080506026) q[0];
sx q[0];
rz(-1.2895786) q[0];
sx q[0];
rz(3.0628231) q[0];
x q[1];
rz(2.8396222) q[2];
sx q[2];
rz(-1.540138) q[2];
sx q[2];
rz(-1.4238525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8164067) q[1];
sx q[1];
rz(-0.46811399) q[1];
sx q[1];
rz(2.8312842) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1377596) q[3];
sx q[3];
rz(-1.4272913) q[3];
sx q[3];
rz(1.4583064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3382807) q[2];
sx q[2];
rz(-0.25140005) q[2];
sx q[2];
rz(2.5982762) q[2];
rz(1.733689) q[3];
sx q[3];
rz(-1.7270154) q[3];
sx q[3];
rz(2.2442815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5498085) q[0];
sx q[0];
rz(-1.1487288) q[0];
sx q[0];
rz(2.5445004) q[0];
rz(1.0626556) q[1];
sx q[1];
rz(-0.63908827) q[1];
sx q[1];
rz(2.7426681) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9358109) q[0];
sx q[0];
rz(-1.1880186) q[0];
sx q[0];
rz(-2.1924928) q[0];
rz(-pi) q[1];
x q[1];
rz(2.127803) q[2];
sx q[2];
rz(-1.5502056) q[2];
sx q[2];
rz(-2.2228732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0860721) q[1];
sx q[1];
rz(-2.2694553) q[1];
sx q[1];
rz(-2.3585206) q[1];
rz(-pi) q[2];
rz(-2.2406897) q[3];
sx q[3];
rz(-0.76053166) q[3];
sx q[3];
rz(2.4128715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4661633) q[2];
sx q[2];
rz(-1.9025981) q[2];
sx q[2];
rz(-2.0002401) q[2];
rz(2.2232248) q[3];
sx q[3];
rz(-2.4310302) q[3];
sx q[3];
rz(0.58951283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68701768) q[0];
sx q[0];
rz(-0.37746754) q[0];
sx q[0];
rz(-2.9662509) q[0];
rz(-1.31458) q[1];
sx q[1];
rz(-1.8914696) q[1];
sx q[1];
rz(-0.66741991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9475896) q[0];
sx q[0];
rz(-2.0663446) q[0];
sx q[0];
rz(0.13107863) q[0];
x q[1];
rz(1.0834789) q[2];
sx q[2];
rz(-1.0333152) q[2];
sx q[2];
rz(-1.4097241) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2167425) q[1];
sx q[1];
rz(-1.5702254) q[1];
sx q[1];
rz(0.3158098) q[1];
x q[2];
rz(2.5273058) q[3];
sx q[3];
rz(-1.7899872) q[3];
sx q[3];
rz(-3.1295071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0257001) q[2];
sx q[2];
rz(-0.95558715) q[2];
sx q[2];
rz(-1.0828177) q[2];
rz(-1.8358021) q[3];
sx q[3];
rz(-2.3256153) q[3];
sx q[3];
rz(-2.344237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0037271518) q[0];
sx q[0];
rz(-2.9481413) q[0];
sx q[0];
rz(0.61099148) q[0];
rz(-0.56426471) q[1];
sx q[1];
rz(-2.1128426) q[1];
sx q[1];
rz(-2.4345523) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5093197) q[0];
sx q[0];
rz(-1.600926) q[0];
sx q[0];
rz(-0.79674652) q[0];
x q[1];
rz(-1.9025175) q[2];
sx q[2];
rz(-1.9608856) q[2];
sx q[2];
rz(1.6029101) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1752211) q[1];
sx q[1];
rz(-0.67357682) q[1];
sx q[1];
rz(2.3994956) q[1];
x q[2];
rz(-0.50102542) q[3];
sx q[3];
rz(-1.3210557) q[3];
sx q[3];
rz(1.4755032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6881312) q[2];
sx q[2];
rz(-0.42669272) q[2];
sx q[2];
rz(2.1235662) q[2];
rz(-0.38797837) q[3];
sx q[3];
rz(-1.4464902) q[3];
sx q[3];
rz(-0.33897266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729048) q[0];
sx q[0];
rz(-2.0141116) q[0];
sx q[0];
rz(-2.6197523) q[0];
rz(2.8537967) q[1];
sx q[1];
rz(-2.5591873) q[1];
sx q[1];
rz(-0.60307455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0780132) q[0];
sx q[0];
rz(-1.619864) q[0];
sx q[0];
rz(0.087875333) q[0];
rz(2.4580946) q[2];
sx q[2];
rz(-2.0216235) q[2];
sx q[2];
rz(1.766654) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.52023639) q[1];
sx q[1];
rz(-1.6545611) q[1];
sx q[1];
rz(-2.7158974) q[1];
rz(0.92080812) q[3];
sx q[3];
rz(-2.7201456) q[3];
sx q[3];
rz(-1.309066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4294943) q[2];
sx q[2];
rz(-2.3488729) q[2];
sx q[2];
rz(-2.402044) q[2];
rz(-1.0150917) q[3];
sx q[3];
rz(-0.52777094) q[3];
sx q[3];
rz(3.1225865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58448142) q[0];
sx q[0];
rz(-0.33753532) q[0];
sx q[0];
rz(0.8648411) q[0];
rz(0.64019126) q[1];
sx q[1];
rz(-2.4914111) q[1];
sx q[1];
rz(-2.6472299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011117293) q[0];
sx q[0];
rz(-2.0164765) q[0];
sx q[0];
rz(2.9167239) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96832595) q[2];
sx q[2];
rz(-1.5366388) q[2];
sx q[2];
rz(0.57847018) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3013005) q[1];
sx q[1];
rz(-1.5363785) q[1];
sx q[1];
rz(3.084206) q[1];
x q[2];
rz(-2.3253489) q[3];
sx q[3];
rz(-0.67749575) q[3];
sx q[3];
rz(1.578581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2690027) q[2];
sx q[2];
rz(-2.5073017) q[2];
sx q[2];
rz(-2.1799901) q[2];
rz(-1.6327935) q[3];
sx q[3];
rz(-1.6274933) q[3];
sx q[3];
rz(3.0646724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2831869) q[0];
sx q[0];
rz(-2.7981693) q[0];
sx q[0];
rz(2.0984233) q[0];
rz(-0.6768325) q[1];
sx q[1];
rz(-2.4969641) q[1];
sx q[1];
rz(0.72789311) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90157824) q[0];
sx q[0];
rz(-1.9838617) q[0];
sx q[0];
rz(-3.0204642) q[0];
rz(2.3995326) q[2];
sx q[2];
rz(-1.5681609) q[2];
sx q[2];
rz(-1.2330841) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54503578) q[1];
sx q[1];
rz(-1.2378197) q[1];
sx q[1];
rz(-0.92988981) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58957995) q[3];
sx q[3];
rz(-1.8581207) q[3];
sx q[3];
rz(-0.048487566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68208838) q[2];
sx q[2];
rz(-0.75152087) q[2];
sx q[2];
rz(-0.53306836) q[2];
rz(-1.3047949) q[3];
sx q[3];
rz(-0.47520906) q[3];
sx q[3];
rz(-2.3615725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.099139) q[0];
sx q[0];
rz(-0.071439698) q[0];
sx q[0];
rz(-3.1169917) q[0];
rz(-3.091231) q[1];
sx q[1];
rz(-2.5854526) q[1];
sx q[1];
rz(-2.3742026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3022223) q[0];
sx q[0];
rz(-1.5451533) q[0];
sx q[0];
rz(0.54021949) q[0];
rz(-0.51961706) q[2];
sx q[2];
rz(-0.37562464) q[2];
sx q[2];
rz(0.56359282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.394262) q[1];
sx q[1];
rz(-2.39191) q[1];
sx q[1];
rz(0.59632991) q[1];
rz(-pi) q[2];
rz(-0.57524753) q[3];
sx q[3];
rz(-1.9048371) q[3];
sx q[3];
rz(-0.0041242139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6699827) q[2];
sx q[2];
rz(-0.47405425) q[2];
sx q[2];
rz(-1.6101884) q[2];
rz(-1.0129741) q[3];
sx q[3];
rz(-1.6771202) q[3];
sx q[3];
rz(0.60232919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7811964) q[0];
sx q[0];
rz(-0.66806) q[0];
sx q[0];
rz(0.49041954) q[0];
rz(-2.7410653) q[1];
sx q[1];
rz(-0.35919765) q[1];
sx q[1];
rz(1.5708956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2119031) q[0];
sx q[0];
rz(-2.4725998) q[0];
sx q[0];
rz(0.59697911) q[0];
rz(-pi) q[1];
rz(-0.8428085) q[2];
sx q[2];
rz(-2.669816) q[2];
sx q[2];
rz(0.17201422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8388673) q[1];
sx q[1];
rz(-1.6714393) q[1];
sx q[1];
rz(-1.9221646) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83028806) q[3];
sx q[3];
rz(-0.92819302) q[3];
sx q[3];
rz(-0.45906367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5639497) q[2];
sx q[2];
rz(-2.9190013) q[2];
sx q[2];
rz(0.27601784) q[2];
rz(-0.66463071) q[3];
sx q[3];
rz(-2.1641927) q[3];
sx q[3];
rz(-0.19653921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4589602) q[0];
sx q[0];
rz(-0.17550547) q[0];
sx q[0];
rz(1.5631787) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-2.1140607) q[1];
sx q[1];
rz(3.077502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16312137) q[0];
sx q[0];
rz(-0.45867094) q[0];
sx q[0];
rz(-1.6244829) q[0];
rz(-0.61707492) q[2];
sx q[2];
rz(-1.1960746) q[2];
sx q[2];
rz(-0.20301698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7193072) q[1];
sx q[1];
rz(-1.6537084) q[1];
sx q[1];
rz(2.0598498) q[1];
x q[2];
rz(-2.552382) q[3];
sx q[3];
rz(-2.2452045) q[3];
sx q[3];
rz(-0.28705825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85612193) q[2];
sx q[2];
rz(-2.8557114) q[2];
sx q[2];
rz(-2.2355283) q[2];
rz(1.9237349) q[3];
sx q[3];
rz(-1.8702312) q[3];
sx q[3];
rz(2.0876032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.333552) q[0];
sx q[0];
rz(-1.522949) q[0];
sx q[0];
rz(-0.20725313) q[0];
rz(-2.8605657) q[1];
sx q[1];
rz(-1.749975) q[1];
sx q[1];
rz(-0.71798807) q[1];
rz(0.40194527) q[2];
sx q[2];
rz(-1.6464169) q[2];
sx q[2];
rz(0.34435549) q[2];
rz(-0.47697502) q[3];
sx q[3];
rz(-1.5834482) q[3];
sx q[3];
rz(-0.96848828) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
