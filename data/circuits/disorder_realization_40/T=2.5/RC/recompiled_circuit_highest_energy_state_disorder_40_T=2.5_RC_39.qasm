OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1542198) q[0];
sx q[0];
rz(-1.5232824) q[0];
sx q[0];
rz(-2.9989938) q[0];
rz(0.59347403) q[1];
sx q[1];
rz(-2.4886517) q[1];
sx q[1];
rz(-1.0452193) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0526177) q[0];
sx q[0];
rz(-1.0427022) q[0];
sx q[0];
rz(-1.0333956) q[0];
x q[1];
rz(-2.2085243) q[2];
sx q[2];
rz(-1.6003813) q[2];
sx q[2];
rz(2.3552908) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3732304) q[1];
sx q[1];
rz(-1.6715896) q[1];
sx q[1];
rz(1.412315) q[1];
rz(-pi) q[2];
rz(-1.0020489) q[3];
sx q[3];
rz(-1.6709329) q[3];
sx q[3];
rz(1.1796234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1656701) q[2];
sx q[2];
rz(-2.133635) q[2];
sx q[2];
rz(-0.21318501) q[2];
rz(2.2255157) q[3];
sx q[3];
rz(-0.35917425) q[3];
sx q[3];
rz(-2.7540414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857472) q[0];
sx q[0];
rz(-0.88749945) q[0];
sx q[0];
rz(-2.6708653) q[0];
rz(-2.0384906) q[1];
sx q[1];
rz(-2.6328937) q[1];
sx q[1];
rz(-0.57977605) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.267201) q[0];
sx q[0];
rz(-1.5109332) q[0];
sx q[0];
rz(0.52459985) q[0];
x q[1];
rz(-0.51885278) q[2];
sx q[2];
rz(-2.5343072) q[2];
sx q[2];
rz(2.3407206) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6439542) q[1];
sx q[1];
rz(-2.4259287) q[1];
sx q[1];
rz(2.1791451) q[1];
rz(1.9409431) q[3];
sx q[3];
rz(-1.8164299) q[3];
sx q[3];
rz(0.64526886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-0.85593587) q[2];
sx q[2];
rz(-0.095005438) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-0.78229457) q[3];
sx q[3];
rz(1.6949722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59548241) q[0];
sx q[0];
rz(-0.70850104) q[0];
sx q[0];
rz(-2.8693759) q[0];
rz(3.0369924) q[1];
sx q[1];
rz(-1.0403386) q[1];
sx q[1];
rz(1.6786172) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3264831) q[0];
sx q[0];
rz(-2.4158951) q[0];
sx q[0];
rz(-1.4952907) q[0];
rz(-pi) q[1];
rz(-0.06426364) q[2];
sx q[2];
rz(-2.3552364) q[2];
sx q[2];
rz(0.83640316) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.926907) q[1];
sx q[1];
rz(-0.72539293) q[1];
sx q[1];
rz(0.43154413) q[1];
x q[2];
rz(0.36896173) q[3];
sx q[3];
rz(-1.6272244) q[3];
sx q[3];
rz(-1.4101613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0766803) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(-2.7317375) q[2];
rz(3.0900433) q[3];
sx q[3];
rz(-2.1487273) q[3];
sx q[3];
rz(-0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16875295) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(2.4421316) q[0];
rz(0.93060023) q[1];
sx q[1];
rz(-1.7761296) q[1];
sx q[1];
rz(1.0994937) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65419009) q[0];
sx q[0];
rz(-2.497597) q[0];
sx q[0];
rz(2.0430768) q[0];
rz(-1.0223939) q[2];
sx q[2];
rz(-1.363593) q[2];
sx q[2];
rz(2.0275379) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2840957) q[1];
sx q[1];
rz(-1.2624143) q[1];
sx q[1];
rz(-0.33607884) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64208018) q[3];
sx q[3];
rz(-1.9028946) q[3];
sx q[3];
rz(1.4043216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3804271) q[2];
sx q[2];
rz(-2.0695504) q[2];
sx q[2];
rz(-2.2011444) q[2];
rz(0.56194168) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(2.565062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1171653) q[0];
sx q[0];
rz(-3.0553525) q[0];
sx q[0];
rz(2.1249791) q[0];
rz(-2.2158465) q[1];
sx q[1];
rz(-0.75030202) q[1];
sx q[1];
rz(-3.064916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6252687) q[0];
sx q[0];
rz(-1.4051354) q[0];
sx q[0];
rz(-1.1722086) q[0];
x q[1];
rz(1.8163278) q[2];
sx q[2];
rz(-1.7707728) q[2];
sx q[2];
rz(0.94975805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.58248427) q[1];
sx q[1];
rz(-0.6000207) q[1];
sx q[1];
rz(1.9514854) q[1];
rz(-pi) q[2];
rz(-1.2604976) q[3];
sx q[3];
rz(-2.3633133) q[3];
sx q[3];
rz(0.33184127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76806796) q[2];
sx q[2];
rz(-0.7879476) q[2];
sx q[2];
rz(-0.12316556) q[2];
rz(0.67433107) q[3];
sx q[3];
rz(-0.47334039) q[3];
sx q[3];
rz(0.82790747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3097565) q[0];
sx q[0];
rz(-2.9530544) q[0];
sx q[0];
rz(-0.48102608) q[0];
rz(-0.24093974) q[1];
sx q[1];
rz(-2.6048581) q[1];
sx q[1];
rz(3.0066838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19770007) q[0];
sx q[0];
rz(-2.0257844) q[0];
sx q[0];
rz(-0.83857341) q[0];
rz(-2.606484) q[2];
sx q[2];
rz(-1.4819179) q[2];
sx q[2];
rz(-0.8366226) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.125306) q[1];
sx q[1];
rz(-1.2963489) q[1];
sx q[1];
rz(-0.80933615) q[1];
rz(0.14860146) q[3];
sx q[3];
rz(-2.6872748) q[3];
sx q[3];
rz(0.065520614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2065108) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(-1.8180465) q[2];
rz(2.8694618) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-0.14565295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226591) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(-0.65852273) q[0];
rz(-1.5497442) q[1];
sx q[1];
rz(-2.1287983) q[1];
sx q[1];
rz(0.50643593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50473266) q[0];
sx q[0];
rz(-1.9397282) q[0];
sx q[0];
rz(-1.4999267) q[0];
rz(-3.1108702) q[2];
sx q[2];
rz(-1.0892158) q[2];
sx q[2];
rz(-0.55847634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5665019) q[1];
sx q[1];
rz(-3.021214) q[1];
sx q[1];
rz(0.68415227) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5793367) q[3];
sx q[3];
rz(-1.5361353) q[3];
sx q[3];
rz(3.0140585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8062313) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(-1.2696421) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-3.1277872) q[3];
sx q[3];
rz(0.55240101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7954471) q[0];
sx q[0];
rz(-0.2427635) q[0];
sx q[0];
rz(-0.41918293) q[0];
rz(3.0508793) q[1];
sx q[1];
rz(-2.2234629) q[1];
sx q[1];
rz(1.027164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61495435) q[0];
sx q[0];
rz(-1.1246343) q[0];
sx q[0];
rz(-0.63132186) q[0];
rz(-pi) q[1];
rz(-1.5628603) q[2];
sx q[2];
rz(-1.5479701) q[2];
sx q[2];
rz(2.9812682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13090573) q[1];
sx q[1];
rz(-1.3584777) q[1];
sx q[1];
rz(-0.10879083) q[1];
rz(-pi) q[2];
x q[2];
rz(0.054612463) q[3];
sx q[3];
rz(-2.255618) q[3];
sx q[3];
rz(2.7582061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98296982) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(2.1515004) q[2];
rz(2.8236735) q[3];
sx q[3];
rz(-2.3293994) q[3];
sx q[3];
rz(-3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7541499) q[0];
sx q[0];
rz(-2.4328572) q[0];
sx q[0];
rz(-2.5842174) q[0];
rz(-2.7793461) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-2.974496) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97826577) q[0];
sx q[0];
rz(-1.8370795) q[0];
sx q[0];
rz(2.8243218) q[0];
rz(-pi) q[1];
rz(2.9252824) q[2];
sx q[2];
rz(-0.71291332) q[2];
sx q[2];
rz(2.0305433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7726248) q[1];
sx q[1];
rz(-2.4030622) q[1];
sx q[1];
rz(2.6767297) q[1];
x q[2];
rz(-1.4792419) q[3];
sx q[3];
rz(-0.6855489) q[3];
sx q[3];
rz(0.26622546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1780213) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(-2.7131405) q[2];
rz(0.7306478) q[3];
sx q[3];
rz(-2.0615536) q[3];
sx q[3];
rz(0.60034269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927602) q[0];
sx q[0];
rz(-1.6535783) q[0];
sx q[0];
rz(1.1195419) q[0];
rz(-2.4730832) q[1];
sx q[1];
rz(-2.5709277) q[1];
sx q[1];
rz(0.25984919) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52212673) q[0];
sx q[0];
rz(-1.3106951) q[0];
sx q[0];
rz(1.0570265) q[0];
rz(2.75801) q[2];
sx q[2];
rz(-1.1337987) q[2];
sx q[2];
rz(-2.008977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47626859) q[1];
sx q[1];
rz(-0.55108738) q[1];
sx q[1];
rz(2.4399806) q[1];
rz(-pi) q[2];
rz(-1.3803452) q[3];
sx q[3];
rz(-0.69882876) q[3];
sx q[3];
rz(2.8052398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4708289) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(-2.0551576) q[2];
rz(0.49232617) q[3];
sx q[3];
rz(-0.84572518) q[3];
sx q[3];
rz(0.72211784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54871854) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(-1.2552352) q[1];
sx q[1];
rz(-1.3175169) q[1];
sx q[1];
rz(-3.0138737) q[1];
rz(-1.3484501) q[2];
sx q[2];
rz(-0.4081453) q[2];
sx q[2];
rz(-1.2884507) q[2];
rz(-0.063379824) q[3];
sx q[3];
rz(-2.2324003) q[3];
sx q[3];
rz(1.5008055) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
