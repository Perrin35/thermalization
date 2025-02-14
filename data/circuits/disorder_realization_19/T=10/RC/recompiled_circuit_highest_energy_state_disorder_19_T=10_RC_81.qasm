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
rz(0.28111464) q[0];
sx q[0];
rz(1.7908362) q[0];
sx q[0];
rz(10.375441) q[0];
rz(0.22275337) q[1];
sx q[1];
rz(-2.8877701) q[1];
sx q[1];
rz(2.4737127) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1199824) q[0];
sx q[0];
rz(-1.9304781) q[0];
sx q[0];
rz(-1.5403454) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6654885) q[2];
sx q[2];
rz(-0.166278) q[2];
sx q[2];
rz(-2.4537697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8409792) q[1];
sx q[1];
rz(-1.4491699) q[1];
sx q[1];
rz(-3.0517314) q[1];
rz(-pi) q[2];
rz(0.018785211) q[3];
sx q[3];
rz(-0.21078706) q[3];
sx q[3];
rz(0.8715521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5369947) q[2];
sx q[2];
rz(-0.94702417) q[2];
sx q[2];
rz(-2.9941518) q[2];
rz(1.3747181) q[3];
sx q[3];
rz(-1.7539897) q[3];
sx q[3];
rz(-1.714777) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0680256) q[0];
sx q[0];
rz(-1.0140714) q[0];
sx q[0];
rz(-2.9378939) q[0];
rz(-3.0286466) q[1];
sx q[1];
rz(-2.4599383) q[1];
sx q[1];
rz(-3.122094) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18960139) q[0];
sx q[0];
rz(-0.87965779) q[0];
sx q[0];
rz(2.7777872) q[0];
rz(1.9505368) q[2];
sx q[2];
rz(-1.4342208) q[2];
sx q[2];
rz(1.6644163) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9128117) q[1];
sx q[1];
rz(-1.8789018) q[1];
sx q[1];
rz(-0.56477115) q[1];
rz(-pi) q[2];
rz(-0.89159151) q[3];
sx q[3];
rz(-2.0304095) q[3];
sx q[3];
rz(-0.63320049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.514275) q[2];
sx q[2];
rz(-1.5760199) q[2];
sx q[2];
rz(0.19228284) q[2];
rz(0.2187885) q[3];
sx q[3];
rz(-1.8407121) q[3];
sx q[3];
rz(1.7910262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.669765) q[0];
sx q[0];
rz(-0.28194031) q[0];
sx q[0];
rz(-2.6523253) q[0];
rz(-1.4504112) q[1];
sx q[1];
rz(-1.9905118) q[1];
sx q[1];
rz(-0.53389186) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.691026) q[0];
sx q[0];
rz(-1.8922784) q[0];
sx q[0];
rz(1.5430928) q[0];
x q[1];
rz(0.30812906) q[2];
sx q[2];
rz(-0.14585431) q[2];
sx q[2];
rz(3.0792091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7545751) q[1];
sx q[1];
rz(-1.1711517) q[1];
sx q[1];
rz(-0.15367457) q[1];
rz(-pi) q[2];
rz(2.4160014) q[3];
sx q[3];
rz(-1.4342566) q[3];
sx q[3];
rz(-1.2927593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.78407946) q[2];
sx q[2];
rz(-2.1046905) q[2];
sx q[2];
rz(2.4927523) q[2];
rz(2.8283548) q[3];
sx q[3];
rz(-2.1514838) q[3];
sx q[3];
rz(1.8296957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.677815) q[0];
sx q[0];
rz(-2.2721993) q[0];
sx q[0];
rz(0.75611269) q[0];
rz(-0.4568049) q[1];
sx q[1];
rz(-2.017338) q[1];
sx q[1];
rz(0.66064984) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9220369) q[0];
sx q[0];
rz(-1.7496373) q[0];
sx q[0];
rz(-0.40750176) q[0];
rz(-pi) q[1];
rz(-0.22539745) q[2];
sx q[2];
rz(-1.9323914) q[2];
sx q[2];
rz(2.0394675) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7380437) q[1];
sx q[1];
rz(-1.7521841) q[1];
sx q[1];
rz(2.2970143) q[1];
rz(-pi) q[2];
rz(-2.0235429) q[3];
sx q[3];
rz(-1.2635316) q[3];
sx q[3];
rz(-1.6045517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75055355) q[2];
sx q[2];
rz(-1.6497352) q[2];
sx q[2];
rz(2.1569815) q[2];
rz(-1.6330968) q[3];
sx q[3];
rz(-1.0909785) q[3];
sx q[3];
rz(-0.35198894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9462747) q[0];
sx q[0];
rz(-1.0349422) q[0];
sx q[0];
rz(0.6629194) q[0];
rz(1.7341057) q[1];
sx q[1];
rz(-0.88277849) q[1];
sx q[1];
rz(0.9443121) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1809826) q[0];
sx q[0];
rz(-2.2274605) q[0];
sx q[0];
rz(1.3707536) q[0];
x q[1];
rz(-0.92534368) q[2];
sx q[2];
rz(-2.1133377) q[2];
sx q[2];
rz(2.0069356) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1270557) q[1];
sx q[1];
rz(-0.98871283) q[1];
sx q[1];
rz(-2.0392557) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7529519) q[3];
sx q[3];
rz(-2.5699617) q[3];
sx q[3];
rz(1.1370259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0974836) q[2];
sx q[2];
rz(-0.58610836) q[2];
sx q[2];
rz(3.0666472) q[2];
rz(-2.1465541) q[3];
sx q[3];
rz(-2.0375662) q[3];
sx q[3];
rz(-1.1546571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5465882) q[0];
sx q[0];
rz(-1.1052479) q[0];
sx q[0];
rz(-2.8357491) q[0];
rz(-2.2587237) q[1];
sx q[1];
rz(-0.63778937) q[1];
sx q[1];
rz(0.84552228) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8354685) q[0];
sx q[0];
rz(-1.7996664) q[0];
sx q[0];
rz(-0.86979515) q[0];
rz(-0.0011618753) q[2];
sx q[2];
rz(-1.2865304) q[2];
sx q[2];
rz(-0.53476108) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.90585044) q[1];
sx q[1];
rz(-1.0585856) q[1];
sx q[1];
rz(2.3925969) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85725089) q[3];
sx q[3];
rz(-1.4298059) q[3];
sx q[3];
rz(-0.14234358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66158295) q[2];
sx q[2];
rz(-2.2955387) q[2];
sx q[2];
rz(1.0711077) q[2];
rz(0.24460159) q[3];
sx q[3];
rz(-0.94541234) q[3];
sx q[3];
rz(1.6036114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59413183) q[0];
sx q[0];
rz(-2.5762711) q[0];
sx q[0];
rz(-2.2070337) q[0];
rz(2.4681828) q[1];
sx q[1];
rz(-2.5118561) q[1];
sx q[1];
rz(-3.0455132) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.729674) q[0];
sx q[0];
rz(-0.11765495) q[0];
sx q[0];
rz(1.3484701) q[0];
rz(1.3911909) q[2];
sx q[2];
rz(-1.7471004) q[2];
sx q[2];
rz(1.7352895) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90111946) q[1];
sx q[1];
rz(-0.34209004) q[1];
sx q[1];
rz(-0.53804496) q[1];
rz(-pi) q[2];
rz(-1.6953498) q[3];
sx q[3];
rz(-2.4807408) q[3];
sx q[3];
rz(-1.648267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38976321) q[2];
sx q[2];
rz(-1.7182173) q[2];
sx q[2];
rz(0.54330379) q[2];
rz(-3.1144888) q[3];
sx q[3];
rz(-0.47309858) q[3];
sx q[3];
rz(-2.134034) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84697023) q[0];
sx q[0];
rz(-0.68055081) q[0];
sx q[0];
rz(2.3934225) q[0];
rz(1.4717884) q[1];
sx q[1];
rz(-1.4133778) q[1];
sx q[1];
rz(-0.69563037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78139898) q[0];
sx q[0];
rz(-1.8273786) q[0];
sx q[0];
rz(3.0220023) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7901445) q[2];
sx q[2];
rz(-2.2168014) q[2];
sx q[2];
rz(-1.1644808) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4804617) q[1];
sx q[1];
rz(-0.5782402) q[1];
sx q[1];
rz(-0.64774143) q[1];
rz(-pi) q[2];
rz(2.3588156) q[3];
sx q[3];
rz(-0.15412384) q[3];
sx q[3];
rz(-1.085404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2601605) q[2];
sx q[2];
rz(-0.83923927) q[2];
sx q[2];
rz(3.0879171) q[2];
rz(-1.3850877) q[3];
sx q[3];
rz(-1.799492) q[3];
sx q[3];
rz(-0.16689859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78228918) q[0];
sx q[0];
rz(-1.259869) q[0];
sx q[0];
rz(-2.2480929) q[0];
rz(-2.9529849) q[1];
sx q[1];
rz(-0.74917561) q[1];
sx q[1];
rz(-0.20208727) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24465626) q[0];
sx q[0];
rz(-0.74743987) q[0];
sx q[0];
rz(-1.5267751) q[0];
rz(2.0941433) q[2];
sx q[2];
rz(-1.6138895) q[2];
sx q[2];
rz(0.76165253) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.34101453) q[1];
sx q[1];
rz(-1.0381925) q[1];
sx q[1];
rz(2.8982603) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7703448) q[3];
sx q[3];
rz(-1.8168257) q[3];
sx q[3];
rz(0.061316874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.14822745) q[2];
sx q[2];
rz(-1.3067641) q[2];
sx q[2];
rz(-2.4723049) q[2];
rz(-2.3547442) q[3];
sx q[3];
rz(-1.9653178) q[3];
sx q[3];
rz(-2.440786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4871019) q[0];
sx q[0];
rz(-2.632532) q[0];
sx q[0];
rz(-1.85602) q[0];
rz(0.14699832) q[1];
sx q[1];
rz(-1.1451984) q[1];
sx q[1];
rz(0.95050341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72785801) q[0];
sx q[0];
rz(-1.04299) q[0];
sx q[0];
rz(-2.1377273) q[0];
rz(-pi) q[1];
rz(-1.4661524) q[2];
sx q[2];
rz(-2.6015835) q[2];
sx q[2];
rz(-1.7162985) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9933789) q[1];
sx q[1];
rz(-1.0328173) q[1];
sx q[1];
rz(2.5003273) q[1];
x q[2];
rz(-1.688429) q[3];
sx q[3];
rz(-1.7087987) q[3];
sx q[3];
rz(1.3641629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46108437) q[2];
sx q[2];
rz(-0.66404873) q[2];
sx q[2];
rz(-0.39592478) q[2];
rz(0.33216533) q[3];
sx q[3];
rz(-1.9931404) q[3];
sx q[3];
rz(-0.90126669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.111515) q[0];
sx q[0];
rz(-2.5319396) q[0];
sx q[0];
rz(1.4421705) q[0];
rz(-0.036399966) q[1];
sx q[1];
rz(-1.2926688) q[1];
sx q[1];
rz(-1.7784437) q[1];
rz(-1.2716952) q[2];
sx q[2];
rz(-0.80812412) q[2];
sx q[2];
rz(-3.0199188) q[2];
rz(-0.069478819) q[3];
sx q[3];
rz(-2.4705349) q[3];
sx q[3];
rz(0.51746838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
