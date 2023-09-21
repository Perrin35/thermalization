OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.84500399) q[0];
sx q[0];
rz(-0.72405976) q[0];
sx q[0];
rz(1.4847423) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(-2.2533921) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9125497) q[0];
sx q[0];
rz(-2.9648844) q[0];
sx q[0];
rz(2.8852709) q[0];
rz(-0.16562478) q[2];
sx q[2];
rz(-1.0399482) q[2];
sx q[2];
rz(2.3373375) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8623212) q[1];
sx q[1];
rz(-2.9746378) q[1];
sx q[1];
rz(-0.0032940666) q[1];
x q[2];
rz(2.4429951) q[3];
sx q[3];
rz(-1.5454096) q[3];
sx q[3];
rz(2.4217525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28329864) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-1.1179914) q[2];
rz(-2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(-3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730826) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(3.0156946) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1032216) q[0];
sx q[0];
rz(-1.201259) q[0];
sx q[0];
rz(-1.0612556) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6217896) q[2];
sx q[2];
rz(-2.1513373) q[2];
sx q[2];
rz(2.7948126) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1344578) q[1];
sx q[1];
rz(-1.4440698) q[1];
sx q[1];
rz(-3.111372) q[1];
rz(-pi) q[2];
rz(-1.5829854) q[3];
sx q[3];
rz(-0.72353957) q[3];
sx q[3];
rz(-0.42750588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(-1.4240501) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(2.5542636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5858784) q[0];
sx q[0];
rz(-0.782574) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(-0.084005984) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(1.9817339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1201598) q[0];
sx q[0];
rz(-2.624357) q[0];
sx q[0];
rz(-1.73818) q[0];
rz(2.1986387) q[2];
sx q[2];
rz(-1.1477594) q[2];
sx q[2];
rz(1.1689651) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68418903) q[1];
sx q[1];
rz(-1.1154798) q[1];
sx q[1];
rz(1.9940358) q[1];
x q[2];
rz(2.1163164) q[3];
sx q[3];
rz(-1.5226411) q[3];
sx q[3];
rz(-0.84850509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(0.1082871) q[2];
rz(0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(-0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9765587) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(-1.6595586) q[0];
rz(0.7011134) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(-0.28465095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6391622) q[0];
sx q[0];
rz(-2.6579755) q[0];
sx q[0];
rz(1.731427) q[0];
rz(-0.95732032) q[2];
sx q[2];
rz(-2.7942604) q[2];
sx q[2];
rz(2.6383102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.804867) q[1];
sx q[1];
rz(-1.6972099) q[1];
sx q[1];
rz(-0.079041914) q[1];
x q[2];
rz(-2.2399726) q[3];
sx q[3];
rz(-1.9603143) q[3];
sx q[3];
rz(0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.231679) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(-2.5740734) q[2];
rz(2.7275758) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48859566) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(1.8150785) q[0];
rz(-1.9891706) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-2.2036536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4612761) q[0];
sx q[0];
rz(-0.094705908) q[0];
sx q[0];
rz(1.8041496) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3100501) q[2];
sx q[2];
rz(-2.4387896) q[2];
sx q[2];
rz(-1.0520832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.93463072) q[1];
sx q[1];
rz(-0.043255581) q[1];
sx q[1];
rz(-0.91911493) q[1];
x q[2];
rz(2.786817) q[3];
sx q[3];
rz(-1.6924904) q[3];
sx q[3];
rz(2.0718758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-1.0413292) q[2];
rz(1.8390309) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(1.8235824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(0.55737108) q[0];
rz(0.56466651) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(-0.55647892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0628478) q[0];
sx q[0];
rz(-1.9468465) q[0];
sx q[0];
rz(-1.3486805) q[0];
x q[1];
rz(0.17369341) q[2];
sx q[2];
rz(-0.86280338) q[2];
sx q[2];
rz(0.049335418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0007799) q[1];
sx q[1];
rz(-1.8218578) q[1];
sx q[1];
rz(-2.954133) q[1];
rz(-pi) q[2];
rz(3.0895124) q[3];
sx q[3];
rz(-0.77643231) q[3];
sx q[3];
rz(-1.6805502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53211987) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790134) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(-0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(2.6409805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4780873) q[0];
sx q[0];
rz(-2.9604719) q[0];
sx q[0];
rz(-1.3785133) q[0];
x q[1];
rz(-1.074965) q[2];
sx q[2];
rz(-2.5410286) q[2];
sx q[2];
rz(-1.1662837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7350406) q[1];
sx q[1];
rz(-0.41932377) q[1];
sx q[1];
rz(-0.29962824) q[1];
rz(-pi) q[2];
x q[2];
rz(2.72243) q[3];
sx q[3];
rz(-2.3401642) q[3];
sx q[3];
rz(-0.49405801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(-2.7006941) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(-2.6202257) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0764517) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(1.0674397) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(-1.6759466) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6976801) q[0];
sx q[0];
rz(-1.7051823) q[0];
sx q[0];
rz(2.0356376) q[0];
rz(-pi) q[1];
rz(0.18657121) q[2];
sx q[2];
rz(-1.6779043) q[2];
sx q[2];
rz(-0.69498108) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2579736) q[1];
sx q[1];
rz(-2.5273364) q[1];
sx q[1];
rz(-0.96354624) q[1];
rz(-pi) q[2];
rz(-2.500324) q[3];
sx q[3];
rz(-2.0740168) q[3];
sx q[3];
rz(-0.36676952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(1.6652997) q[2];
rz(-2.2680797) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2575689) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(-2.4043758) q[0];
rz(0.018741477) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-0.92528701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6773274) q[0];
sx q[0];
rz(-2.263501) q[0];
sx q[0];
rz(0.6662743) q[0];
rz(-pi) q[1];
rz(0.16889062) q[2];
sx q[2];
rz(-2.4714111) q[2];
sx q[2];
rz(2.9549213) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9463897) q[1];
sx q[1];
rz(-2.9159947) q[1];
sx q[1];
rz(2.5220847) q[1];
rz(-1.0989283) q[3];
sx q[3];
rz(-0.77583757) q[3];
sx q[3];
rz(-0.27234205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(1.0037237) q[2];
rz(-0.090099661) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(-1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1019679) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(-0.75751799) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19125464) q[0];
sx q[0];
rz(-1.7636824) q[0];
sx q[0];
rz(0.14001503) q[0];
rz(-pi) q[1];
rz(-2.1610356) q[2];
sx q[2];
rz(-1.7014116) q[2];
sx q[2];
rz(-0.57934258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47093686) q[1];
sx q[1];
rz(-0.48479143) q[1];
sx q[1];
rz(2.0623341) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71000368) q[3];
sx q[3];
rz(-1.5492951) q[3];
sx q[3];
rz(1.8332675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99047986) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(-2.347836) q[2];
rz(-0.63888597) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(-1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6486075) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-2.5095148) q[2];
sx q[2];
rz(-0.52543228) q[2];
sx q[2];
rz(0.57731522) q[2];
rz(-0.32140857) q[3];
sx q[3];
rz(-2.1264429) q[3];
sx q[3];
rz(0.85254729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];