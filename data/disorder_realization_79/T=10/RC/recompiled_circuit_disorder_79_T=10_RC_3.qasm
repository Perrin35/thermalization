OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(0.19302364) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8773168) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(-0.067406128) q[0];
rz(-0.24129759) q[2];
sx q[2];
rz(-1.9120875) q[2];
sx q[2];
rz(1.2933921) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.623917) q[1];
sx q[1];
rz(-2.3530934) q[1];
sx q[1];
rz(-0.6448402) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0055553) q[3];
sx q[3];
rz(-1.2473277) q[3];
sx q[3];
rz(-2.2068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(-1.0985628) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(1.1014972) q[3];
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
rz(-1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(-0.20794491) q[0];
rz(0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.6764486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6786878) q[0];
sx q[0];
rz(-1.5683187) q[0];
sx q[0];
rz(-2.4173827) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1128088) q[2];
sx q[2];
rz(-2.6070242) q[2];
sx q[2];
rz(2.2528258) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84491731) q[1];
sx q[1];
rz(-1.0122609) q[1];
sx q[1];
rz(-2.1293473) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37104718) q[3];
sx q[3];
rz(-0.78073946) q[3];
sx q[3];
rz(-2.3649529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1229822) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(3.0318276) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96317545) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(0.57488817) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(2.3410472) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.752906) q[0];
sx q[0];
rz(-1.1645002) q[0];
sx q[0];
rz(2.9857062) q[0];
rz(-pi) q[1];
rz(0.66652253) q[2];
sx q[2];
rz(-2.1608673) q[2];
sx q[2];
rz(-0.18596622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.264818) q[1];
sx q[1];
rz(-2.7285828) q[1];
sx q[1];
rz(-3.1289711) q[1];
rz(-pi) q[2];
rz(-1.21739) q[3];
sx q[3];
rz(-0.73261515) q[3];
sx q[3];
rz(-2.4917847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(2.676679) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(1.0850614) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61468716) q[0];
sx q[0];
rz(-0.58262107) q[0];
sx q[0];
rz(2.7132062) q[0];
rz(-pi) q[1];
rz(-2.1539139) q[2];
sx q[2];
rz(-2.7104125) q[2];
sx q[2];
rz(-2.0476066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.19238732) q[1];
sx q[1];
rz(-2.9678223) q[1];
sx q[1];
rz(2.5237571) q[1];
x q[2];
rz(-1.9784847) q[3];
sx q[3];
rz(-0.99916047) q[3];
sx q[3];
rz(0.15743263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-0.26369035) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8624449) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(2.7092343) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(2.1599105) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0504793) q[0];
sx q[0];
rz(-1.0266773) q[0];
sx q[0];
rz(-1.0247466) q[0];
x q[1];
rz(1.8129187) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(2.9850933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9202068) q[1];
sx q[1];
rz(-1.1077987) q[1];
sx q[1];
rz(1.7358001) q[1];
x q[2];
rz(-0.97335191) q[3];
sx q[3];
rz(-1.5095599) q[3];
sx q[3];
rz(-0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.6550487) q[2];
sx q[2];
rz(2.6521818) q[2];
rz(-2.126746) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(-1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(-1.3308446) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(-2.9763124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3344791) q[0];
sx q[0];
rz(-1.665859) q[0];
sx q[0];
rz(-1.6001742) q[0];
x q[1];
rz(-1.500962) q[2];
sx q[2];
rz(-0.88431057) q[2];
sx q[2];
rz(-1.7718466) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0300671) q[1];
sx q[1];
rz(-1.0377874) q[1];
sx q[1];
rz(2.0433321) q[1];
rz(0.23541707) q[3];
sx q[3];
rz(-2.0093007) q[3];
sx q[3];
rz(2.8188843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91339397) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(-1.6983263) q[2];
rz(0.36744395) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(0.34926397) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(2.4051037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9518785) q[0];
sx q[0];
rz(-1.6197546) q[0];
sx q[0];
rz(3.1244856) q[0];
x q[1];
rz(0.24918208) q[2];
sx q[2];
rz(-1.2417925) q[2];
sx q[2];
rz(-0.56275425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6474825) q[1];
sx q[1];
rz(-0.78299114) q[1];
sx q[1];
rz(-2.8857735) q[1];
rz(0.50076671) q[3];
sx q[3];
rz(-1.1874677) q[3];
sx q[3];
rz(0.70481833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(-1.8435562) q[0];
rz(0.80728665) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(2.2198026) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73496504) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(-2.448003) q[0];
rz(2.8554101) q[2];
sx q[2];
rz(-2.2879062) q[2];
sx q[2];
rz(1.66695) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39107716) q[1];
sx q[1];
rz(-0.32954307) q[1];
sx q[1];
rz(0.68117546) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7796302) q[3];
sx q[3];
rz(-2.4224671) q[3];
sx q[3];
rz(-1.5428839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(-1.3575859) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(0.67970651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98201671) q[0];
sx q[0];
rz(-1.7734818) q[0];
sx q[0];
rz(2.6152339) q[0];
rz(-pi) q[1];
rz(-2.26608) q[2];
sx q[2];
rz(-1.7098223) q[2];
sx q[2];
rz(-1.8347486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3956086) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(-0.70396522) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71984843) q[3];
sx q[3];
rz(-2.3206629) q[3];
sx q[3];
rz(-1.6380701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
rz(0.05242059) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(0.7146548) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(-1.9649327) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3424073) q[0];
sx q[0];
rz(-2.9710794) q[0];
sx q[0];
rz(1.8605581) q[0];
rz(-pi) q[1];
rz(-0.93675905) q[2];
sx q[2];
rz(-2.0271795) q[2];
sx q[2];
rz(-0.87181834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51917045) q[1];
sx q[1];
rz(-1.4960438) q[1];
sx q[1];
rz(-2.1291817) q[1];
rz(-0.15575274) q[3];
sx q[3];
rz(-0.51625801) q[3];
sx q[3];
rz(1.932365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8979793) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(2.0991142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(0.30504967) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(-1.5626004) q[3];
sx q[3];
rz(-1.9215487) q[3];
sx q[3];
rz(2.4196845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
