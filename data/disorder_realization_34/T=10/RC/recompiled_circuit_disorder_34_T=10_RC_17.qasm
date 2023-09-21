OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(-0.43963471) q[0];
sx q[0];
rz(3.0602732) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(4.4250017) q[1];
sx q[1];
rz(11.783574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62277943) q[0];
sx q[0];
rz(-2.8840027) q[0];
sx q[0];
rz(1.326484) q[0];
x q[1];
rz(2.9847203) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(-3.0770609) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1374955) q[1];
sx q[1];
rz(-1.6988519) q[1];
sx q[1];
rz(-1.1346243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4814032) q[3];
sx q[3];
rz(-1.8854183) q[3];
sx q[3];
rz(-2.6320626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2471182) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(0.11581126) q[2];
rz(-1.5995021) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2540934) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(0.19533531) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(0.23981747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0175184) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(-0.56970861) q[0];
x q[1];
rz(1.2698783) q[2];
sx q[2];
rz(-1.3687203) q[2];
sx q[2];
rz(2.2965455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5259243) q[1];
sx q[1];
rz(-2.3500372) q[1];
sx q[1];
rz(-0.073992373) q[1];
rz(-pi) q[2];
rz(-2.2858743) q[3];
sx q[3];
rz(-0.021589605) q[3];
sx q[3];
rz(-1.74687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6372765) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(1.3416393) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(-1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2770237) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(0.66816107) q[0];
rz(-1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(1.0659165) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.802793) q[0];
sx q[0];
rz(-0.25164139) q[0];
sx q[0];
rz(1.7358857) q[0];
x q[1];
rz(0.5553603) q[2];
sx q[2];
rz(-1.4484324) q[2];
sx q[2];
rz(0.33081474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6560116) q[1];
sx q[1];
rz(-2.1335568) q[1];
sx q[1];
rz(2.807711) q[1];
rz(-3.021574) q[3];
sx q[3];
rz(-2.2285322) q[3];
sx q[3];
rz(-1.1413871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.187591) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(2.4348863) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(-1.6548086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0006205) q[0];
sx q[0];
rz(-1.9019433) q[0];
sx q[0];
rz(-1.8676057) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9183667) q[2];
sx q[2];
rz(-2.1516557) q[2];
sx q[2];
rz(-0.52085224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6390037) q[1];
sx q[1];
rz(-1.2856312) q[1];
sx q[1];
rz(1.4147948) q[1];
x q[2];
rz(-2.3151822) q[3];
sx q[3];
rz(-1.5658169) q[3];
sx q[3];
rz(-1.5116364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(2.13307) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0356045) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.460357) q[0];
rz(-1.5530855) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(-3.1255186) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8668629) q[0];
sx q[0];
rz(-2.1511973) q[0];
sx q[0];
rz(2.4972563) q[0];
x q[1];
rz(1.5526505) q[2];
sx q[2];
rz(-0.43076736) q[2];
sx q[2];
rz(-2.919163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9863723) q[1];
sx q[1];
rz(-1.9934137) q[1];
sx q[1];
rz(-2.1834454) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5752951) q[3];
sx q[3];
rz(-1.9330977) q[3];
sx q[3];
rz(-0.11210657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(-2.0444929) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(-2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2621433) q[0];
sx q[0];
rz(-0.95603846) q[0];
sx q[0];
rz(0.26457796) q[0];
rz(-2.9333026) q[2];
sx q[2];
rz(-2.0049094) q[2];
sx q[2];
rz(1.2693894) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.743222) q[1];
sx q[1];
rz(-2.1757158) q[1];
sx q[1];
rz(-0.066992316) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25572689) q[3];
sx q[3];
rz(-2.7284107) q[3];
sx q[3];
rz(-0.43858389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9796824) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(2.5612223) q[0];
rz(-1.0549818) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(-2.4408128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.472815) q[0];
sx q[0];
rz(-2.6758709) q[0];
sx q[0];
rz(-1.8229683) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8220196) q[2];
sx q[2];
rz(-1.6817037) q[2];
sx q[2];
rz(0.061352913) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5138445) q[1];
sx q[1];
rz(-2.5690837) q[1];
sx q[1];
rz(1.4036914) q[1];
x q[2];
rz(-0.084214597) q[3];
sx q[3];
rz(-2.299752) q[3];
sx q[3];
rz(-2.6654411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78684029) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(2.2498806) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.007908) q[0];
sx q[0];
rz(-0.79151756) q[0];
sx q[0];
rz(-1.7931116) q[0];
x q[1];
rz(0.58015577) q[2];
sx q[2];
rz(-1.9556502) q[2];
sx q[2];
rz(2.5600524) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6139639) q[1];
sx q[1];
rz(-0.91273897) q[1];
sx q[1];
rz(-1.0440473) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31657747) q[3];
sx q[3];
rz(-2.205875) q[3];
sx q[3];
rz(2.3659335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6531758) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(-0.50576058) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(-0.23323664) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-0.72203565) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(-1.3649712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9997172) q[0];
sx q[0];
rz(-1.5293855) q[0];
sx q[0];
rz(-0.74383225) q[0];
rz(-0.083655595) q[2];
sx q[2];
rz(-1.9803279) q[2];
sx q[2];
rz(1.6905897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2024467) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(1.1254315) q[1];
rz(-pi) q[2];
rz(-0.98718231) q[3];
sx q[3];
rz(-0.35161388) q[3];
sx q[3];
rz(-2.3021063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1086796) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(0.96819425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6598845) q[0];
sx q[0];
rz(-1.9511001) q[0];
sx q[0];
rz(-0.36614059) q[0];
rz(2.7231611) q[2];
sx q[2];
rz(-1.1088587) q[2];
sx q[2];
rz(-0.07428169) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89322461) q[1];
sx q[1];
rz(-1.7176504) q[1];
sx q[1];
rz(2.8180772) q[1];
rz(-0.00023437436) q[3];
sx q[3];
rz(-1.2486542) q[3];
sx q[3];
rz(2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.7989981) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447727) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(1.8832072) q[2];
sx q[2];
rz(-0.84597833) q[2];
sx q[2];
rz(-0.99920263) q[2];
rz(2.1063741) q[3];
sx q[3];
rz(-0.62789161) q[3];
sx q[3];
rz(-0.37554489) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
