OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(-3.0602732) q[0];
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
rz(-1.9570219) q[0];
sx q[0];
rz(-1.6324568) q[0];
sx q[0];
rz(-1.3205359) q[0];
rz(2.9847203) q[2];
sx q[2];
rz(-1.1638767) q[2];
sx q[2];
rz(-3.0770609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3738457) q[1];
sx q[1];
rz(-2.0031553) q[1];
sx q[1];
rz(-3.0004764) q[1];
rz(0.26773914) q[3];
sx q[3];
rz(-2.8149238) q[3];
sx q[3];
rz(-0.79145811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2471182) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88749921) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(-0.37503606) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(0.23981747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0175184) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(-2.571884) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2698783) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(-2.2965455) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0071734) q[1];
sx q[1];
rz(-1.6234142) q[1];
sx q[1];
rz(0.79018553) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2858743) q[3];
sx q[3];
rz(-3.120003) q[3];
sx q[3];
rz(-1.3947226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2770237) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(-2.4734316) q[0];
rz(1.4913303) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(-1.0659165) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0720172) q[0];
sx q[0];
rz(-1.5298651) q[0];
sx q[0];
rz(1.3224365) q[0];
rz(-2.5862323) q[2];
sx q[2];
rz(-1.6931603) q[2];
sx q[2];
rz(2.8107779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6560116) q[1];
sx q[1];
rz(-2.1335568) q[1];
sx q[1];
rz(-2.807711) q[1];
rz(-pi) q[2];
x q[2];
rz(3.021574) q[3];
sx q[3];
rz(-2.2285322) q[3];
sx q[3];
rz(1.1413871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(-1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.4867841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1409722) q[0];
sx q[0];
rz(-1.9019433) q[0];
sx q[0];
rz(-1.8676057) q[0];
rz(-0.47866486) q[2];
sx q[2];
rz(-0.66648167) q[2];
sx q[2];
rz(2.037231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6390037) q[1];
sx q[1];
rz(-1.2856312) q[1];
sx q[1];
rz(1.7267978) q[1];
rz(-pi) q[2];
rz(-0.0067700245) q[3];
sx q[3];
rz(-2.3151708) q[3];
sx q[3];
rz(0.063746728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(1.0085227) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.460357) q[0];
rz(-1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(3.1255186) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8668629) q[0];
sx q[0];
rz(-0.99039536) q[0];
sx q[0];
rz(-0.6443364) q[0];
rz(-pi) q[1];
rz(-1.5526505) q[2];
sx q[2];
rz(-2.7108253) q[2];
sx q[2];
rz(0.22242966) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2542418) q[1];
sx q[1];
rz(-0.72853959) q[1];
sx q[1];
rz(-0.90708797) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(0.90676701) q[0];
rz(0.81470195) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(-1.9168568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46366102) q[0];
sx q[0];
rz(-1.3555129) q[0];
sx q[0];
rz(-0.93925516) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9903509) q[2];
sx q[2];
rz(-0.47861368) q[2];
sx q[2];
rz(1.4065557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21056255) q[1];
sx q[1];
rz(-1.5157053) q[1];
sx q[1];
rz(0.964826) q[1];
rz(-0.25572689) q[3];
sx q[3];
rz(-2.7284107) q[3];
sx q[3];
rz(0.43858389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99888745) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619103) q[0];
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
rz(0.66877767) q[0];
sx q[0];
rz(-0.46572177) q[0];
sx q[0];
rz(1.8229683) q[0];
rz(-pi) q[1];
rz(2.8220196) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(3.0802397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62774819) q[1];
sx q[1];
rz(-0.57250896) q[1];
sx q[1];
rz(1.7379012) q[1];
x q[2];
rz(3.0573781) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(2.6654411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3646399) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(0.022162612) q[2];
rz(2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(2.7409592) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3547524) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(-1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(-2.2498806) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4211593) q[0];
sx q[0];
rz(-1.4132858) q[0];
sx q[0];
rz(-0.79173761) q[0];
rz(-pi) q[1];
rz(-2.5052091) q[2];
sx q[2];
rz(-0.6837662) q[2];
sx q[2];
rz(-0.46905876) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8430431) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(0.72960735) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1710839) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(-1.8613929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(0.78424224) q[2];
rz(0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(2.8083535) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(1.3649712) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7576335) q[0];
sx q[0];
rz(-0.74476349) q[0];
sx q[0];
rz(-0.061116771) q[0];
x q[1];
rz(-1.3806254) q[2];
sx q[2];
rz(-0.41751465) q[2];
sx q[2];
rz(1.6585569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.078799876) q[1];
sx q[1];
rz(-1.2275057) q[1];
sx q[1];
rz(-2.417056) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8678719) q[3];
sx q[3];
rz(-1.3798514) q[3];
sx q[3];
rz(1.2862658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0329131) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(0.03406295) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(1.4779133) q[0];
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
rz(1.6598845) q[0];
sx q[0];
rz(-1.1904926) q[0];
sx q[0];
rz(2.7754521) q[0];
x q[1];
rz(0.41843157) q[2];
sx q[2];
rz(-1.1088587) q[2];
sx q[2];
rz(0.07428169) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8755175) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(-0.43550272) q[1];
rz(-pi) q[2];
rz(-3.1413583) q[3];
sx q[3];
rz(-1.8929385) q[3];
sx q[3];
rz(2.5664267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.39682) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(0.33404074) q[2];
sx q[2];
rz(-2.3636849) q[2];
sx q[2];
rz(-1.452527) q[2];
rz(0.35477521) q[3];
sx q[3];
rz(-2.1003902) q[3];
sx q[3];
rz(-1.0082705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
