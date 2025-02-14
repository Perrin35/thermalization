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
rz(-0.24124423) q[0];
sx q[0];
rz(-0.21114199) q[0];
sx q[0];
rz(-2.733732) q[0];
rz(1.3762228) q[1];
sx q[1];
rz(-1.9150182) q[1];
sx q[1];
rz(3.0768375) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0000251) q[0];
sx q[0];
rz(-1.5760826) q[0];
sx q[0];
rz(-1.686211) q[0];
rz(-pi) q[1];
rz(1.1062293) q[2];
sx q[2];
rz(-0.63815123) q[2];
sx q[2];
rz(-1.3582071) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8246461) q[1];
sx q[1];
rz(-1.7120518) q[1];
sx q[1];
rz(-1.5819751) q[1];
rz(2.9969941) q[3];
sx q[3];
rz(-0.62718348) q[3];
sx q[3];
rz(-2.7759564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.18225081) q[2];
sx q[2];
rz(-2.0472517) q[2];
sx q[2];
rz(-1.135929) q[2];
rz(-0.7302537) q[3];
sx q[3];
rz(-1.7466931) q[3];
sx q[3];
rz(1.5203389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.092904329) q[0];
sx q[0];
rz(-2.1850259) q[0];
sx q[0];
rz(-0.055305716) q[0];
rz(2.7703908) q[1];
sx q[1];
rz(-0.35938811) q[1];
sx q[1];
rz(0.41195437) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5980976) q[0];
sx q[0];
rz(-0.5306705) q[0];
sx q[0];
rz(0.12904386) q[0];
x q[1];
rz(-2.8438004) q[2];
sx q[2];
rz(-1.2504559) q[2];
sx q[2];
rz(3.0646119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.69393278) q[1];
sx q[1];
rz(-0.74697633) q[1];
sx q[1];
rz(-0.20836094) q[1];
rz(1.0874416) q[3];
sx q[3];
rz(-2.4341741) q[3];
sx q[3];
rz(-0.055295769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61341316) q[2];
sx q[2];
rz(-1.6772905) q[2];
sx q[2];
rz(2.2589653) q[2];
rz(-1.708185) q[3];
sx q[3];
rz(-1.2127168) q[3];
sx q[3];
rz(2.9679756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16647896) q[0];
sx q[0];
rz(-0.011703165) q[0];
sx q[0];
rz(-2.574918) q[0];
rz(-2.7239679) q[1];
sx q[1];
rz(-1.5387225) q[1];
sx q[1];
rz(-1.8142987) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1786363) q[0];
sx q[0];
rz(-1.6554852) q[0];
sx q[0];
rz(0.87015193) q[0];
x q[1];
rz(0.94608091) q[2];
sx q[2];
rz(-1.6951121) q[2];
sx q[2];
rz(0.92783463) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3440888) q[1];
sx q[1];
rz(-0.3080536) q[1];
sx q[1];
rz(-1.5750242) q[1];
rz(-2.684115) q[3];
sx q[3];
rz(-2.2443534) q[3];
sx q[3];
rz(2.9411773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7495482) q[2];
sx q[2];
rz(-2.0508524) q[2];
sx q[2];
rz(-0.90847477) q[2];
rz(-0.79489094) q[3];
sx q[3];
rz(-2.0386233) q[3];
sx q[3];
rz(-3.095043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6598776) q[0];
sx q[0];
rz(-2.2292723) q[0];
sx q[0];
rz(0.6657486) q[0];
rz(-0.84716973) q[1];
sx q[1];
rz(-0.24996346) q[1];
sx q[1];
rz(-0.4096823) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2100135) q[0];
sx q[0];
rz(-2.0635475) q[0];
sx q[0];
rz(1.7442305) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3607765) q[2];
sx q[2];
rz(-2.2998126) q[2];
sx q[2];
rz(1.8052124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0627111) q[1];
sx q[1];
rz(-2.1797175) q[1];
sx q[1];
rz(-2.8669861) q[1];
x q[2];
rz(1.5104896) q[3];
sx q[3];
rz(-2.6968287) q[3];
sx q[3];
rz(0.30032762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72982558) q[2];
sx q[2];
rz(-2.4335786) q[2];
sx q[2];
rz(-1.7626308) q[2];
rz(-1.1369368) q[3];
sx q[3];
rz(-1.7863019) q[3];
sx q[3];
rz(0.8374477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1788504) q[0];
sx q[0];
rz(-2.0806291) q[0];
sx q[0];
rz(0.41819292) q[0];
rz(-1.7182982) q[1];
sx q[1];
rz(-1.1885252) q[1];
sx q[1];
rz(-1.4994015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20883501) q[0];
sx q[0];
rz(-1.3262188) q[0];
sx q[0];
rz(2.6840058) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8756494) q[2];
sx q[2];
rz(-1.8902167) q[2];
sx q[2];
rz(-1.4679421) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9019431) q[1];
sx q[1];
rz(-0.22615438) q[1];
sx q[1];
rz(-1.8312901) q[1];
rz(-pi) q[2];
rz(3.0762227) q[3];
sx q[3];
rz(-1.4413068) q[3];
sx q[3];
rz(-2.0355899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2749918) q[2];
sx q[2];
rz(-1.9860705) q[2];
sx q[2];
rz(-2.9633813) q[2];
rz(-0.37402672) q[3];
sx q[3];
rz(-0.97291294) q[3];
sx q[3];
rz(1.9988352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.48443925) q[0];
sx q[0];
rz(-1.8171808) q[0];
sx q[0];
rz(1.4763747) q[0];
rz(0.58552512) q[1];
sx q[1];
rz(-0.89591566) q[1];
sx q[1];
rz(-2.4077328) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39291362) q[0];
sx q[0];
rz(-2.3322912) q[0];
sx q[0];
rz(-0.63776638) q[0];
rz(-0.81686498) q[2];
sx q[2];
rz(-2.7079963) q[2];
sx q[2];
rz(-0.24092937) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8283353) q[1];
sx q[1];
rz(-2.0861107) q[1];
sx q[1];
rz(-0.22586598) q[1];
rz(-pi) q[2];
rz(-2.0229983) q[3];
sx q[3];
rz(-1.7590932) q[3];
sx q[3];
rz(-1.888243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0840941) q[2];
sx q[2];
rz(-1.7666662) q[2];
sx q[2];
rz(-1.6424087) q[2];
rz(1.5205787) q[3];
sx q[3];
rz(-0.64783827) q[3];
sx q[3];
rz(-0.15629855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3868189) q[0];
sx q[0];
rz(-0.91803011) q[0];
sx q[0];
rz(-1.3936438) q[0];
rz(0.58760324) q[1];
sx q[1];
rz(-1.5005451) q[1];
sx q[1];
rz(0.29744068) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5818243) q[0];
sx q[0];
rz(-2.3399669) q[0];
sx q[0];
rz(-1.0645435) q[0];
x q[1];
rz(1.8413421) q[2];
sx q[2];
rz(-1.5890117) q[2];
sx q[2];
rz(-1.2684938) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.45876452) q[1];
sx q[1];
rz(-1.5611083) q[1];
sx q[1];
rz(-0.89114706) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92771156) q[3];
sx q[3];
rz(-1.7684019) q[3];
sx q[3];
rz(1.9324945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8224767) q[2];
sx q[2];
rz(-1.3603223) q[2];
sx q[2];
rz(2.759867) q[2];
rz(0.60025674) q[3];
sx q[3];
rz(-2.468942) q[3];
sx q[3];
rz(0.265358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6622019) q[0];
sx q[0];
rz(-1.5150161) q[0];
sx q[0];
rz(1.8040682) q[0];
rz(3.1321101) q[1];
sx q[1];
rz(-0.9597221) q[1];
sx q[1];
rz(1.7705852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66390304) q[0];
sx q[0];
rz(-2.7743917) q[0];
sx q[0];
rz(-0.3911256) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1751733) q[2];
sx q[2];
rz(-2.3766209) q[2];
sx q[2];
rz(-1.6993831) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0095894) q[1];
sx q[1];
rz(-1.6958964) q[1];
sx q[1];
rz(-2.3522207) q[1];
rz(2.9389006) q[3];
sx q[3];
rz(-0.42169844) q[3];
sx q[3];
rz(-1.7360753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.81801549) q[2];
sx q[2];
rz(-2.24021) q[2];
sx q[2];
rz(-2.9991007) q[2];
rz(-0.88322181) q[3];
sx q[3];
rz(-1.3775237) q[3];
sx q[3];
rz(1.6787329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0505117) q[0];
sx q[0];
rz(-3.0497146) q[0];
sx q[0];
rz(1.2485414) q[0];
rz(-2.7342791) q[1];
sx q[1];
rz(-1.7946323) q[1];
sx q[1];
rz(1.1936845) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7418967) q[0];
sx q[0];
rz(-2.7177161) q[0];
sx q[0];
rz(0.90452832) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2676554) q[2];
sx q[2];
rz(-2.1997582) q[2];
sx q[2];
rz(0.12128092) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15577573) q[1];
sx q[1];
rz(-2.0627506) q[1];
sx q[1];
rz(-1.1883177) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8305998) q[3];
sx q[3];
rz(-1.570829) q[3];
sx q[3];
rz(2.3268229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0672037) q[2];
sx q[2];
rz(-2.1838102) q[2];
sx q[2];
rz(2.0286782) q[2];
rz(3.0985966) q[3];
sx q[3];
rz(-1.4837416) q[3];
sx q[3];
rz(-2.903741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.064903108) q[0];
sx q[0];
rz(-1.7273128) q[0];
sx q[0];
rz(2.9458556) q[0];
rz(-1.2204569) q[1];
sx q[1];
rz(-0.38206044) q[1];
sx q[1];
rz(-0.97741309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5223434) q[0];
sx q[0];
rz(-0.27740208) q[0];
sx q[0];
rz(2.7707556) q[0];
x q[1];
rz(1.270411) q[2];
sx q[2];
rz(-1.7160048) q[2];
sx q[2];
rz(1.5086439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5557438) q[1];
sx q[1];
rz(-0.80744707) q[1];
sx q[1];
rz(0.84485742) q[1];
x q[2];
rz(-1.8582088) q[3];
sx q[3];
rz(-1.1927146) q[3];
sx q[3];
rz(1.9668129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0105373) q[2];
sx q[2];
rz(-2.8752893) q[2];
sx q[2];
rz(1.4492501) q[2];
rz(2.2629755) q[3];
sx q[3];
rz(-1.0613469) q[3];
sx q[3];
rz(-1.7884458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286248) q[0];
sx q[0];
rz(-1.9504564) q[0];
sx q[0];
rz(1.7497028) q[0];
rz(1.3799008) q[1];
sx q[1];
rz(-1.6086144) q[1];
sx q[1];
rz(-0.10133941) q[1];
rz(2.352667) q[2];
sx q[2];
rz(-0.25392214) q[2];
sx q[2];
rz(1.8352652) q[2];
rz(1.8203406) q[3];
sx q[3];
rz(-0.93946447) q[3];
sx q[3];
rz(2.3238121) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
