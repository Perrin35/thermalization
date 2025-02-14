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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(-2.4315779) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(-1.8097872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84146554) q[0];
sx q[0];
rz(-1.3464768) q[0];
sx q[0];
rz(-2.8594144) q[0];
rz(0.34320404) q[2];
sx q[2];
rz(-1.3319974) q[2];
sx q[2];
rz(-2.3607139) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.097447473) q[1];
sx q[1];
rz(-2.4821979) q[1];
sx q[1];
rz(-1.9800746) q[1];
rz(-2.2043742) q[3];
sx q[3];
rz(-1.9301747) q[3];
sx q[3];
rz(2.0546953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4189202) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(-0.81532064) q[2];
rz(2.3319862) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.078449) q[0];
sx q[0];
rz(-1.0493295) q[0];
sx q[0];
rz(3.0310042) q[0];
rz(-2.1965006) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5623215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5961869) q[0];
sx q[0];
rz(-2.2431787) q[0];
sx q[0];
rz(2.5491712) q[0];
x q[1];
rz(-2.7209372) q[2];
sx q[2];
rz(-1.1569945) q[2];
sx q[2];
rz(2.5977573) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64792127) q[1];
sx q[1];
rz(-2.646138) q[1];
sx q[1];
rz(0.22101553) q[1];
x q[2];
rz(2.3195295) q[3];
sx q[3];
rz(-1.2044319) q[3];
sx q[3];
rz(1.0613393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.14757806) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(2.6914524) q[2];
rz(-2.0077997) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(0.97761893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.5728773) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(-2.8602588) q[0];
rz(-1.7549134) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(0.6001572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2413797) q[0];
sx q[0];
rz(-0.29415392) q[0];
sx q[0];
rz(0.16071975) q[0];
rz(-pi) q[1];
rz(-3.0855028) q[2];
sx q[2];
rz(-0.96470913) q[2];
sx q[2];
rz(-0.99585271) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8474947) q[1];
sx q[1];
rz(-1.0007684) q[1];
sx q[1];
rz(-1.1815726) q[1];
rz(-0.5731606) q[3];
sx q[3];
rz(-1.7148682) q[3];
sx q[3];
rz(-2.746563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0967789) q[2];
sx q[2];
rz(-1.1324977) q[2];
sx q[2];
rz(-1.0463932) q[2];
rz(0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(-0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9699049) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(-2.0890253) q[0];
rz(-1.196208) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(2.722091) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013166817) q[0];
sx q[0];
rz(-1.3072398) q[0];
sx q[0];
rz(0.58851425) q[0];
rz(-pi) q[1];
rz(2.2050836) q[2];
sx q[2];
rz(-1.3872996) q[2];
sx q[2];
rz(-1.8734135) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12246015) q[1];
sx q[1];
rz(-1.9472709) q[1];
sx q[1];
rz(-0.89577641) q[1];
rz(-pi) q[2];
rz(-2.6330248) q[3];
sx q[3];
rz(-1.3150104) q[3];
sx q[3];
rz(-0.05318197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9516051) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(-2.0513963) q[2];
rz(2.6103141) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200318) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.3979727) q[0];
rz(0.13294237) q[1];
sx q[1];
rz(-1.3016737) q[1];
sx q[1];
rz(-3.1403819) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0069790445) q[0];
sx q[0];
rz(-1.9513357) q[0];
sx q[0];
rz(1.5881722) q[0];
rz(-pi) q[1];
rz(-1.0693477) q[2];
sx q[2];
rz(-2.8184888) q[2];
sx q[2];
rz(2.1954775) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7442786) q[1];
sx q[1];
rz(-0.71956735) q[1];
sx q[1];
rz(-0.47081635) q[1];
x q[2];
rz(1.7204956) q[3];
sx q[3];
rz(-1.5970917) q[3];
sx q[3];
rz(0.61248518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2289537) q[2];
sx q[2];
rz(-1.3372083) q[2];
sx q[2];
rz(0.21052989) q[2];
rz(-2.1052965) q[3];
sx q[3];
rz(-0.784289) q[3];
sx q[3];
rz(1.9265296) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482816) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(-2.7358828) q[0];
rz(-1.7871208) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(-0.92232651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94033891) q[0];
sx q[0];
rz(-2.8329599) q[0];
sx q[0];
rz(-1.2135394) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10293691) q[2];
sx q[2];
rz(-2.2349226) q[2];
sx q[2];
rz(0.31598202) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4757918) q[1];
sx q[1];
rz(-2.5195055) q[1];
sx q[1];
rz(1.3740963) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9590553) q[3];
sx q[3];
rz(-1.045533) q[3];
sx q[3];
rz(2.1059259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4313844) q[2];
sx q[2];
rz(-1.9269383) q[2];
sx q[2];
rz(2.7396438) q[2];
rz(1.2881783) q[3];
sx q[3];
rz(-0.86789075) q[3];
sx q[3];
rz(-1.8624381) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61889082) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(-0.064362137) q[0];
rz(-2.3035658) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(-1.3305957) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143599) q[0];
sx q[0];
rz(-1.4454675) q[0];
sx q[0];
rz(2.5663239) q[0];
rz(-pi) q[1];
rz(0.60637252) q[2];
sx q[2];
rz(-1.1506216) q[2];
sx q[2];
rz(-0.5768896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3804662) q[1];
sx q[1];
rz(-2.5674324) q[1];
sx q[1];
rz(-2.7284184) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17223151) q[3];
sx q[3];
rz(-0.96559286) q[3];
sx q[3];
rz(-1.3394522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4055206) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(-0.68354496) q[2];
rz(0.73379597) q[3];
sx q[3];
rz(-1.4132696) q[3];
sx q[3];
rz(0.84764135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34506327) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(-1.1731359) q[0];
rz(-2.7660811) q[1];
sx q[1];
rz(-0.48986062) q[1];
sx q[1];
rz(-1.0367905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4597367) q[0];
sx q[0];
rz(-1.5899183) q[0];
sx q[0];
rz(-3.1179423) q[0];
rz(0.16536819) q[2];
sx q[2];
rz(-1.8435119) q[2];
sx q[2];
rz(-2.337526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4040096) q[1];
sx q[1];
rz(-2.3669457) q[1];
sx q[1];
rz(0.20139106) q[1];
x q[2];
rz(-2.6508207) q[3];
sx q[3];
rz(-1.0939071) q[3];
sx q[3];
rz(-0.27974883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5658687) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(-2.4110528) q[2];
rz(2.9324487) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(-1.8714347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65649477) q[0];
sx q[0];
rz(-2.7182343) q[0];
sx q[0];
rz(0.04743162) q[0];
rz(2.9196396) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(-0.91112959) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.52235) q[0];
sx q[0];
rz(-1.5372835) q[0];
sx q[0];
rz(-1.5977809) q[0];
rz(0.92911559) q[2];
sx q[2];
rz(-1.5826028) q[2];
sx q[2];
rz(0.75071834) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.92420134) q[1];
sx q[1];
rz(-1.5834406) q[1];
sx q[1];
rz(-0.36751698) q[1];
rz(-pi) q[2];
rz(-0.20424517) q[3];
sx q[3];
rz(-1.0139272) q[3];
sx q[3];
rz(2.000252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86965108) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(0.16858777) q[2];
rz(0.21909675) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99437) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(0.45743531) q[0];
rz(1.9153197) q[1];
sx q[1];
rz(-0.33718449) q[1];
sx q[1];
rz(1.7785243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6493312) q[0];
sx q[0];
rz(-2.1773585) q[0];
sx q[0];
rz(1.6442293) q[0];
rz(-2.5117433) q[2];
sx q[2];
rz(-1.7694663) q[2];
sx q[2];
rz(2.3185286) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41238989) q[1];
sx q[1];
rz(-2.7250184) q[1];
sx q[1];
rz(0.42804407) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2533204) q[3];
sx q[3];
rz(-1.5049045) q[3];
sx q[3];
rz(0.94277387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7676131) q[2];
sx q[2];
rz(-2.5869936) q[2];
sx q[2];
rz(-0.45200959) q[2];
rz(0.40397817) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(-2.6192464) q[1];
sx q[1];
rz(-2.6112687) q[1];
sx q[1];
rz(-1.912259) q[1];
rz(-0.73595388) q[2];
sx q[2];
rz(-2.7685168) q[2];
sx q[2];
rz(-1.4442486) q[2];
rz(-0.26351874) q[3];
sx q[3];
rz(-1.7150039) q[3];
sx q[3];
rz(1.1506469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
