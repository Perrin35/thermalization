OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(-0.50146377) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(-2.745435) q[1];
sx q[1];
rz(-0.3224386) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0599521) q[0];
sx q[0];
rz(-1.2498706) q[0];
sx q[0];
rz(3.0517464) q[0];
rz(2.2519977) q[2];
sx q[2];
rz(-1.0729562) q[2];
sx q[2];
rz(1.3537458) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8391708) q[1];
sx q[1];
rz(-1.9170554) q[1];
sx q[1];
rz(0.94567169) q[1];
x q[2];
rz(-2.206771) q[3];
sx q[3];
rz(-1.1879731) q[3];
sx q[3];
rz(0.32360199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(2.5855529) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(-0.94579831) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(-3.0325586) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5599247) q[0];
sx q[0];
rz(-1.7413057) q[0];
sx q[0];
rz(-1.3397564) q[0];
rz(-pi) q[1];
rz(-0.63983812) q[2];
sx q[2];
rz(-2.8176762) q[2];
sx q[2];
rz(1.4878291) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1034531) q[1];
sx q[1];
rz(-2.1557169) q[1];
sx q[1];
rz(2.1409047) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7440967) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(1.7931995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(-1.2634574) q[2];
rz(-0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-2.6599191) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7307229) q[0];
sx q[0];
rz(-1.0856837) q[0];
sx q[0];
rz(2.7184125) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0273676) q[2];
sx q[2];
rz(-0.66514981) q[2];
sx q[2];
rz(0.84012023) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8183221) q[1];
sx q[1];
rz(-2.1774877) q[1];
sx q[1];
rz(0.89241772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2504911) q[3];
sx q[3];
rz(-0.84935969) q[3];
sx q[3];
rz(-0.47618714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(-2.6417007) q[2];
rz(-0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.19701476) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(0.552185) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(-1.2447371) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23403215) q[0];
sx q[0];
rz(-1.2313156) q[0];
sx q[0];
rz(-1.5690804) q[0];
rz(-2.142698) q[2];
sx q[2];
rz(-2.9036387) q[2];
sx q[2];
rz(1.1513125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4501805) q[1];
sx q[1];
rz(-2.1069063) q[1];
sx q[1];
rz(0.24713534) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3454868) q[3];
sx q[3];
rz(-1.011519) q[3];
sx q[3];
rz(2.6252281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2923979) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(1.9909031) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0634336) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.6385471) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016184729) q[0];
sx q[0];
rz(-0.59690079) q[0];
sx q[0];
rz(1.7993268) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6790752) q[2];
sx q[2];
rz(-1.6533972) q[2];
sx q[2];
rz(-0.29512197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.45080966) q[1];
sx q[1];
rz(-2.8125617) q[1];
sx q[1];
rz(1.8915218) q[1];
rz(-pi) q[2];
rz(-0.50932192) q[3];
sx q[3];
rz(-1.29042) q[3];
sx q[3];
rz(-2.544739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6333255) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(-1.2379237) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5313107) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-0.26671985) q[0];
rz(-2.5807014) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(0.7985324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0621588) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(0.14607231) q[0];
rz(1.919826) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(-2.7407321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8171463) q[1];
sx q[1];
rz(-1.518461) q[1];
sx q[1];
rz(-0.99919341) q[1];
rz(-pi) q[2];
rz(0.800662) q[3];
sx q[3];
rz(-1.1122397) q[3];
sx q[3];
rz(0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9810527) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(-2.7748761) q[2];
rz(-1.2612873) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(-3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40903184) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(0.60638705) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(0.46404776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2460829) q[0];
sx q[0];
rz(-1.1469139) q[0];
sx q[0];
rz(2.258582) q[0];
rz(-1.8115933) q[2];
sx q[2];
rz(-1.5511998) q[2];
sx q[2];
rz(0.400825) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5210919) q[1];
sx q[1];
rz(-1.2695841) q[1];
sx q[1];
rz(2.0786933) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4887605) q[3];
sx q[3];
rz(-1.335252) q[3];
sx q[3];
rz(1.9627278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93418926) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(2.8835473) q[2];
rz(1.1856273) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19514062) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(-2.7602957) q[0];
rz(-3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(1.7000748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7521116) q[0];
sx q[0];
rz(-1.6342667) q[0];
sx q[0];
rz(-1.585929) q[0];
rz(2.5289815) q[2];
sx q[2];
rz(-1.5393886) q[2];
sx q[2];
rz(-0.14234662) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2802906) q[1];
sx q[1];
rz(-2.960627) q[1];
sx q[1];
rz(0.7678395) q[1];
rz(-0.33393319) q[3];
sx q[3];
rz(-1.6113558) q[3];
sx q[3];
rz(-0.12727748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1429446) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(0.22658919) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(-0.8297689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7609693) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(1.3990078) q[0];
rz(0.31708583) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6331659) q[0];
sx q[0];
rz(-2.0619443) q[0];
sx q[0];
rz(1.7500061) q[0];
x q[1];
rz(2.1836906) q[2];
sx q[2];
rz(-0.79682486) q[2];
sx q[2];
rz(0.56607841) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7762737) q[1];
sx q[1];
rz(-2.1734108) q[1];
sx q[1];
rz(1.7425294) q[1];
rz(-pi) q[2];
rz(2.968623) q[3];
sx q[3];
rz(-1.9860455) q[3];
sx q[3];
rz(-0.017661974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2150779) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(2.731936) q[2];
rz(-2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(-1.4051399) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(1.7260889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1499407) q[0];
sx q[0];
rz(-1.9016001) q[0];
sx q[0];
rz(1.8206235) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24869463) q[2];
sx q[2];
rz(-1.9242312) q[2];
sx q[2];
rz(-0.24662185) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.532383) q[1];
sx q[1];
rz(-1.7222952) q[1];
sx q[1];
rz(1.3551559) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9726082) q[3];
sx q[3];
rz(-1.0645234) q[3];
sx q[3];
rz(-0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(-2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60349764) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(-1.7514501) q[2];
sx q[2];
rz(-1.3144819) q[2];
sx q[2];
rz(-1.1983295) q[2];
rz(-0.55660558) q[3];
sx q[3];
rz(-1.442933) q[3];
sx q[3];
rz(1.8419151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];