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
rz(6.364967) q[0];
sx q[0];
rz(9.9262417) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(3.5377503) q[1];
sx q[1];
rz(9.1023394) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0599521) q[0];
sx q[0];
rz(-1.8917221) q[0];
sx q[0];
rz(3.0517464) q[0];
rz(-pi) q[1];
rz(0.85876043) q[2];
sx q[2];
rz(-0.81958629) q[2];
sx q[2];
rz(2.8263) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0281801) q[1];
sx q[1];
rz(-0.98787381) q[1];
sx q[1];
rz(-2.7229573) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93482165) q[3];
sx q[3];
rz(-1.9536195) q[3];
sx q[3];
rz(2.8179907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(2.5855529) q[2];
rz(-0.83267823) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-2.9843176) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(3.0325586) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5599247) q[0];
sx q[0];
rz(-1.7413057) q[0];
sx q[0];
rz(1.8018363) q[0];
rz(-pi) q[1];
rz(1.7686339) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(-0.82222647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1034531) q[1];
sx q[1];
rz(-2.1557169) q[1];
sx q[1];
rz(-1.000688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7440967) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(-1.3483931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2382425) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(-1.2634574) q[2];
rz(2.8144828) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(1.3431312) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-2.6599191) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96268481) q[0];
sx q[0];
rz(-0.6324397) q[0];
sx q[0];
rz(2.2326367) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9573612) q[2];
sx q[2];
rz(-1.295225) q[2];
sx q[2];
rz(-0.36188175) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.958365) q[1];
sx q[1];
rz(-2.1124358) q[1];
sx q[1];
rz(0.72802131) q[1];
rz(-pi) q[2];
rz(-2.5211469) q[3];
sx q[3];
rz(-0.9471604) q[3];
sx q[3];
rz(0.40943957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(-0.49989191) q[2];
rz(-0.56097427) q[3];
sx q[3];
rz(-1.2596954) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9445779) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(-1.5532956) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.8968556) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9024076) q[0];
sx q[0];
rz(-0.33948487) q[0];
sx q[0];
rz(0.0048588077) q[0];
rz(-1.7719901) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(-0.97845562) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9089531) q[1];
sx q[1];
rz(-2.5563572) q[1];
sx q[1];
rz(1.1802243) q[1];
rz(-0.57085412) q[3];
sx q[3];
rz(-1.3802765) q[3];
sx q[3];
rz(-1.1754456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84919471) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(1.1506895) q[2];
rz(1.4771279) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(-2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(-1.6385471) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8836356) q[0];
sx q[0];
rz(-2.1501318) q[0];
sx q[0];
rz(0.15276076) q[0];
x q[1];
rz(-1.4625174) q[2];
sx q[2];
rz(-1.4881954) q[2];
sx q[2];
rz(2.8464707) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4245783) q[1];
sx q[1];
rz(-1.4687521) q[1];
sx q[1];
rz(1.8841519) q[1];
rz(-pi) q[2];
rz(1.2522069) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(-2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5082671) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(-1.903669) q[2];
rz(2.0189019) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(0.56089127) q[1];
sx q[1];
rz(-1.2979049) q[1];
sx q[1];
rz(-0.7985324) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.417056) q[0];
sx q[0];
rz(-0.31053156) q[0];
sx q[0];
rz(1.0945555) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49662923) q[2];
sx q[2];
rz(-1.8804669) q[2];
sx q[2];
rz(-1.3336099) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9763899) q[1];
sx q[1];
rz(-0.573728) q[1];
sx q[1];
rz(1.6673253) q[1];
rz(2.1875601) q[3];
sx q[3];
rz(-0.87152374) q[3];
sx q[3];
rz(-1.9642533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(2.7748761) q[2];
rz(-1.8803053) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(0.19730332) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(-2.6775449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34940091) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(2.6130136) q[0];
x q[1];
rz(-1.6527962) q[2];
sx q[2];
rz(-2.900015) q[2];
sx q[2];
rz(1.2496003) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0276427) q[1];
sx q[1];
rz(-1.0877891) q[1];
sx q[1];
rz(2.7999858) q[1];
x q[2];
rz(2.9052827) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(1.1856273) q[3];
sx q[3];
rz(-1.6262755) q[3];
sx q[3];
rz(0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-0.19514062) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(-3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.4415178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38948108) q[0];
sx q[0];
rz(-1.6342667) q[0];
sx q[0];
rz(-1.585929) q[0];
rz(-pi) q[1];
rz(1.6091789) q[2];
sx q[2];
rz(-0.95853171) q[2];
sx q[2];
rz(1.4505163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.86130202) q[1];
sx q[1];
rz(-2.960627) q[1];
sx q[1];
rz(0.7678395) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12318792) q[3];
sx q[3];
rz(-0.33629575) q[3];
sx q[3];
rz(1.5598701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7609693) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(-2.8245068) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(-2.1549966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50842677) q[0];
sx q[0];
rz(-2.0619443) q[0];
sx q[0];
rz(1.3915865) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87395845) q[2];
sx q[2];
rz(-1.9947589) q[2];
sx q[2];
rz(1.4615814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4793195) q[1];
sx q[1];
rz(-2.5179177) q[1];
sx q[1];
rz(-0.24346607) q[1];
x q[2];
rz(-1.9916233) q[3];
sx q[3];
rz(-1.7289366) q[3];
sx q[3];
rz(-1.4827673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(2.731936) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(-2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7423994) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(-1.7364527) q[0];
rz(0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.7260889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6380499) q[0];
sx q[0];
rz(-1.8068131) q[0];
sx q[0];
rz(-0.34061265) q[0];
rz(-pi) q[1];
rz(2.892898) q[2];
sx q[2];
rz(-1.9242312) q[2];
sx q[2];
rz(2.8949708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6092097) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(1.7864368) q[1];
rz(-pi) q[2];
rz(2.9726082) q[3];
sx q[3];
rz(-1.0645234) q[3];
sx q[3];
rz(2.3224725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71904174) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-3.0174875) q[2];
rz(0.96578807) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(-2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.538095) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(2.8812257) q[2];
sx q[2];
rz(-1.7454864) q[2];
sx q[2];
rz(-2.7228552) q[2];
rz(-2.9028805) q[3];
sx q[3];
rz(-2.5720027) q[3];
sx q[3];
rz(-3.0726074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
