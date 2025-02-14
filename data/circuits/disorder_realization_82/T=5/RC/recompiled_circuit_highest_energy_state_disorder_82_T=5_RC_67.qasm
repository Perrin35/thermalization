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
rz(1.0225811) q[0];
sx q[0];
rz(5.2978088) q[0];
sx q[0];
rz(10.537416) q[0];
rz(0.39482173) q[1];
sx q[1];
rz(-1.1345175) q[1];
sx q[1];
rz(-1.1910103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1853088) q[0];
sx q[0];
rz(-0.6260159) q[0];
sx q[0];
rz(-0.52395384) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.27287) q[2];
sx q[2];
rz(-2.7162144) q[2];
sx q[2];
rz(-0.59051248) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8419876) q[1];
sx q[1];
rz(-1.5041495) q[1];
sx q[1];
rz(-1.503809) q[1];
x q[2];
rz(-0.43893473) q[3];
sx q[3];
rz(-1.1263444) q[3];
sx q[3];
rz(1.4107454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4224008) q[2];
sx q[2];
rz(-2.9802608) q[2];
sx q[2];
rz(-2.7154229) q[2];
rz(2.4424477) q[3];
sx q[3];
rz(-1.9482502) q[3];
sx q[3];
rz(-2.7307811) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8461269) q[0];
sx q[0];
rz(-2.182425) q[0];
sx q[0];
rz(2.3058983) q[0];
rz(-0.524638) q[1];
sx q[1];
rz(-1.6181889) q[1];
sx q[1];
rz(-1.5118648) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98868167) q[0];
sx q[0];
rz(-0.67711867) q[0];
sx q[0];
rz(-1.0557014) q[0];
rz(-2.2316547) q[2];
sx q[2];
rz(-0.99915394) q[2];
sx q[2];
rz(0.89905587) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0084739) q[1];
sx q[1];
rz(-0.50280675) q[1];
sx q[1];
rz(2.7333906) q[1];
rz(2.5901661) q[3];
sx q[3];
rz(-1.2987483) q[3];
sx q[3];
rz(0.37032933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52874804) q[2];
sx q[2];
rz(-1.5809487) q[2];
sx q[2];
rz(-0.35231248) q[2];
rz(2.6490372) q[3];
sx q[3];
rz(-1.1198606) q[3];
sx q[3];
rz(-0.39920863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6414129) q[0];
sx q[0];
rz(-0.68920511) q[0];
sx q[0];
rz(2.8885762) q[0];
rz(-0.5223271) q[1];
sx q[1];
rz(-0.42799196) q[1];
sx q[1];
rz(-0.36756757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2151486) q[0];
sx q[0];
rz(-0.068000168) q[0];
sx q[0];
rz(1.1626194) q[0];
rz(-pi) q[1];
rz(2.1001294) q[2];
sx q[2];
rz(-0.69165889) q[2];
sx q[2];
rz(-1.4931607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2985017) q[1];
sx q[1];
rz(-1.6502145) q[1];
sx q[1];
rz(-1.3629713) q[1];
rz(1.713836) q[3];
sx q[3];
rz(-1.4554001) q[3];
sx q[3];
rz(0.6108329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0005409) q[2];
sx q[2];
rz(-0.29033574) q[2];
sx q[2];
rz(-0.42624897) q[2];
rz(-0.78921562) q[3];
sx q[3];
rz(-2.0246678) q[3];
sx q[3];
rz(-1.5392019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5696647) q[0];
sx q[0];
rz(-2.8043788) q[0];
sx q[0];
rz(-2.9544882) q[0];
rz(-1.7722173) q[1];
sx q[1];
rz(-1.2792842) q[1];
sx q[1];
rz(-0.61417907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2046004) q[0];
sx q[0];
rz(-1.0796756) q[0];
sx q[0];
rz(3.0762927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34030045) q[2];
sx q[2];
rz(-2.0893761) q[2];
sx q[2];
rz(-0.22835635) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5803197) q[1];
sx q[1];
rz(-2.4497708) q[1];
sx q[1];
rz(-2.0136334) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37482609) q[3];
sx q[3];
rz(-1.4901461) q[3];
sx q[3];
rz(-1.4750236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6807009) q[2];
sx q[2];
rz(-2.2368175) q[2];
sx q[2];
rz(1.2659849) q[2];
rz(2.8547817) q[3];
sx q[3];
rz(-2.3293142) q[3];
sx q[3];
rz(0.12107818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97499651) q[0];
sx q[0];
rz(-2.2280405) q[0];
sx q[0];
rz(1.9972557) q[0];
rz(2.3727349) q[1];
sx q[1];
rz(-0.91173333) q[1];
sx q[1];
rz(1.2164046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.394687) q[0];
sx q[0];
rz(-1.1228859) q[0];
sx q[0];
rz(-0.42952092) q[0];
rz(-pi) q[1];
rz(0.94444176) q[2];
sx q[2];
rz(-2.9778096) q[2];
sx q[2];
rz(0.46609571) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7256355) q[1];
sx q[1];
rz(-2.8548988) q[1];
sx q[1];
rz(1.8724724) q[1];
rz(-pi) q[2];
rz(0.36190108) q[3];
sx q[3];
rz(-1.6503157) q[3];
sx q[3];
rz(-1.9504296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6151109) q[2];
sx q[2];
rz(-0.59978008) q[2];
sx q[2];
rz(-2.6748924) q[2];
rz(-0.99272054) q[3];
sx q[3];
rz(-1.4027184) q[3];
sx q[3];
rz(1.7464975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863131) q[0];
sx q[0];
rz(-1.7236973) q[0];
sx q[0];
rz(0.41743761) q[0];
rz(-2.7159269) q[1];
sx q[1];
rz(-2.3047431) q[1];
sx q[1];
rz(-2.4136037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6052743) q[0];
sx q[0];
rz(-1.124255) q[0];
sx q[0];
rz(-0.98636084) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9404782) q[2];
sx q[2];
rz(-2.3916247) q[2];
sx q[2];
rz(0.88378831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5983819) q[1];
sx q[1];
rz(-1.9459459) q[1];
sx q[1];
rz(0.30483706) q[1];
rz(-2.1393439) q[3];
sx q[3];
rz(-1.2268775) q[3];
sx q[3];
rz(-3.0483766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.723168) q[2];
sx q[2];
rz(-1.6614513) q[2];
sx q[2];
rz(-1.1600912) q[2];
rz(-2.5892819) q[3];
sx q[3];
rz(-0.6902802) q[3];
sx q[3];
rz(2.6753329) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1138678) q[0];
sx q[0];
rz(-0.92806569) q[0];
sx q[0];
rz(3.0357251) q[0];
rz(-1.8766807) q[1];
sx q[1];
rz(-2.5906339) q[1];
sx q[1];
rz(1.6884621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3336847) q[0];
sx q[0];
rz(-0.8561058) q[0];
sx q[0];
rz(2.7175275) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12281863) q[2];
sx q[2];
rz(-1.747532) q[2];
sx q[2];
rz(-2.3251102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4857551) q[1];
sx q[1];
rz(-1.7110516) q[1];
sx q[1];
rz(1.2698238) q[1];
rz(-pi) q[2];
rz(-2.9960521) q[3];
sx q[3];
rz(-2.2284433) q[3];
sx q[3];
rz(-1.6931134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.037131) q[2];
sx q[2];
rz(-2.0115325) q[2];
sx q[2];
rz(1.5360443) q[2];
rz(-1.7724841) q[3];
sx q[3];
rz(-0.57309279) q[3];
sx q[3];
rz(-1.3376741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8417514) q[0];
sx q[0];
rz(-1.890269) q[0];
sx q[0];
rz(-1.1773671) q[0];
rz(1.8384701) q[1];
sx q[1];
rz(-1.6603371) q[1];
sx q[1];
rz(-0.67828137) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2934351) q[0];
sx q[0];
rz(-0.57852902) q[0];
sx q[0];
rz(2.047477) q[0];
rz(-pi) q[1];
rz(0.2770284) q[2];
sx q[2];
rz(-1.4619383) q[2];
sx q[2];
rz(-1.9831374) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.57215124) q[1];
sx q[1];
rz(-1.2673693) q[1];
sx q[1];
rz(-1.3511168) q[1];
rz(2.5362064) q[3];
sx q[3];
rz(-1.6035214) q[3];
sx q[3];
rz(0.42057366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9699041) q[2];
sx q[2];
rz(-2.6320612) q[2];
sx q[2];
rz(-0.92863885) q[2];
rz(1.5019794) q[3];
sx q[3];
rz(-2.00311) q[3];
sx q[3];
rz(1.5665215) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.11012) q[0];
sx q[0];
rz(-2.0125772) q[0];
sx q[0];
rz(1.8054777) q[0];
rz(-2.8307092) q[1];
sx q[1];
rz(-1.5300405) q[1];
sx q[1];
rz(1.8119887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7696849) q[0];
sx q[0];
rz(-1.2331672) q[0];
sx q[0];
rz(2.4680074) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63116341) q[2];
sx q[2];
rz(-2.1710325) q[2];
sx q[2];
rz(-1.7886666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8927129) q[1];
sx q[1];
rz(-1.5096438) q[1];
sx q[1];
rz(0.11492954) q[1];
rz(-pi) q[2];
rz(-0.90985591) q[3];
sx q[3];
rz(-1.6279216) q[3];
sx q[3];
rz(-1.2965681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0335994) q[2];
sx q[2];
rz(-0.42069837) q[2];
sx q[2];
rz(-1.1619953) q[2];
rz(-1.1437931) q[3];
sx q[3];
rz(-1.9547209) q[3];
sx q[3];
rz(-2.7403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4556295) q[0];
sx q[0];
rz(-2.3271053) q[0];
sx q[0];
rz(-0.90675768) q[0];
rz(-2.7470159) q[1];
sx q[1];
rz(-2.5386609) q[1];
sx q[1];
rz(-2.0049863) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2897091) q[0];
sx q[0];
rz(-2.7628081) q[0];
sx q[0];
rz(2.3251371) q[0];
rz(-pi) q[1];
x q[1];
rz(0.074176057) q[2];
sx q[2];
rz(-1.2741421) q[2];
sx q[2];
rz(-2.3919472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1882595) q[1];
sx q[1];
rz(-1.7112977) q[1];
sx q[1];
rz(-0.85660117) q[1];
rz(-pi) q[2];
rz(-1.8618989) q[3];
sx q[3];
rz(-1.2754692) q[3];
sx q[3];
rz(2.3307523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11791983) q[2];
sx q[2];
rz(-1.6912141) q[2];
sx q[2];
rz(-2.0898021) q[2];
rz(-0.037083179) q[3];
sx q[3];
rz(-0.6784234) q[3];
sx q[3];
rz(-1.3858494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7294075) q[0];
sx q[0];
rz(-0.89898983) q[0];
sx q[0];
rz(0.40869024) q[0];
rz(-2.9858934) q[1];
sx q[1];
rz(-1.0380048) q[1];
sx q[1];
rz(0.17539594) q[1];
rz(-0.88798203) q[2];
sx q[2];
rz(-1.8941634) q[2];
sx q[2];
rz(-1.8695199) q[2];
rz(1.8694405) q[3];
sx q[3];
rz(-2.1437672) q[3];
sx q[3];
rz(-1.7199316) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
