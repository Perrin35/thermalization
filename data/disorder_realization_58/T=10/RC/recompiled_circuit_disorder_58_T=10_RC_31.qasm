OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9392202) q[0];
sx q[0];
rz(-0.4063172) q[0];
sx q[0];
rz(0.82011861) q[0];
rz(-0.36110538) q[1];
sx q[1];
rz(0.63280025) q[1];
sx q[1];
rz(11.735698) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93847371) q[0];
sx q[0];
rz(-1.1955368) q[0];
sx q[0];
rz(-1.24036) q[0];
rz(-pi) q[1];
rz(2.4029998) q[2];
sx q[2];
rz(-1.2054218) q[2];
sx q[2];
rz(1.5473168) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4315011) q[1];
sx q[1];
rz(-0.7845062) q[1];
sx q[1];
rz(1.1169408) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44736638) q[3];
sx q[3];
rz(-1.2036185) q[3];
sx q[3];
rz(-0.51608738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7044907) q[2];
sx q[2];
rz(-2.7097242) q[2];
sx q[2];
rz(3.0214156) q[2];
rz(-1.9834571) q[3];
sx q[3];
rz(-1.7430256) q[3];
sx q[3];
rz(-2.2226298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08081089) q[0];
sx q[0];
rz(-1.3588384) q[0];
sx q[0];
rz(2.2297915) q[0];
rz(-2.3520825) q[1];
sx q[1];
rz(-2.1531838) q[1];
sx q[1];
rz(0.3266913) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.806843) q[0];
sx q[0];
rz(-0.93870367) q[0];
sx q[0];
rz(-2.500781) q[0];
x q[1];
rz(-1.8961043) q[2];
sx q[2];
rz(-0.46519687) q[2];
sx q[2];
rz(-2.543769) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5888728) q[1];
sx q[1];
rz(-0.96541222) q[1];
sx q[1];
rz(-2.6436716) q[1];
rz(-pi) q[2];
rz(-0.88300206) q[3];
sx q[3];
rz(-1.6371173) q[3];
sx q[3];
rz(1.6460713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5102753) q[2];
sx q[2];
rz(-0.77284208) q[2];
sx q[2];
rz(1.2871683) q[2];
rz(3.0316947) q[3];
sx q[3];
rz(-1.4108312) q[3];
sx q[3];
rz(1.7597642) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54863769) q[0];
sx q[0];
rz(-2.4023963) q[0];
sx q[0];
rz(0.32989311) q[0];
rz(0.27711162) q[1];
sx q[1];
rz(-1.3169293) q[1];
sx q[1];
rz(1.057391) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7085416) q[0];
sx q[0];
rz(-1.2305224) q[0];
sx q[0];
rz(3.1197381) q[0];
x q[1];
rz(-0.34122841) q[2];
sx q[2];
rz(-1.382302) q[2];
sx q[2];
rz(-0.54268062) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96506572) q[1];
sx q[1];
rz(-0.9170734) q[1];
sx q[1];
rz(2.7214126) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3789165) q[3];
sx q[3];
rz(-1.1550316) q[3];
sx q[3];
rz(-2.2172745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7827591) q[2];
sx q[2];
rz(-2.0262572) q[2];
sx q[2];
rz(-1.3624297) q[2];
rz(0.6247012) q[3];
sx q[3];
rz(-2.0420859) q[3];
sx q[3];
rz(-0.81956285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
rz(1.3574922) q[0];
sx q[0];
rz(-0.52629137) q[0];
sx q[0];
rz(2.6065361) q[0];
rz(1.1401945) q[1];
sx q[1];
rz(-1.3277206) q[1];
sx q[1];
rz(0.16539703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7968984) q[0];
sx q[0];
rz(-1.3860774) q[0];
sx q[0];
rz(0.12634191) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16343127) q[2];
sx q[2];
rz(-1.7497517) q[2];
sx q[2];
rz(0.93726678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4277028) q[1];
sx q[1];
rz(-1.0298924) q[1];
sx q[1];
rz(1.9546024) q[1];
x q[2];
rz(1.9966647) q[3];
sx q[3];
rz(-0.59755675) q[3];
sx q[3];
rz(1.4262517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7455204) q[2];
sx q[2];
rz(-1.3362276) q[2];
sx q[2];
rz(-2.9361434) q[2];
rz(1.127634) q[3];
sx q[3];
rz(-1.9879568) q[3];
sx q[3];
rz(-1.8849461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66185343) q[0];
sx q[0];
rz(-2.2275708) q[0];
sx q[0];
rz(-0.35476312) q[0];
rz(1.154249) q[1];
sx q[1];
rz(-2.216979) q[1];
sx q[1];
rz(2.9096471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7888768) q[0];
sx q[0];
rz(-2.0577601) q[0];
sx q[0];
rz(2.6608174) q[0];
rz(-pi) q[1];
rz(2.4460692) q[2];
sx q[2];
rz(-1.6791653) q[2];
sx q[2];
rz(0.34851375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0143118) q[1];
sx q[1];
rz(-2.8297272) q[1];
sx q[1];
rz(0.77906268) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36851818) q[3];
sx q[3];
rz(-1.5526062) q[3];
sx q[3];
rz(-2.8858678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0014687) q[2];
sx q[2];
rz(-1.8073558) q[2];
sx q[2];
rz(-1.2403963) q[2];
rz(-2.5455348) q[3];
sx q[3];
rz(-1.8362935) q[3];
sx q[3];
rz(-1.8381455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0971138) q[0];
sx q[0];
rz(-3.0713186) q[0];
sx q[0];
rz(-0.20275673) q[0];
rz(0.98908201) q[1];
sx q[1];
rz(-1.6977856) q[1];
sx q[1];
rz(-2.1441377) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0925804) q[0];
sx q[0];
rz(-0.870734) q[0];
sx q[0];
rz(-2.9655365) q[0];
rz(-pi) q[1];
rz(-0.41001292) q[2];
sx q[2];
rz(-2.3668681) q[2];
sx q[2];
rz(1.1020401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7194697) q[1];
sx q[1];
rz(-2.649253) q[1];
sx q[1];
rz(-2.0960397) q[1];
rz(-2.0632083) q[3];
sx q[3];
rz(-0.78740722) q[3];
sx q[3];
rz(1.9099964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22770195) q[2];
sx q[2];
rz(-1.9753549) q[2];
sx q[2];
rz(0.4513936) q[2];
rz(-2.732892) q[3];
sx q[3];
rz(-1.5390076) q[3];
sx q[3];
rz(1.2020948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8354427) q[0];
sx q[0];
rz(-0.4040443) q[0];
sx q[0];
rz(-2.5174482) q[0];
rz(-1.6250601) q[1];
sx q[1];
rz(-2.8846278) q[1];
sx q[1];
rz(-2.5278032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543024) q[0];
sx q[0];
rz(-2.3809732) q[0];
sx q[0];
rz(-1.8421696) q[0];
x q[1];
rz(-1.6348398) q[2];
sx q[2];
rz(-1.2876858) q[2];
sx q[2];
rz(-1.3751021) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2144199) q[1];
sx q[1];
rz(-1.629429) q[1];
sx q[1];
rz(-3.0374132) q[1];
rz(1.1055787) q[3];
sx q[3];
rz(-1.675616) q[3];
sx q[3];
rz(3.1302111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7081786) q[2];
sx q[2];
rz(-2.2257979) q[2];
sx q[2];
rz(1.9667352) q[2];
rz(-0.60837778) q[3];
sx q[3];
rz(-1.404168) q[3];
sx q[3];
rz(1.4233937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2414918) q[0];
sx q[0];
rz(-1.2986203) q[0];
sx q[0];
rz(-2.6532145) q[0];
rz(1.5178559) q[1];
sx q[1];
rz(-1.3987712) q[1];
sx q[1];
rz(0.98446313) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280633) q[0];
sx q[0];
rz(-2.5941879) q[0];
sx q[0];
rz(3.0696966) q[0];
x q[1];
rz(0.77163561) q[2];
sx q[2];
rz(-1.8688335) q[2];
sx q[2];
rz(2.3723797) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24482803) q[1];
sx q[1];
rz(-0.36537376) q[1];
sx q[1];
rz(-3.1246964) q[1];
x q[2];
rz(-0.55926178) q[3];
sx q[3];
rz(-2.5945633) q[3];
sx q[3];
rz(-0.83838851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70665923) q[2];
sx q[2];
rz(-1.443053) q[2];
sx q[2];
rz(0.53517503) q[2];
rz(1.0501856) q[3];
sx q[3];
rz(-1.8015367) q[3];
sx q[3];
rz(1.0092658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64131367) q[0];
sx q[0];
rz(-2.2135493) q[0];
sx q[0];
rz(-2.4556659) q[0];
rz(-0.39086875) q[1];
sx q[1];
rz(-0.94376826) q[1];
sx q[1];
rz(2.2156782) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029322421) q[0];
sx q[0];
rz(-1.7824714) q[0];
sx q[0];
rz(-0.14365833) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7943031) q[2];
sx q[2];
rz(-2.1199193) q[2];
sx q[2];
rz(-1.4505475) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.029286413) q[1];
sx q[1];
rz(-0.46488133) q[1];
sx q[1];
rz(-1.2390562) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52243201) q[3];
sx q[3];
rz(-0.2345095) q[3];
sx q[3];
rz(1.0684551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9459076) q[2];
sx q[2];
rz(-0.47161272) q[2];
sx q[2];
rz(-0.53608981) q[2];
rz(-2.7219971) q[3];
sx q[3];
rz(-0.59046888) q[3];
sx q[3];
rz(1.4887811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560028) q[0];
sx q[0];
rz(-0.35878006) q[0];
sx q[0];
rz(-0.39500239) q[0];
rz(1.5123873) q[1];
sx q[1];
rz(-1.869166) q[1];
sx q[1];
rz(2.1283456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7618338) q[0];
sx q[0];
rz(-1.8031617) q[0];
sx q[0];
rz(-2.9985715) q[0];
rz(-2.7878694) q[2];
sx q[2];
rz(-1.612066) q[2];
sx q[2];
rz(-1.8281787) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4229065) q[1];
sx q[1];
rz(-1.4740853) q[1];
sx q[1];
rz(0.42214091) q[1];
rz(2.4471531) q[3];
sx q[3];
rz(-0.89530066) q[3];
sx q[3];
rz(-1.465786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44832486) q[2];
sx q[2];
rz(-1.5193181) q[2];
sx q[2];
rz(0.75941336) q[2];
rz(1.365472) q[3];
sx q[3];
rz(-1.1080192) q[3];
sx q[3];
rz(0.56994462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7286745) q[0];
sx q[0];
rz(-1.2644132) q[0];
sx q[0];
rz(1.2046474) q[0];
rz(0.57327523) q[1];
sx q[1];
rz(-1.9075867) q[1];
sx q[1];
rz(1.8214068) q[1];
rz(2.0157094) q[2];
sx q[2];
rz(-1.5013668) q[2];
sx q[2];
rz(-2.8253386) q[2];
rz(1.6462973) q[3];
sx q[3];
rz(-2.4423238) q[3];
sx q[3];
rz(-2.8764976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
