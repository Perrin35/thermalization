OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6717186) q[0];
sx q[0];
rz(-1.1550386) q[0];
sx q[0];
rz(1.4299097) q[0];
rz(-2.6861796) q[1];
sx q[1];
rz(-2.5502584) q[1];
sx q[1];
rz(1.1270181) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4658964) q[0];
sx q[0];
rz(-1.5592335) q[0];
sx q[0];
rz(-1.8013062) q[0];
x q[1];
rz(1.3052218) q[2];
sx q[2];
rz(-1.9023393) q[2];
sx q[2];
rz(1.3217373) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0309567) q[1];
sx q[1];
rz(-2.1345761) q[1];
sx q[1];
rz(0.33302506) q[1];
x q[2];
rz(-1.5886878) q[3];
sx q[3];
rz(-1.761529) q[3];
sx q[3];
rz(0.73303661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0226125) q[2];
sx q[2];
rz(-0.70465124) q[2];
sx q[2];
rz(1.4878147) q[2];
rz(-0.37114272) q[3];
sx q[3];
rz(-1.1266339) q[3];
sx q[3];
rz(2.3625372) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30598518) q[0];
sx q[0];
rz(-1.3779409) q[0];
sx q[0];
rz(0.66407472) q[0];
rz(-0.56655073) q[1];
sx q[1];
rz(-1.605875) q[1];
sx q[1];
rz(2.9552592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31748397) q[0];
sx q[0];
rz(-0.58360277) q[0];
sx q[0];
rz(-1.9822796) q[0];
rz(0.9695942) q[2];
sx q[2];
rz(-0.45192138) q[2];
sx q[2];
rz(2.4659803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9606321) q[1];
sx q[1];
rz(-0.69968984) q[1];
sx q[1];
rz(-1.8262499) q[1];
rz(1.4545031) q[3];
sx q[3];
rz(-2.6175559) q[3];
sx q[3];
rz(-0.91479036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8191007) q[2];
sx q[2];
rz(-1.842061) q[2];
sx q[2];
rz(1.5796278) q[2];
rz(-1.1473848) q[3];
sx q[3];
rz(-2.8473144) q[3];
sx q[3];
rz(1.3636205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086562432) q[0];
sx q[0];
rz(-1.4921621) q[0];
sx q[0];
rz(-0.30558875) q[0];
rz(-1.2809523) q[1];
sx q[1];
rz(-2.8260904) q[1];
sx q[1];
rz(-1.5234647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9010503) q[0];
sx q[0];
rz(-0.43473703) q[0];
sx q[0];
rz(-0.57603277) q[0];
x q[1];
rz(1.5554713) q[2];
sx q[2];
rz(-1.6159247) q[2];
sx q[2];
rz(-0.082782291) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.87464) q[1];
sx q[1];
rz(-2.2494509) q[1];
sx q[1];
rz(1.1043504) q[1];
x q[2];
rz(-2.9039246) q[3];
sx q[3];
rz(-2.1690627) q[3];
sx q[3];
rz(3.0868343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4898701) q[2];
sx q[2];
rz(-1.8946596) q[2];
sx q[2];
rz(-2.9215802) q[2];
rz(-3.0618727) q[3];
sx q[3];
rz(-2.1646175) q[3];
sx q[3];
rz(0.93130934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7109969) q[0];
sx q[0];
rz(-1.3871223) q[0];
sx q[0];
rz(3.1331449) q[0];
rz(-1.9795817) q[1];
sx q[1];
rz(-1.153667) q[1];
sx q[1];
rz(1.1147503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8793068) q[0];
sx q[0];
rz(-1.382785) q[0];
sx q[0];
rz(-1.8858869) q[0];
rz(-pi) q[1];
rz(-2.5774245) q[2];
sx q[2];
rz(-2.0551066) q[2];
sx q[2];
rz(2.6251305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47059688) q[1];
sx q[1];
rz(-2.6866815) q[1];
sx q[1];
rz(-0.90683337) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9825806) q[3];
sx q[3];
rz(-1.9596385) q[3];
sx q[3];
rz(1.4847311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.70725012) q[2];
sx q[2];
rz(-1.0184526) q[2];
sx q[2];
rz(-0.84558359) q[2];
rz(-2.8905458) q[3];
sx q[3];
rz(-1.2716525) q[3];
sx q[3];
rz(-0.74560753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71395981) q[0];
sx q[0];
rz(-2.7290955) q[0];
sx q[0];
rz(-2.1115671) q[0];
rz(-0.59459844) q[1];
sx q[1];
rz(-1.4232891) q[1];
sx q[1];
rz(1.3535708) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0145688) q[0];
sx q[0];
rz(-1.5480124) q[0];
sx q[0];
rz(1.5915074) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7556023) q[2];
sx q[2];
rz(-1.2535742) q[2];
sx q[2];
rz(2.3776907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7723778) q[1];
sx q[1];
rz(-0.90638834) q[1];
sx q[1];
rz(-0.64488761) q[1];
x q[2];
rz(1.0566125) q[3];
sx q[3];
rz(-2.100889) q[3];
sx q[3];
rz(2.7831654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9109965) q[2];
sx q[2];
rz(-1.6350919) q[2];
sx q[2];
rz(0.050749151) q[2];
rz(-2.0283608) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(-2.7509403) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0473061) q[0];
sx q[0];
rz(-2.081649) q[0];
sx q[0];
rz(2.768709) q[0];
rz(1.8376384) q[1];
sx q[1];
rz(-1.2645489) q[1];
sx q[1];
rz(-1.0252999) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7367497) q[0];
sx q[0];
rz(-1.4944634) q[0];
sx q[0];
rz(-0.45053225) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16496678) q[2];
sx q[2];
rz(-1.8968762) q[2];
sx q[2];
rz(3.0086089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2220425) q[1];
sx q[1];
rz(-1.371438) q[1];
sx q[1];
rz(0.87032236) q[1];
x q[2];
rz(2.9168374) q[3];
sx q[3];
rz(-1.156575) q[3];
sx q[3];
rz(2.1638768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.23124) q[2];
sx q[2];
rz(-2.9515036) q[2];
sx q[2];
rz(-1.4300038) q[2];
rz(1.6010239) q[3];
sx q[3];
rz(-2.4547596) q[3];
sx q[3];
rz(-1.6704667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98825276) q[0];
sx q[0];
rz(-1.3000458) q[0];
sx q[0];
rz(2.7100995) q[0];
rz(1.2639812) q[1];
sx q[1];
rz(-1.5411721) q[1];
sx q[1];
rz(-2.4914609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1558485) q[0];
sx q[0];
rz(-1.2835644) q[0];
sx q[0];
rz(1.656507) q[0];
rz(-pi) q[1];
rz(1.6986786) q[2];
sx q[2];
rz(-2.4543216) q[2];
sx q[2];
rz(-0.56035794) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0228867) q[1];
sx q[1];
rz(-0.71057767) q[1];
sx q[1];
rz(-1.3936682) q[1];
rz(-pi) q[2];
rz(2.0286719) q[3];
sx q[3];
rz(-1.5639389) q[3];
sx q[3];
rz(-0.53660652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55107895) q[2];
sx q[2];
rz(-1.3912018) q[2];
sx q[2];
rz(3.1371269) q[2];
rz(3.1320069) q[3];
sx q[3];
rz(-1.3349345) q[3];
sx q[3];
rz(2.7981304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9850605) q[0];
sx q[0];
rz(-2.8592906) q[0];
sx q[0];
rz(0.16047934) q[0];
rz(-1.7566682) q[1];
sx q[1];
rz(-0.79912186) q[1];
sx q[1];
rz(2.8327732) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64381448) q[0];
sx q[0];
rz(-0.22499946) q[0];
sx q[0];
rz(-2.1589222) q[0];
rz(-0.71134348) q[2];
sx q[2];
rz(-1.2738196) q[2];
sx q[2];
rz(-2.1893196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9659979) q[1];
sx q[1];
rz(-1.9456989) q[1];
sx q[1];
rz(2.7213827) q[1];
rz(-1.2587955) q[3];
sx q[3];
rz(-1.9083169) q[3];
sx q[3];
rz(-0.58956335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34362346) q[2];
sx q[2];
rz(-0.60911959) q[2];
sx q[2];
rz(-1.5416175) q[2];
rz(-2.4566417) q[3];
sx q[3];
rz(-1.5870321) q[3];
sx q[3];
rz(1.5233013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110382) q[0];
sx q[0];
rz(-1.7031952) q[0];
sx q[0];
rz(2.5901929) q[0];
rz(-2.664387) q[1];
sx q[1];
rz(-0.91608945) q[1];
sx q[1];
rz(0.83470693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4906368) q[0];
sx q[0];
rz(-2.0623364) q[0];
sx q[0];
rz(1.0440582) q[0];
rz(-pi) q[1];
rz(-2.4407225) q[2];
sx q[2];
rz(-1.5972023) q[2];
sx q[2];
rz(-0.7757265) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44025035) q[1];
sx q[1];
rz(-1.6052488) q[1];
sx q[1];
rz(-2.1828096) q[1];
rz(-2.1797997) q[3];
sx q[3];
rz(-1.0474221) q[3];
sx q[3];
rz(2.9030622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7029552) q[2];
sx q[2];
rz(-2.1795887) q[2];
sx q[2];
rz(-0.23923242) q[2];
rz(0.3244102) q[3];
sx q[3];
rz(-1.7041465) q[3];
sx q[3];
rz(-1.602406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6259554) q[0];
sx q[0];
rz(-0.13787585) q[0];
sx q[0];
rz(2.4249518) q[0];
rz(-0.58249885) q[1];
sx q[1];
rz(-1.3526252) q[1];
sx q[1];
rz(-0.31002054) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4449883) q[0];
sx q[0];
rz(-1.9194845) q[0];
sx q[0];
rz(-0.83072386) q[0];
rz(2.7298195) q[2];
sx q[2];
rz(-1.9343107) q[2];
sx q[2];
rz(-2.1902254) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2523732) q[1];
sx q[1];
rz(-0.82102005) q[1];
sx q[1];
rz(2.0682425) q[1];
rz(-pi) q[2];
rz(1.0729372) q[3];
sx q[3];
rz(-2.0611576) q[3];
sx q[3];
rz(-0.89807248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50905716) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(-1.4198111) q[2];
rz(-1.696473) q[3];
sx q[3];
rz(-0.71075478) q[3];
sx q[3];
rz(2.3783309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2390908) q[0];
sx q[0];
rz(-2.1949174) q[0];
sx q[0];
rz(1.3876023) q[0];
rz(1.6613962) q[1];
sx q[1];
rz(-1.7330685) q[1];
sx q[1];
rz(0.20191244) q[1];
rz(-1.6461862) q[2];
sx q[2];
rz(-1.8373377) q[2];
sx q[2];
rz(-1.1129825) q[2];
rz(-1.0140513) q[3];
sx q[3];
rz(-0.76273668) q[3];
sx q[3];
rz(-0.098165011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
