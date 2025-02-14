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
rz(-1.3938067) q[0];
sx q[0];
rz(-0.90016794) q[0];
sx q[0];
rz(1.7322487) q[0];
rz(0.70977587) q[1];
sx q[1];
rz(-2.8083399) q[1];
sx q[1];
rz(0.38212734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6674675) q[0];
sx q[0];
rz(-2.29974) q[0];
sx q[0];
rz(-1.390662) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3163739) q[2];
sx q[2];
rz(-1.7013936) q[2];
sx q[2];
rz(0.25564627) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5772878) q[1];
sx q[1];
rz(-1.9365942) q[1];
sx q[1];
rz(-1.1759961) q[1];
rz(-pi) q[2];
rz(3.0218533) q[3];
sx q[3];
rz(-2.2594707) q[3];
sx q[3];
rz(2.7514091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0673151) q[2];
sx q[2];
rz(-1.5019608) q[2];
sx q[2];
rz(-0.074946694) q[2];
rz(1.8986757) q[3];
sx q[3];
rz(-0.26151812) q[3];
sx q[3];
rz(-2.7459131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13260929) q[0];
sx q[0];
rz(-0.41985303) q[0];
sx q[0];
rz(-2.5123151) q[0];
rz(2.3043326) q[1];
sx q[1];
rz(-2.3910797) q[1];
sx q[1];
rz(-0.57964051) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4052359) q[0];
sx q[0];
rz(-1.7742298) q[0];
sx q[0];
rz(0.74958165) q[0];
rz(0.39548611) q[2];
sx q[2];
rz(-1.8137534) q[2];
sx q[2];
rz(1.1370575) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.9379679) q[1];
sx q[1];
rz(-1.1826828) q[1];
sx q[1];
rz(-1.2118503) q[1];
x q[2];
rz(1.5906672) q[3];
sx q[3];
rz(-1.7345034) q[3];
sx q[3];
rz(0.57396561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89977449) q[2];
sx q[2];
rz(-2.9782229) q[2];
sx q[2];
rz(-2.3972798) q[2];
rz(0.36738473) q[3];
sx q[3];
rz(-1.2920047) q[3];
sx q[3];
rz(-1.415409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431417) q[0];
sx q[0];
rz(-0.95040584) q[0];
sx q[0];
rz(-1.0071734) q[0];
rz(-3.1305283) q[1];
sx q[1];
rz(-2.8304351) q[1];
sx q[1];
rz(0.99753582) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0031834) q[0];
sx q[0];
rz(-1.5965726) q[0];
sx q[0];
rz(1.7083712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80834016) q[2];
sx q[2];
rz(-1.6534717) q[2];
sx q[2];
rz(2.9660564) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.0094781965) q[1];
sx q[1];
rz(-1.0524228) q[1];
sx q[1];
rz(-1.5794157) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5854534) q[3];
sx q[3];
rz(-1.0590226) q[3];
sx q[3];
rz(-2.0320333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.02115383) q[2];
sx q[2];
rz(-2.8595371) q[2];
sx q[2];
rz(-0.8417449) q[2];
rz(0.65819955) q[3];
sx q[3];
rz(-0.87116146) q[3];
sx q[3];
rz(0.021520821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72092527) q[0];
sx q[0];
rz(-3.0166716) q[0];
sx q[0];
rz(-2.4125873) q[0];
rz(-2.3782102) q[1];
sx q[1];
rz(-2.5114676) q[1];
sx q[1];
rz(2.7785832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1358534) q[0];
sx q[0];
rz(-0.38977888) q[0];
sx q[0];
rz(-0.71758349) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3173787) q[2];
sx q[2];
rz(-1.1477648) q[2];
sx q[2];
rz(1.8366448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2263086) q[1];
sx q[1];
rz(-1.7123509) q[1];
sx q[1];
rz(2.841921) q[1];
rz(-pi) q[2];
rz(-1.6239802) q[3];
sx q[3];
rz(-1.4442354) q[3];
sx q[3];
rz(-2.3317331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9558941) q[2];
sx q[2];
rz(-2.318435) q[2];
sx q[2];
rz(1.4996747) q[2];
rz(-2.5540292) q[3];
sx q[3];
rz(-2.1524119) q[3];
sx q[3];
rz(0.51830083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8374306) q[0];
sx q[0];
rz(-0.70697933) q[0];
sx q[0];
rz(-2.8564603) q[0];
rz(0.25310165) q[1];
sx q[1];
rz(-2.1123835) q[1];
sx q[1];
rz(2.0957799) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3110712) q[0];
sx q[0];
rz(-1.9609299) q[0];
sx q[0];
rz(2.5907787) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1509456) q[2];
sx q[2];
rz(-1.1641181) q[2];
sx q[2];
rz(0.57957725) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1872028) q[1];
sx q[1];
rz(-2.7437401) q[1];
sx q[1];
rz(-0.14601645) q[1];
x q[2];
rz(0.67774421) q[3];
sx q[3];
rz(-1.7278921) q[3];
sx q[3];
rz(-1.398063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73695856) q[2];
sx q[2];
rz(-2.0817231) q[2];
sx q[2];
rz(2.5671379) q[2];
rz(-2.5308841) q[3];
sx q[3];
rz(-0.52053958) q[3];
sx q[3];
rz(-2.0846539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33559281) q[0];
sx q[0];
rz(-1.7570423) q[0];
sx q[0];
rz(-0.35032508) q[0];
rz(-2.8490745) q[1];
sx q[1];
rz(-3.0151093) q[1];
sx q[1];
rz(2.3622321) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8838046) q[0];
sx q[0];
rz(-1.6291926) q[0];
sx q[0];
rz(-1.4517054) q[0];
rz(1.3931403) q[2];
sx q[2];
rz(-2.0442932) q[2];
sx q[2];
rz(2.908978) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7747969) q[1];
sx q[1];
rz(-1.8908943) q[1];
sx q[1];
rz(1.7843549) q[1];
rz(-1.4102226) q[3];
sx q[3];
rz(-2.0856557) q[3];
sx q[3];
rz(-0.60100473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.902035) q[2];
sx q[2];
rz(-1.6263447) q[2];
sx q[2];
rz(-2.1774192) q[2];
rz(2.9686109) q[3];
sx q[3];
rz(-1.0123342) q[3];
sx q[3];
rz(-0.13121901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59703374) q[0];
sx q[0];
rz(-1.9439161) q[0];
sx q[0];
rz(-0.069393754) q[0];
rz(-1.767905) q[1];
sx q[1];
rz(-1.8927788) q[1];
sx q[1];
rz(0.51838851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64053464) q[0];
sx q[0];
rz(-2.448521) q[0];
sx q[0];
rz(1.3035167) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0807627) q[2];
sx q[2];
rz(-1.9980717) q[2];
sx q[2];
rz(-3.0025122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9897108) q[1];
sx q[1];
rz(-3.0009785) q[1];
sx q[1];
rz(0.18528823) q[1];
rz(2.9767738) q[3];
sx q[3];
rz(-1.6281152) q[3];
sx q[3];
rz(-2.8652193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9531276) q[2];
sx q[2];
rz(-1.1981107) q[2];
sx q[2];
rz(0.24448621) q[2];
rz(1.9237579) q[3];
sx q[3];
rz(-3.0471424) q[3];
sx q[3];
rz(-0.34734669) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0511047) q[0];
sx q[0];
rz(-2.8489887) q[0];
sx q[0];
rz(-2.4421413) q[0];
rz(0.72169101) q[1];
sx q[1];
rz(-1.3612008) q[1];
sx q[1];
rz(-1.9500505) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7724567) q[0];
sx q[0];
rz(-1.4459608) q[0];
sx q[0];
rz(-2.9278276) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6416816) q[2];
sx q[2];
rz(-2.8770718) q[2];
sx q[2];
rz(-1.7873639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8766159) q[1];
sx q[1];
rz(-0.54538762) q[1];
sx q[1];
rz(1.7343069) q[1];
x q[2];
rz(1.9189214) q[3];
sx q[3];
rz(-1.1468107) q[3];
sx q[3];
rz(-0.67397398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2193853) q[2];
sx q[2];
rz(-0.81874138) q[2];
sx q[2];
rz(-1.7259664) q[2];
rz(2.7042232) q[3];
sx q[3];
rz(-2.8919817) q[3];
sx q[3];
rz(2.6387446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9472083) q[0];
sx q[0];
rz(-1.2057065) q[0];
sx q[0];
rz(-2.6956287) q[0];
rz(1.3392316) q[1];
sx q[1];
rz(-2.4868592) q[1];
sx q[1];
rz(-1.9816678) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8561607) q[0];
sx q[0];
rz(-2.2158379) q[0];
sx q[0];
rz(-2.0061532) q[0];
x q[1];
rz(3.0662698) q[2];
sx q[2];
rz(-1.3100071) q[2];
sx q[2];
rz(-2.0465849) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36148188) q[1];
sx q[1];
rz(-1.0358397) q[1];
sx q[1];
rz(-1.8423716) q[1];
rz(-0.58034133) q[3];
sx q[3];
rz(-0.99513061) q[3];
sx q[3];
rz(0.87596303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7356073) q[2];
sx q[2];
rz(-2.1129825) q[2];
sx q[2];
rz(2.410991) q[2];
rz(-0.90100151) q[3];
sx q[3];
rz(-0.42947072) q[3];
sx q[3];
rz(-3.1072531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.63448298) q[0];
sx q[0];
rz(-2.6431838) q[0];
sx q[0];
rz(-0.46257567) q[0];
rz(-1.0628465) q[1];
sx q[1];
rz(-1.7741508) q[1];
sx q[1];
rz(-0.06632334) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6505867) q[0];
sx q[0];
rz(-2.1438476) q[0];
sx q[0];
rz(-1.2038403) q[0];
rz(0.1758258) q[2];
sx q[2];
rz(-1.9011902) q[2];
sx q[2];
rz(-1.2820455) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8005261) q[1];
sx q[1];
rz(-0.99565804) q[1];
sx q[1];
rz(-2.2209206) q[1];
rz(-pi) q[2];
rz(1.9948694) q[3];
sx q[3];
rz(-1.3641285) q[3];
sx q[3];
rz(0.28462946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.037420951) q[2];
sx q[2];
rz(-2.6080242) q[2];
sx q[2];
rz(1.9083692) q[2];
rz(2.6836266) q[3];
sx q[3];
rz(-2.870324) q[3];
sx q[3];
rz(2.7428194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11252277) q[0];
sx q[0];
rz(-1.4334913) q[0];
sx q[0];
rz(1.7423472) q[0];
rz(-0.15432547) q[1];
sx q[1];
rz(-1.4230774) q[1];
sx q[1];
rz(1.9565061) q[1];
rz(1.1823282) q[2];
sx q[2];
rz(-1.2797838) q[2];
sx q[2];
rz(-1.01576) q[2];
rz(0.88944351) q[3];
sx q[3];
rz(-2.4460197) q[3];
sx q[3];
rz(-1.0424436) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
