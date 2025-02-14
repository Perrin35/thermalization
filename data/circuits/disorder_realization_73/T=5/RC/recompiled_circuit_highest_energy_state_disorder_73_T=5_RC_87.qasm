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
rz(-2.7849164) q[0];
sx q[0];
rz(1.6950322) q[0];
sx q[0];
rz(8.5230081) q[0];
rz(-1.8094485) q[1];
sx q[1];
rz(-0.44337115) q[1];
sx q[1];
rz(-2.3763357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63432307) q[0];
sx q[0];
rz(-1.028484) q[0];
sx q[0];
rz(-0.42748283) q[0];
rz(-pi) q[1];
rz(-1.9087431) q[2];
sx q[2];
rz(-1.3409919) q[2];
sx q[2];
rz(0.99529642) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74435788) q[1];
sx q[1];
rz(-1.2325979) q[1];
sx q[1];
rz(-1.5282791) q[1];
rz(1.1422994) q[3];
sx q[3];
rz(-2.836578) q[3];
sx q[3];
rz(-1.7465137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8454933) q[2];
sx q[2];
rz(-1.0881492) q[2];
sx q[2];
rz(1.4720526) q[2];
rz(1.4269525) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(-0.53372395) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1832837) q[0];
sx q[0];
rz(-1.5077718) q[0];
sx q[0];
rz(-2.2928152) q[0];
rz(0.035004184) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(0.95357198) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35915056) q[0];
sx q[0];
rz(-1.1138565) q[0];
sx q[0];
rz(-0.63553973) q[0];
x q[1];
rz(2.8368901) q[2];
sx q[2];
rz(-1.1977473) q[2];
sx q[2];
rz(2.3359131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1998905) q[1];
sx q[1];
rz(-2.7459956) q[1];
sx q[1];
rz(-2.0989749) q[1];
rz(0.66324948) q[3];
sx q[3];
rz(-2.2683072) q[3];
sx q[3];
rz(2.8467941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9780875) q[2];
sx q[2];
rz(-1.7814025) q[2];
sx q[2];
rz(2.2920442) q[2];
rz(0.16544011) q[3];
sx q[3];
rz(-1.6980419) q[3];
sx q[3];
rz(-2.3918242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7716832) q[0];
sx q[0];
rz(-0.92215466) q[0];
sx q[0];
rz(-0.40192303) q[0];
rz(2.9882714) q[1];
sx q[1];
rz(-2.9647277) q[1];
sx q[1];
rz(2.0053999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72969681) q[0];
sx q[0];
rz(-1.171291) q[0];
sx q[0];
rz(3.1285498) q[0];
rz(-pi) q[1];
rz(0.49164518) q[2];
sx q[2];
rz(-2.6043713) q[2];
sx q[2];
rz(-1.1342837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3553149) q[1];
sx q[1];
rz(-1.8031562) q[1];
sx q[1];
rz(1.3241598) q[1];
rz(0.69092423) q[3];
sx q[3];
rz(-2.5657585) q[3];
sx q[3];
rz(1.3802647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0574657) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(-3.086536) q[2];
rz(-1.7940686) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(2.1373035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36443001) q[0];
sx q[0];
rz(-2.936852) q[0];
sx q[0];
rz(-1.5675911) q[0];
rz(-2.4500627) q[1];
sx q[1];
rz(-0.66929308) q[1];
sx q[1];
rz(-2.9561668) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3464934) q[0];
sx q[0];
rz(-2.7072664) q[0];
sx q[0];
rz(0.40652911) q[0];
x q[1];
rz(-1.9965956) q[2];
sx q[2];
rz(-1.123872) q[2];
sx q[2];
rz(1.1584692) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37764964) q[1];
sx q[1];
rz(-2.7353706) q[1];
sx q[1];
rz(-0.88546125) q[1];
rz(-0.33230761) q[3];
sx q[3];
rz(-2.5704567) q[3];
sx q[3];
rz(-2.6337998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4039679) q[2];
sx q[2];
rz(-1.4289958) q[2];
sx q[2];
rz(-0.0052304012) q[2];
rz(-0.38749203) q[3];
sx q[3];
rz(-0.5046851) q[3];
sx q[3];
rz(-1.8385878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4150647) q[0];
sx q[0];
rz(-0.98353493) q[0];
sx q[0];
rz(2.5256185) q[0];
rz(1.1888602) q[1];
sx q[1];
rz(-2.523246) q[1];
sx q[1];
rz(0.51330769) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0831871) q[0];
sx q[0];
rz(-1.289201) q[0];
sx q[0];
rz(-0.0076936184) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31536021) q[2];
sx q[2];
rz(-0.7509481) q[2];
sx q[2];
rz(2.2436004) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9037902) q[1];
sx q[1];
rz(-0.93751838) q[1];
sx q[1];
rz(-1.6602181) q[1];
rz(3.0232459) q[3];
sx q[3];
rz(-2.0353298) q[3];
sx q[3];
rz(1.9637513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1335699) q[2];
sx q[2];
rz(-1.8654537) q[2];
sx q[2];
rz(0.7424773) q[2];
rz(1.4437458) q[3];
sx q[3];
rz(-1.4092813) q[3];
sx q[3];
rz(2.0402563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3105069) q[0];
sx q[0];
rz(-1.0467014) q[0];
sx q[0];
rz(2.781784) q[0];
rz(2.2911435) q[1];
sx q[1];
rz(-1.9603739) q[1];
sx q[1];
rz(2.6757619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1547928) q[0];
sx q[0];
rz(-1.354711) q[0];
sx q[0];
rz(-0.84439069) q[0];
rz(-pi) q[1];
rz(-1.7066782) q[2];
sx q[2];
rz(-2.3294221) q[2];
sx q[2];
rz(0.093531724) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40738525) q[1];
sx q[1];
rz(-1.8115329) q[1];
sx q[1];
rz(-1.6233141) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5683018) q[3];
sx q[3];
rz(-2.0985641) q[3];
sx q[3];
rz(2.9506369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74238527) q[2];
sx q[2];
rz(-1.9611605) q[2];
sx q[2];
rz(0.36435374) q[2];
rz(-2.897701) q[3];
sx q[3];
rz(-3.0920005) q[3];
sx q[3];
rz(-1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57109443) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(-1.97557) q[0];
rz(-1.0150389) q[1];
sx q[1];
rz(-2.6045585) q[1];
sx q[1];
rz(0.38645116) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84595118) q[0];
sx q[0];
rz(-1.6014454) q[0];
sx q[0];
rz(1.7804734) q[0];
rz(-pi) q[1];
rz(-1.3141733) q[2];
sx q[2];
rz(-1.2701747) q[2];
sx q[2];
rz(-1.8150309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4402953) q[1];
sx q[1];
rz(-3.0616425) q[1];
sx q[1];
rz(1.8953807) q[1];
rz(1.3644045) q[3];
sx q[3];
rz(-2.5897621) q[3];
sx q[3];
rz(-1.3760374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40547907) q[2];
sx q[2];
rz(-1.4961286) q[2];
sx q[2];
rz(-0.22769895) q[2];
rz(-2.0685711) q[3];
sx q[3];
rz(-2.3927972) q[3];
sx q[3];
rz(-2.8480215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88184083) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(-2.7135799) q[0];
rz(1.0182861) q[1];
sx q[1];
rz(-2.7052453) q[1];
sx q[1];
rz(-2.2705618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4081602) q[0];
sx q[0];
rz(-1.0116315) q[0];
sx q[0];
rz(0.59691043) q[0];
rz(2.0823048) q[2];
sx q[2];
rz(-2.7370484) q[2];
sx q[2];
rz(-1.900857) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.58295137) q[1];
sx q[1];
rz(-1.4498596) q[1];
sx q[1];
rz(-1.4287097) q[1];
rz(-2.4610041) q[3];
sx q[3];
rz(-1.1201123) q[3];
sx q[3];
rz(-0.87928538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13228664) q[2];
sx q[2];
rz(-2.3460178) q[2];
sx q[2];
rz(-1.2786678) q[2];
rz(2.5631185) q[3];
sx q[3];
rz(-2.5090802) q[3];
sx q[3];
rz(1.9875897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4671675) q[0];
sx q[0];
rz(-0.217087) q[0];
sx q[0];
rz(0.59980741) q[0];
rz(1.2990052) q[1];
sx q[1];
rz(-1.5312342) q[1];
sx q[1];
rz(-2.4997247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80461795) q[0];
sx q[0];
rz(-1.0952767) q[0];
sx q[0];
rz(-1.1853976) q[0];
x q[1];
rz(-1.4207877) q[2];
sx q[2];
rz(-2.370655) q[2];
sx q[2];
rz(1.2074167) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.228926) q[1];
sx q[1];
rz(-0.70799457) q[1];
sx q[1];
rz(-0.028381746) q[1];
rz(-pi) q[2];
rz(1.0439375) q[3];
sx q[3];
rz(-2.2552285) q[3];
sx q[3];
rz(1.2795283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.32144) q[2];
sx q[2];
rz(-1.4250616) q[2];
sx q[2];
rz(-0.610262) q[2];
rz(-2.6214456) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(0.49351969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34058061) q[0];
sx q[0];
rz(-1.8980674) q[0];
sx q[0];
rz(2.2035759) q[0];
rz(0.34171379) q[1];
sx q[1];
rz(-1.7594756) q[1];
sx q[1];
rz(-2.9843073) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2153335) q[0];
sx q[0];
rz(-0.65500998) q[0];
sx q[0];
rz(2.8192855) q[0];
rz(-pi) q[1];
rz(-0.75802676) q[2];
sx q[2];
rz(-1.6322513) q[2];
sx q[2];
rz(-3.0830992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54522911) q[1];
sx q[1];
rz(-0.62178722) q[1];
sx q[1];
rz(0.53485628) q[1];
rz(-1.9795762) q[3];
sx q[3];
rz(-1.805154) q[3];
sx q[3];
rz(0.63692585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9226795) q[2];
sx q[2];
rz(-0.34505406) q[2];
sx q[2];
rz(1.2738796) q[2];
rz(3.1359361) q[3];
sx q[3];
rz(-1.7001245) q[3];
sx q[3];
rz(-2.1605632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8511178) q[0];
sx q[0];
rz(-1.1707476) q[0];
sx q[0];
rz(0.049402417) q[0];
rz(-0.46650096) q[1];
sx q[1];
rz(-1.7955045) q[1];
sx q[1];
rz(-2.4972965) q[1];
rz(-0.51085761) q[2];
sx q[2];
rz(-0.96561531) q[2];
sx q[2];
rz(-2.2462788) q[2];
rz(-1.6531108) q[3];
sx q[3];
rz(-2.001279) q[3];
sx q[3];
rz(1.6829987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
