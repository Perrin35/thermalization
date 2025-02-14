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
rz(1.3321441) q[1];
sx q[1];
rz(-2.6982215) q[1];
sx q[1];
rz(2.3763357) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0527406) q[0];
sx q[0];
rz(-2.46457) q[0];
sx q[0];
rz(2.1734326) q[0];
rz(-pi) q[1];
rz(-1.2328495) q[2];
sx q[2];
rz(-1.8006007) q[2];
sx q[2];
rz(0.99529642) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87188497) q[1];
sx q[1];
rz(-0.34075865) q[1];
sx q[1];
rz(-3.0213256) q[1];
x q[2];
rz(1.2918858) q[3];
sx q[3];
rz(-1.6959012) q[3];
sx q[3];
rz(-0.58663128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29609933) q[2];
sx q[2];
rz(-2.0534434) q[2];
sx q[2];
rz(-1.6695401) q[2];
rz(1.4269525) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(2.6078687) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1832837) q[0];
sx q[0];
rz(-1.5077718) q[0];
sx q[0];
rz(-0.84877745) q[0];
rz(0.035004184) q[1];
sx q[1];
rz(-1.9947546) q[1];
sx q[1];
rz(-0.95357198) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6153107) q[0];
sx q[0];
rz(-2.1325975) q[0];
sx q[0];
rz(1.022382) q[0];
rz(-pi) q[1];
rz(-2.8368901) q[2];
sx q[2];
rz(-1.1977473) q[2];
sx q[2];
rz(-2.3359131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1998905) q[1];
sx q[1];
rz(-2.7459956) q[1];
sx q[1];
rz(-1.0426177) q[1];
rz(2.4783432) q[3];
sx q[3];
rz(-0.8732855) q[3];
sx q[3];
rz(-0.29479858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9780875) q[2];
sx q[2];
rz(-1.3601902) q[2];
sx q[2];
rz(2.2920442) q[2];
rz(-2.9761525) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(2.3918242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7716832) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(0.40192303) q[0];
rz(0.15332128) q[1];
sx q[1];
rz(-0.17686495) q[1];
sx q[1];
rz(2.0053999) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72969681) q[0];
sx q[0];
rz(-1.9703016) q[0];
sx q[0];
rz(-0.013042838) q[0];
x q[1];
rz(1.2966782) q[2];
sx q[2];
rz(-2.038836) q[2];
sx q[2];
rz(2.5646891) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0979251) q[1];
sx q[1];
rz(-2.8043724) q[1];
sx q[1];
rz(2.340576) q[1];
rz(-pi) q[2];
rz(2.6776776) q[3];
sx q[3];
rz(-1.216421) q[3];
sx q[3];
rz(-2.7257435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0574657) q[2];
sx q[2];
rz(-2.2630313) q[2];
sx q[2];
rz(3.086536) q[2];
rz(-1.7940686) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(-1.0042892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.36443001) q[0];
sx q[0];
rz(-0.20474064) q[0];
sx q[0];
rz(1.5740016) q[0];
rz(-0.69152999) q[1];
sx q[1];
rz(-0.66929308) q[1];
sx q[1];
rz(2.9561668) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3464934) q[0];
sx q[0];
rz(-0.43432626) q[0];
sx q[0];
rz(0.40652911) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9965956) q[2];
sx q[2];
rz(-1.123872) q[2];
sx q[2];
rz(-1.1584692) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5925827) q[1];
sx q[1];
rz(-1.8235778) q[1];
sx q[1];
rz(-1.2493253) q[1];
rz(-2.5957362) q[3];
sx q[3];
rz(-1.3935157) q[3];
sx q[3];
rz(-1.3455678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7376248) q[2];
sx q[2];
rz(-1.7125968) q[2];
sx q[2];
rz(-0.0052304012) q[2];
rz(-0.38749203) q[3];
sx q[3];
rz(-2.6369075) q[3];
sx q[3];
rz(-1.3030049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-2.4150647) q[0];
sx q[0];
rz(-2.1580577) q[0];
sx q[0];
rz(-2.5256185) q[0];
rz(1.9527324) q[1];
sx q[1];
rz(-0.6183466) q[1];
sx q[1];
rz(0.51330769) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0831871) q[0];
sx q[0];
rz(-1.8523916) q[0];
sx q[0];
rz(3.133899) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31536021) q[2];
sx q[2];
rz(-0.7509481) q[2];
sx q[2];
rz(2.2436004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9037902) q[1];
sx q[1];
rz(-2.2040743) q[1];
sx q[1];
rz(1.6602181) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3393976) q[3];
sx q[3];
rz(-2.6632892) q[3];
sx q[3];
rz(-1.4372642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1335699) q[2];
sx q[2];
rz(-1.2761389) q[2];
sx q[2];
rz(-2.3991154) q[2];
rz(1.4437458) q[3];
sx q[3];
rz(-1.4092813) q[3];
sx q[3];
rz(2.0402563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3105069) q[0];
sx q[0];
rz(-1.0467014) q[0];
sx q[0];
rz(2.781784) q[0];
rz(-2.2911435) q[1];
sx q[1];
rz(-1.9603739) q[1];
sx q[1];
rz(-2.6757619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3208) q[0];
sx q[0];
rz(-2.3893836) q[0];
sx q[0];
rz(1.8899931) q[0];
x q[1];
rz(-0.1419576) q[2];
sx q[2];
rz(-0.76830155) q[2];
sx q[2];
rz(-0.10266081) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9512144) q[1];
sx q[1];
rz(-2.8953027) q[1];
sx q[1];
rz(0.21065335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57329081) q[3];
sx q[3];
rz(-1.0430286) q[3];
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
rz(-0.049592169) q[3];
sx q[3];
rz(1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5704982) q[0];
sx q[0];
rz(-0.97935337) q[0];
sx q[0];
rz(1.1660227) q[0];
rz(2.1265538) q[1];
sx q[1];
rz(-0.53703419) q[1];
sx q[1];
rz(2.7551415) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7313663) q[0];
sx q[0];
rz(-1.7803734) q[0];
sx q[0];
rz(-3.1102577) q[0];
x q[1];
rz(-0.3101686) q[2];
sx q[2];
rz(-1.3259238) q[2];
sx q[2];
rz(2.819811) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1931175) q[1];
sx q[1];
rz(-1.5962692) q[1];
sx q[1];
rz(1.6465882) q[1];
rz(1.0284958) q[3];
sx q[3];
rz(-1.6784378) q[3];
sx q[3];
rz(0.018317761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7361136) q[2];
sx q[2];
rz(-1.4961286) q[2];
sx q[2];
rz(2.9138937) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88184083) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(-0.42801273) q[0];
rz(1.0182861) q[1];
sx q[1];
rz(-0.4363474) q[1];
sx q[1];
rz(2.2705618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329417) q[0];
sx q[0];
rz(-2.067446) q[0];
sx q[0];
rz(-0.92306851) q[0];
rz(-2.9350014) q[2];
sx q[2];
rz(-1.9211413) q[2];
sx q[2];
rz(-1.7889345) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6884675) q[1];
sx q[1];
rz(-0.18632132) q[1];
sx q[1];
rz(-0.86155714) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0139129) q[3];
sx q[3];
rz(-0.96864163) q[3];
sx q[3];
rz(-0.35246655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.009306) q[2];
sx q[2];
rz(-2.3460178) q[2];
sx q[2];
rz(1.2786678) q[2];
rz(-2.5631185) q[3];
sx q[3];
rz(-0.63251248) q[3];
sx q[3];
rz(-1.1540029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.4671675) q[0];
sx q[0];
rz(-0.217087) q[0];
sx q[0];
rz(-0.59980741) q[0];
rz(-1.2990052) q[1];
sx q[1];
rz(-1.5312342) q[1];
sx q[1];
rz(2.4997247) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079499809) q[0];
sx q[0];
rz(-0.60264093) q[0];
sx q[0];
rz(-0.63061611) q[0];
rz(-2.3360905) q[2];
sx q[2];
rz(-1.4664716) q[2];
sx q[2];
rz(-0.25539216) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.95002) q[1];
sx q[1];
rz(-0.86314647) q[1];
sx q[1];
rz(1.5950844) q[1];
x q[2];
rz(0.75662002) q[3];
sx q[3];
rz(-1.9709658) q[3];
sx q[3];
rz(2.4979046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.32144) q[2];
sx q[2];
rz(-1.4250616) q[2];
sx q[2];
rz(0.610262) q[2];
rz(-2.6214456) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(0.49351969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34058061) q[0];
sx q[0];
rz(-1.8980674) q[0];
sx q[0];
rz(-2.2035759) q[0];
rz(2.7998789) q[1];
sx q[1];
rz(-1.382117) q[1];
sx q[1];
rz(0.15728532) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92625916) q[0];
sx q[0];
rz(-2.4865827) q[0];
sx q[0];
rz(2.8192855) q[0];
rz(-pi) q[1];
rz(-3.0523275) q[2];
sx q[2];
rz(-2.3815739) q[2];
sx q[2];
rz(1.4475198) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5648236) q[1];
sx q[1];
rz(-1.2693468) q[1];
sx q[1];
rz(2.5891073) q[1];
rz(1.1620164) q[3];
sx q[3];
rz(-1.805154) q[3];
sx q[3];
rz(-2.5046668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9226795) q[2];
sx q[2];
rz(-2.7965386) q[2];
sx q[2];
rz(-1.867713) q[2];
rz(3.1359361) q[3];
sx q[3];
rz(-1.4414682) q[3];
sx q[3];
rz(-0.98102942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8511178) q[0];
sx q[0];
rz(-1.9708451) q[0];
sx q[0];
rz(-3.0921902) q[0];
rz(0.46650096) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(0.51085761) q[2];
sx q[2];
rz(-2.1759773) q[2];
sx q[2];
rz(0.89531384) q[2];
rz(2.9644184) q[3];
sx q[3];
rz(-0.43779793) q[3];
sx q[3];
rz(1.4878185) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
