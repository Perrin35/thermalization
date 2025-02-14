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
rz(0.27999347) q[0];
sx q[0];
rz(2.1092829) q[0];
sx q[0];
rz(10.584925) q[0];
rz(2.2634444) q[1];
sx q[1];
rz(-0.4385837) q[1];
sx q[1];
rz(-2.261472) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297247) q[0];
sx q[0];
rz(-2.1639185) q[0];
sx q[0];
rz(-0.2082227) q[0];
rz(-pi) q[1];
rz(3.0612602) q[2];
sx q[2];
rz(-1.972562) q[2];
sx q[2];
rz(0.062594819) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82950006) q[1];
sx q[1];
rz(-2.1662427) q[1];
sx q[1];
rz(0.19705806) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9008242) q[3];
sx q[3];
rz(-0.46564049) q[3];
sx q[3];
rz(-0.2351528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1151513) q[2];
sx q[2];
rz(-2.5796964) q[2];
sx q[2];
rz(-1.119841) q[2];
rz(-1.0037054) q[3];
sx q[3];
rz(-2.7916838) q[3];
sx q[3];
rz(-2.2899535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72964662) q[0];
sx q[0];
rz(-0.8929407) q[0];
sx q[0];
rz(-2.6935284) q[0];
rz(-1.0753151) q[1];
sx q[1];
rz(-1.0766462) q[1];
sx q[1];
rz(-0.040464673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4141636) q[0];
sx q[0];
rz(-2.232904) q[0];
sx q[0];
rz(0.24358769) q[0];
x q[1];
rz(1.7732271) q[2];
sx q[2];
rz(-0.42095612) q[2];
sx q[2];
rz(-1.6594002) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3738714) q[1];
sx q[1];
rz(-1.5295856) q[1];
sx q[1];
rz(2.8727864) q[1];
rz(-pi) q[2];
rz(-3.1149666) q[3];
sx q[3];
rz(-2.8288719) q[3];
sx q[3];
rz(1.3440162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.86541092) q[2];
sx q[2];
rz(-2.2409596) q[2];
sx q[2];
rz(0.5033699) q[2];
rz(1.6940176) q[3];
sx q[3];
rz(-1.9083128) q[3];
sx q[3];
rz(-2.2342822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68509787) q[0];
sx q[0];
rz(-2.1879897) q[0];
sx q[0];
rz(-0.29846919) q[0];
rz(1.586277) q[1];
sx q[1];
rz(-1.7319873) q[1];
sx q[1];
rz(1.635199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28373805) q[0];
sx q[0];
rz(-1.8687792) q[0];
sx q[0];
rz(1.4922754) q[0];
rz(-1.9048018) q[2];
sx q[2];
rz(-2.4270505) q[2];
sx q[2];
rz(1.3151907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8019234) q[1];
sx q[1];
rz(-0.69782062) q[1];
sx q[1];
rz(0.20734792) q[1];
rz(-pi) q[2];
rz(-1.3812341) q[3];
sx q[3];
rz(-1.5175356) q[3];
sx q[3];
rz(1.9587751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7895268) q[2];
sx q[2];
rz(-2.3730998) q[2];
sx q[2];
rz(-2.3699769) q[2];
rz(-1.7302892) q[3];
sx q[3];
rz(-1.227042) q[3];
sx q[3];
rz(2.3119149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4378433) q[0];
sx q[0];
rz(-2.4898536) q[0];
sx q[0];
rz(-1.334345) q[0];
rz(-1.5127381) q[1];
sx q[1];
rz(-1.6213497) q[1];
sx q[1];
rz(3.0049862) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2481987) q[0];
sx q[0];
rz(-1.8780439) q[0];
sx q[0];
rz(-2.224438) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0380546) q[2];
sx q[2];
rz(-1.8788726) q[2];
sx q[2];
rz(-3.0335226) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97164153) q[1];
sx q[1];
rz(-1.2593237) q[1];
sx q[1];
rz(0.71572742) q[1];
rz(-pi) q[2];
rz(1.9196687) q[3];
sx q[3];
rz(-1.2563946) q[3];
sx q[3];
rz(-1.507198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85285774) q[2];
sx q[2];
rz(-1.4790269) q[2];
sx q[2];
rz(-0.23957254) q[2];
rz(2.1824956) q[3];
sx q[3];
rz(-1.9775016) q[3];
sx q[3];
rz(-1.0378999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.1921444) q[0];
sx q[0];
rz(-1.0553772) q[0];
sx q[0];
rz(1.5978093) q[0];
rz(-1.1100618) q[1];
sx q[1];
rz(-1.203457) q[1];
sx q[1];
rz(-1.4470626) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8657762) q[0];
sx q[0];
rz(-1.5779243) q[0];
sx q[0];
rz(-2.843186) q[0];
rz(-pi) q[1];
rz(-1.6844774) q[2];
sx q[2];
rz(-2.7964513) q[2];
sx q[2];
rz(-0.98649401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6441783) q[1];
sx q[1];
rz(-1.7458054) q[1];
sx q[1];
rz(2.464556) q[1];
x q[2];
rz(0.54357678) q[3];
sx q[3];
rz(-1.936562) q[3];
sx q[3];
rz(0.32777946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99628228) q[2];
sx q[2];
rz(-2.0029533) q[2];
sx q[2];
rz(-0.81926695) q[2];
rz(-0.94775689) q[3];
sx q[3];
rz(-2.3534333) q[3];
sx q[3];
rz(1.88571) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6325697) q[0];
sx q[0];
rz(-0.13795723) q[0];
sx q[0];
rz(-2.5914958) q[0];
rz(0.24712786) q[1];
sx q[1];
rz(-1.1526266) q[1];
sx q[1];
rz(-2.1955042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9086015) q[0];
sx q[0];
rz(-2.6938022) q[0];
sx q[0];
rz(-0.93771387) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3576415) q[2];
sx q[2];
rz(-1.6317131) q[2];
sx q[2];
rz(2.3062381) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.94080776) q[1];
sx q[1];
rz(-2.0403123) q[1];
sx q[1];
rz(-1.6108617) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4660353) q[3];
sx q[3];
rz(-1.3332061) q[3];
sx q[3];
rz(-1.4723569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3863824) q[2];
sx q[2];
rz(-0.55042616) q[2];
sx q[2];
rz(-2.4605301) q[2];
rz(-2.3303473) q[3];
sx q[3];
rz(-1.0439876) q[3];
sx q[3];
rz(-0.55027858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022245971) q[0];
sx q[0];
rz(-2.5211054) q[0];
sx q[0];
rz(2.3959809) q[0];
rz(1.9781808) q[1];
sx q[1];
rz(-2.5201576) q[1];
sx q[1];
rz(-0.17722873) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39800571) q[0];
sx q[0];
rz(-1.9015802) q[0];
sx q[0];
rz(0.28689204) q[0];
rz(-pi) q[1];
rz(3.0783098) q[2];
sx q[2];
rz(-1.6613792) q[2];
sx q[2];
rz(-0.16646756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2043269) q[1];
sx q[1];
rz(-0.64548739) q[1];
sx q[1];
rz(2.2068404) q[1];
x q[2];
rz(-1.528622) q[3];
sx q[3];
rz(-0.87648731) q[3];
sx q[3];
rz(0.75533849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54997286) q[2];
sx q[2];
rz(-2.1592996) q[2];
sx q[2];
rz(-1.6922692) q[2];
rz(-0.6226522) q[3];
sx q[3];
rz(-2.6144274) q[3];
sx q[3];
rz(1.3194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016668884) q[0];
sx q[0];
rz(-1.7542087) q[0];
sx q[0];
rz(1.6296052) q[0];
rz(0.92388693) q[1];
sx q[1];
rz(-1.4199363) q[1];
sx q[1];
rz(-0.61663827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23417191) q[0];
sx q[0];
rz(-1.4185832) q[0];
sx q[0];
rz(2.175075) q[0];
x q[1];
rz(-0.086313649) q[2];
sx q[2];
rz(-0.46937916) q[2];
sx q[2];
rz(0.98938194) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0860446) q[1];
sx q[1];
rz(-1.6564072) q[1];
sx q[1];
rz(2.5216091) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8455335) q[3];
sx q[3];
rz(-1.0062075) q[3];
sx q[3];
rz(-2.4170411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3672678) q[2];
sx q[2];
rz(-0.35949817) q[2];
sx q[2];
rz(-2.9765339) q[2];
rz(0.74602357) q[3];
sx q[3];
rz(-1.6227928) q[3];
sx q[3];
rz(-3.0987958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54166722) q[0];
sx q[0];
rz(-0.2921108) q[0];
sx q[0];
rz(2.7274729) q[0];
rz(-2.7060624) q[1];
sx q[1];
rz(-2.6256517) q[1];
sx q[1];
rz(-1.863265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0062388) q[0];
sx q[0];
rz(-1.4246294) q[0];
sx q[0];
rz(1.8144905) q[0];
x q[1];
rz(2.8299471) q[2];
sx q[2];
rz(-1.7362744) q[2];
sx q[2];
rz(2.8143425) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8929249) q[1];
sx q[1];
rz(-1.5589412) q[1];
sx q[1];
rz(-0.067889386) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0429562) q[3];
sx q[3];
rz(-1.5164638) q[3];
sx q[3];
rz(2.899037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4182959) q[2];
sx q[2];
rz(-1.7593242) q[2];
sx q[2];
rz(0.74271512) q[2];
rz(2.3501979) q[3];
sx q[3];
rz(-0.68358889) q[3];
sx q[3];
rz(1.7323823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.047569711) q[0];
sx q[0];
rz(-0.3293193) q[0];
sx q[0];
rz(0.91162115) q[0];
rz(2.1564663) q[1];
sx q[1];
rz(-2.7148425) q[1];
sx q[1];
rz(-1.8034579) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1250004) q[0];
sx q[0];
rz(-0.58412281) q[0];
sx q[0];
rz(-2.16433) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9495487) q[2];
sx q[2];
rz(-1.2745924) q[2];
sx q[2];
rz(2.6282994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9836925) q[1];
sx q[1];
rz(-2.165598) q[1];
sx q[1];
rz(0.48911323) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84566859) q[3];
sx q[3];
rz(-0.4825497) q[3];
sx q[3];
rz(-3.1180814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7868339) q[2];
sx q[2];
rz(-1.6666731) q[2];
sx q[2];
rz(-1.1653398) q[2];
rz(1.745584) q[3];
sx q[3];
rz(-1.952012) q[3];
sx q[3];
rz(1.5497807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.5042481) q[0];
sx q[0];
rz(-1.5865542) q[0];
sx q[0];
rz(-2.4028548) q[0];
rz(-0.63738102) q[1];
sx q[1];
rz(-1.5370054) q[1];
sx q[1];
rz(-0.76269033) q[1];
rz(1.1470096) q[2];
sx q[2];
rz(-1.4374566) q[2];
sx q[2];
rz(-0.79604436) q[2];
rz(-2.4828667) q[3];
sx q[3];
rz(-2.4729508) q[3];
sx q[3];
rz(-2.0153468) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
