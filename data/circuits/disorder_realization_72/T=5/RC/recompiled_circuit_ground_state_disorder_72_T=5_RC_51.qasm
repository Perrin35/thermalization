OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44719625) q[0];
sx q[0];
rz(1.6262654) q[0];
sx q[0];
rz(11.716421) q[0];
rz(1.3810459) q[1];
sx q[1];
rz(-0.84264207) q[1];
sx q[1];
rz(0.050328644) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6724191) q[0];
sx q[0];
rz(-0.19939961) q[0];
sx q[0];
rz(-0.43740793) q[0];
x q[1];
rz(2.1904616) q[2];
sx q[2];
rz(-1.4339087) q[2];
sx q[2];
rz(2.3093683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.58568629) q[1];
sx q[1];
rz(-1.3765322) q[1];
sx q[1];
rz(-1.4997086) q[1];
rz(-pi) q[2];
rz(1.5123472) q[3];
sx q[3];
rz(-2.2104467) q[3];
sx q[3];
rz(2.9906038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7488659) q[2];
sx q[2];
rz(-1.8957081) q[2];
sx q[2];
rz(-2.5457814) q[2];
rz(2.3303604) q[3];
sx q[3];
rz(-1.2657284) q[3];
sx q[3];
rz(0.59734145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3101462) q[0];
sx q[0];
rz(-2.6409288) q[0];
sx q[0];
rz(-1.798382) q[0];
rz(-1.385618) q[1];
sx q[1];
rz(-2.0090397) q[1];
sx q[1];
rz(-1.064942) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32414054) q[0];
sx q[0];
rz(-1.5498501) q[0];
sx q[0];
rz(-3.0684242) q[0];
rz(-pi) q[1];
rz(-1.0491514) q[2];
sx q[2];
rz(-1.8966833) q[2];
sx q[2];
rz(-0.87363315) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8044195) q[1];
sx q[1];
rz(-1.7111756) q[1];
sx q[1];
rz(-1.5650723) q[1];
x q[2];
rz(-0.029964793) q[3];
sx q[3];
rz(-1.6707722) q[3];
sx q[3];
rz(-2.4520017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.500835) q[2];
sx q[2];
rz(-1.0229599) q[2];
sx q[2];
rz(-0.41066059) q[2];
rz(1.9985577) q[3];
sx q[3];
rz(-0.7468907) q[3];
sx q[3];
rz(0.66543287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2086585) q[0];
sx q[0];
rz(-1.7993131) q[0];
sx q[0];
rz(0.68516532) q[0];
rz(-0.67600983) q[1];
sx q[1];
rz(-1.6513377) q[1];
sx q[1];
rz(0.73838678) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3078416) q[0];
sx q[0];
rz(-2.0594119) q[0];
sx q[0];
rz(1.9355544) q[0];
rz(-pi) q[1];
rz(-1.5398847) q[2];
sx q[2];
rz(-2.3031903) q[2];
sx q[2];
rz(0.67275999) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8505212) q[1];
sx q[1];
rz(-1.5187361) q[1];
sx q[1];
rz(1.4223767) q[1];
rz(-pi) q[2];
rz(0.66716828) q[3];
sx q[3];
rz(-0.83339429) q[3];
sx q[3];
rz(-1.7478314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43397778) q[2];
sx q[2];
rz(-1.4716163) q[2];
sx q[2];
rz(0.99226704) q[2];
rz(-2.9662002) q[3];
sx q[3];
rz(-2.6447191) q[3];
sx q[3];
rz(2.8309256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6076412) q[0];
sx q[0];
rz(-2.0168004) q[0];
sx q[0];
rz(-0.24096179) q[0];
rz(0.85722771) q[1];
sx q[1];
rz(-2.3576184) q[1];
sx q[1];
rz(1.8704174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64528148) q[0];
sx q[0];
rz(-1.4110312) q[0];
sx q[0];
rz(2.0156572) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0849146) q[2];
sx q[2];
rz(-2.0419697) q[2];
sx q[2];
rz(-0.65907329) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1041455) q[1];
sx q[1];
rz(-2.0097783) q[1];
sx q[1];
rz(0.46260117) q[1];
rz(-pi) q[2];
rz(-1.1587093) q[3];
sx q[3];
rz(-2.2818687) q[3];
sx q[3];
rz(2.3726316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6406389) q[2];
sx q[2];
rz(-2.1827938) q[2];
sx q[2];
rz(-1.0183838) q[2];
rz(-1.9765249) q[3];
sx q[3];
rz(-2.1924721) q[3];
sx q[3];
rz(-0.44786662) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1027706) q[0];
sx q[0];
rz(-2.3508681) q[0];
sx q[0];
rz(-3.0550509) q[0];
rz(1.4718919) q[1];
sx q[1];
rz(-2.4457928) q[1];
sx q[1];
rz(2.5873628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3329268) q[0];
sx q[0];
rz(-0.73360591) q[0];
sx q[0];
rz(0.33302078) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7090899) q[2];
sx q[2];
rz(-0.90511647) q[2];
sx q[2];
rz(-2.4709357) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2947301) q[1];
sx q[1];
rz(-2.1536441) q[1];
sx q[1];
rz(1.0148878) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1301) q[3];
sx q[3];
rz(-1.8365068) q[3];
sx q[3];
rz(1.2540224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8467796) q[2];
sx q[2];
rz(-0.3598992) q[2];
sx q[2];
rz(2.2100718) q[2];
rz(0.31342634) q[3];
sx q[3];
rz(-0.64160186) q[3];
sx q[3];
rz(-2.6581367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.038092) q[0];
sx q[0];
rz(-2.3374538) q[0];
sx q[0];
rz(-1.7622129) q[0];
rz(-2.3992959) q[1];
sx q[1];
rz(-1.392044) q[1];
sx q[1];
rz(-1.0662063) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48165769) q[0];
sx q[0];
rz(-1.5391506) q[0];
sx q[0];
rz(3.0040347) q[0];
x q[1];
rz(3.0020038) q[2];
sx q[2];
rz(-1.0537647) q[2];
sx q[2];
rz(1.6986183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6247793) q[1];
sx q[1];
rz(-1.32059) q[1];
sx q[1];
rz(-2.183319) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4857852) q[3];
sx q[3];
rz(-1.8411311) q[3];
sx q[3];
rz(0.47209376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1833056) q[2];
sx q[2];
rz(-0.97465193) q[2];
sx q[2];
rz(-2.8948114) q[2];
rz(2.0731481) q[3];
sx q[3];
rz(-1.9717792) q[3];
sx q[3];
rz(1.0698414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088155) q[0];
sx q[0];
rz(-2.1822073) q[0];
sx q[0];
rz(2.7799613) q[0];
rz(1.5879141) q[1];
sx q[1];
rz(-1.5767153) q[1];
sx q[1];
rz(1.9379001) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7772352) q[0];
sx q[0];
rz(-0.89221749) q[0];
sx q[0];
rz(2.4410072) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33233099) q[2];
sx q[2];
rz(-2.6164101) q[2];
sx q[2];
rz(-3.0572403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1275529) q[1];
sx q[1];
rz(-1.3689835) q[1];
sx q[1];
rz(2.1743808) q[1];
rz(-1.2343812) q[3];
sx q[3];
rz(-1.9877476) q[3];
sx q[3];
rz(0.41301504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92538539) q[2];
sx q[2];
rz(-1.6758827) q[2];
sx q[2];
rz(-2.1290131) q[2];
rz(-2.8818164) q[3];
sx q[3];
rz(-2.8949819) q[3];
sx q[3];
rz(-0.40464211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15553661) q[0];
sx q[0];
rz(-1.8590834) q[0];
sx q[0];
rz(0.78722659) q[0];
rz(2.2975445) q[1];
sx q[1];
rz(-2.7828352) q[1];
sx q[1];
rz(1.2740096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9106261) q[0];
sx q[0];
rz(-1.2342802) q[0];
sx q[0];
rz(-1.6331768) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73972685) q[2];
sx q[2];
rz(-1.8056895) q[2];
sx q[2];
rz(-1.7559831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2383409) q[1];
sx q[1];
rz(-2.1834186) q[1];
sx q[1];
rz(-2.5765221) q[1];
x q[2];
rz(-2.9926821) q[3];
sx q[3];
rz(-2.732548) q[3];
sx q[3];
rz(0.27924505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0969703) q[2];
sx q[2];
rz(-1.5644093) q[2];
sx q[2];
rz(-0.38193646) q[2];
rz(2.0197798) q[3];
sx q[3];
rz(-2.8048281) q[3];
sx q[3];
rz(0.063966123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065351039) q[0];
sx q[0];
rz(-0.92472804) q[0];
sx q[0];
rz(1.4724154) q[0];
rz(0.30638254) q[1];
sx q[1];
rz(-0.52426052) q[1];
sx q[1];
rz(-2.870097) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58518325) q[0];
sx q[0];
rz(-2.5228254) q[0];
sx q[0];
rz(-2.2947001) q[0];
rz(-pi) q[1];
rz(-2.3112663) q[2];
sx q[2];
rz(-1.8486075) q[2];
sx q[2];
rz(-0.72208126) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4541536) q[1];
sx q[1];
rz(-1.0311797) q[1];
sx q[1];
rz(-2.1381049) q[1];
x q[2];
rz(-0.90320194) q[3];
sx q[3];
rz(-2.4861504) q[3];
sx q[3];
rz(0.02070657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.86758119) q[2];
sx q[2];
rz(-2.0143955) q[2];
sx q[2];
rz(2.4809044) q[2];
rz(0.88179669) q[3];
sx q[3];
rz(-2.5776358) q[3];
sx q[3];
rz(2.2753415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26447403) q[0];
sx q[0];
rz(-1.6082123) q[0];
sx q[0];
rz(-1.3382753) q[0];
rz(1.0461461) q[1];
sx q[1];
rz(-2.1144861) q[1];
sx q[1];
rz(-0.17790067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.698211) q[0];
sx q[0];
rz(-0.92635768) q[0];
sx q[0];
rz(0.81581913) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98921241) q[2];
sx q[2];
rz(-1.4228369) q[2];
sx q[2];
rz(-1.4799679) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5452176) q[1];
sx q[1];
rz(-1.8918953) q[1];
sx q[1];
rz(-1.9209747) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37306656) q[3];
sx q[3];
rz(-1.698768) q[3];
sx q[3];
rz(-1.7440918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.88343128) q[2];
sx q[2];
rz(-1.0940172) q[2];
sx q[2];
rz(2.2056313) q[2];
rz(-0.59518138) q[3];
sx q[3];
rz(-1.429957) q[3];
sx q[3];
rz(2.2619251) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7137939) q[0];
sx q[0];
rz(-1.6811163) q[0];
sx q[0];
rz(1.2291193) q[0];
rz(-1.3270558) q[1];
sx q[1];
rz(-1.5057024) q[1];
sx q[1];
rz(1.1485752) q[1];
rz(1.003895) q[2];
sx q[2];
rz(-0.2093506) q[2];
sx q[2];
rz(1.9471526) q[2];
rz(-0.86092864) q[3];
sx q[3];
rz(-2.7459675) q[3];
sx q[3];
rz(-0.31386123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
