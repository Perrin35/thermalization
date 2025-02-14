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
rz(-0.65547216) q[0];
sx q[0];
rz(-2.1894426) q[0];
sx q[0];
rz(0.093753554) q[0];
rz(1.3682415) q[1];
sx q[1];
rz(-1.5803087) q[1];
sx q[1];
rz(-0.66722792) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038644636) q[0];
sx q[0];
rz(-1.1955452) q[0];
sx q[0];
rz(1.0374336) q[0];
x q[1];
rz(2.480443) q[2];
sx q[2];
rz(-0.3232269) q[2];
sx q[2];
rz(1.8913392) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88244438) q[1];
sx q[1];
rz(-1.8516774) q[1];
sx q[1];
rz(-0.85390635) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62105234) q[3];
sx q[3];
rz(-0.9844616) q[3];
sx q[3];
rz(2.1042657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5253456) q[2];
sx q[2];
rz(-1.6051925) q[2];
sx q[2];
rz(3.1180535) q[2];
rz(-0.25447887) q[3];
sx q[3];
rz(-1.1602217) q[3];
sx q[3];
rz(-2.3833073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1940521) q[0];
sx q[0];
rz(-1.3964615) q[0];
sx q[0];
rz(2.5260455) q[0];
rz(2.5415892) q[1];
sx q[1];
rz(-1.2702076) q[1];
sx q[1];
rz(-0.34872762) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48744142) q[0];
sx q[0];
rz(-2.2117227) q[0];
sx q[0];
rz(-1.6665672) q[0];
x q[1];
rz(-1.2288413) q[2];
sx q[2];
rz(-1.8445024) q[2];
sx q[2];
rz(3.0897763) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.80417577) q[1];
sx q[1];
rz(-1.8693718) q[1];
sx q[1];
rz(1.0945303) q[1];
rz(-pi) q[2];
rz(2.8564151) q[3];
sx q[3];
rz(-2.0582407) q[3];
sx q[3];
rz(-2.2219258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49150026) q[2];
sx q[2];
rz(-1.3050175) q[2];
sx q[2];
rz(2.4864062) q[2];
rz(0.16264597) q[3];
sx q[3];
rz(-2.197367) q[3];
sx q[3];
rz(2.9338525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96838897) q[0];
sx q[0];
rz(-2.0252616) q[0];
sx q[0];
rz(-0.094245687) q[0];
rz(-1.3000129) q[1];
sx q[1];
rz(-2.2424707) q[1];
sx q[1];
rz(-1.947044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.081095) q[0];
sx q[0];
rz(-2.0682242) q[0];
sx q[0];
rz(1.2242076) q[0];
x q[1];
rz(-1.7060227) q[2];
sx q[2];
rz(-0.36721729) q[2];
sx q[2];
rz(1.4892422) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49385168) q[1];
sx q[1];
rz(-1.5268699) q[1];
sx q[1];
rz(-1.3141339) q[1];
x q[2];
rz(2.7500509) q[3];
sx q[3];
rz(-2.3394008) q[3];
sx q[3];
rz(1.4456309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.82075787) q[2];
sx q[2];
rz(-1.6944378) q[2];
sx q[2];
rz(0.40795946) q[2];
rz(-1.7631433) q[3];
sx q[3];
rz(-0.89865509) q[3];
sx q[3];
rz(0.11817008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5648062) q[0];
sx q[0];
rz(-1.7498359) q[0];
sx q[0];
rz(0.80192649) q[0];
rz(2.7469514) q[1];
sx q[1];
rz(-1.0486187) q[1];
sx q[1];
rz(-0.035042979) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5749767) q[0];
sx q[0];
rz(-2.1144951) q[0];
sx q[0];
rz(-2.7046142) q[0];
rz(-pi) q[1];
rz(1.3801244) q[2];
sx q[2];
rz(-1.9004656) q[2];
sx q[2];
rz(-0.53527385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1815419) q[1];
sx q[1];
rz(-1.6142577) q[1];
sx q[1];
rz(-0.04219136) q[1];
x q[2];
rz(0.48864103) q[3];
sx q[3];
rz(-2.7436769) q[3];
sx q[3];
rz(-0.15109381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0798215) q[2];
sx q[2];
rz(-1.062919) q[2];
sx q[2];
rz(1.7064095) q[2];
rz(-1.7363413) q[3];
sx q[3];
rz(-2.0515714) q[3];
sx q[3];
rz(2.0315571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6768796) q[0];
sx q[0];
rz(-1.1486624) q[0];
sx q[0];
rz(-1.1418463) q[0];
rz(-2.6308718) q[1];
sx q[1];
rz(-1.2341713) q[1];
sx q[1];
rz(1.4265192) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1349484) q[0];
sx q[0];
rz(-2.0922959) q[0];
sx q[0];
rz(0.64158507) q[0];
x q[1];
rz(-1.2586589) q[2];
sx q[2];
rz(-1.3545879) q[2];
sx q[2];
rz(2.5463111) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99014501) q[1];
sx q[1];
rz(-1.6832608) q[1];
sx q[1];
rz(-0.97655762) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1633888) q[3];
sx q[3];
rz(-1.4168882) q[3];
sx q[3];
rz(-0.28095442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4003754) q[2];
sx q[2];
rz(-1.8343238) q[2];
sx q[2];
rz(-0.2571787) q[2];
rz(2.2687965) q[3];
sx q[3];
rz(-1.3215439) q[3];
sx q[3];
rz(1.1375554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208825) q[0];
sx q[0];
rz(-2.6873984) q[0];
sx q[0];
rz(-2.0632451) q[0];
rz(0.9067761) q[1];
sx q[1];
rz(-0.6858784) q[1];
sx q[1];
rz(0.041898601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8804733) q[0];
sx q[0];
rz(-0.47396892) q[0];
sx q[0];
rz(1.649545) q[0];
x q[1];
rz(2.2380377) q[2];
sx q[2];
rz(-1.4741401) q[2];
sx q[2];
rz(0.54723155) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4087569) q[1];
sx q[1];
rz(-2.6799767) q[1];
sx q[1];
rz(-2.2240586) q[1];
rz(-pi) q[2];
rz(-0.53421212) q[3];
sx q[3];
rz(-0.90972661) q[3];
sx q[3];
rz(-0.40701696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28356734) q[2];
sx q[2];
rz(-2.4391386) q[2];
sx q[2];
rz(-2.2502327) q[2];
rz(0.71409613) q[3];
sx q[3];
rz(-2.4024506) q[3];
sx q[3];
rz(2.0785296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5080268) q[0];
sx q[0];
rz(-1.897568) q[0];
sx q[0];
rz(2.4666393) q[0];
rz(-1.630111) q[1];
sx q[1];
rz(-2.5210896) q[1];
sx q[1];
rz(-0.57704467) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70449146) q[0];
sx q[0];
rz(-0.54973212) q[0];
sx q[0];
rz(2.0389407) q[0];
rz(0.43319019) q[2];
sx q[2];
rz(-1.4775039) q[2];
sx q[2];
rz(0.46822883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8220209) q[1];
sx q[1];
rz(-2.4720925) q[1];
sx q[1];
rz(1.4386402) q[1];
rz(-0.87793276) q[3];
sx q[3];
rz(-2.0725277) q[3];
sx q[3];
rz(0.89661613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6881037) q[2];
sx q[2];
rz(-1.1223531) q[2];
sx q[2];
rz(0.92448676) q[2];
rz(1.9214572) q[3];
sx q[3];
rz(-2.852738) q[3];
sx q[3];
rz(-0.060062241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39153758) q[0];
sx q[0];
rz(-2.4704762) q[0];
sx q[0];
rz(1.7975988) q[0];
rz(1.9186107) q[1];
sx q[1];
rz(-1.5822625) q[1];
sx q[1];
rz(1.5917684) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61448288) q[0];
sx q[0];
rz(-1.7507995) q[0];
sx q[0];
rz(-1.2935782) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28002589) q[2];
sx q[2];
rz(-0.38869263) q[2];
sx q[2];
rz(-0.9946781) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0613411) q[1];
sx q[1];
rz(-2.3480573) q[1];
sx q[1];
rz(2.6166051) q[1];
rz(-pi) q[2];
rz(0.76390169) q[3];
sx q[3];
rz(-1.6391067) q[3];
sx q[3];
rz(2.1602283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0595155) q[2];
sx q[2];
rz(-1.0249219) q[2];
sx q[2];
rz(2.9151741) q[2];
rz(-0.039479937) q[3];
sx q[3];
rz(-1.6624781) q[3];
sx q[3];
rz(1.1994908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74462849) q[0];
sx q[0];
rz(-1.9292984) q[0];
sx q[0];
rz(-0.82897559) q[0];
rz(-3.0204311) q[1];
sx q[1];
rz(-2.1794901) q[1];
sx q[1];
rz(1.4519579) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4273116) q[0];
sx q[0];
rz(-1.7085008) q[0];
sx q[0];
rz(-2.4192823) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0214861) q[2];
sx q[2];
rz(-0.7996489) q[2];
sx q[2];
rz(-2.9422408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.87178225) q[1];
sx q[1];
rz(-0.93913684) q[1];
sx q[1];
rz(0.58806432) q[1];
rz(-pi) q[2];
rz(-2.1690458) q[3];
sx q[3];
rz(-0.64278162) q[3];
sx q[3];
rz(1.9721667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5539603) q[2];
sx q[2];
rz(-2.3162737) q[2];
sx q[2];
rz(0.86453214) q[2];
rz(-0.65166059) q[3];
sx q[3];
rz(-1.0385907) q[3];
sx q[3];
rz(0.017875044) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.637735) q[0];
sx q[0];
rz(-2.2561769) q[0];
sx q[0];
rz(0.46022415) q[0];
rz(1.3453311) q[1];
sx q[1];
rz(-2.5049152) q[1];
sx q[1];
rz(-1.771079) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84708285) q[0];
sx q[0];
rz(-2.6381603) q[0];
sx q[0];
rz(2.2398021) q[0];
rz(-2.7659594) q[2];
sx q[2];
rz(-2.4690095) q[2];
sx q[2];
rz(-3.1194558) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3440086) q[1];
sx q[1];
rz(-1.3329778) q[1];
sx q[1];
rz(1.306755) q[1];
x q[2];
rz(-1.7779782) q[3];
sx q[3];
rz(-1.6231114) q[3];
sx q[3];
rz(-2.4567043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7728077) q[2];
sx q[2];
rz(-1.7926755) q[2];
sx q[2];
rz(0.14300145) q[2];
rz(2.9755106) q[3];
sx q[3];
rz(-0.69286418) q[3];
sx q[3];
rz(-2.3486229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-2.3731257) q[0];
sx q[0];
rz(-1.65092) q[0];
sx q[0];
rz(-1.7478818) q[0];
rz(1.985818) q[1];
sx q[1];
rz(-2.0230237) q[1];
sx q[1];
rz(0.3442234) q[1];
rz(-2.5603603) q[2];
sx q[2];
rz(-2.317133) q[2];
sx q[2];
rz(1.9036507) q[2];
rz(-0.96859453) q[3];
sx q[3];
rz(-1.852949) q[3];
sx q[3];
rz(-0.52976086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
