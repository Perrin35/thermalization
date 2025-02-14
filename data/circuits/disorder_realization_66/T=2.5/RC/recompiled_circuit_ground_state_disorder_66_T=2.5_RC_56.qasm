OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3341137) q[0];
sx q[0];
rz(-2.4456094) q[0];
sx q[0];
rz(0.88710436) q[0];
rz(-0.71169418) q[1];
sx q[1];
rz(-1.0625755) q[1];
sx q[1];
rz(-0.33023155) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6905211) q[0];
sx q[0];
rz(-2.1981648) q[0];
sx q[0];
rz(2.5894707) q[0];
rz(-pi) q[1];
rz(-0.55203931) q[2];
sx q[2];
rz(-0.62969724) q[2];
sx q[2];
rz(2.3384467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8793927) q[1];
sx q[1];
rz(-2.6385251) q[1];
sx q[1];
rz(-2.9285953) q[1];
rz(-pi) q[2];
rz(2.0014181) q[3];
sx q[3];
rz(-1.4994687) q[3];
sx q[3];
rz(0.50702107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7547138) q[2];
sx q[2];
rz(-2.7140996) q[2];
sx q[2];
rz(-1.0270366) q[2];
rz(2.2586281) q[3];
sx q[3];
rz(-1.3436147) q[3];
sx q[3];
rz(0.55356717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.3565732) q[0];
sx q[0];
rz(-3.0632601) q[0];
sx q[0];
rz(-2.8739492) q[0];
rz(1.74125) q[1];
sx q[1];
rz(-2.5282271) q[1];
sx q[1];
rz(-0.5563446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69795495) q[0];
sx q[0];
rz(-2.58444) q[0];
sx q[0];
rz(1.0865666) q[0];
rz(-pi) q[1];
rz(0.70425561) q[2];
sx q[2];
rz(-1.8374763) q[2];
sx q[2];
rz(-1.8530588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9832335) q[1];
sx q[1];
rz(-0.69770798) q[1];
sx q[1];
rz(-2.1844728) q[1];
rz(-pi) q[2];
rz(2.2675927) q[3];
sx q[3];
rz(-0.79588875) q[3];
sx q[3];
rz(0.1800285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23655015) q[2];
sx q[2];
rz(-1.6908129) q[2];
sx q[2];
rz(1.8897918) q[2];
rz(1.8768138) q[3];
sx q[3];
rz(-2.4301961) q[3];
sx q[3];
rz(1.0714162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65997893) q[0];
sx q[0];
rz(-0.55584207) q[0];
sx q[0];
rz(2.3023093) q[0];
rz(-1.5045769) q[1];
sx q[1];
rz(-1.7103651) q[1];
sx q[1];
rz(1.7617216) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8913075) q[0];
sx q[0];
rz(-1.6943185) q[0];
sx q[0];
rz(3.059292) q[0];
rz(1.9154869) q[2];
sx q[2];
rz(-0.57248964) q[2];
sx q[2];
rz(0.16333157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8306337) q[1];
sx q[1];
rz(-1.1006025) q[1];
sx q[1];
rz(0.36010919) q[1];
rz(-pi) q[2];
rz(-2.5702137) q[3];
sx q[3];
rz(-1.0927754) q[3];
sx q[3];
rz(0.92438053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.69924361) q[2];
sx q[2];
rz(-1.9244497) q[2];
sx q[2];
rz(-0.71845636) q[2];
rz(-2.5213304) q[3];
sx q[3];
rz(-0.93554997) q[3];
sx q[3];
rz(-0.81234318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45850596) q[0];
sx q[0];
rz(-1.9106671) q[0];
sx q[0];
rz(1.7339535) q[0];
rz(0.55157026) q[1];
sx q[1];
rz(-3.0307814) q[1];
sx q[1];
rz(-1.3161906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5410446) q[0];
sx q[0];
rz(-1.2886314) q[0];
sx q[0];
rz(-0.94775424) q[0];
x q[1];
rz(-1.8502683) q[2];
sx q[2];
rz(-1.5827465) q[2];
sx q[2];
rz(-0.751647) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0312457) q[1];
sx q[1];
rz(-1.7025885) q[1];
sx q[1];
rz(0.40459569) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4689581) q[3];
sx q[3];
rz(-0.97390538) q[3];
sx q[3];
rz(1.2850645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59336415) q[2];
sx q[2];
rz(-0.43539384) q[2];
sx q[2];
rz(2.1335874) q[2];
rz(0.27901444) q[3];
sx q[3];
rz(-2.3271826) q[3];
sx q[3];
rz(-2.6830955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27610436) q[0];
sx q[0];
rz(-1.8115598) q[0];
sx q[0];
rz(-2.8651067) q[0];
rz(-2.8522988) q[1];
sx q[1];
rz(-1.1951059) q[1];
sx q[1];
rz(0.2624661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2967178) q[0];
sx q[0];
rz(-2.8709559) q[0];
sx q[0];
rz(2.4225745) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5868411) q[2];
sx q[2];
rz(-0.59140059) q[2];
sx q[2];
rz(0.41447869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.561132) q[1];
sx q[1];
rz(-1.9101686) q[1];
sx q[1];
rz(-2.5575693) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5632202) q[3];
sx q[3];
rz(-1.6993853) q[3];
sx q[3];
rz(-1.8687539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7485973) q[2];
sx q[2];
rz(-2.5454919) q[2];
sx q[2];
rz(1.7945012) q[2];
rz(-0.10284452) q[3];
sx q[3];
rz(-1.9072073) q[3];
sx q[3];
rz(1.6245406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1384001) q[0];
sx q[0];
rz(-1.5035368) q[0];
sx q[0];
rz(0.060977161) q[0];
rz(1.2334476) q[1];
sx q[1];
rz(-1.3950709) q[1];
sx q[1];
rz(1.5256418) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73404445) q[0];
sx q[0];
rz(-1.5015232) q[0];
sx q[0];
rz(0.53014) q[0];
x q[1];
rz(1.5360918) q[2];
sx q[2];
rz(-0.74932171) q[2];
sx q[2];
rz(-2.9108436) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1970343) q[1];
sx q[1];
rz(-1.6813633) q[1];
sx q[1];
rz(-1.8227897) q[1];
rz(-2.9750729) q[3];
sx q[3];
rz(-0.92574471) q[3];
sx q[3];
rz(1.695961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6700651) q[2];
sx q[2];
rz(-1.3482886) q[2];
sx q[2];
rz(-1.8111551) q[2];
rz(1.1250251) q[3];
sx q[3];
rz(-2.2671813) q[3];
sx q[3];
rz(-0.037954656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9576981) q[0];
sx q[0];
rz(-1.4426458) q[0];
sx q[0];
rz(2.884602) q[0];
rz(-0.43486241) q[1];
sx q[1];
rz(-0.75006524) q[1];
sx q[1];
rz(0.90352568) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5241476) q[0];
sx q[0];
rz(-3.0492231) q[0];
sx q[0];
rz(1.6156107) q[0];
rz(0.74600474) q[2];
sx q[2];
rz(-1.6705319) q[2];
sx q[2];
rz(-1.6670645) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82680399) q[1];
sx q[1];
rz(-1.6182403) q[1];
sx q[1];
rz(-3.080009) q[1];
rz(-2.1210455) q[3];
sx q[3];
rz(-1.9686832) q[3];
sx q[3];
rz(-2.3031421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8187108) q[2];
sx q[2];
rz(-1.5078397) q[2];
sx q[2];
rz(1.4212849) q[2];
rz(2.4317702) q[3];
sx q[3];
rz(-2.0028508) q[3];
sx q[3];
rz(-0.53409725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1555971) q[0];
sx q[0];
rz(-0.43102145) q[0];
sx q[0];
rz(2.0565597) q[0];
rz(2.494508) q[1];
sx q[1];
rz(-1.9661463) q[1];
sx q[1];
rz(-0.064402493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2875117) q[0];
sx q[0];
rz(-0.72923822) q[0];
sx q[0];
rz(0.50042787) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4738381) q[2];
sx q[2];
rz(-0.96558324) q[2];
sx q[2];
rz(-2.3507694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6820087) q[1];
sx q[1];
rz(-1.1136076) q[1];
sx q[1];
rz(2.3383635) q[1];
rz(-pi) q[2];
rz(2.2256718) q[3];
sx q[3];
rz(-2.7034211) q[3];
sx q[3];
rz(2.8000591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.881968) q[2];
sx q[2];
rz(-2.7590064) q[2];
sx q[2];
rz(0.043370334) q[2];
rz(2.4671593) q[3];
sx q[3];
rz(-1.784227) q[3];
sx q[3];
rz(1.9969214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44716537) q[0];
sx q[0];
rz(-2.5612216) q[0];
sx q[0];
rz(0.32661435) q[0];
rz(-2.0475552) q[1];
sx q[1];
rz(-2.2973165) q[1];
sx q[1];
rz(-0.48386595) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9897468) q[0];
sx q[0];
rz(-0.92597472) q[0];
sx q[0];
rz(0.65412866) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6695374) q[2];
sx q[2];
rz(-2.6723571) q[2];
sx q[2];
rz(3.0735441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8344675) q[1];
sx q[1];
rz(-0.79453429) q[1];
sx q[1];
rz(1.267872) q[1];
rz(-pi) q[2];
rz(2.8362379) q[3];
sx q[3];
rz(-2.476446) q[3];
sx q[3];
rz(2.0572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6924374) q[2];
sx q[2];
rz(-1.038895) q[2];
sx q[2];
rz(0.078744002) q[2];
rz(-0.76830831) q[3];
sx q[3];
rz(-1.2945622) q[3];
sx q[3];
rz(0.19101846) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67097265) q[0];
sx q[0];
rz(-2.5363531) q[0];
sx q[0];
rz(-0.11669267) q[0];
rz(2.281588) q[1];
sx q[1];
rz(-1.8000894) q[1];
sx q[1];
rz(2.2681627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78055843) q[0];
sx q[0];
rz(-1.3300899) q[0];
sx q[0];
rz(-0.81674256) q[0];
rz(1.8928746) q[2];
sx q[2];
rz(-2.0587741) q[2];
sx q[2];
rz(0.26103324) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5705986) q[1];
sx q[1];
rz(-2.6254404) q[1];
sx q[1];
rz(0.25324778) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6386581) q[3];
sx q[3];
rz(-2.181119) q[3];
sx q[3];
rz(0.44059702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42150911) q[2];
sx q[2];
rz(-0.24622157) q[2];
sx q[2];
rz(2.1923547) q[2];
rz(1.2468437) q[3];
sx q[3];
rz(-1.4693762) q[3];
sx q[3];
rz(-2.5505572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0134037) q[0];
sx q[0];
rz(-1.0008151) q[0];
sx q[0];
rz(-1.6304954) q[0];
rz(-0.36698256) q[1];
sx q[1];
rz(-1.783168) q[1];
sx q[1];
rz(2.5617243) q[1];
rz(-1.0147711) q[2];
sx q[2];
rz(-0.53511878) q[2];
sx q[2];
rz(-2.6352885) q[2];
rz(1.1506769) q[3];
sx q[3];
rz(-2.2169866) q[3];
sx q[3];
rz(0.22056072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
