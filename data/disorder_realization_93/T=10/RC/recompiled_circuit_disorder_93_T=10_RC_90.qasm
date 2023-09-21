OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(4.5259024) q[0];
sx q[0];
rz(10.685267) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7605654) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(-2.7678124) q[0];
rz(-pi) q[1];
x q[1];
rz(1.107723) q[2];
sx q[2];
rz(-2.827364) q[2];
sx q[2];
rz(0.83480922) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0337861) q[1];
sx q[1];
rz(-2.4596655) q[1];
sx q[1];
rz(-2.4297907) q[1];
rz(0.49433319) q[3];
sx q[3];
rz(-1.9088073) q[3];
sx q[3];
rz(0.61275834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(-2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0682003) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(0.67990047) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(1.4555567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50549492) q[0];
sx q[0];
rz(-0.60244766) q[0];
sx q[0];
rz(-1.2732182) q[0];
x q[1];
rz(-2.935264) q[2];
sx q[2];
rz(-2.7596139) q[2];
sx q[2];
rz(2.1260335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.204139) q[1];
sx q[1];
rz(-1.0480282) q[1];
sx q[1];
rz(3.1259895) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7299269) q[3];
sx q[3];
rz(-1.5951612) q[3];
sx q[3];
rz(1.0893694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7130647) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-0.18243608) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(0.3119719) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75333726) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.172539) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7045672) q[0];
sx q[0];
rz(-0.66172681) q[0];
sx q[0];
rz(-0.5163124) q[0];
rz(2.3861888) q[2];
sx q[2];
rz(-1.3284151) q[2];
sx q[2];
rz(0.6005477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8923924) q[1];
sx q[1];
rz(-1.3274267) q[1];
sx q[1];
rz(-2.1057486) q[1];
x q[2];
rz(-0.058733744) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(0.55610031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0744434) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(-2.6180843) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9273705) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(-0.58332304) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24387118) q[2];
sx q[2];
rz(-0.3393617) q[2];
sx q[2];
rz(-1.7983758) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6701811) q[1];
sx q[1];
rz(-1.3172611) q[1];
sx q[1];
rz(0.86004911) q[1];
rz(0.13682271) q[3];
sx q[3];
rz(-0.98383437) q[3];
sx q[3];
rz(2.6329991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52544242) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(0.564044) q[2];
rz(-0.28856746) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(2.585876) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-0.99194828) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6447727) q[0];
sx q[0];
rz(-2.9368375) q[0];
sx q[0];
rz(2.5909008) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8736585) q[2];
sx q[2];
rz(-1.5622683) q[2];
sx q[2];
rz(-2.0451343) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76584133) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(0.30893107) q[1];
x q[2];
rz(1.4335853) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(-1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1296967) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(2.9233542) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(-1.7165002) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544918) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(1.3160734) q[0];
x q[1];
rz(3.0503057) q[2];
sx q[2];
rz(-1.316615) q[2];
sx q[2];
rz(-2.9448178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1849991) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(-2.0208298) q[1];
rz(0.85961996) q[3];
sx q[3];
rz(-0.73577995) q[3];
sx q[3];
rz(-1.3501292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(1.5527027) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(2.3197876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215295) q[0];
sx q[0];
rz(-0.95525817) q[0];
sx q[0];
rz(-2.2277742) q[0];
x q[1];
rz(0.49634883) q[2];
sx q[2];
rz(-0.90663547) q[2];
sx q[2];
rz(-1.4585444) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51470876) q[1];
sx q[1];
rz(-0.85299546) q[1];
sx q[1];
rz(-1.8741329) q[1];
x q[2];
rz(1.9415226) q[3];
sx q[3];
rz(-2.126174) q[3];
sx q[3];
rz(2.6851418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3372779) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(-0.24469963) q[2];
rz(-3.0120567) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(-1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.8364505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7628521) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(0.61600323) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4370286) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(-1.2003843) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9040363) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(1.1035641) q[1];
rz(3.1180624) q[3];
sx q[3];
rz(-1.9100034) q[3];
sx q[3];
rz(-2.6587405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.7548783) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(1.6053258) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(0.70294356) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84530572) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(-1.9592497) q[0];
rz(-1.6417575) q[2];
sx q[2];
rz(-0.5103726) q[2];
sx q[2];
rz(0.77529782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.326509) q[1];
sx q[1];
rz(-1.0352967) q[1];
sx q[1];
rz(0.19170796) q[1];
rz(0.33396696) q[3];
sx q[3];
rz(-0.98200646) q[3];
sx q[3];
rz(-0.86563084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90074173) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(-1.9158069) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43384957) q[0];
sx q[0];
rz(-1.187547) q[0];
sx q[0];
rz(-1.3462523) q[0];
rz(2.8154545) q[2];
sx q[2];
rz(-1.9378127) q[2];
sx q[2];
rz(2.0740167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84506449) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(1.3045842) q[1];
x q[2];
rz(-1.1115132) q[3];
sx q[3];
rz(-1.1772403) q[3];
sx q[3];
rz(-2.4126088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(1.998385) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-2.2767699) q[2];
sx q[2];
rz(-1.0100126) q[2];
sx q[2];
rz(1.0457912) q[2];
rz(-0.91924304) q[3];
sx q[3];
rz(-1.5083434) q[3];
sx q[3];
rz(-1.2283243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];