OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3888336) q[0];
sx q[0];
rz(-1.692481) q[0];
sx q[0];
rz(-1.4504855) q[0];
rz(0.9736355) q[1];
sx q[1];
rz(-1.7042301) q[1];
sx q[1];
rz(-0.91926423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5333119) q[0];
sx q[0];
rz(-1.5402147) q[0];
sx q[0];
rz(-1.5536867) q[0];
rz(-pi) q[1];
rz(3.0212901) q[2];
sx q[2];
rz(-1.4805111) q[2];
sx q[2];
rz(1.9419958) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1303006) q[1];
sx q[1];
rz(-0.71349547) q[1];
sx q[1];
rz(1.43047) q[1];
rz(0.84897016) q[3];
sx q[3];
rz(-0.59576407) q[3];
sx q[3];
rz(1.9879544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8093402) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.7285041) q[2];
rz(-2.938802) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(-0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8973812) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(-0.71075034) q[0];
rz(-0.76849014) q[1];
sx q[1];
rz(-2.0748506) q[1];
sx q[1];
rz(1.01952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7026414) q[0];
sx q[0];
rz(-2.0176689) q[0];
sx q[0];
rz(1.5269482) q[0];
x q[1];
rz(3.0654991) q[2];
sx q[2];
rz(-2.00092) q[2];
sx q[2];
rz(-0.96904749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4162035) q[1];
sx q[1];
rz(-2.7692911) q[1];
sx q[1];
rz(-1.6561693) q[1];
x q[2];
rz(-1.4367661) q[3];
sx q[3];
rz(-0.55826) q[3];
sx q[3];
rz(-2.9001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76356137) q[2];
sx q[2];
rz(-1.8969994) q[2];
sx q[2];
rz(-1.2223318) q[2];
rz(1.9289121) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(-0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0843622) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(1.2179751) q[0];
rz(-2.6257264) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(-2.2191494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9410011) q[0];
sx q[0];
rz(-0.070644826) q[0];
sx q[0];
rz(-2.5289422) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0751245) q[2];
sx q[2];
rz(-2.2309003) q[2];
sx q[2];
rz(-0.65048993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53239142) q[1];
sx q[1];
rz(-2.2431886) q[1];
sx q[1];
rz(2.3228541) q[1];
rz(-pi) q[2];
rz(1.4018784) q[3];
sx q[3];
rz(-1.2460105) q[3];
sx q[3];
rz(1.5320154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.018365232) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(-1.8017192) q[2];
rz(-0.54723048) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(-0.71119285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0860586) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(0.15922971) q[0];
rz(-3.1314462) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(-2.8841282) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.793593) q[0];
sx q[0];
rz(-1.9257015) q[0];
sx q[0];
rz(-0.76323842) q[0];
rz(2.3472559) q[2];
sx q[2];
rz(-1.6172501) q[2];
sx q[2];
rz(0.72512324) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3066669) q[1];
sx q[1];
rz(-1.4713227) q[1];
sx q[1];
rz(2.5132781) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80735029) q[3];
sx q[3];
rz(-2.0424543) q[3];
sx q[3];
rz(-2.014267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(-0.47284687) q[2];
rz(2.4371448) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(-2.4956467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5064297) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(0.065486431) q[0];
rz(-2.7323515) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(1.4261036) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5201841) q[0];
sx q[0];
rz(-1.4778293) q[0];
sx q[0];
rz(2.9220102) q[0];
rz(-pi) q[1];
x q[1];
rz(1.481856) q[2];
sx q[2];
rz(-1.9581902) q[2];
sx q[2];
rz(0.26814869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3196484) q[1];
sx q[1];
rz(-1.8488171) q[1];
sx q[1];
rz(1.8026428) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17284837) q[3];
sx q[3];
rz(-1.8573055) q[3];
sx q[3];
rz(2.2423173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9466729) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(-1.8348414) q[2];
rz(-1.9715747) q[3];
sx q[3];
rz(-0.40478671) q[3];
sx q[3];
rz(-3.034333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829247) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(-2.7401127) q[0];
rz(-0.67277706) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(0.85404095) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9977048) q[0];
sx q[0];
rz(-1.524462) q[0];
sx q[0];
rz(1.3530988) q[0];
rz(-2.9561192) q[2];
sx q[2];
rz(-1.4234241) q[2];
sx q[2];
rz(-0.81306785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8484162) q[1];
sx q[1];
rz(-2.5418315) q[1];
sx q[1];
rz(-2.9648215) q[1];
rz(-pi) q[2];
rz(1.1229418) q[3];
sx q[3];
rz(-2.1338226) q[3];
sx q[3];
rz(-2.8729968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.65471571) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(-1.1775449) q[2];
rz(-0.55142895) q[3];
sx q[3];
rz(-1.5588372) q[3];
sx q[3];
rz(-2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42171445) q[0];
sx q[0];
rz(-0.28110176) q[0];
sx q[0];
rz(-1.9792492) q[0];
rz(-1.1516736) q[1];
sx q[1];
rz(-1.78777) q[1];
sx q[1];
rz(-2.2231359) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3716482) q[0];
sx q[0];
rz(-1.3027096) q[0];
sx q[0];
rz(-1.9140585) q[0];
rz(2.8807441) q[2];
sx q[2];
rz(-2.2493304) q[2];
sx q[2];
rz(-1.4229753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9388401) q[1];
sx q[1];
rz(-1.563464) q[1];
sx q[1];
rz(1.5775024) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1226467) q[3];
sx q[3];
rz(-0.3780685) q[3];
sx q[3];
rz(0.44071769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1071757) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(2.6386063) q[2];
rz(-0.77477396) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(-1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4485432) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.5013129) q[0];
rz(-2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(-2.0223845) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6014746) q[0];
sx q[0];
rz(-1.0154503) q[0];
sx q[0];
rz(-1.1309654) q[0];
rz(-pi) q[1];
rz(2.0430492) q[2];
sx q[2];
rz(-2.611428) q[2];
sx q[2];
rz(1.4250172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70485605) q[1];
sx q[1];
rz(-0.67187998) q[1];
sx q[1];
rz(-0.19499548) q[1];
rz(-2.7469278) q[3];
sx q[3];
rz(-1.695444) q[3];
sx q[3];
rz(-0.11388557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1559653) q[2];
sx q[2];
rz(-0.95981821) q[2];
sx q[2];
rz(0.36925527) q[2];
rz(2.8333832) q[3];
sx q[3];
rz(-2.1741368) q[3];
sx q[3];
rz(-1.6109899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592598) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(0.34307137) q[0];
rz(1.169091) q[1];
sx q[1];
rz(-1.3645376) q[1];
sx q[1];
rz(1.7128568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0853303) q[0];
sx q[0];
rz(-1.3199727) q[0];
sx q[0];
rz(-0.40928264) q[0];
x q[1];
rz(-0.83447225) q[2];
sx q[2];
rz(-1.021046) q[2];
sx q[2];
rz(2.2817734) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6862168) q[1];
sx q[1];
rz(-1.1955402) q[1];
sx q[1];
rz(-0.065836716) q[1];
x q[2];
rz(2.6713624) q[3];
sx q[3];
rz(-2.2080126) q[3];
sx q[3];
rz(1.5772082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(-1.7500056) q[2];
rz(2.7348147) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(-0.87882915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10130356) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(-1.783675) q[0];
rz(1.2376002) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(-2.3416669) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185266) q[0];
sx q[0];
rz(-2.3481821) q[0];
sx q[0];
rz(1.2454459) q[0];
x q[1];
rz(-0.55028693) q[2];
sx q[2];
rz(-1.3105416) q[2];
sx q[2];
rz(1.257892) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0837896) q[1];
sx q[1];
rz(-0.68843319) q[1];
sx q[1];
rz(0.24569421) q[1];
rz(-pi) q[2];
rz(2.8915845) q[3];
sx q[3];
rz(-2.0431314) q[3];
sx q[3];
rz(-1.3046622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(3.0685032) q[2];
rz(-2.8857005) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(1.6528486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8681317) q[0];
sx q[0];
rz(-1.4556226) q[0];
sx q[0];
rz(1.8709394) q[0];
rz(-1.3399667) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(-3.1392787) q[2];
sx q[2];
rz(-2.4321767) q[2];
sx q[2];
rz(-1.3117758) q[2];
rz(2.1382016) q[3];
sx q[3];
rz(-1.0091253) q[3];
sx q[3];
rz(-0.71970018) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
