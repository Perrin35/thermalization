OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(-3.1414519) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1228468) q[0];
sx q[0];
rz(-1.393264) q[0];
sx q[0];
rz(-1.8621423) q[0];
rz(0.5483746) q[2];
sx q[2];
rz(-1.8273241) q[2];
sx q[2];
rz(2.246849) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27543435) q[1];
sx q[1];
rz(-0.94238102) q[1];
sx q[1];
rz(-0.98316146) q[1];
rz(1.3532577) q[3];
sx q[3];
rz(-1.4635524) q[3];
sx q[3];
rz(1.6579962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.45941916) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(1.9127282) q[2];
rz(1.7284283) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6035778) q[0];
sx q[0];
rz(-1.5025654) q[0];
sx q[0];
rz(-2.1287825) q[0];
rz(-0.027659841) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-2.0181296) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9119308) q[0];
sx q[0];
rz(-0.46470416) q[0];
sx q[0];
rz(1.4421411) q[0];
rz(1.6337109) q[2];
sx q[2];
rz(-2.3496369) q[2];
sx q[2];
rz(2.6346249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9065735) q[1];
sx q[1];
rz(-2.0779013) q[1];
sx q[1];
rz(-0.97096918) q[1];
x q[2];
rz(2.1334322) q[3];
sx q[3];
rz(-0.99346549) q[3];
sx q[3];
rz(-1.2853704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3479487) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(0.91903764) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27750257) q[0];
sx q[0];
rz(-2.9798177) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-0.69349849) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(-2.0085874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95110675) q[0];
sx q[0];
rz(-1.5740984) q[0];
sx q[0];
rz(-1.7130873) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79046952) q[2];
sx q[2];
rz(-1.1970453) q[2];
sx q[2];
rz(-2.8369396) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6269835) q[1];
sx q[1];
rz(-1.2830462) q[1];
sx q[1];
rz(-2.1327553) q[1];
rz(2.3660925) q[3];
sx q[3];
rz(-0.79458487) q[3];
sx q[3];
rz(-2.9426108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8901849) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(-0.039316468) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(1.1608634) q[0];
rz(-0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(0.13555759) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4281222) q[0];
sx q[0];
rz(-0.8849511) q[0];
sx q[0];
rz(-1.1629521) q[0];
x q[1];
rz(-3.0669078) q[2];
sx q[2];
rz(-1.1876145) q[2];
sx q[2];
rz(-2.243724) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5038824) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(-1.7590982) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9837708) q[3];
sx q[3];
rz(-1.2994248) q[3];
sx q[3];
rz(1.8872758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23665145) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(0.87990749) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(2.1039466) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-1.0132382) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(1.0838881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5686544) q[0];
sx q[0];
rz(-1.7540115) q[0];
sx q[0];
rz(-1.3230447) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8439699) q[2];
sx q[2];
rz(-1.8130842) q[2];
sx q[2];
rz(0.66852409) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0005972) q[1];
sx q[1];
rz(-2.6037569) q[1];
sx q[1];
rz(1.7829893) q[1];
x q[2];
rz(2.6277222) q[3];
sx q[3];
rz(-1.5268945) q[3];
sx q[3];
rz(-1.9572789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(-0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-0.50306815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(-2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(2.24618) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0953513) q[0];
sx q[0];
rz(-1.6086676) q[0];
sx q[0];
rz(-0.34237679) q[0];
x q[1];
rz(0.32161153) q[2];
sx q[2];
rz(-2.4827637) q[2];
sx q[2];
rz(-2.4668601) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4737282) q[1];
sx q[1];
rz(-1.6147991) q[1];
sx q[1];
rz(1.7394789) q[1];
rz(-pi) q[2];
rz(-1.6020681) q[3];
sx q[3];
rz(-3.0668689) q[3];
sx q[3];
rz(0.63262313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0078997) q[2];
sx q[2];
rz(-2.7329625) q[2];
sx q[2];
rz(2.3383979) q[2];
rz(1.9512272) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068709277) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(-0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(-1.8849467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6709267) q[0];
sx q[0];
rz(-1.6111272) q[0];
sx q[0];
rz(-1.3192024) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0357382) q[2];
sx q[2];
rz(-1.3824116) q[2];
sx q[2];
rz(2.9525194) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.41374846) q[1];
sx q[1];
rz(-1.2894948) q[1];
sx q[1];
rz(-1.6378228) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4323746) q[3];
sx q[3];
rz(-1.3243444) q[3];
sx q[3];
rz(-0.49530464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98823035) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.777565) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(-1.6114657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664739) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(-2.8334154) q[0];
rz(-3.0691052) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(-0.38696188) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.455084) q[0];
sx q[0];
rz(-2.5065656) q[0];
sx q[0];
rz(-1.1839266) q[0];
rz(-pi) q[1];
rz(-2.1808124) q[2];
sx q[2];
rz(-0.79723251) q[2];
sx q[2];
rz(2.0955992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19883979) q[1];
sx q[1];
rz(-1.6457335) q[1];
sx q[1];
rz(-0.67237512) q[1];
x q[2];
rz(-1.3513226) q[3];
sx q[3];
rz(-2.05728) q[3];
sx q[3];
rz(-2.6284077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5179634) q[2];
sx q[2];
rz(-1.7771746) q[2];
sx q[2];
rz(-0.44000885) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(-2.0429042) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7002174) q[0];
sx q[0];
rz(-2.2451631) q[0];
sx q[0];
rz(0.81654878) q[0];
rz(-pi) q[1];
rz(-2.6463037) q[2];
sx q[2];
rz(-0.9501179) q[2];
sx q[2];
rz(0.6492614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2652492) q[1];
sx q[1];
rz(-2.0298404) q[1];
sx q[1];
rz(2.6190119) q[1];
rz(-pi) q[2];
rz(-2.6155869) q[3];
sx q[3];
rz(-2.3283539) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.07842841) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(-2.4196529) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-0.75073457) q[3];
sx q[3];
rz(-0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14770517) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(-1.7396897) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83090529) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(-3.0124245) q[0];
rz(-pi) q[1];
rz(-0.14342587) q[2];
sx q[2];
rz(-2.9529245) q[2];
sx q[2];
rz(-0.27486899) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5857196) q[1];
sx q[1];
rz(-1.266529) q[1];
sx q[1];
rz(2.2833706) q[1];
rz(-0.63125061) q[3];
sx q[3];
rz(-1.4025941) q[3];
sx q[3];
rz(1.9025308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6587276) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.520291) q[2];
rz(-2.5907717) q[3];
sx q[3];
rz(-0.80533022) q[3];
sx q[3];
rz(2.4441392) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(2.2254754) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-0.73268391) q[2];
sx q[2];
rz(-0.97326836) q[2];
sx q[2];
rz(1.2748981) q[2];
rz(2.0879073) q[3];
sx q[3];
rz(-1.5862982) q[3];
sx q[3];
rz(1.3261212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
