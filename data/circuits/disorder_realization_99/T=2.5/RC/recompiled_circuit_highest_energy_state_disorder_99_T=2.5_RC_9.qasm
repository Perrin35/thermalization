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
rz(-1.2054297) q[0];
sx q[0];
rz(-3.087145) q[0];
sx q[0];
rz(2.2757538) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(-0.058252637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0308273) q[0];
sx q[0];
rz(-1.7351307) q[0];
sx q[0];
rz(3.0940233) q[0];
rz(2.0664332) q[2];
sx q[2];
rz(-1.505739) q[2];
sx q[2];
rz(-0.046147271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4076947) q[1];
sx q[1];
rz(-2.0996033) q[1];
sx q[1];
rz(-0.66441345) q[1];
rz(-pi) q[2];
rz(-0.32367651) q[3];
sx q[3];
rz(-1.2403245) q[3];
sx q[3];
rz(-3.0728473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3129348) q[2];
sx q[2];
rz(-0.88064319) q[2];
sx q[2];
rz(-0.9642967) q[2];
rz(3.0621081) q[3];
sx q[3];
rz(-1.8389523) q[3];
sx q[3];
rz(-2.3704119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(3.1281779) q[0];
sx q[0];
rz(-0.34341136) q[0];
sx q[0];
rz(1.6012023) q[0];
rz(-0.84114289) q[1];
sx q[1];
rz(-1.7336188) q[1];
sx q[1];
rz(-2.2841563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4682789) q[0];
sx q[0];
rz(-2.7029917) q[0];
sx q[0];
rz(1.5579786) q[0];
x q[1];
rz(-0.47274752) q[2];
sx q[2];
rz(-0.98927151) q[2];
sx q[2];
rz(-0.66752226) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7405058) q[1];
sx q[1];
rz(-1.1336278) q[1];
sx q[1];
rz(2.565041) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8624865) q[3];
sx q[3];
rz(-1.4785023) q[3];
sx q[3];
rz(-1.5971605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88605827) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(1.2096171) q[2];
rz(1.7602734) q[3];
sx q[3];
rz(-1.8807024) q[3];
sx q[3];
rz(-0.59022006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7968314) q[0];
sx q[0];
rz(-1.826257) q[0];
sx q[0];
rz(-1.3321846) q[0];
rz(-1.6650797) q[1];
sx q[1];
rz(-1.8107199) q[1];
sx q[1];
rz(0.61980334) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7999511) q[0];
sx q[0];
rz(-1.5864306) q[0];
sx q[0];
rz(-0.077341103) q[0];
rz(-1.8263426) q[2];
sx q[2];
rz(-1.7339306) q[2];
sx q[2];
rz(-2.852885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9528954) q[1];
sx q[1];
rz(-2.8159375) q[1];
sx q[1];
rz(-0.22518994) q[1];
x q[2];
rz(-2.8768646) q[3];
sx q[3];
rz(-0.94784289) q[3];
sx q[3];
rz(2.5271551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1406113) q[2];
sx q[2];
rz(-0.32537127) q[2];
sx q[2];
rz(0.78835431) q[2];
rz(-1.2761448) q[3];
sx q[3];
rz(-1.2272464) q[3];
sx q[3];
rz(0.86281323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9901504) q[0];
sx q[0];
rz(-1.3862415) q[0];
sx q[0];
rz(2.0748806) q[0];
rz(1.7881296) q[1];
sx q[1];
rz(-0.86694327) q[1];
sx q[1];
rz(-1.1563168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0321372) q[0];
sx q[0];
rz(-2.9988213) q[0];
sx q[0];
rz(-1.4383063) q[0];
rz(-pi) q[1];
rz(-2.1422191) q[2];
sx q[2];
rz(-3.093155) q[2];
sx q[2];
rz(0.20666176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9741648) q[1];
sx q[1];
rz(-1.8596223) q[1];
sx q[1];
rz(-1.1812737) q[1];
rz(-pi) q[2];
rz(-2.4772024) q[3];
sx q[3];
rz(-1.5186465) q[3];
sx q[3];
rz(-0.96168226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68916965) q[2];
sx q[2];
rz(-2.4781879) q[2];
sx q[2];
rz(-0.058825292) q[2];
rz(0.63932738) q[3];
sx q[3];
rz(-1.2213734) q[3];
sx q[3];
rz(-2.2842469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0462129) q[0];
sx q[0];
rz(-1.4154499) q[0];
sx q[0];
rz(0.010490622) q[0];
rz(1.5785716) q[1];
sx q[1];
rz(-2.450921) q[1];
sx q[1];
rz(-0.48935997) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6883514) q[0];
sx q[0];
rz(-0.64435378) q[0];
sx q[0];
rz(-2.8768566) q[0];
rz(-pi) q[1];
rz(-0.8801062) q[2];
sx q[2];
rz(-1.890213) q[2];
sx q[2];
rz(-1.3605786) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6076679) q[1];
sx q[1];
rz(-1.9235538) q[1];
sx q[1];
rz(-0.69490786) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33133026) q[3];
sx q[3];
rz(-3.1090571) q[3];
sx q[3];
rz(1.0573385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2995149) q[2];
sx q[2];
rz(-1.0966417) q[2];
sx q[2];
rz(-0.00046029885) q[2];
rz(0.67356235) q[3];
sx q[3];
rz(-0.898415) q[3];
sx q[3];
rz(-2.4323997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.8507268) q[0];
sx q[0];
rz(-1.8274266) q[0];
sx q[0];
rz(-1.3036183) q[0];
rz(0.29779008) q[1];
sx q[1];
rz(-0.89556634) q[1];
sx q[1];
rz(-0.29388014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0336908) q[0];
sx q[0];
rz(-1.0009888) q[0];
sx q[0];
rz(-0.66500647) q[0];
rz(-pi) q[1];
rz(-0.89048902) q[2];
sx q[2];
rz(-1.2273623) q[2];
sx q[2];
rz(0.78989093) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9781295) q[1];
sx q[1];
rz(-0.67369622) q[1];
sx q[1];
rz(-2.3232295) q[1];
rz(-pi) q[2];
rz(1.3354662) q[3];
sx q[3];
rz(-1.8782065) q[3];
sx q[3];
rz(-0.056005489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6728354) q[2];
sx q[2];
rz(-1.7054649) q[2];
sx q[2];
rz(-2.920816) q[2];
rz(1.8686434) q[3];
sx q[3];
rz(-1.3947398) q[3];
sx q[3];
rz(1.9934191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4918168) q[0];
sx q[0];
rz(-2.5204372) q[0];
sx q[0];
rz(-2.5671) q[0];
rz(-2.5993787) q[1];
sx q[1];
rz(-1.503399) q[1];
sx q[1];
rz(2.7508459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0076535) q[0];
sx q[0];
rz(-1.4263194) q[0];
sx q[0];
rz(2.5424315) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1152045) q[2];
sx q[2];
rz(-1.5261391) q[2];
sx q[2];
rz(0.80151973) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4225715) q[1];
sx q[1];
rz(-0.41677058) q[1];
sx q[1];
rz(-0.013643459) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7758632) q[3];
sx q[3];
rz(-1.5363294) q[3];
sx q[3];
rz(1.7548387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6812402) q[2];
sx q[2];
rz(-1.9443024) q[2];
sx q[2];
rz(-0.61275068) q[2];
rz(-2.1593306) q[3];
sx q[3];
rz(-1.3145612) q[3];
sx q[3];
rz(2.6764892) q[3];
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
rz(-pi/2) q[3];
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
rz(2.9923582) q[0];
sx q[0];
rz(-1.5748011) q[0];
sx q[0];
rz(-3.0862578) q[0];
rz(1.5845567) q[1];
sx q[1];
rz(-1.0880071) q[1];
sx q[1];
rz(0.43992821) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2038368) q[0];
sx q[0];
rz(-1.8169662) q[0];
sx q[0];
rz(0.47292821) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49100809) q[2];
sx q[2];
rz(-1.3898808) q[2];
sx q[2];
rz(2.343585) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7288558) q[1];
sx q[1];
rz(-1.2555033) q[1];
sx q[1];
rz(0.58087491) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8596419) q[3];
sx q[3];
rz(-2.1190756) q[3];
sx q[3];
rz(-1.5632526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.7679193) q[2];
sx q[2];
rz(-1.7934711) q[2];
sx q[2];
rz(0.31201735) q[2];
rz(-2.2219374) q[3];
sx q[3];
rz(-1.3684401) q[3];
sx q[3];
rz(-2.9254204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8298518) q[0];
sx q[0];
rz(-0.98167247) q[0];
sx q[0];
rz(-0.028129015) q[0];
rz(-1.6955388) q[1];
sx q[1];
rz(-1.1808993) q[1];
sx q[1];
rz(-0.45949724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579679) q[0];
sx q[0];
rz(-3.0959385) q[0];
sx q[0];
rz(-2.245957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3192603) q[2];
sx q[2];
rz(-1.8206545) q[2];
sx q[2];
rz(0.21669924) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5002437) q[1];
sx q[1];
rz(-0.8952221) q[1];
sx q[1];
rz(-0.18570034) q[1];
rz(-pi) q[2];
rz(-0.13925456) q[3];
sx q[3];
rz(-1.1920658) q[3];
sx q[3];
rz(-1.9851353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0405937) q[2];
sx q[2];
rz(-1.7969635) q[2];
sx q[2];
rz(-0.25516587) q[2];
rz(0.98662871) q[3];
sx q[3];
rz(-2.7268703) q[3];
sx q[3];
rz(1.423098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3808909) q[0];
sx q[0];
rz(-2.07169) q[0];
sx q[0];
rz(1.799452) q[0];
rz(0.0019207151) q[1];
sx q[1];
rz(-2.0989959) q[1];
sx q[1];
rz(-1.049918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494019) q[0];
sx q[0];
rz(-1.7190542) q[0];
sx q[0];
rz(-1.3561983) q[0];
x q[1];
rz(1.5776063) q[2];
sx q[2];
rz(-0.17596888) q[2];
sx q[2];
rz(-0.74534742) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.109595) q[1];
sx q[1];
rz(-1.9383311) q[1];
sx q[1];
rz(1.8046741) q[1];
rz(-pi) q[2];
rz(0.2157507) q[3];
sx q[3];
rz(-0.85236824) q[3];
sx q[3];
rz(-1.6529447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2125825) q[2];
sx q[2];
rz(-1.2139823) q[2];
sx q[2];
rz(-0.49599656) q[2];
rz(1.6811194) q[3];
sx q[3];
rz(-1.7998327) q[3];
sx q[3];
rz(-1.9546485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.7620508) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(-2.5024391) q[1];
sx q[1];
rz(-1.879138) q[1];
sx q[1];
rz(-2.383147) q[1];
rz(-0.46001804) q[2];
sx q[2];
rz(-0.5207396) q[2];
sx q[2];
rz(0.034551721) q[2];
rz(-2.4927327) q[3];
sx q[3];
rz(-0.93001233) q[3];
sx q[3];
rz(-0.42540023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
