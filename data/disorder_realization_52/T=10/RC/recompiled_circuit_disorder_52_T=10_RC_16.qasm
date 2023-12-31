OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(-2.2572416) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9031154) q[0];
sx q[0];
rz(-1.7503386) q[0];
sx q[0];
rz(-1.9360916) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1332364) q[2];
sx q[2];
rz(-0.5651606) q[2];
sx q[2];
rz(-1.1430119) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4554162) q[1];
sx q[1];
rz(-1.505291) q[1];
sx q[1];
rz(-2.8951365) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2803335) q[3];
sx q[3];
rz(-1.2951295) q[3];
sx q[3];
rz(-1.1112801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-0.65594977) q[2];
rz(-1.2077228) q[3];
sx q[3];
rz(-1.175712) q[3];
sx q[3];
rz(-0.99457994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2475964) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(-0.43352747) q[0];
rz(-2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-3.1343592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1996961) q[0];
sx q[0];
rz(-1.4775839) q[0];
sx q[0];
rz(1.1774363) q[0];
rz(0.8231926) q[2];
sx q[2];
rz(-2.7232812) q[2];
sx q[2];
rz(0.38564607) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4870287) q[1];
sx q[1];
rz(-1.9332758) q[1];
sx q[1];
rz(-1.3351721) q[1];
rz(-pi) q[2];
rz(0.98874436) q[3];
sx q[3];
rz(-2.3855004) q[3];
sx q[3];
rz(2.3088629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-0.80792892) q[2];
sx q[2];
rz(2.4439404) q[2];
rz(-0.12156045) q[3];
sx q[3];
rz(-1.9024885) q[3];
sx q[3];
rz(-2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6737297) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(1.3695705) q[0];
rz(1.9000152) q[1];
sx q[1];
rz(-1.7280271) q[1];
sx q[1];
rz(2.8799768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152552) q[0];
sx q[0];
rz(-1.5887504) q[0];
sx q[0];
rz(3.1319588) q[0];
rz(-pi) q[1];
rz(-0.82654731) q[2];
sx q[2];
rz(-1.6625704) q[2];
sx q[2];
rz(1.2395791) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9237823) q[1];
sx q[1];
rz(-0.79499704) q[1];
sx q[1];
rz(1.4008821) q[1];
rz(-pi) q[2];
rz(-1.5393799) q[3];
sx q[3];
rz(-2.0835702) q[3];
sx q[3];
rz(0.98785066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6802784) q[2];
sx q[2];
rz(-1.5051944) q[2];
sx q[2];
rz(0.20351163) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(-2.8620499) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97776425) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(-1.571636) q[0];
rz(2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(-1.2483695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6682537) q[0];
sx q[0];
rz(-0.11675294) q[0];
sx q[0];
rz(2.1658685) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4291971) q[2];
sx q[2];
rz(-0.27370307) q[2];
sx q[2];
rz(1.4272387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9256546) q[1];
sx q[1];
rz(-0.92080322) q[1];
sx q[1];
rz(-1.710612) q[1];
rz(-pi) q[2];
rz(-0.72327153) q[3];
sx q[3];
rz(-1.7687106) q[3];
sx q[3];
rz(-0.42641446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1127597) q[2];
sx q[2];
rz(-1.3112105) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(1.2083496) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(-0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355302) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(2.3663882) q[0];
rz(-0.40183055) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(-0.90243375) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8309801) q[0];
sx q[0];
rz(-1.6164301) q[0];
sx q[0];
rz(2.5020585) q[0];
x q[1];
rz(-2.6685153) q[2];
sx q[2];
rz(-0.50191754) q[2];
sx q[2];
rz(2.9830473) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.24689281) q[1];
sx q[1];
rz(-1.0122074) q[1];
sx q[1];
rz(2.488399) q[1];
x q[2];
rz(-1.8875214) q[3];
sx q[3];
rz(-1.5734908) q[3];
sx q[3];
rz(0.56841422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9465785) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(-1.9449332) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8114132) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(-1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(2.9398289) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.266298) q[0];
sx q[0];
rz(-1.504997) q[0];
sx q[0];
rz(-1.6732897) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74886518) q[2];
sx q[2];
rz(-2.5267234) q[2];
sx q[2];
rz(-2.2392337) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4060496) q[1];
sx q[1];
rz(-2.4821401) q[1];
sx q[1];
rz(-0.72592782) q[1];
rz(-pi) q[2];
rz(-0.50815177) q[3];
sx q[3];
rz(-2.3414632) q[3];
sx q[3];
rz(-1.8964163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9266944) q[2];
sx q[2];
rz(-1.5025257) q[2];
sx q[2];
rz(0.26724896) q[2];
rz(2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(2.2657623) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-0.057549495) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(0.94271359) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1377624) q[0];
sx q[0];
rz(-0.98235213) q[0];
sx q[0];
rz(-1.0914735) q[0];
rz(-0.96037453) q[2];
sx q[2];
rz(-2.865961) q[2];
sx q[2];
rz(2.9396217) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.645694) q[1];
sx q[1];
rz(-1.8897448) q[1];
sx q[1];
rz(-1.1369399) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0122635) q[3];
sx q[3];
rz(-0.55674508) q[3];
sx q[3];
rz(-0.41551513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(-2.7589202) q[2];
rz(-1.0391957) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(-2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7169749) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(-0.27012816) q[0];
rz(-0.62942901) q[1];
sx q[1];
rz(-2.4286178) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2487508) q[0];
sx q[0];
rz(-0.98309702) q[0];
sx q[0];
rz(-2.6177004) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7173041) q[2];
sx q[2];
rz(-0.97587817) q[2];
sx q[2];
rz(-2.811424) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.64360196) q[1];
sx q[1];
rz(-2.0703348) q[1];
sx q[1];
rz(1.040578) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5408526) q[3];
sx q[3];
rz(-0.91120126) q[3];
sx q[3];
rz(-0.75920446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(-1.930687) q[2];
rz(2.8816913) q[3];
sx q[3];
rz(-2.5172958) q[3];
sx q[3];
rz(2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07638409) q[0];
sx q[0];
rz(-0.56088352) q[0];
sx q[0];
rz(2.912345) q[0];
rz(2.8385838) q[1];
sx q[1];
rz(-1.7508933) q[1];
sx q[1];
rz(1.4607666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47855908) q[0];
sx q[0];
rz(-1.2872818) q[0];
sx q[0];
rz(-1.589033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3912348) q[2];
sx q[2];
rz(-1.944724) q[2];
sx q[2];
rz(-1.1073081) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2158828) q[1];
sx q[1];
rz(-0.52535086) q[1];
sx q[1];
rz(0.044563091) q[1];
rz(-0.38735729) q[3];
sx q[3];
rz(-2.318317) q[3];
sx q[3];
rz(-0.92126095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(-0.6955859) q[2];
rz(0.43186489) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(0.7152043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5678976) q[0];
sx q[0];
rz(-1.2117813) q[0];
sx q[0];
rz(2.8046872) q[0];
rz(0.20740549) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(0.38063231) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363627) q[0];
sx q[0];
rz(-2.9281033) q[0];
sx q[0];
rz(1.2345033) q[0];
rz(-pi) q[1];
rz(0.62217525) q[2];
sx q[2];
rz(-2.2969349) q[2];
sx q[2];
rz(-0.13571339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.759728) q[1];
sx q[1];
rz(-1.7654395) q[1];
sx q[1];
rz(1.0641644) q[1];
rz(-0.57772824) q[3];
sx q[3];
rz(-1.2168222) q[3];
sx q[3];
rz(-0.83565328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(-1.1364737) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-0.80871964) q[3];
sx q[3];
rz(1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(2.1144755) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(-2.8935511) q[2];
sx q[2];
rz(-1.3477865) q[2];
sx q[2];
rz(-1.8852521) q[2];
rz(1.1626701) q[3];
sx q[3];
rz(-2.6682731) q[3];
sx q[3];
rz(-0.70262739) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
