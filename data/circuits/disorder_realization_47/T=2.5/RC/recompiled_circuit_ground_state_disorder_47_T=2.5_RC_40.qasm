OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.2831777) q[0];
sx q[0];
rz(-1.5812961) q[0];
sx q[0];
rz(-0.66891447) q[0];
rz(0.34612292) q[1];
sx q[1];
rz(-2.3200413) q[1];
sx q[1];
rz(-0.93924826) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72059435) q[0];
sx q[0];
rz(-1.2828553) q[0];
sx q[0];
rz(1.9701411) q[0];
x q[1];
rz(-1.16667) q[2];
sx q[2];
rz(-1.6138645) q[2];
sx q[2];
rz(1.7801746) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5377944) q[1];
sx q[1];
rz(-1.1883573) q[1];
sx q[1];
rz(-2.8549744) q[1];
x q[2];
rz(0.40832728) q[3];
sx q[3];
rz(-0.75330594) q[3];
sx q[3];
rz(1.7024794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5657438) q[2];
sx q[2];
rz(-1.0172903) q[2];
sx q[2];
rz(-3.1080833) q[2];
rz(-1.2014028) q[3];
sx q[3];
rz(-1.7824495) q[3];
sx q[3];
rz(-1.0777773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031438436) q[0];
sx q[0];
rz(-2.8903676) q[0];
sx q[0];
rz(-0.95397368) q[0];
rz(1.6961478) q[1];
sx q[1];
rz(-1.0547538) q[1];
sx q[1];
rz(1.9704069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508043) q[0];
sx q[0];
rz(-0.91676869) q[0];
sx q[0];
rz(-0.79933856) q[0];
rz(1.2759802) q[2];
sx q[2];
rz(-2.1871242) q[2];
sx q[2];
rz(0.64586879) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3747159) q[1];
sx q[1];
rz(-1.3424113) q[1];
sx q[1];
rz(-1.5414052) q[1];
rz(2.8071142) q[3];
sx q[3];
rz(-1.1568767) q[3];
sx q[3];
rz(-2.778307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.235432) q[2];
sx q[2];
rz(-1.4305328) q[2];
sx q[2];
rz(1.1492427) q[2];
rz(0.033128459) q[3];
sx q[3];
rz(-1.6413611) q[3];
sx q[3];
rz(-2.5455425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0843435) q[0];
sx q[0];
rz(-2.3488022) q[0];
sx q[0];
rz(-1.0302011) q[0];
rz(1.239981) q[1];
sx q[1];
rz(-1.3571309) q[1];
sx q[1];
rz(0.99303594) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6246374) q[0];
sx q[0];
rz(-2.0258396) q[0];
sx q[0];
rz(1.5173453) q[0];
rz(0.15056653) q[2];
sx q[2];
rz(-2.4685301) q[2];
sx q[2];
rz(1.7889495) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60336411) q[1];
sx q[1];
rz(-2.1346655) q[1];
sx q[1];
rz(-2.4634874) q[1];
x q[2];
rz(2.7168324) q[3];
sx q[3];
rz(-0.83213193) q[3];
sx q[3];
rz(-2.5855541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21732907) q[2];
sx q[2];
rz(-2.0254878) q[2];
sx q[2];
rz(-0.35476157) q[2];
rz(-0.99700704) q[3];
sx q[3];
rz(-1.6544147) q[3];
sx q[3];
rz(-2.2009489) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16911258) q[0];
sx q[0];
rz(-1.7153772) q[0];
sx q[0];
rz(1.1068363) q[0];
rz(0.86889443) q[1];
sx q[1];
rz(-0.51742253) q[1];
sx q[1];
rz(2.4443764) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91928673) q[0];
sx q[0];
rz(-2.5450052) q[0];
sx q[0];
rz(-0.35937341) q[0];
rz(1.4585232) q[2];
sx q[2];
rz(-0.85190433) q[2];
sx q[2];
rz(-2.7996705) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.47095151) q[1];
sx q[1];
rz(-1.1837675) q[1];
sx q[1];
rz(2.7426038) q[1];
rz(-3.117048) q[3];
sx q[3];
rz(-1.5460137) q[3];
sx q[3];
rz(-1.9056232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2775468) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(-2.8475658) q[2];
rz(2.0781519) q[3];
sx q[3];
rz(-1.682351) q[3];
sx q[3];
rz(-2.3991876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.71642891) q[0];
sx q[0];
rz(-2.9603781) q[0];
sx q[0];
rz(2.678405) q[0];
rz(2.7376392) q[1];
sx q[1];
rz(-1.633176) q[1];
sx q[1];
rz(-0.24872669) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8549268) q[0];
sx q[0];
rz(-1.5724941) q[0];
sx q[0];
rz(1.7683074) q[0];
x q[1];
rz(-2.6004148) q[2];
sx q[2];
rz(-1.3480617) q[2];
sx q[2];
rz(2.7991852) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7686004) q[1];
sx q[1];
rz(-1.2443719) q[1];
sx q[1];
rz(-1.1821163) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90051265) q[3];
sx q[3];
rz(-0.95275646) q[3];
sx q[3];
rz(1.6992118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2948461) q[2];
sx q[2];
rz(-1.8241901) q[2];
sx q[2];
rz(1.741629) q[2];
rz(-1.2830265) q[3];
sx q[3];
rz(-1.3834407) q[3];
sx q[3];
rz(0.2259026) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95752174) q[0];
sx q[0];
rz(-0.83860832) q[0];
sx q[0];
rz(0.34570178) q[0];
rz(-3.0905981) q[1];
sx q[1];
rz(-1.2268927) q[1];
sx q[1];
rz(-0.7720224) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8851903) q[0];
sx q[0];
rz(-1.2202383) q[0];
sx q[0];
rz(-1.3023443) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4441024) q[2];
sx q[2];
rz(-0.65595731) q[2];
sx q[2];
rz(-2.3672589) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1353957) q[1];
sx q[1];
rz(-2.4675698) q[1];
sx q[1];
rz(-2.1042906) q[1];
rz(-pi) q[2];
rz(2.8528522) q[3];
sx q[3];
rz(-1.0613228) q[3];
sx q[3];
rz(-2.4276707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8649851) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(0.15743206) q[2];
rz(0.43846798) q[3];
sx q[3];
rz(-2.4003568) q[3];
sx q[3];
rz(-0.95611519) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601783) q[0];
sx q[0];
rz(-1.0670476) q[0];
sx q[0];
rz(2.881158) q[0];
rz(-2.5632437) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(-0.15301212) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3653202) q[0];
sx q[0];
rz(-2.0663102) q[0];
sx q[0];
rz(-0.97712626) q[0];
rz(-pi) q[1];
rz(0.084916755) q[2];
sx q[2];
rz(-1.6667637) q[2];
sx q[2];
rz(-1.7964448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2305829) q[1];
sx q[1];
rz(-2.1248159) q[1];
sx q[1];
rz(-1.0826712) q[1];
rz(-pi) q[2];
rz(0.32185359) q[3];
sx q[3];
rz(-1.1328837) q[3];
sx q[3];
rz(1.9066325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1795307) q[2];
sx q[2];
rz(-0.099055812) q[2];
sx q[2];
rz(-2.8774101) q[2];
rz(-1.8386486) q[3];
sx q[3];
rz(-1.7641726) q[3];
sx q[3];
rz(-1.1072268) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95689479) q[0];
sx q[0];
rz(-2.1883924) q[0];
sx q[0];
rz(-1.4554998) q[0];
rz(0.9616583) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(2.2419194) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5035928) q[0];
sx q[0];
rz(-1.3595306) q[0];
sx q[0];
rz(2.6799623) q[0];
rz(-pi) q[1];
rz(1.3191965) q[2];
sx q[2];
rz(-2.4008022) q[2];
sx q[2];
rz(2.0513926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3916495) q[1];
sx q[1];
rz(-2.1575621) q[1];
sx q[1];
rz(3.090078) q[1];
x q[2];
rz(-0.03735383) q[3];
sx q[3];
rz(-1.3539697) q[3];
sx q[3];
rz(1.9268394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4453033) q[2];
sx q[2];
rz(-2.5407007) q[2];
sx q[2];
rz(0.75330934) q[2];
rz(0.2937915) q[3];
sx q[3];
rz(-1.1283504) q[3];
sx q[3];
rz(2.0297348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927354) q[0];
sx q[0];
rz(-0.98965544) q[0];
sx q[0];
rz(2.2897172) q[0];
rz(1.7333671) q[1];
sx q[1];
rz(-1.6586761) q[1];
sx q[1];
rz(-0.3826938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9155884) q[0];
sx q[0];
rz(-1.3289641) q[0];
sx q[0];
rz(-0.42892021) q[0];
rz(0.76994728) q[2];
sx q[2];
rz(-1.4048409) q[2];
sx q[2];
rz(2.6024352) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.574461) q[1];
sx q[1];
rz(-1.8715011) q[1];
sx q[1];
rz(1.6733132) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48608853) q[3];
sx q[3];
rz(-1.4975784) q[3];
sx q[3];
rz(-1.40772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7644299) q[2];
sx q[2];
rz(-0.96401507) q[2];
sx q[2];
rz(0.60297472) q[2];
rz(-1.0682586) q[3];
sx q[3];
rz(-1.4385185) q[3];
sx q[3];
rz(2.300613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0076440796) q[0];
sx q[0];
rz(-1.6772062) q[0];
sx q[0];
rz(2.7879047) q[0];
rz(2.0091281) q[1];
sx q[1];
rz(-2.1734889) q[1];
sx q[1];
rz(0.38945928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9326453) q[0];
sx q[0];
rz(-2.0660127) q[0];
sx q[0];
rz(1.3071278) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.816941) q[2];
sx q[2];
rz(-1.6707641) q[2];
sx q[2];
rz(2.113518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88122565) q[1];
sx q[1];
rz(-2.0282451) q[1];
sx q[1];
rz(-1.264099) q[1];
rz(0.090062304) q[3];
sx q[3];
rz(-2.6115993) q[3];
sx q[3];
rz(1.8217721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13122095) q[2];
sx q[2];
rz(-1.6348569) q[2];
sx q[2];
rz(-1.1466675) q[2];
rz(-1.5366588) q[3];
sx q[3];
rz(-0.48303548) q[3];
sx q[3];
rz(-0.35227942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1491886) q[0];
sx q[0];
rz(-2.438899) q[0];
sx q[0];
rz(-2.6371523) q[0];
rz(3.0814677) q[1];
sx q[1];
rz(-0.45777121) q[1];
sx q[1];
rz(-0.4578185) q[1];
rz(1.4627152) q[2];
sx q[2];
rz(-1.6440132) q[2];
sx q[2];
rz(-3.0655412) q[2];
rz(-1.6144013) q[3];
sx q[3];
rz(-2.4666967) q[3];
sx q[3];
rz(-1.1793292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
