OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(3.1363857) q[0];
rz(1.7967254) q[1];
sx q[1];
rz(4.2978573) q[1];
sx q[1];
rz(8.2351091) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1253163) q[0];
sx q[0];
rz(-1.9128886) q[0];
sx q[0];
rz(-2.5972511) q[0];
x q[1];
rz(0.18007937) q[2];
sx q[2];
rz(-1.4305978) q[2];
sx q[2];
rz(2.7498498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3559239) q[1];
sx q[1];
rz(-2.7969116) q[1];
sx q[1];
rz(2.0211401) q[1];
rz(-1.0359997) q[3];
sx q[3];
rz(-0.70123226) q[3];
sx q[3];
rz(2.0555156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4253915) q[2];
sx q[2];
rz(-0.74906817) q[2];
sx q[2];
rz(0.5973967) q[2];
rz(-1.776009) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(1.2805773) q[3];
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
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7146724) q[0];
sx q[0];
rz(-2.5898114) q[0];
sx q[0];
rz(0.33357099) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(3.0283668) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9202168) q[0];
sx q[0];
rz(-2.7144103) q[0];
sx q[0];
rz(-1.4772052) q[0];
rz(1.5510773) q[2];
sx q[2];
rz(-1.644051) q[2];
sx q[2];
rz(0.8578701) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.37564056) q[1];
sx q[1];
rz(-1.1728219) q[1];
sx q[1];
rz(-2.6938733) q[1];
rz(-1.232997) q[3];
sx q[3];
rz(-2.5099953) q[3];
sx q[3];
rz(0.88463569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-0.48015067) q[2];
sx q[2];
rz(2.9193027) q[2];
rz(-0.22953454) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(3.0139794) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022284431) q[0];
sx q[0];
rz(-0.89953178) q[0];
sx q[0];
rz(-0.22247252) q[0];
rz(-pi) q[1];
rz(0.87330841) q[2];
sx q[2];
rz(-1.0190939) q[2];
sx q[2];
rz(-1.1683977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.322791) q[1];
sx q[1];
rz(-1.4885508) q[1];
sx q[1];
rz(-1.3311885) q[1];
rz(1.913194) q[3];
sx q[3];
rz(-2.4857593) q[3];
sx q[3];
rz(-2.0369903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7814653) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(-2.8313417) q[2];
rz(-2.7919853) q[3];
sx q[3];
rz(-2.4571556) q[3];
sx q[3];
rz(1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.67868245) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-2.983685) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(0.37240949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83852864) q[0];
sx q[0];
rz(-1.6922564) q[0];
sx q[0];
rz(-2.2785506) q[0];
rz(-pi) q[1];
rz(-0.33505586) q[2];
sx q[2];
rz(-2.1261566) q[2];
sx q[2];
rz(-0.43713883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.47632521) q[1];
sx q[1];
rz(-1.7204666) q[1];
sx q[1];
rz(1.1900351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4140671) q[3];
sx q[3];
rz(-1.993506) q[3];
sx q[3];
rz(-0.62733516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(-2.8692029) q[2];
rz(0.90879905) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(-2.6866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-2.5407476) q[0];
sx q[0];
rz(2.0102665) q[0];
rz(-0.5979901) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(0.67684832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9899321) q[0];
sx q[0];
rz(-2.6940072) q[0];
sx q[0];
rz(2.9582892) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9025181) q[2];
sx q[2];
rz(-0.74547807) q[2];
sx q[2];
rz(-2.4649232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10830282) q[1];
sx q[1];
rz(-1.3160333) q[1];
sx q[1];
rz(-1.5441896) q[1];
rz(-pi) q[2];
rz(-1.5951717) q[3];
sx q[3];
rz(-1.9642324) q[3];
sx q[3];
rz(0.76373053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1375492) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(-1.1452902) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1260219) q[0];
sx q[0];
rz(-2.4937622) q[0];
sx q[0];
rz(1.4591249) q[0];
rz(-2.3964264) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-0.93313342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6720848) q[0];
sx q[0];
rz(-0.55969319) q[0];
sx q[0];
rz(-0.99225386) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51322333) q[2];
sx q[2];
rz(-2.4729118) q[2];
sx q[2];
rz(0.19323397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3999403) q[1];
sx q[1];
rz(-1.5685023) q[1];
sx q[1];
rz(1.2029349) q[1];
rz(-0.91306367) q[3];
sx q[3];
rz(-1.4300031) q[3];
sx q[3];
rz(-1.6933683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-2.9170673) q[2];
sx q[2];
rz(-2.8549426) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(-2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0934802) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(2.1215227) q[1];
sx q[1];
rz(-1.4694829) q[1];
sx q[1];
rz(2.1405623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0953656) q[0];
sx q[0];
rz(-1.2903178) q[0];
sx q[0];
rz(0.038897184) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83058968) q[2];
sx q[2];
rz(-0.48362728) q[2];
sx q[2];
rz(-0.73053503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8965473) q[1];
sx q[1];
rz(-2.7198615) q[1];
sx q[1];
rz(-1.1827724) q[1];
rz(0.6506728) q[3];
sx q[3];
rz(-1.7400029) q[3];
sx q[3];
rz(0.31864877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.004185685) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(-3.0380847) q[2];
rz(2.5497656) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25933927) q[0];
sx q[0];
rz(-1.0162202) q[0];
sx q[0];
rz(-2.5296339) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-0.38988316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1471257) q[0];
sx q[0];
rz(-1.9305221) q[0];
sx q[0];
rz(-0.61867449) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6340918) q[2];
sx q[2];
rz(-1.2050873) q[2];
sx q[2];
rz(-2.4909004) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.08323616) q[1];
sx q[1];
rz(-2.5857158) q[1];
sx q[1];
rz(0.91169375) q[1];
rz(0.26384683) q[3];
sx q[3];
rz(-2.2328394) q[3];
sx q[3];
rz(2.2784233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3020246) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-2.4273382) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(-0.57202655) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-0.51819003) q[1];
sx q[1];
rz(-0.27639595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24647507) q[0];
sx q[0];
rz(-2.7013489) q[0];
sx q[0];
rz(-2.759139) q[0];
rz(0.70119621) q[2];
sx q[2];
rz(-0.29507911) q[2];
sx q[2];
rz(-1.8104749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.285187) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(1.3473947) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60864457) q[3];
sx q[3];
rz(-1.3666271) q[3];
sx q[3];
rz(-0.58231402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(0.12750553) q[2];
rz(0.036711983) q[3];
sx q[3];
rz(-0.82009411) q[3];
sx q[3];
rz(-1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(0.35274831) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(-1.030285) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47980967) q[0];
sx q[0];
rz(-0.74736887) q[0];
sx q[0];
rz(-2.8519147) q[0];
rz(-pi) q[1];
rz(-2.7416517) q[2];
sx q[2];
rz(-2.415014) q[2];
sx q[2];
rz(-2.4253997) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5926338) q[1];
sx q[1];
rz(-2.0083798) q[1];
sx q[1];
rz(-1.9220819) q[1];
x q[2];
rz(0.45566166) q[3];
sx q[3];
rz(-0.66484287) q[3];
sx q[3];
rz(-1.1841707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0742005) q[2];
sx q[2];
rz(-1.8204764) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(2.4009005) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-0.091879524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-3.0108903) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(0.66979349) q[1];
sx q[1];
rz(-0.61129807) q[1];
sx q[1];
rz(1.5652464) q[1];
rz(0.43351602) q[2];
sx q[2];
rz(-1.3939861) q[2];
sx q[2];
rz(1.3409333) q[2];
rz(-0.99707281) q[3];
sx q[3];
rz(-0.96521796) q[3];
sx q[3];
rz(-0.60346606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];