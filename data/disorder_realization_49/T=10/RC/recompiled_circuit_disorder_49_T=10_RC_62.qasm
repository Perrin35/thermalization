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
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(1.1896689) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1253163) q[0];
sx q[0];
rz(-1.2287041) q[0];
sx q[0];
rz(-2.5972511) q[0];
x q[1];
rz(2.9615133) q[2];
sx q[2];
rz(-1.7109949) q[2];
sx q[2];
rz(2.7498498) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7835044) q[1];
sx q[1];
rz(-1.4231829) q[1];
sx q[1];
rz(-1.258177) q[1];
x q[2];
rz(0.40640229) q[3];
sx q[3];
rz(-0.98234017) q[3];
sx q[3];
rz(1.7455846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71620119) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(2.544196) q[2];
rz(1.776009) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(-1.2805773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7146724) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(-2.0479653) q[1];
sx q[1];
rz(-2.239614) q[1];
sx q[1];
rz(-3.0283668) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118606) q[0];
sx q[0];
rz(-1.1456053) q[0];
sx q[0];
rz(-3.0990764) q[0];
rz(0.26248652) q[2];
sx q[2];
rz(-0.075857698) q[2];
sx q[2];
rz(2.5469317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51616878) q[1];
sx q[1];
rz(-0.58991573) q[1];
sx q[1];
rz(-2.3708458) q[1];
x q[2];
rz(-2.9037335) q[3];
sx q[3];
rz(-0.97994643) q[3];
sx q[3];
rz(-1.8464586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-2.9193027) q[2];
rz(-0.22953454) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-3.0632609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753321) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(-0.88071841) q[0];
rz(2.0942028) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(-3.0139794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1193082) q[0];
sx q[0];
rz(-2.2420609) q[0];
sx q[0];
rz(0.22247252) q[0];
x q[1];
rz(-2.2682842) q[2];
sx q[2];
rz(-2.1224988) q[2];
sx q[2];
rz(-1.973195) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4096654) q[1];
sx q[1];
rz(-1.3320142) q[1];
sx q[1];
rz(-3.0569397) q[1];
rz(-2.8887799) q[3];
sx q[3];
rz(-0.95889927) q[3];
sx q[3];
rz(-1.6143527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7814653) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(2.8313417) q[2];
rz(-2.7919853) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(1.7278956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629102) q[0];
sx q[0];
rz(-1.8197729) q[0];
sx q[0];
rz(-0.15790766) q[0];
rz(0.36610106) q[1];
sx q[1];
rz(-2.3524275) q[1];
sx q[1];
rz(-2.7691832) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.303064) q[0];
sx q[0];
rz(-1.6922564) q[0];
sx q[0];
rz(2.2785506) q[0];
rz(0.33505586) q[2];
sx q[2];
rz(-2.1261566) q[2];
sx q[2];
rz(0.43713883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73788961) q[1];
sx q[1];
rz(-0.40778128) q[1];
sx q[1];
rz(-1.9562734) q[1];
rz(2.4140671) q[3];
sx q[3];
rz(-1.993506) q[3];
sx q[3];
rz(-2.5142575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41775122) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(-0.27238971) q[2];
rz(-0.90879905) q[3];
sx q[3];
rz(-0.89509982) q[3];
sx q[3];
rz(-0.4549543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(-2.0102665) q[0];
rz(0.5979901) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(2.4647443) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9899321) q[0];
sx q[0];
rz(-2.6940072) q[0];
sx q[0];
rz(-0.18330343) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73111515) q[2];
sx q[2];
rz(-1.4094681) q[2];
sx q[2];
rz(-0.71691712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9280793) q[1];
sx q[1];
rz(-2.885474) q[1];
sx q[1];
rz(3.03979) q[1];
x q[2];
rz(0.39354126) q[3];
sx q[3];
rz(-1.5482836) q[3];
sx q[3];
rz(-2.3251806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0040434917) q[2];
sx q[2];
rz(-1.9876336) q[2];
sx q[2];
rz(1.9963025) q[2];
rz(-2.9925313) q[3];
sx q[3];
rz(-2.5429433) q[3];
sx q[3];
rz(-1.0173652) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1260219) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(-1.6824678) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-0.93313342) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81290302) q[0];
sx q[0];
rz(-1.110154) q[0];
sx q[0];
rz(-0.33005379) q[0];
rz(-0.51322333) q[2];
sx q[2];
rz(-2.4729118) q[2];
sx q[2];
rz(0.19323397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9698525) q[1];
sx q[1];
rz(-1.9386567) q[1];
sx q[1];
rz(0.0024585558) q[1];
rz(2.9643781) q[3];
sx q[3];
rz(-0.92068499) q[3];
sx q[3];
rz(-0.23055102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(2.8549426) q[2];
rz(2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(0.4683032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0934802) q[0];
sx q[0];
rz(-1.1756287) q[0];
sx q[0];
rz(1.4666784) q[0];
rz(2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(-2.1405623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769343) q[0];
sx q[0];
rz(-1.6081728) q[0];
sx q[0];
rz(-1.2901165) q[0];
x q[1];
rz(1.940735) q[2];
sx q[2];
rz(-1.2518034) q[2];
sx q[2];
rz(1.621304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66623464) q[1];
sx q[1];
rz(-1.959414) q[1];
sx q[1];
rz(-0.16814853) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4909199) q[3];
sx q[3];
rz(-1.4015897) q[3];
sx q[3];
rz(-2.8229439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.137407) q[2];
sx q[2];
rz(-1.2265394) q[2];
sx q[2];
rz(0.10350791) q[2];
rz(-0.59182709) q[3];
sx q[3];
rz(-1.8154362) q[3];
sx q[3];
rz(-2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25933927) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(-0.6119588) q[0];
rz(1.3423963) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(-0.38988316) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99446699) q[0];
sx q[0];
rz(-1.2110706) q[0];
sx q[0];
rz(-0.61867449) q[0];
x q[1];
rz(-1.9837603) q[2];
sx q[2];
rz(-1.0997084) q[2];
sx q[2];
rz(-2.4177891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.08323616) q[1];
sx q[1];
rz(-0.55587686) q[1];
sx q[1];
rz(0.91169375) q[1];
rz(1.893703) q[3];
sx q[3];
rz(-0.70526037) q[3];
sx q[3];
rz(-0.44912072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8395681) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(2.4273382) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-2.0102746) q[3];
sx q[3];
rz(-2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488688) q[0];
sx q[0];
rz(-0.435193) q[0];
sx q[0];
rz(-0.77734787) q[0];
rz(-0.84689394) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(0.27639595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6733288) q[0];
sx q[0];
rz(-1.7305166) q[0];
sx q[0];
rz(2.7295652) q[0];
x q[1];
rz(1.3771636) q[2];
sx q[2];
rz(-1.7948705) q[2];
sx q[2];
rz(-0.60806882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8917577) q[1];
sx q[1];
rz(-1.7913622) q[1];
sx q[1];
rz(2.9794429) q[1];
rz(-1.3235839) q[3];
sx q[3];
rz(-0.97655481) q[3];
sx q[3];
rz(-2.0127399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85236621) q[2];
sx q[2];
rz(-1.5091395) q[2];
sx q[2];
rz(-0.12750553) q[2];
rz(0.036711983) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(1.2344853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4002832) q[0];
sx q[0];
rz(-0.4929587) q[0];
sx q[0];
rz(-2.7888443) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-2.1513217) q[1];
sx q[1];
rz(-2.1113077) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47980967) q[0];
sx q[0];
rz(-0.74736887) q[0];
sx q[0];
rz(-0.28967793) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2376386) q[2];
sx q[2];
rz(-0.91234708) q[2];
sx q[2];
rz(1.2308987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5489588) q[1];
sx q[1];
rz(-1.1332129) q[1];
sx q[1];
rz(-1.2195107) q[1];
rz(1.9029721) q[3];
sx q[3];
rz(-2.1579451) q[3];
sx q[3];
rz(1.4004933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0673922) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(-2.8004048) q[2];
rz(-0.7406922) q[3];
sx q[3];
rz(-0.98978981) q[3];
sx q[3];
rz(-3.0497131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
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
rz(-0.40217051) q[2];
sx q[2];
rz(-2.675534) q[2];
sx q[2];
rz(0.13327394) q[2];
rz(0.66486799) q[3];
sx q[3];
rz(-2.3330199) q[3];
sx q[3];
rz(-1.4521269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
