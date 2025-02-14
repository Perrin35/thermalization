OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.858415) q[0];
sx q[0];
rz(-1.5602966) q[0];
sx q[0];
rz(-2.4726782) q[0];
rz(-2.7954697) q[1];
sx q[1];
rz(-0.82155138) q[1];
sx q[1];
rz(0.93924826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72059435) q[0];
sx q[0];
rz(-1.8587374) q[0];
sx q[0];
rz(1.1714515) q[0];
x q[1];
rz(-3.0947565) q[2];
sx q[2];
rz(-1.1670665) q[2];
sx q[2];
rz(-2.9138034) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60379825) q[1];
sx q[1];
rz(-1.9532353) q[1];
sx q[1];
rz(-2.8549744) q[1];
rz(0.71066831) q[3];
sx q[3];
rz(-1.2957199) q[3];
sx q[3];
rz(2.9675837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5657438) q[2];
sx q[2];
rz(-2.1243024) q[2];
sx q[2];
rz(0.033509342) q[2];
rz(-1.2014028) q[3];
sx q[3];
rz(-1.3591432) q[3];
sx q[3];
rz(-2.0638154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.031438436) q[0];
sx q[0];
rz(-0.25122508) q[0];
sx q[0];
rz(2.187619) q[0];
rz(1.6961478) q[1];
sx q[1];
rz(-1.0547538) q[1];
sx q[1];
rz(1.9704069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508043) q[0];
sx q[0];
rz(-0.91676869) q[0];
sx q[0];
rz(-0.79933856) q[0];
x q[1];
rz(0.38925217) q[2];
sx q[2];
rz(-2.4667423) q[2];
sx q[2];
rz(-2.9794326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5038472) q[1];
sx q[1];
rz(-0.23023573) q[1];
sx q[1];
rz(0.12576018) q[1];
x q[2];
rz(-0.33447845) q[3];
sx q[3];
rz(-1.984716) q[3];
sx q[3];
rz(2.778307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.235432) q[2];
sx q[2];
rz(-1.4305328) q[2];
sx q[2];
rz(-1.99235) q[2];
rz(3.1084642) q[3];
sx q[3];
rz(-1.6413611) q[3];
sx q[3];
rz(2.5455425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0572492) q[0];
sx q[0];
rz(-2.3488022) q[0];
sx q[0];
rz(-2.1113915) q[0];
rz(-1.239981) q[1];
sx q[1];
rz(-1.3571309) q[1];
sx q[1];
rz(2.1485567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030331197) q[0];
sx q[0];
rz(-1.6188038) q[0];
sx q[0];
rz(0.45560776) q[0];
rz(-pi) q[1];
rz(-2.4740691) q[2];
sx q[2];
rz(-1.6644396) q[2];
sx q[2];
rz(0.33622959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5807016) q[1];
sx q[1];
rz(-2.1296394) q[1];
sx q[1];
rz(2.2528095) q[1];
rz(-pi) q[2];
rz(-1.1458323) q[3];
sx q[3];
rz(-2.309955) q[3];
sx q[3];
rz(-1.1475565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.21732907) q[2];
sx q[2];
rz(-1.1161048) q[2];
sx q[2];
rz(0.35476157) q[2];
rz(0.99700704) q[3];
sx q[3];
rz(-1.6544147) q[3];
sx q[3];
rz(-0.94064373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724801) q[0];
sx q[0];
rz(-1.4262154) q[0];
sx q[0];
rz(-1.1068363) q[0];
rz(2.2726982) q[1];
sx q[1];
rz(-2.6241701) q[1];
sx q[1];
rz(-0.69721627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91928673) q[0];
sx q[0];
rz(-0.59658748) q[0];
sx q[0];
rz(-0.35937341) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6830695) q[2];
sx q[2];
rz(-2.2896883) q[2];
sx q[2];
rz(-2.7996705) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8839421) q[1];
sx q[1];
rz(-1.9387987) q[1];
sx q[1];
rz(-1.1543399) q[1];
x q[2];
rz(2.3512164) q[3];
sx q[3];
rz(-3.1067143) q[3];
sx q[3];
rz(0.45524516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8640459) q[2];
sx q[2];
rz(-0.20865455) q[2];
sx q[2];
rz(-0.29402688) q[2];
rz(-1.0634408) q[3];
sx q[3];
rz(-1.682351) q[3];
sx q[3];
rz(0.74240509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251637) q[0];
sx q[0];
rz(-0.18121457) q[0];
sx q[0];
rz(-2.678405) q[0];
rz(-2.7376392) q[1];
sx q[1];
rz(-1.633176) q[1];
sx q[1];
rz(0.24872669) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8578019) q[0];
sx q[0];
rz(-1.3732855) q[0];
sx q[0];
rz(-3.1398612) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8291446) q[2];
sx q[2];
rz(-2.097192) q[2];
sx q[2];
rz(-2.0451982) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.675093) q[1];
sx q[1];
rz(-0.50216253) q[1];
sx q[1];
rz(0.84169047) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4047923) q[3];
sx q[3];
rz(-1.0399138) q[3];
sx q[3];
rz(-2.8394222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8467466) q[2];
sx q[2];
rz(-1.8241901) q[2];
sx q[2];
rz(1.3999636) q[2];
rz(1.2830265) q[3];
sx q[3];
rz(-1.7581519) q[3];
sx q[3];
rz(-2.9156901) q[3];
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
rz(2.1840709) q[0];
sx q[0];
rz(-0.83860832) q[0];
sx q[0];
rz(-0.34570178) q[0];
rz(0.050994571) q[1];
sx q[1];
rz(-1.2268927) q[1];
sx q[1];
rz(2.3695703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2098198) q[0];
sx q[0];
rz(-2.7034524) q[0];
sx q[0];
rz(2.5140425) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2228681) q[2];
sx q[2];
rz(-1.6479392) q[2];
sx q[2];
rz(2.2445298) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.009338) q[1];
sx q[1];
rz(-1.2478095) q[1];
sx q[1];
rz(2.1732974) q[1];
rz(-pi) q[2];
rz(-1.0431248) q[3];
sx q[3];
rz(-1.3195795) q[3];
sx q[3];
rz(-0.71300292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2766075) q[2];
sx q[2];
rz(-1.425068) q[2];
sx q[2];
rz(-2.9841606) q[2];
rz(0.43846798) q[3];
sx q[3];
rz(-0.74123588) q[3];
sx q[3];
rz(0.95611519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081414374) q[0];
sx q[0];
rz(-1.0670476) q[0];
sx q[0];
rz(-2.881158) q[0];
rz(2.5632437) q[1];
sx q[1];
rz(-1.3713505) q[1];
sx q[1];
rz(-2.9885805) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4082913) q[0];
sx q[0];
rz(-2.387945) q[0];
sx q[0];
rz(2.338999) q[0];
x q[1];
rz(1.4744841) q[2];
sx q[2];
rz(-1.6553214) q[2];
sx q[2];
rz(-0.21749228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93219705) q[1];
sx q[1];
rz(-1.160566) q[1];
sx q[1];
rz(0.61102976) q[1];
rz(2.1649394) q[3];
sx q[3];
rz(-0.53721957) q[3];
sx q[3];
rz(-1.2402759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.962062) q[2];
sx q[2];
rz(-3.0425368) q[2];
sx q[2];
rz(0.2641826) q[2];
rz(-1.3029441) q[3];
sx q[3];
rz(-1.3774201) q[3];
sx q[3];
rz(-1.1072268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-2.1846979) q[0];
sx q[0];
rz(-2.1883924) q[0];
sx q[0];
rz(1.4554998) q[0];
rz(2.1799344) q[1];
sx q[1];
rz(-1.7627629) q[1];
sx q[1];
rz(0.89967322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5035928) q[0];
sx q[0];
rz(-1.7820621) q[0];
sx q[0];
rz(-2.6799623) q[0];
rz(-pi) q[1];
rz(-1.8223962) q[2];
sx q[2];
rz(-2.4008022) q[2];
sx q[2];
rz(2.0513926) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4845061) q[1];
sx q[1];
rz(-0.58875798) q[1];
sx q[1];
rz(-1.4935054) q[1];
x q[2];
rz(3.1042388) q[3];
sx q[3];
rz(-1.3539697) q[3];
sx q[3];
rz(-1.2147533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6962894) q[2];
sx q[2];
rz(-2.5407007) q[2];
sx q[2];
rz(0.75330934) q[2];
rz(-2.8478012) q[3];
sx q[3];
rz(-2.0132422) q[3];
sx q[3];
rz(1.1118579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488572) q[0];
sx q[0];
rz(-0.98965544) q[0];
sx q[0];
rz(-0.85187546) q[0];
rz(-1.4082255) q[1];
sx q[1];
rz(-1.4829166) q[1];
sx q[1];
rz(-2.7588989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9155884) q[0];
sx q[0];
rz(-1.3289641) q[0];
sx q[0];
rz(-2.7126724) q[0];
x q[1];
rz(0.23613249) q[2];
sx q[2];
rz(-2.3575767) q[2];
sx q[2];
rz(-1.2003984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.574461) q[1];
sx q[1];
rz(-1.2700915) q[1];
sx q[1];
rz(-1.4682795) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(1.7644299) q[2];
sx q[2];
rz(-2.1775776) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1339486) q[0];
sx q[0];
rz(-1.6772062) q[0];
sx q[0];
rz(-2.7879047) q[0];
rz(1.1324646) q[1];
sx q[1];
rz(-0.9681038) q[1];
sx q[1];
rz(0.38945928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4492466) q[0];
sx q[0];
rz(-0.55583891) q[0];
sx q[0];
rz(-2.6920094) q[0];
x q[1];
rz(-1.3246516) q[2];
sx q[2];
rz(-1.4708286) q[2];
sx q[2];
rz(-1.0280746) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88122565) q[1];
sx q[1];
rz(-1.1133476) q[1];
sx q[1];
rz(1.264099) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5181481) q[3];
sx q[3];
rz(-1.0431759) q[3];
sx q[3];
rz(1.2155346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13122095) q[2];
sx q[2];
rz(-1.6348569) q[2];
sx q[2];
rz(-1.1466675) q[2];
rz(-1.6049339) q[3];
sx q[3];
rz(-0.48303548) q[3];
sx q[3];
rz(-2.7893132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99240408) q[0];
sx q[0];
rz(-0.70269361) q[0];
sx q[0];
rz(0.50444034) q[0];
rz(3.0814677) q[1];
sx q[1];
rz(-0.45777121) q[1];
sx q[1];
rz(-0.4578185) q[1];
rz(3.0679476) q[2];
sx q[2];
rz(-1.4630058) q[2];
sx q[2];
rz(-1.502682) q[2];
rz(1.5271913) q[3];
sx q[3];
rz(-2.4666967) q[3];
sx q[3];
rz(-1.1793292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
