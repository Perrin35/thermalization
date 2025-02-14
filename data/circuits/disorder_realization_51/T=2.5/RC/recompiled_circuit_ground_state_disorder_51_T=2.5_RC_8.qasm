OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2188323) q[0];
sx q[0];
rz(-0.063475944) q[0];
sx q[0];
rz(1.7712964) q[0];
rz(0.81075794) q[1];
sx q[1];
rz(-1.4501362) q[1];
sx q[1];
rz(0.22051799) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761326) q[0];
sx q[0];
rz(-1.5031843) q[0];
sx q[0];
rz(-2.3071482) q[0];
x q[1];
rz(2.0283255) q[2];
sx q[2];
rz(-1.2731247) q[2];
sx q[2];
rz(-3.0106017) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.045956417) q[1];
sx q[1];
rz(-1.7929309) q[1];
sx q[1];
rz(1.5595116) q[1];
rz(-pi) q[2];
rz(-1.0706717) q[3];
sx q[3];
rz(-0.4005188) q[3];
sx q[3];
rz(2.5409215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.084879547) q[2];
sx q[2];
rz(-3.1352545) q[2];
sx q[2];
rz(2.32178) q[2];
rz(0.020717185) q[3];
sx q[3];
rz(-2.8262704) q[3];
sx q[3];
rz(-1.0516385) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3725975) q[0];
sx q[0];
rz(-0.1854493) q[0];
sx q[0];
rz(0.73082596) q[0];
rz(-1.4251047) q[1];
sx q[1];
rz(-2.4162879) q[1];
sx q[1];
rz(-2.6740429) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0529405) q[0];
sx q[0];
rz(-3.0346617) q[0];
sx q[0];
rz(1.8631975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7317762) q[2];
sx q[2];
rz(-1.5404319) q[2];
sx q[2];
rz(-1.8259468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7035571) q[1];
sx q[1];
rz(-2.272637) q[1];
sx q[1];
rz(-1.3753533) q[1];
rz(0.83722026) q[3];
sx q[3];
rz(-1.7927815) q[3];
sx q[3];
rz(2.8573621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9759489) q[2];
sx q[2];
rz(-3.1303945) q[2];
sx q[2];
rz(-0.33640081) q[2];
rz(0.30763704) q[3];
sx q[3];
rz(-0.010713723) q[3];
sx q[3];
rz(2.352534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0129608) q[0];
sx q[0];
rz(-0.19879453) q[0];
sx q[0];
rz(-0.13191731) q[0];
rz(-1.7031274) q[1];
sx q[1];
rz(-1.4142282) q[1];
sx q[1];
rz(1.6579113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0795201) q[0];
sx q[0];
rz(-1.2925081) q[0];
sx q[0];
rz(-2.5132283) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5882115) q[2];
sx q[2];
rz(-1.546374) q[2];
sx q[2];
rz(-1.2877854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36101389) q[1];
sx q[1];
rz(-0.11615651) q[1];
sx q[1];
rz(-2.0347982) q[1];
rz(0.68861945) q[3];
sx q[3];
rz(-3.0303749) q[3];
sx q[3];
rz(-2.2764549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5710473) q[2];
sx q[2];
rz(-1.3344301) q[2];
sx q[2];
rz(1.4578311) q[2];
rz(-0.7782065) q[3];
sx q[3];
rz(-3.1359378) q[3];
sx q[3];
rz(-0.2125423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0197765) q[0];
sx q[0];
rz(-1.3874522) q[0];
sx q[0];
rz(-0.29086599) q[0];
rz(-1.5930755) q[1];
sx q[1];
rz(-0.80268186) q[1];
sx q[1];
rz(-3.1150418) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.247581) q[0];
sx q[0];
rz(-2.6373626) q[0];
sx q[0];
rz(2.3419041) q[0];
x q[1];
rz(-1.5174116) q[2];
sx q[2];
rz(-2.0211271) q[2];
sx q[2];
rz(0.93102294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3624653) q[1];
sx q[1];
rz(-2.315265) q[1];
sx q[1];
rz(-2.4646548) q[1];
x q[2];
rz(3.13843) q[3];
sx q[3];
rz(-1.5517276) q[3];
sx q[3];
rz(2.9657488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2619541) q[2];
sx q[2];
rz(-3.1231572) q[2];
sx q[2];
rz(-0.43799841) q[2];
rz(-2.0336464) q[3];
sx q[3];
rz(-3.1236533) q[3];
sx q[3];
rz(-1.7855135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22439013) q[0];
sx q[0];
rz(-0.49773911) q[0];
sx q[0];
rz(0.14601953) q[0];
rz(-3.0085425) q[1];
sx q[1];
rz(-1.0597884) q[1];
sx q[1];
rz(2.6859247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5600223) q[0];
sx q[0];
rz(-3.0455493) q[0];
sx q[0];
rz(0.21012975) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5495569) q[2];
sx q[2];
rz(-1.3130273) q[2];
sx q[2];
rz(1.3883615) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4670438) q[1];
sx q[1];
rz(-1.70948) q[1];
sx q[1];
rz(-1.8013181) q[1];
rz(-0.066448575) q[3];
sx q[3];
rz(-0.22804582) q[3];
sx q[3];
rz(-1.4888637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8749775) q[2];
sx q[2];
rz(-0.019860331) q[2];
sx q[2];
rz(1.2561426) q[2];
rz(-2.6202294) q[3];
sx q[3];
rz(-0.25786906) q[3];
sx q[3];
rz(2.7969587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94619036) q[0];
sx q[0];
rz(-2.7552216) q[0];
sx q[0];
rz(0.56889164) q[0];
rz(-2.9733114) q[1];
sx q[1];
rz(-1.556309) q[1];
sx q[1];
rz(0.010244244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83547384) q[0];
sx q[0];
rz(-0.17806192) q[0];
sx q[0];
rz(2.8639069) q[0];
x q[1];
rz(3.1222759) q[2];
sx q[2];
rz(-1.7512808) q[2];
sx q[2];
rz(-1.3767124) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6737187) q[1];
sx q[1];
rz(-1.5281902) q[1];
sx q[1];
rz(-3.1070479) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63416173) q[3];
sx q[3];
rz(-1.595257) q[3];
sx q[3];
rz(2.0195877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16113082) q[2];
sx q[2];
rz(-2.9176596) q[2];
sx q[2];
rz(2.0437608) q[2];
rz(0.3499507) q[3];
sx q[3];
rz(-1.8643458) q[3];
sx q[3];
rz(2.2866975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2646645) q[0];
sx q[0];
rz(-2.9783037) q[0];
sx q[0];
rz(0.40633416) q[0];
rz(-1.3111275) q[1];
sx q[1];
rz(-2.6268112) q[1];
sx q[1];
rz(-0.11022551) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0144008) q[0];
sx q[0];
rz(-1.6699635) q[0];
sx q[0];
rz(-1.5983461) q[0];
rz(-0.0037721088) q[2];
sx q[2];
rz(-1.5654148) q[2];
sx q[2];
rz(0.37056915) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69198845) q[1];
sx q[1];
rz(-1.9983728) q[1];
sx q[1];
rz(3.1319992) q[1];
x q[2];
rz(-0.87894999) q[3];
sx q[3];
rz(-0.96262299) q[3];
sx q[3];
rz(1.5076306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.280764) q[2];
sx q[2];
rz(-3.132756) q[2];
sx q[2];
rz(-2.3542118) q[2];
rz(0.27960882) q[3];
sx q[3];
rz(-0.044737261) q[3];
sx q[3];
rz(2.4217822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0019919458) q[0];
sx q[0];
rz(-0.75374341) q[0];
sx q[0];
rz(1.4308223) q[0];
rz(2.9665973) q[1];
sx q[1];
rz(-0.7376968) q[1];
sx q[1];
rz(1.4281248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7044107) q[0];
sx q[0];
rz(-1.5927654) q[0];
sx q[0];
rz(-3.1385413) q[0];
rz(1.578184) q[2];
sx q[2];
rz(-1.566717) q[2];
sx q[2];
rz(2.5830327) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1154394) q[1];
sx q[1];
rz(-2.8255531) q[1];
sx q[1];
rz(2.1137966) q[1];
rz(1.346502) q[3];
sx q[3];
rz(-1.845775) q[3];
sx q[3];
rz(-1.1077653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.78508198) q[2];
sx q[2];
rz(-0.76866895) q[2];
sx q[2];
rz(1.53995) q[2];
rz(-2.8421863) q[3];
sx q[3];
rz(-1.6573903) q[3];
sx q[3];
rz(0.33777344) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47141075) q[0];
sx q[0];
rz(-0.010936745) q[0];
sx q[0];
rz(0.48728824) q[0];
rz(1.4393073) q[1];
sx q[1];
rz(-2.3479562) q[1];
sx q[1];
rz(-0.38584858) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13284616) q[0];
sx q[0];
rz(-1.9632959) q[0];
sx q[0];
rz(-1.5061492) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11211498) q[2];
sx q[2];
rz(-1.9939853) q[2];
sx q[2];
rz(-1.8433169) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1362366) q[1];
sx q[1];
rz(-0.87974977) q[1];
sx q[1];
rz(-1.6137531) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8063784) q[3];
sx q[3];
rz(-1.6044558) q[3];
sx q[3];
rz(2.0884643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7172598) q[2];
sx q[2];
rz(-0.029742664) q[2];
sx q[2];
rz(2.7359803) q[2];
rz(0.82617104) q[3];
sx q[3];
rz(-0.043353733) q[3];
sx q[3];
rz(-0.9894754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9721603) q[0];
sx q[0];
rz(-0.5652453) q[0];
sx q[0];
rz(1.3954847) q[0];
rz(1.5211498) q[1];
sx q[1];
rz(-0.65137678) q[1];
sx q[1];
rz(-1.7924538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5842218) q[0];
sx q[0];
rz(-0.85600805) q[0];
sx q[0];
rz(-0.28066158) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1824532) q[2];
sx q[2];
rz(-2.8504319) q[2];
sx q[2];
rz(-2.3860402) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5384465) q[1];
sx q[1];
rz(-1.6308404) q[1];
sx q[1];
rz(-3.1272496) q[1];
rz(-pi) q[2];
rz(0.030142763) q[3];
sx q[3];
rz(-1.5557914) q[3];
sx q[3];
rz(2.4001013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.32538432) q[2];
sx q[2];
rz(-0.00486972) q[2];
sx q[2];
rz(1.7083141) q[2];
rz(-1.8520744) q[3];
sx q[3];
rz(-0.09434814) q[3];
sx q[3];
rz(1.2963699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077944003) q[0];
sx q[0];
rz(-2.1049121) q[0];
sx q[0];
rz(-2.5459469) q[0];
rz(0.044364914) q[1];
sx q[1];
rz(-0.27635074) q[1];
sx q[1];
rz(1.9377294) q[1];
rz(-1.1677713) q[2];
sx q[2];
rz(-0.51196212) q[2];
sx q[2];
rz(-2.5434824) q[2];
rz(-1.1607004) q[3];
sx q[3];
rz(-1.5103419) q[3];
sx q[3];
rz(1.7132669) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
