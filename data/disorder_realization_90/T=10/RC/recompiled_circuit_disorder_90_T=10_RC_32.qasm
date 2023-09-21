OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(0.13983146) q[0];
sx q[0];
rz(10.034372) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(-0.087021526) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60458175) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(2.1405311) q[0];
x q[1];
rz(3.0084228) q[2];
sx q[2];
rz(-2.3699017) q[2];
sx q[2];
rz(-1.1222249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7056071) q[1];
sx q[1];
rz(-1.7090624) q[1];
sx q[1];
rz(-2.0558753) q[1];
rz(0.10591412) q[3];
sx q[3];
rz(-2.7273791) q[3];
sx q[3];
rz(-0.29990444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0027851) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(-0.3368245) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8841298) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.194838) q[0];
rz(-0.082611235) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9449687) q[0];
sx q[0];
rz(-1.9163016) q[0];
sx q[0];
rz(1.1254805) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2194901) q[2];
sx q[2];
rz(-0.32425913) q[2];
sx q[2];
rz(-2.1622554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5217168) q[1];
sx q[1];
rz(-0.65345018) q[1];
sx q[1];
rz(0.39342777) q[1];
x q[2];
rz(3.0307816) q[3];
sx q[3];
rz(-2.6786945) q[3];
sx q[3];
rz(-0.09679951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3327545) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(-2.9648798) q[2];
rz(0.79408944) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(-0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2394543) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(2.7368271) q[0];
rz(-1.8602712) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(1.3084897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3151911) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(2.1000923) q[0];
rz(-2.7576315) q[2];
sx q[2];
rz(-0.57120354) q[2];
sx q[2];
rz(-2.2307894) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82061003) q[1];
sx q[1];
rz(-0.87249407) q[1];
sx q[1];
rz(0.33336158) q[1];
rz(-pi) q[2];
rz(-1.4806467) q[3];
sx q[3];
rz(-1.1336859) q[3];
sx q[3];
rz(2.3428832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.219316) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(-3.1211839) q[2];
rz(-2.0698047) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0650478) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(1.5989499) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82614233) q[2];
sx q[2];
rz(-0.75781265) q[2];
sx q[2];
rz(-0.10105029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6009439) q[1];
sx q[1];
rz(-1.800866) q[1];
sx q[1];
rz(0.01407108) q[1];
x q[2];
rz(1.3435059) q[3];
sx q[3];
rz(-0.82633457) q[3];
sx q[3];
rz(-0.28971653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2477734) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(-0.34269732) q[2];
rz(-1.6977067) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(2.4263583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(3.0551531) q[0];
rz(1.7516288) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(2.6729029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2738004) q[0];
sx q[0];
rz(-0.81875728) q[0];
sx q[0];
rz(1.0206945) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7519978) q[2];
sx q[2];
rz(-1.8743519) q[2];
sx q[2];
rz(0.62523491) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3846181) q[1];
sx q[1];
rz(-0.24090919) q[1];
sx q[1];
rz(-0.5730281) q[1];
rz(-pi) q[2];
rz(-0.080805578) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9995352) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-0.86432499) q[2];
rz(2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(2.8488081) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(-0.19700225) q[0];
rz(1.3621832) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(-0.33624712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.907619) q[0];
sx q[0];
rz(-0.95171463) q[0];
sx q[0];
rz(-1.0236077) q[0];
x q[1];
rz(0.98724987) q[2];
sx q[2];
rz(-1.0580214) q[2];
sx q[2];
rz(3.0836881) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0834004) q[1];
sx q[1];
rz(-0.5914878) q[1];
sx q[1];
rz(3.0565492) q[1];
rz(-pi) q[2];
rz(2.9365262) q[3];
sx q[3];
rz(-2.4347217) q[3];
sx q[3];
rz(-2.6170078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53529915) q[2];
sx q[2];
rz(-1.665411) q[2];
sx q[2];
rz(-0.84632787) q[2];
rz(1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7476615) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(-0.25892648) q[0];
rz(-1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-2.0475725) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8261995) q[0];
sx q[0];
rz(-2.0059735) q[0];
sx q[0];
rz(0.3776334) q[0];
rz(-pi) q[1];
rz(1.6254243) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(-0.54578997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.104407) q[1];
sx q[1];
rz(-1.6652602) q[1];
sx q[1];
rz(0.027992804) q[1];
rz(-pi) q[2];
x q[2];
rz(1.592698) q[3];
sx q[3];
rz(-1.6688445) q[3];
sx q[3];
rz(-2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(0.32067498) q[2];
rz(-0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(-2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
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
rz(2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(-2.136769) q[0];
rz(-0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(-0.80387962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893387) q[0];
sx q[0];
rz(-2.5544871) q[0];
sx q[0];
rz(-0.071285204) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3857533) q[2];
sx q[2];
rz(-1.0244601) q[2];
sx q[2];
rz(-1.9667369) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56402928) q[1];
sx q[1];
rz(-1.4821577) q[1];
sx q[1];
rz(-0.010239756) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6635267) q[3];
sx q[3];
rz(-1.6420206) q[3];
sx q[3];
rz(1.0458667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(1.0104898) q[2];
rz(-0.60339749) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(-2.5914014) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5378961) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(0.4075152) q[0];
rz(-2.852476) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(2.3908652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2931965) q[0];
sx q[0];
rz(-1.2788532) q[0];
sx q[0];
rz(1.2743203) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27955555) q[2];
sx q[2];
rz(-1.6492372) q[2];
sx q[2];
rz(-0.58141764) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1322861) q[1];
sx q[1];
rz(-2.7704151) q[1];
sx q[1];
rz(0.98307857) q[1];
x q[2];
rz(-0.30131807) q[3];
sx q[3];
rz(-1.5958061) q[3];
sx q[3];
rz(-3.0400288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0687381) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(-1.3646305) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5831379) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.3903842) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(1.6814544) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5064142) q[0];
sx q[0];
rz(-1.1513452) q[0];
sx q[0];
rz(1.0181396) q[0];
rz(-1.5286469) q[2];
sx q[2];
rz(-0.50439207) q[2];
sx q[2];
rz(-1.991589) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3223443) q[1];
sx q[1];
rz(-1.7548314) q[1];
sx q[1];
rz(-0.40477246) q[1];
rz(-pi) q[2];
rz(-0.7695997) q[3];
sx q[3];
rz(-2.1160612) q[3];
sx q[3];
rz(1.256497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.6798518) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-2.6859443) q[3];
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
rz(-2.5158952) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.3810146) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(1.21576) q[2];
sx q[2];
rz(-2.3264865) q[2];
sx q[2];
rz(-2.9441499) q[2];
rz(-1.1827042) q[3];
sx q[3];
rz(-0.93508616) q[3];
sx q[3];
rz(-1.3133776) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];