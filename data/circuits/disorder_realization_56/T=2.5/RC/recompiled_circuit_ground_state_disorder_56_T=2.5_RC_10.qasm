OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96707764) q[0];
sx q[0];
rz(-1.4580026) q[0];
sx q[0];
rz(-2.8414677) q[0];
rz(1.1974273) q[1];
sx q[1];
rz(4.6680968) q[1];
sx q[1];
rz(11.065281) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34632698) q[0];
sx q[0];
rz(-3.1259968) q[0];
sx q[0];
rz(-2.8898284) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.757171) q[2];
sx q[2];
rz(-0.84058207) q[2];
sx q[2];
rz(-0.79193774) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0847124) q[1];
sx q[1];
rz(-2.4507629) q[1];
sx q[1];
rz(-2.7954196) q[1];
x q[2];
rz(-3.1367125) q[3];
sx q[3];
rz(-2.3327851) q[3];
sx q[3];
rz(0.67833662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.40452051) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(-2.4386151) q[2];
rz(2.5718555) q[3];
sx q[3];
rz(-0.8786141) q[3];
sx q[3];
rz(-1.5574633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.1341781) q[0];
sx q[0];
rz(-2.5226722) q[0];
sx q[0];
rz(1.2331569) q[0];
rz(-0.59457072) q[1];
sx q[1];
rz(-1.0392799) q[1];
sx q[1];
rz(-1.0979244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63609517) q[0];
sx q[0];
rz(-1.5042802) q[0];
sx q[0];
rz(-1.7719222) q[0];
rz(1.1851083) q[2];
sx q[2];
rz(-2.0314221) q[2];
sx q[2];
rz(-0.88049358) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6095396) q[1];
sx q[1];
rz(-1.5960448) q[1];
sx q[1];
rz(-1.6482501) q[1];
rz(-pi) q[2];
rz(0.79165801) q[3];
sx q[3];
rz(-1.6362305) q[3];
sx q[3];
rz(1.017788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42830959) q[2];
sx q[2];
rz(-1.8417336) q[2];
sx q[2];
rz(-1.5353047) q[2];
rz(-2.4226268) q[3];
sx q[3];
rz(-0.9459559) q[3];
sx q[3];
rz(-2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3035901) q[0];
sx q[0];
rz(-0.8135697) q[0];
sx q[0];
rz(0.59233061) q[0];
rz(1.2447641) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(2.9959784) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.515265) q[0];
sx q[0];
rz(-1.087731) q[0];
sx q[0];
rz(-1.7240216) q[0];
rz(-3.0119409) q[2];
sx q[2];
rz(-2.2551422) q[2];
sx q[2];
rz(0.42829542) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.023886746) q[1];
sx q[1];
rz(-0.63359208) q[1];
sx q[1];
rz(2.0633477) q[1];
rz(-2.0592732) q[3];
sx q[3];
rz(-1.3311738) q[3];
sx q[3];
rz(-2.3828196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40272063) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(1.9090451) q[2];
rz(-0.23972073) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(-2.6565552) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0204912) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(0.030601587) q[0];
rz(-2.5864511) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(-1.0928924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3557778) q[0];
sx q[0];
rz(-2.2270538) q[0];
sx q[0];
rz(2.2758765) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76538122) q[2];
sx q[2];
rz(-1.6134521) q[2];
sx q[2];
rz(-0.19848196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2177201) q[1];
sx q[1];
rz(-1.3663379) q[1];
sx q[1];
rz(-0.93074284) q[1];
rz(-pi) q[2];
rz(0.83715688) q[3];
sx q[3];
rz(-0.87292307) q[3];
sx q[3];
rz(-0.18750377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0838202) q[2];
sx q[2];
rz(-1.3670992) q[2];
sx q[2];
rz(-2.8687381) q[2];
rz(-1.3911635) q[3];
sx q[3];
rz(-2.302156) q[3];
sx q[3];
rz(-2.1896037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.4680173) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(-1.8433628) q[0];
rz(0.6908373) q[1];
sx q[1];
rz(-1.9419443) q[1];
sx q[1];
rz(-1.5265436) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1896601) q[0];
sx q[0];
rz(-2.2554923) q[0];
sx q[0];
rz(-0.71345274) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.853302) q[2];
sx q[2];
rz(-1.9108678) q[2];
sx q[2];
rz(1.9486537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5954819) q[1];
sx q[1];
rz(-2.271436) q[1];
sx q[1];
rz(-0.033820669) q[1];
rz(-0.68256179) q[3];
sx q[3];
rz(-0.45387156) q[3];
sx q[3];
rz(2.0477432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7319298) q[2];
sx q[2];
rz(-2.3356428) q[2];
sx q[2];
rz(0.16732495) q[2];
rz(1.3903728) q[3];
sx q[3];
rz(-0.53283397) q[3];
sx q[3];
rz(0.65756857) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996027) q[0];
sx q[0];
rz(-2.7775192) q[0];
sx q[0];
rz(1.2784736) q[0];
rz(-1.4356042) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(-0.37567589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091051) q[0];
sx q[0];
rz(-1.300749) q[0];
sx q[0];
rz(0.23464111) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4009052) q[2];
sx q[2];
rz(-1.7180016) q[2];
sx q[2];
rz(-0.81913713) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76498079) q[1];
sx q[1];
rz(-1.3479479) q[1];
sx q[1];
rz(-1.8972023) q[1];
x q[2];
rz(1.814498) q[3];
sx q[3];
rz(-2.3117495) q[3];
sx q[3];
rz(-2.0654701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0659236) q[2];
sx q[2];
rz(-2.0772987) q[2];
sx q[2];
rz(2.2806878) q[2];
rz(-1.9979477) q[3];
sx q[3];
rz(-0.30503169) q[3];
sx q[3];
rz(2.8633269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1160195) q[0];
sx q[0];
rz(-1.8208068) q[0];
sx q[0];
rz(-2.386911) q[0];
rz(0.93983752) q[1];
sx q[1];
rz(-2.266423) q[1];
sx q[1];
rz(1.2831203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1484598) q[0];
sx q[0];
rz(-2.8348036) q[0];
sx q[0];
rz(1.5932142) q[0];
rz(-1.1129598) q[2];
sx q[2];
rz(-1.2414059) q[2];
sx q[2];
rz(3.0731346) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2055473) q[1];
sx q[1];
rz(-2.6523771) q[1];
sx q[1];
rz(0.65566109) q[1];
rz(-0.8831034) q[3];
sx q[3];
rz(-1.1576443) q[3];
sx q[3];
rz(2.2108111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54520404) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(2.0580573) q[2];
rz(0.74294535) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(-1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899406) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(2.3866744) q[0];
rz(0.49916357) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(0.69854936) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7592418) q[0];
sx q[0];
rz(-2.1672591) q[0];
sx q[0];
rz(-1.9475219) q[0];
x q[1];
rz(2.8590747) q[2];
sx q[2];
rz(-0.48136863) q[2];
sx q[2];
rz(3.0309739) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6900151) q[1];
sx q[1];
rz(-1.5908341) q[1];
sx q[1];
rz(-1.0385787) q[1];
x q[2];
rz(1.9898765) q[3];
sx q[3];
rz(-1.5718096) q[3];
sx q[3];
rz(-0.59194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90427202) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(2.3392056) q[2];
rz(-2.5896416) q[3];
sx q[3];
rz(-2.3410083) q[3];
sx q[3];
rz(-0.62057453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9629843) q[0];
sx q[0];
rz(-1.7010138) q[0];
sx q[0];
rz(-2.0340023) q[0];
rz(-1.7604609) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(-1.507087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3920306) q[0];
sx q[0];
rz(-1.1774585) q[0];
sx q[0];
rz(-1.7225207) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92169833) q[2];
sx q[2];
rz(-0.49921152) q[2];
sx q[2];
rz(1.8550903) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17738755) q[1];
sx q[1];
rz(-1.2086165) q[1];
sx q[1];
rz(-1.4153348) q[1];
rz(0.082066925) q[3];
sx q[3];
rz(-2.3130199) q[3];
sx q[3];
rz(-2.9794793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0570809) q[2];
sx q[2];
rz(-0.19597404) q[2];
sx q[2];
rz(-1.2723119) q[2];
rz(-1.48014) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(-0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24499527) q[0];
sx q[0];
rz(-2.5801881) q[0];
sx q[0];
rz(-1.8580612) q[0];
rz(3.1237579) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(3.0160115) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9367919) q[0];
sx q[0];
rz(-1.556634) q[0];
sx q[0];
rz(-1.0165748) q[0];
rz(-pi) q[1];
rz(0.56350817) q[2];
sx q[2];
rz(-0.76125188) q[2];
sx q[2];
rz(0.13680563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6148147) q[1];
sx q[1];
rz(-1.1452951) q[1];
sx q[1];
rz(-2.6761901) q[1];
rz(-2.2108881) q[3];
sx q[3];
rz(-0.32880201) q[3];
sx q[3];
rz(2.9931675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77832001) q[2];
sx q[2];
rz(-1.2691701) q[2];
sx q[2];
rz(-0.8030836) q[2];
rz(0.17656365) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(-2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9618027) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(-1.8190307) q[1];
sx q[1];
rz(-1.3643199) q[1];
sx q[1];
rz(1.551569) q[1];
rz(0.24257913) q[2];
sx q[2];
rz(-1.6318113) q[2];
sx q[2];
rz(3.0616374) q[2];
rz(-1.6311036) q[3];
sx q[3];
rz(-1.2338439) q[3];
sx q[3];
rz(-3.0046786) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
