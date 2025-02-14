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
rz(-1.9441654) q[1];
sx q[1];
rz(-1.5265042) q[1];
sx q[1];
rz(-1.6405029) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59812058) q[0];
sx q[0];
rz(-1.5859005) q[0];
sx q[0];
rz(-1.5669109) q[0];
rz(-pi) q[1];
rz(-2.757171) q[2];
sx q[2];
rz(-2.3010106) q[2];
sx q[2];
rz(0.79193774) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3565271) q[1];
sx q[1];
rz(-1.7887113) q[1];
sx q[1];
rz(2.480605) q[1];
x q[2];
rz(1.5656823) q[3];
sx q[3];
rz(-2.3795914) q[3];
sx q[3];
rz(0.68540547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40452051) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(2.4386151) q[2];
rz(2.5718555) q[3];
sx q[3];
rz(-0.8786141) q[3];
sx q[3];
rz(-1.5574633) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1341781) q[0];
sx q[0];
rz(-2.5226722) q[0];
sx q[0];
rz(-1.9084357) q[0];
rz(-0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(1.0979244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5220338) q[0];
sx q[0];
rz(-0.21169835) q[0];
sx q[0];
rz(-1.2489399) q[0];
x q[1];
rz(-2.4929201) q[2];
sx q[2];
rz(-0.59174109) q[2];
sx q[2];
rz(1.6206738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2757551) q[1];
sx q[1];
rz(-3.0601353) q[1];
sx q[1];
rz(1.8862731) q[1];
x q[2];
rz(-3.0497562) q[3];
sx q[3];
rz(-0.79376924) q[3];
sx q[3];
rz(-2.5240999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42830959) q[2];
sx q[2];
rz(-1.8417336) q[2];
sx q[2];
rz(1.5353047) q[2];
rz(-0.71896583) q[3];
sx q[3];
rz(-0.9459559) q[3];
sx q[3];
rz(2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83800256) q[0];
sx q[0];
rz(-2.328023) q[0];
sx q[0];
rz(0.59233061) q[0];
rz(-1.8968286) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(2.9959784) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94731047) q[0];
sx q[0];
rz(-2.6366451) q[0];
sx q[0];
rz(-0.28316747) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12965173) q[2];
sx q[2];
rz(-0.88645049) q[2];
sx q[2];
rz(-2.7132972) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.023886746) q[1];
sx q[1];
rz(-2.5080006) q[1];
sx q[1];
rz(-1.078245) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0907902) q[3];
sx q[3];
rz(-2.601805) q[3];
sx q[3];
rz(-1.2320428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40272063) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(1.2325475) q[2];
rz(2.9018719) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(-2.6565552) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0204912) q[0];
sx q[0];
rz(-2.4395269) q[0];
sx q[0];
rz(3.1109911) q[0];
rz(0.55514151) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(-1.0928924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3557778) q[0];
sx q[0];
rz(-2.2270538) q[0];
sx q[0];
rz(-2.2758765) q[0];
rz(-1.6299155) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(1.3313683) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92387256) q[1];
sx q[1];
rz(-1.3663379) q[1];
sx q[1];
rz(-2.2108498) q[1];
rz(-0.84597702) q[3];
sx q[3];
rz(-2.1095157) q[3];
sx q[3];
rz(-1.2332476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0577724) q[2];
sx q[2];
rz(-1.3670992) q[2];
sx q[2];
rz(-2.8687381) q[2];
rz(-1.3911635) q[3];
sx q[3];
rz(-0.83943668) q[3];
sx q[3];
rz(-0.951989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4680173) q[0];
sx q[0];
rz(-1.3382358) q[0];
sx q[0];
rz(-1.8433628) q[0];
rz(-2.4507554) q[1];
sx q[1];
rz(-1.1996484) q[1];
sx q[1];
rz(1.5265436) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1896601) q[0];
sx q[0];
rz(-2.2554923) q[0];
sx q[0];
rz(0.71345274) q[0];
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
rz(-2.4936798) q[1];
sx q[1];
rz(-2.4402752) q[1];
sx q[1];
rz(-1.6108684) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36187343) q[3];
sx q[3];
rz(-1.2905777) q[3];
sx q[3];
rz(1.1079196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7319298) q[2];
sx q[2];
rz(-0.80594984) q[2];
sx q[2];
rz(0.16732495) q[2];
rz(-1.7512199) q[3];
sx q[3];
rz(-0.53283397) q[3];
sx q[3];
rz(-2.4840241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.14199) q[0];
sx q[0];
rz(-0.36407343) q[0];
sx q[0];
rz(1.863119) q[0];
rz(1.4356042) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(-2.7659168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77462353) q[0];
sx q[0];
rz(-1.7967829) q[0];
sx q[0];
rz(1.8480728) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2907545) q[2];
sx q[2];
rz(-2.9172643) q[2];
sx q[2];
rz(-1.6825324) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.73114016) q[1];
sx q[1];
rz(-1.8888432) q[1];
sx q[1];
rz(0.2348301) q[1];
rz(-pi) q[2];
rz(-2.385684) q[3];
sx q[3];
rz(-1.7497853) q[3];
sx q[3];
rz(-2.8132015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0756691) q[2];
sx q[2];
rz(-2.0772987) q[2];
sx q[2];
rz(2.2806878) q[2];
rz(1.9979477) q[3];
sx q[3];
rz(-0.30503169) q[3];
sx q[3];
rz(0.27826571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025573108) q[0];
sx q[0];
rz(-1.8208068) q[0];
sx q[0];
rz(0.75468165) q[0];
rz(0.93983752) q[1];
sx q[1];
rz(-2.266423) q[1];
sx q[1];
rz(1.2831203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5425576) q[0];
sx q[0];
rz(-1.5640266) q[0];
sx q[0];
rz(1.2640796) q[0];
rz(-pi) q[1];
rz(2.2290984) q[2];
sx q[2];
rz(-0.55710885) q[2];
sx q[2];
rz(-1.0583641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2312113) q[1];
sx q[1];
rz(-1.2802135) q[1];
sx q[1];
rz(-0.39931764) q[1];
rz(-pi) q[2];
rz(-0.96638443) q[3];
sx q[3];
rz(-0.78456402) q[3];
sx q[3];
rz(-2.0469637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54520404) q[2];
sx q[2];
rz(-1.117492) q[2];
sx q[2];
rz(1.0835353) q[2];
rz(0.74294535) q[3];
sx q[3];
rz(-1.6094306) q[3];
sx q[3];
rz(1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899406) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(2.3866744) q[0];
rz(2.6424291) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(2.4430433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23125739) q[0];
sx q[0];
rz(-0.6930348) q[0];
sx q[0];
rz(0.49654754) q[0];
rz(-0.28251799) q[2];
sx q[2];
rz(-0.48136863) q[2];
sx q[2];
rz(3.0309739) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9883635) q[1];
sx q[1];
rz(-2.6090342) q[1];
sx q[1];
rz(-1.531324) q[1];
rz(-1.9898765) q[3];
sx q[3];
rz(-1.5718096) q[3];
sx q[3];
rz(-2.5496497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90427202) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(2.3392056) q[2];
rz(2.5896416) q[3];
sx q[3];
rz(-0.80058432) q[3];
sx q[3];
rz(-0.62057453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17860831) q[0];
sx q[0];
rz(-1.7010138) q[0];
sx q[0];
rz(-1.1075903) q[0];
rz(-1.3811318) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(1.507087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0213623) q[0];
sx q[0];
rz(-1.430738) q[0];
sx q[0];
rz(-2.7441478) q[0];
x q[1];
rz(-0.31839026) q[2];
sx q[2];
rz(-1.9620775) q[2];
sx q[2];
rz(1.9992145) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.448882) q[1];
sx q[1];
rz(-1.4254942) q[1];
sx q[1];
rz(2.7753745) q[1];
rz(-pi) q[2];
rz(-2.3146995) q[3];
sx q[3];
rz(-1.6312459) q[3];
sx q[3];
rz(1.7884487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0845118) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(-1.2723119) q[2];
rz(-1.48014) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(2.541466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24499527) q[0];
sx q[0];
rz(-0.56140459) q[0];
sx q[0];
rz(1.2835314) q[0];
rz(3.1237579) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(3.0160115) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3572306) q[0];
sx q[0];
rz(-1.0166369) q[0];
sx q[0];
rz(-0.016655075) q[0];
x q[1];
rz(0.56350817) q[2];
sx q[2];
rz(-2.3803408) q[2];
sx q[2];
rz(-0.13680563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5267779) q[1];
sx q[1];
rz(-1.1452951) q[1];
sx q[1];
rz(-2.6761901) q[1];
x q[2];
rz(2.2108881) q[3];
sx q[3];
rz(-2.8127906) q[3];
sx q[3];
rz(2.9931675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77832001) q[2];
sx q[2];
rz(-1.8724226) q[2];
sx q[2];
rz(0.8030836) q[2];
rz(0.17656365) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(-2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9618027) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(1.322562) q[1];
sx q[1];
rz(-1.3643199) q[1];
sx q[1];
rz(1.551569) q[1];
rz(-0.24904574) q[2];
sx q[2];
rz(-0.24998827) q[2];
sx q[2];
rz(-1.892358) q[2];
rz(1.5104891) q[3];
sx q[3];
rz(-1.2338439) q[3];
sx q[3];
rz(-3.0046786) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
