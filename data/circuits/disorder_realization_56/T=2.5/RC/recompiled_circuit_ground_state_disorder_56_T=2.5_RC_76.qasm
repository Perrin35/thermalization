OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(-1.6835901) q[0];
sx q[0];
rz(-0.30012497) q[0];
rz(1.1974273) q[1];
sx q[1];
rz(-1.6150885) q[1];
sx q[1];
rz(1.6405029) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34632698) q[0];
sx q[0];
rz(-3.1259968) q[0];
sx q[0];
rz(-0.25176425) q[0];
rz(1.1741224) q[2];
sx q[2];
rz(-0.80840092) q[2];
sx q[2];
rz(-1.3371181) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7607494) q[1];
sx q[1];
rz(-0.92807209) q[1];
sx q[1];
rz(1.2973143) q[1];
rz(-pi) q[2];
x q[2];
rz(2.332791) q[3];
sx q[3];
rz(-1.5743269) q[3];
sx q[3];
rz(0.88909066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40452051) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(-0.70297757) q[2];
rz(0.56973714) q[3];
sx q[3];
rz(-2.2629786) q[3];
sx q[3];
rz(1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
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
rz(0.59457072) q[1];
sx q[1];
rz(-1.0392799) q[1];
sx q[1];
rz(-2.0436683) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6195589) q[0];
sx q[0];
rz(-0.21169835) q[0];
sx q[0];
rz(-1.2489399) q[0];
rz(-2.4929201) q[2];
sx q[2];
rz(-0.59174109) q[2];
sx q[2];
rz(1.6206738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2757551) q[1];
sx q[1];
rz(-3.0601353) q[1];
sx q[1];
rz(-1.8862731) q[1];
rz(2.3499346) q[3];
sx q[3];
rz(-1.6362305) q[3];
sx q[3];
rz(-1.017788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42830959) q[2];
sx q[2];
rz(-1.8417336) q[2];
sx q[2];
rz(1.606288) q[2];
rz(-0.71896583) q[3];
sx q[3];
rz(-2.1956367) q[3];
sx q[3];
rz(-2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.14561428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.515265) q[0];
sx q[0];
rz(-2.0538616) q[0];
sx q[0];
rz(1.7240216) q[0];
x q[1];
rz(-1.4136366) q[2];
sx q[2];
rz(-2.4470235) q[2];
sx q[2];
rz(2.5099011) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1864724) q[1];
sx q[1];
rz(-1.8545517) q[1];
sx q[1];
rz(-2.1452745) q[1];
rz(-2.0508025) q[3];
sx q[3];
rz(-2.601805) q[3];
sx q[3];
rz(1.2320428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.738872) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(-1.2325475) q[2];
rz(-2.9018719) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(2.6565552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0204912) q[0];
sx q[0];
rz(-2.4395269) q[0];
sx q[0];
rz(3.1109911) q[0];
rz(2.5864511) q[1];
sx q[1];
rz(-1.4786913) q[1];
sx q[1];
rz(2.0487002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83704889) q[0];
sx q[0];
rz(-2.2185159) q[0];
sx q[0];
rz(-2.442028) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76538122) q[2];
sx q[2];
rz(-1.6134521) q[2];
sx q[2];
rz(0.19848196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2177201) q[1];
sx q[1];
rz(-1.3663379) q[1];
sx q[1];
rz(0.93074284) q[1];
x q[2];
rz(0.83715688) q[3];
sx q[3];
rz(-2.2686696) q[3];
sx q[3];
rz(0.18750377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0577724) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(0.27285451) q[2];
rz(1.3911635) q[3];
sx q[3];
rz(-2.302156) q[3];
sx q[3];
rz(-0.951989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4680173) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(-1.8433628) q[0];
rz(-0.6908373) q[1];
sx q[1];
rz(-1.1996484) q[1];
sx q[1];
rz(1.615049) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0218791) q[0];
sx q[0];
rz(-2.102431) q[0];
sx q[0];
rz(-0.74703947) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2882907) q[2];
sx q[2];
rz(-1.9108678) q[2];
sx q[2];
rz(1.9486537) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5954819) q[1];
sx q[1];
rz(-2.271436) q[1];
sx q[1];
rz(3.107772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7797192) q[3];
sx q[3];
rz(-1.851015) q[3];
sx q[3];
rz(-2.0336731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4096628) q[2];
sx q[2];
rz(-0.80594984) q[2];
sx q[2];
rz(-2.9742677) q[2];
rz(1.3903728) q[3];
sx q[3];
rz(-0.53283397) q[3];
sx q[3];
rz(-2.4840241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.14199) q[0];
sx q[0];
rz(-0.36407343) q[0];
sx q[0];
rz(1.2784736) q[0];
rz(1.7059884) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(2.7659168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73248752) q[0];
sx q[0];
rz(-1.300749) q[0];
sx q[0];
rz(0.23464111) q[0];
rz(-0.85083811) q[2];
sx q[2];
rz(-0.22432835) q[2];
sx q[2];
rz(1.4590603) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7573812) q[1];
sx q[1];
rz(-2.7486292) q[1];
sx q[1];
rz(-2.1860366) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.814498) q[3];
sx q[3];
rz(-2.3117495) q[3];
sx q[3];
rz(2.0654701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0659236) q[2];
sx q[2];
rz(-2.0772987) q[2];
sx q[2];
rz(-2.2806878) q[2];
rz(-1.143645) q[3];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.1160195) q[0];
sx q[0];
rz(-1.3207859) q[0];
sx q[0];
rz(0.75468165) q[0];
rz(0.93983752) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(1.8584724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5425576) q[0];
sx q[0];
rz(-1.5640266) q[0];
sx q[0];
rz(-1.2640796) q[0];
rz(-2.0286328) q[2];
sx q[2];
rz(-1.9001868) q[2];
sx q[2];
rz(-0.068458099) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2055473) q[1];
sx q[1];
rz(-2.6523771) q[1];
sx q[1];
rz(2.4859316) q[1];
rz(-pi) q[2];
rz(-0.8831034) q[3];
sx q[3];
rz(-1.1576443) q[3];
sx q[3];
rz(2.2108111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54520404) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(2.0580573) q[2];
rz(-2.3986473) q[3];
sx q[3];
rz(-1.6094306) q[3];
sx q[3];
rz(1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899406) q[0];
sx q[0];
rz(-2.3567663) q[0];
sx q[0];
rz(0.75491828) q[0];
rz(-0.49916357) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(2.4430433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23125739) q[0];
sx q[0];
rz(-0.6930348) q[0];
sx q[0];
rz(-0.49654754) q[0];
x q[1];
rz(2.6766308) q[2];
sx q[2];
rz(-1.7002281) q[2];
sx q[2];
rz(-1.7120106) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15322916) q[1];
sx q[1];
rz(-2.6090342) q[1];
sx q[1];
rz(1.6102686) q[1];
x q[2];
rz(1.9898765) q[3];
sx q[3];
rz(-1.5718096) q[3];
sx q[3];
rz(2.5496497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2373206) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(0.80238706) q[2];
rz(-2.5896416) q[3];
sx q[3];
rz(-0.80058432) q[3];
sx q[3];
rz(0.62057453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
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
rz(1.6345056) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7716146) q[0];
sx q[0];
rz(-2.7214337) q[0];
sx q[0];
rz(2.7922947) q[0];
rz(-0.31839026) q[2];
sx q[2];
rz(-1.1795151) q[2];
sx q[2];
rz(1.1423781) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17738755) q[1];
sx q[1];
rz(-1.9329762) q[1];
sx q[1];
rz(1.7262579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3146995) q[3];
sx q[3];
rz(-1.5103467) q[3];
sx q[3];
rz(-1.7884487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0845118) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(-1.8692807) q[2];
rz(1.48014) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8965974) q[0];
sx q[0];
rz(-2.5801881) q[0];
sx q[0];
rz(1.2835314) q[0];
rz(-3.1237579) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(0.12558118) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3572306) q[0];
sx q[0];
rz(-1.0166369) q[0];
sx q[0];
rz(0.016655075) q[0];
rz(-pi) q[1];
rz(-0.67809503) q[2];
sx q[2];
rz(-1.1934278) q[2];
sx q[2];
rz(-1.2785778) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.24841845) q[1];
sx q[1];
rz(-1.1496953) q[1];
sx q[1];
rz(-2.0401272) q[1];
rz(-2.2108881) q[3];
sx q[3];
rz(-2.8127906) q[3];
sx q[3];
rz(-2.9931675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77832001) q[2];
sx q[2];
rz(-1.8724226) q[2];
sx q[2];
rz(-2.3385091) q[2];
rz(2.965029) q[3];
sx q[3];
rz(-1.3253788) q[3];
sx q[3];
rz(0.77809063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17978996) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(1.8190307) q[1];
sx q[1];
rz(-1.7772728) q[1];
sx q[1];
rz(-1.5900236) q[1];
rz(-2.8990135) q[2];
sx q[2];
rz(-1.6318113) q[2];
sx q[2];
rz(3.0616374) q[2];
rz(-1.5104891) q[3];
sx q[3];
rz(-1.9077488) q[3];
sx q[3];
rz(0.13691402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
