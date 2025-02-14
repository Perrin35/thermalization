OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9118948) q[0];
sx q[0];
rz(-1.2461646) q[0];
sx q[0];
rz(-1.5204313) q[0];
rz(0.35960943) q[1];
sx q[1];
rz(-2.8728027) q[1];
sx q[1];
rz(-0.51528817) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20977192) q[0];
sx q[0];
rz(-1.2270667) q[0];
sx q[0];
rz(-1.3552257) q[0];
x q[1];
rz(1.2688387) q[2];
sx q[2];
rz(-1.0450604) q[2];
sx q[2];
rz(0.56509226) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2559515) q[1];
sx q[1];
rz(-0.4547387) q[1];
sx q[1];
rz(0.498147) q[1];
rz(0.73909594) q[3];
sx q[3];
rz(-2.1198273) q[3];
sx q[3];
rz(-2.9573553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9109853) q[2];
sx q[2];
rz(-1.5829986) q[2];
sx q[2];
rz(0.67414635) q[2];
rz(0.19787431) q[3];
sx q[3];
rz(-1.9015692) q[3];
sx q[3];
rz(1.6363293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.823371) q[0];
sx q[0];
rz(-2.8405393) q[0];
sx q[0];
rz(0.2163042) q[0];
rz(-3.0220616) q[1];
sx q[1];
rz(-1.1263589) q[1];
sx q[1];
rz(1.4215887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57935539) q[0];
sx q[0];
rz(-2.6174712) q[0];
sx q[0];
rz(-1.4116794) q[0];
rz(-1.3590712) q[2];
sx q[2];
rz(-1.1781409) q[2];
sx q[2];
rz(-0.77048877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2188637) q[1];
sx q[1];
rz(-2.0170171) q[1];
sx q[1];
rz(0.34901825) q[1];
rz(2.1561863) q[3];
sx q[3];
rz(-0.93275654) q[3];
sx q[3];
rz(0.93321484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8372832) q[2];
sx q[2];
rz(-1.195636) q[2];
sx q[2];
rz(-0.99011123) q[2];
rz(-0.75634161) q[3];
sx q[3];
rz(-0.6476616) q[3];
sx q[3];
rz(-0.17624779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15683098) q[0];
sx q[0];
rz(-2.3909843) q[0];
sx q[0];
rz(-1.1358776) q[0];
rz(-2.1319977) q[1];
sx q[1];
rz(-1.986809) q[1];
sx q[1];
rz(-0.44581595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71190055) q[0];
sx q[0];
rz(-1.4662379) q[0];
sx q[0];
rz(-3.0719396) q[0];
x q[1];
rz(0.03568825) q[2];
sx q[2];
rz(-0.70377398) q[2];
sx q[2];
rz(-1.698871) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.43128428) q[1];
sx q[1];
rz(-2.388622) q[1];
sx q[1];
rz(0.94628872) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6942676) q[3];
sx q[3];
rz(-1.619391) q[3];
sx q[3];
rz(0.20125869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77057114) q[2];
sx q[2];
rz(-2.1393445) q[2];
sx q[2];
rz(2.501343) q[2];
rz(0.69784969) q[3];
sx q[3];
rz(-1.4638487) q[3];
sx q[3];
rz(-2.8016134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3014389) q[0];
sx q[0];
rz(-2.2541663) q[0];
sx q[0];
rz(1.812717) q[0];
rz(1.9245194) q[1];
sx q[1];
rz(-0.71958676) q[1];
sx q[1];
rz(-2.365239) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6049801) q[0];
sx q[0];
rz(-2.1663499) q[0];
sx q[0];
rz(2.0230629) q[0];
rz(-0.66993454) q[2];
sx q[2];
rz(-0.70730722) q[2];
sx q[2];
rz(-0.41367999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94521071) q[1];
sx q[1];
rz(-2.8627765) q[1];
sx q[1];
rz(1.6120738) q[1];
rz(0.0018423041) q[3];
sx q[3];
rz(-2.4569965) q[3];
sx q[3];
rz(2.1948333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2886469) q[2];
sx q[2];
rz(-2.8198346) q[2];
sx q[2];
rz(1.9256437) q[2];
rz(-1.1207885) q[3];
sx q[3];
rz(-1.8782764) q[3];
sx q[3];
rz(0.88302511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585606) q[0];
sx q[0];
rz(-2.164542) q[0];
sx q[0];
rz(-0.46347722) q[0];
rz(1.0889168) q[1];
sx q[1];
rz(-1.8810279) q[1];
sx q[1];
rz(1.8225398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4396297) q[0];
sx q[0];
rz(-1.8314529) q[0];
sx q[0];
rz(-1.319665) q[0];
x q[1];
rz(0.16541188) q[2];
sx q[2];
rz(-2.3369419) q[2];
sx q[2];
rz(-1.9555447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0214349) q[1];
sx q[1];
rz(-1.48367) q[1];
sx q[1];
rz(-2.1539262) q[1];
rz(-pi) q[2];
rz(-1.1518258) q[3];
sx q[3];
rz(-1.4052466) q[3];
sx q[3];
rz(-0.74932428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9702381) q[2];
sx q[2];
rz(-0.78533185) q[2];
sx q[2];
rz(-0.63924092) q[2];
rz(2.3302737) q[3];
sx q[3];
rz(-0.48894426) q[3];
sx q[3];
rz(-1.5343687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92523471) q[0];
sx q[0];
rz(-2.7282867) q[0];
sx q[0];
rz(-0.21892029) q[0];
rz(1.2159329) q[1];
sx q[1];
rz(-0.464012) q[1];
sx q[1];
rz(-1.6485515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0279044) q[0];
sx q[0];
rz(-1.1104212) q[0];
sx q[0];
rz(1.7229863) q[0];
x q[1];
rz(0.0068917787) q[2];
sx q[2];
rz(-2.7240319) q[2];
sx q[2];
rz(1.895895) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1472125) q[1];
sx q[1];
rz(-1.5541422) q[1];
sx q[1];
rz(2.4092595) q[1];
x q[2];
rz(3.0424439) q[3];
sx q[3];
rz(-1.5739023) q[3];
sx q[3];
rz(0.24156027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1047989) q[2];
sx q[2];
rz(-1.7918189) q[2];
sx q[2];
rz(-1.3118504) q[2];
rz(-1.5051684) q[3];
sx q[3];
rz(-1.8421831) q[3];
sx q[3];
rz(2.6817491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5055607) q[0];
sx q[0];
rz(-0.4902896) q[0];
sx q[0];
rz(-1.7946515) q[0];
rz(-1.8147644) q[1];
sx q[1];
rz(-0.83419269) q[1];
sx q[1];
rz(-0.38536513) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037212278) q[0];
sx q[0];
rz(-1.2118441) q[0];
sx q[0];
rz(1.6045553) q[0];
x q[1];
rz(1.0406251) q[2];
sx q[2];
rz(-1.1405924) q[2];
sx q[2];
rz(-2.8456147) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3431541) q[1];
sx q[1];
rz(-1.1249518) q[1];
sx q[1];
rz(-1.6002161) q[1];
x q[2];
rz(-2.6680384) q[3];
sx q[3];
rz(-1.9723168) q[3];
sx q[3];
rz(-2.7724289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4555326) q[2];
sx q[2];
rz(-0.81642381) q[2];
sx q[2];
rz(-1.4875937) q[2];
rz(-2.9758596) q[3];
sx q[3];
rz(-2.3159852) q[3];
sx q[3];
rz(-0.17416557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1626749) q[0];
sx q[0];
rz(-0.62541494) q[0];
sx q[0];
rz(-0.22853525) q[0];
rz(2.8288815) q[1];
sx q[1];
rz(-0.87930185) q[1];
sx q[1];
rz(1.3794587) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40423597) q[0];
sx q[0];
rz(-1.9725058) q[0];
sx q[0];
rz(-1.6549003) q[0];
x q[1];
rz(-1.3993456) q[2];
sx q[2];
rz(-1.8905235) q[2];
sx q[2];
rz(-0.4188183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4743275) q[1];
sx q[1];
rz(-1.9196471) q[1];
sx q[1];
rz(-0.036880924) q[1];
x q[2];
rz(-2.5614235) q[3];
sx q[3];
rz(-2.3456367) q[3];
sx q[3];
rz(2.7027276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.032893) q[2];
sx q[2];
rz(-2.097082) q[2];
sx q[2];
rz(2.1006987) q[2];
rz(-0.4852455) q[3];
sx q[3];
rz(-1.8294561) q[3];
sx q[3];
rz(0.41788873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.8682206) q[0];
sx q[0];
rz(-0.98216787) q[0];
sx q[0];
rz(0.81106538) q[0];
rz(-1.2366933) q[1];
sx q[1];
rz(-0.86634723) q[1];
sx q[1];
rz(2.9467357) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92835125) q[0];
sx q[0];
rz(-1.8272864) q[0];
sx q[0];
rz(-0.24730206) q[0];
rz(-pi) q[1];
rz(-2.2213908) q[2];
sx q[2];
rz(-2.0363931) q[2];
sx q[2];
rz(1.0126142) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5611546) q[1];
sx q[1];
rz(-1.913979) q[1];
sx q[1];
rz(1.0943867) q[1];
rz(-pi) q[2];
rz(2.2236597) q[3];
sx q[3];
rz(-2.0313162) q[3];
sx q[3];
rz(-1.9916818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.98086944) q[2];
sx q[2];
rz(-1.8941433) q[2];
sx q[2];
rz(2.5998739) q[2];
rz(-1.8187652) q[3];
sx q[3];
rz(-1.3397237) q[3];
sx q[3];
rz(-0.26255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498891) q[0];
sx q[0];
rz(-1.5174958) q[0];
sx q[0];
rz(-0.24895915) q[0];
rz(1.2173563) q[1];
sx q[1];
rz(-1.8739506) q[1];
sx q[1];
rz(-0.41044661) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3623724) q[0];
sx q[0];
rz(-1.1817314) q[0];
sx q[0];
rz(-0.44045191) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4538953) q[2];
sx q[2];
rz(-3.0499978) q[2];
sx q[2];
rz(2.8815038) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.16479334) q[1];
sx q[1];
rz(-1.8135957) q[1];
sx q[1];
rz(-1.1006646) q[1];
rz(-0.66907042) q[3];
sx q[3];
rz(-0.72815547) q[3];
sx q[3];
rz(1.6403584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88811389) q[2];
sx q[2];
rz(-1.9733182) q[2];
sx q[2];
rz(-1.1178364) q[2];
rz(0.11876336) q[3];
sx q[3];
rz(-2.4517086) q[3];
sx q[3];
rz(-1.9411055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73040199) q[0];
sx q[0];
rz(-1.406519) q[0];
sx q[0];
rz(2.7997959) q[0];
rz(3.0939915) q[1];
sx q[1];
rz(-1.8879415) q[1];
sx q[1];
rz(-1.3157848) q[1];
rz(-0.81472266) q[2];
sx q[2];
rz(-1.2263032) q[2];
sx q[2];
rz(-2.6487614) q[2];
rz(2.1988395) q[3];
sx q[3];
rz(-2.2136064) q[3];
sx q[3];
rz(2.8436713) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
