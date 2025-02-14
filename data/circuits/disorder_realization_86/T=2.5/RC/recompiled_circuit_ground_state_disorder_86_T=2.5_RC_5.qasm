OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.13088432) q[0];
sx q[0];
rz(-0.67357981) q[0];
sx q[0];
rz(2.3186865) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(-0.6757285) q[1];
sx q[1];
rz(-2.9339209) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0820706) q[0];
sx q[0];
rz(-0.86468177) q[0];
sx q[0];
rz(0.45990277) q[0];
rz(-pi) q[1];
rz(1.0124341) q[2];
sx q[2];
rz(-1.2744546) q[2];
sx q[2];
rz(1.3989965) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.339141) q[1];
sx q[1];
rz(-1.4467518) q[1];
sx q[1];
rz(-1.5088827) q[1];
rz(1.200233) q[3];
sx q[3];
rz(-1.5041122) q[3];
sx q[3];
rz(-2.5197864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1755918) q[2];
sx q[2];
rz(-1.5660183) q[2];
sx q[2];
rz(-1.2292181) q[2];
rz(1.2980596) q[3];
sx q[3];
rz(-1.7652054) q[3];
sx q[3];
rz(1.957533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9640279) q[0];
sx q[0];
rz(-0.25968817) q[0];
sx q[0];
rz(-2.2143256) q[0];
rz(0.19993965) q[1];
sx q[1];
rz(-1.6688469) q[1];
sx q[1];
rz(-0.98958579) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35868506) q[0];
sx q[0];
rz(-1.2099504) q[0];
sx q[0];
rz(3.0263958) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7071955) q[2];
sx q[2];
rz(-1.6100307) q[2];
sx q[2];
rz(-2.5920282) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.687508) q[1];
sx q[1];
rz(-1.6191779) q[1];
sx q[1];
rz(-1.6475186) q[1];
x q[2];
rz(2.5603676) q[3];
sx q[3];
rz(-0.94454256) q[3];
sx q[3];
rz(-3.0581829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9280615) q[2];
sx q[2];
rz(-1.4175043) q[2];
sx q[2];
rz(-0.35378635) q[2];
rz(-0.95156041) q[3];
sx q[3];
rz(-2.3944201) q[3];
sx q[3];
rz(-2.814754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3554409) q[0];
sx q[0];
rz(-2.1936301) q[0];
sx q[0];
rz(2.1308664) q[0];
rz(-3.082869) q[1];
sx q[1];
rz(-1.5208236) q[1];
sx q[1];
rz(2.7322863) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362013) q[0];
sx q[0];
rz(-0.86544468) q[0];
sx q[0];
rz(2.5558997) q[0];
rz(-2.9761613) q[2];
sx q[2];
rz(-1.9768052) q[2];
sx q[2];
rz(2.5819786) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.066591) q[1];
sx q[1];
rz(-1.7670872) q[1];
sx q[1];
rz(2.58954) q[1];
rz(0.037267123) q[3];
sx q[3];
rz(-2.6323595) q[3];
sx q[3];
rz(-1.6515712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0851486) q[2];
sx q[2];
rz(-1.8935545) q[2];
sx q[2];
rz(-0.15963456) q[2];
rz(-1.8917482) q[3];
sx q[3];
rz(-1.7227453) q[3];
sx q[3];
rz(-2.0554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98642629) q[0];
sx q[0];
rz(-0.43743375) q[0];
sx q[0];
rz(-1.0945818) q[0];
rz(1.9704341) q[1];
sx q[1];
rz(-1.1181701) q[1];
sx q[1];
rz(1.0135244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5399234) q[0];
sx q[0];
rz(-3.0832096) q[0];
sx q[0];
rz(2.4235382) q[0];
x q[1];
rz(1.5775575) q[2];
sx q[2];
rz(-0.49443118) q[2];
sx q[2];
rz(-1.0552366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3588874) q[1];
sx q[1];
rz(-0.26391477) q[1];
sx q[1];
rz(2.9857319) q[1];
x q[2];
rz(-1.3427686) q[3];
sx q[3];
rz(-2.6963391) q[3];
sx q[3];
rz(-0.54603117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1363498) q[2];
sx q[2];
rz(-1.3815657) q[2];
sx q[2];
rz(-1.8278149) q[2];
rz(-0.93787307) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(-0.098793678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2336642) q[0];
sx q[0];
rz(-1.6133244) q[0];
sx q[0];
rz(-1.044957) q[0];
rz(-0.16432556) q[1];
sx q[1];
rz(-1.798809) q[1];
sx q[1];
rz(2.9716861) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5066077) q[0];
sx q[0];
rz(-1.7568551) q[0];
sx q[0];
rz(-0.072441262) q[0];
x q[1];
rz(-1.2163148) q[2];
sx q[2];
rz(-1.2574242) q[2];
sx q[2];
rz(-2.2872668) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2881238) q[1];
sx q[1];
rz(-0.92902029) q[1];
sx q[1];
rz(-0.30499129) q[1];
x q[2];
rz(1.4395797) q[3];
sx q[3];
rz(-2.2926712) q[3];
sx q[3];
rz(-3.1027681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5728411) q[2];
sx q[2];
rz(-2.4390287) q[2];
sx q[2];
rz(-2.11002) q[2];
rz(-1.0860363) q[3];
sx q[3];
rz(-0.7998172) q[3];
sx q[3];
rz(-1.4097376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.80424911) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(1.7403437) q[0];
rz(2.7119472) q[1];
sx q[1];
rz(-0.88637561) q[1];
sx q[1];
rz(1.8362129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84112924) q[0];
sx q[0];
rz(-0.80781898) q[0];
sx q[0];
rz(-0.43231583) q[0];
x q[1];
rz(-3.1370509) q[2];
sx q[2];
rz(-2.0097187) q[2];
sx q[2];
rz(-2.2786841) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3682121) q[1];
sx q[1];
rz(-1.5399333) q[1];
sx q[1];
rz(-2.8931151) q[1];
x q[2];
rz(2.1277818) q[3];
sx q[3];
rz(-0.28388043) q[3];
sx q[3];
rz(2.551183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0171011) q[2];
sx q[2];
rz(-1.93511) q[2];
sx q[2];
rz(0.98480946) q[2];
rz(1.5752327) q[3];
sx q[3];
rz(-1.5519578) q[3];
sx q[3];
rz(-0.28520939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8816836) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(-1.786422) q[0];
rz(-0.024070865) q[1];
sx q[1];
rz(-1.6203974) q[1];
sx q[1];
rz(2.7640061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474038) q[0];
sx q[0];
rz(-2.4703175) q[0];
sx q[0];
rz(-2.8888974) q[0];
rz(-pi) q[1];
rz(-0.17654769) q[2];
sx q[2];
rz(-1.7804885) q[2];
sx q[2];
rz(-3.049946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8319893) q[1];
sx q[1];
rz(-1.9506694) q[1];
sx q[1];
rz(1.989945) q[1];
x q[2];
rz(-2.2144775) q[3];
sx q[3];
rz(-2.6171631) q[3];
sx q[3];
rz(2.7853109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.54921237) q[2];
sx q[2];
rz(-2.8040631) q[2];
sx q[2];
rz(-0.40360061) q[2];
rz(0.18925439) q[3];
sx q[3];
rz(-2.0892102) q[3];
sx q[3];
rz(-0.45690593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.6117578) q[0];
sx q[0];
rz(-2.5921322) q[0];
sx q[0];
rz(-0.31103617) q[0];
rz(-0.95651904) q[1];
sx q[1];
rz(-1.6638959) q[1];
sx q[1];
rz(-0.32435736) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3295591) q[0];
sx q[0];
rz(-0.5428462) q[0];
sx q[0];
rz(2.7663648) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5115405) q[2];
sx q[2];
rz(-2.0803335) q[2];
sx q[2];
rz(1.6409724) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2203103) q[1];
sx q[1];
rz(-1.1571572) q[1];
sx q[1];
rz(2.2181554) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1127541) q[3];
sx q[3];
rz(-0.73087091) q[3];
sx q[3];
rz(-2.3422086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.73584622) q[2];
sx q[2];
rz(-0.79068557) q[2];
sx q[2];
rz(-0.22845593) q[2];
rz(-2.6089846) q[3];
sx q[3];
rz(-0.76627982) q[3];
sx q[3];
rz(-2.2920091) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.547895) q[0];
sx q[0];
rz(-0.51012796) q[0];
sx q[0];
rz(-1.8713895) q[0];
rz(0.58865976) q[1];
sx q[1];
rz(-1.207186) q[1];
sx q[1];
rz(0.38280815) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13148334) q[0];
sx q[0];
rz(-1.3416222) q[0];
sx q[0];
rz(3.1160627) q[0];
rz(-2.7076376) q[2];
sx q[2];
rz(-1.7979017) q[2];
sx q[2];
rz(-1.3004829) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8296903) q[1];
sx q[1];
rz(-1.6964579) q[1];
sx q[1];
rz(-0.71255334) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.270003) q[3];
sx q[3];
rz(-2.0473891) q[3];
sx q[3];
rz(1.4780413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1560912) q[2];
sx q[2];
rz(-1.568855) q[2];
sx q[2];
rz(0.4001948) q[2];
rz(0.24724809) q[3];
sx q[3];
rz(-1.6797545) q[3];
sx q[3];
rz(-2.7483773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8409214) q[0];
sx q[0];
rz(-1.8074169) q[0];
sx q[0];
rz(-2.1260496) q[0];
rz(2.4726942) q[1];
sx q[1];
rz(-1.6872419) q[1];
sx q[1];
rz(-1.5060172) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6164613) q[0];
sx q[0];
rz(-0.86989738) q[0];
sx q[0];
rz(-2.4674795) q[0];
x q[1];
rz(2.0066543) q[2];
sx q[2];
rz(-2.0761937) q[2];
sx q[2];
rz(-1.7541898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2336894) q[1];
sx q[1];
rz(-2.0504867) q[1];
sx q[1];
rz(3.0976899) q[1];
x q[2];
rz(-2.2478836) q[3];
sx q[3];
rz(-2.0314616) q[3];
sx q[3];
rz(-3.0627444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16613913) q[2];
sx q[2];
rz(-1.0756805) q[2];
sx q[2];
rz(1.6157185) q[2];
rz(0.29159355) q[3];
sx q[3];
rz(-2.6988131) q[3];
sx q[3];
rz(2.7899138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46353729) q[0];
sx q[0];
rz(-2.1778477) q[0];
sx q[0];
rz(0.5974593) q[0];
rz(0.28868227) q[1];
sx q[1];
rz(-0.79816993) q[1];
sx q[1];
rz(0.023963902) q[1];
rz(0.63715061) q[2];
sx q[2];
rz(-2.646614) q[2];
sx q[2];
rz(2.6029233) q[2];
rz(-1.8221832) q[3];
sx q[3];
rz(-2.0181927) q[3];
sx q[3];
rz(-0.041139091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
