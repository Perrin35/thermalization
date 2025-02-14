OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6588563) q[0];
sx q[0];
rz(-0.98134494) q[0];
sx q[0];
rz(1.097783) q[0];
rz(2.3535347) q[1];
sx q[1];
rz(-1.0507974) q[1];
sx q[1];
rz(0.0069590574) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6774782) q[0];
sx q[0];
rz(-0.65271806) q[0];
sx q[0];
rz(-1.4783036) q[0];
rz(-pi) q[1];
x q[1];
rz(1.263515) q[2];
sx q[2];
rz(-0.6699962) q[2];
sx q[2];
rz(-0.43596632) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6485868) q[1];
sx q[1];
rz(-1.3351591) q[1];
sx q[1];
rz(-0.034502397) q[1];
rz(-pi) q[2];
rz(-0.29421774) q[3];
sx q[3];
rz(-0.3539043) q[3];
sx q[3];
rz(-0.82415165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2385345) q[2];
sx q[2];
rz(-1.3471341) q[2];
sx q[2];
rz(1.6713589) q[2];
rz(-3.0155731) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(-2.457705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11494342) q[0];
sx q[0];
rz(-2.6925955) q[0];
sx q[0];
rz(2.5057416) q[0];
rz(-0.77330971) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(-1.6962475) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.792004) q[0];
sx q[0];
rz(-0.82539648) q[0];
sx q[0];
rz(2.1950153) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47451203) q[2];
sx q[2];
rz(-0.84174985) q[2];
sx q[2];
rz(3.0345033) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0750097) q[1];
sx q[1];
rz(-0.72117469) q[1];
sx q[1];
rz(-0.23951342) q[1];
rz(-pi) q[2];
rz(-1.7023193) q[3];
sx q[3];
rz(-1.6386541) q[3];
sx q[3];
rz(2.0236286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14629743) q[2];
sx q[2];
rz(-1.7701021) q[2];
sx q[2];
rz(-0.13776097) q[2];
rz(2.6855101) q[3];
sx q[3];
rz(-0.59499732) q[3];
sx q[3];
rz(2.6661787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.544203) q[0];
sx q[0];
rz(-1.5837639) q[0];
sx q[0];
rz(0.57918817) q[0];
rz(3.0384565) q[1];
sx q[1];
rz(-1.3214279) q[1];
sx q[1];
rz(2.3613222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0294757) q[0];
sx q[0];
rz(-1.5519841) q[0];
sx q[0];
rz(0.095186724) q[0];
rz(0.68124007) q[2];
sx q[2];
rz(-1.6954172) q[2];
sx q[2];
rz(1.7523426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5835598) q[1];
sx q[1];
rz(-1.7906902) q[1];
sx q[1];
rz(-1.1660485) q[1];
rz(0.87522755) q[3];
sx q[3];
rz(-1.7415294) q[3];
sx q[3];
rz(3.1306992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66037336) q[2];
sx q[2];
rz(-1.4619091) q[2];
sx q[2];
rz(0.090864651) q[2];
rz(1.6615435) q[3];
sx q[3];
rz(-1.3250947) q[3];
sx q[3];
rz(-0.81502325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.018983) q[0];
sx q[0];
rz(-0.18773395) q[0];
sx q[0];
rz(-2.6222099) q[0];
rz(-1.9620365) q[1];
sx q[1];
rz(-0.68125454) q[1];
sx q[1];
rz(-3.093241) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40019401) q[0];
sx q[0];
rz(-2.1693008) q[0];
sx q[0];
rz(0.6299751) q[0];
rz(-2.7555356) q[2];
sx q[2];
rz(-1.1163082) q[2];
sx q[2];
rz(0.41815652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6801694) q[1];
sx q[1];
rz(-1.5986414) q[1];
sx q[1];
rz(-0.48770406) q[1];
x q[2];
rz(0.79508852) q[3];
sx q[3];
rz(-1.3833191) q[3];
sx q[3];
rz(0.80611967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68253303) q[2];
sx q[2];
rz(-0.30632633) q[2];
sx q[2];
rz(-1.691386) q[2];
rz(1.1059443) q[3];
sx q[3];
rz(-1.9927497) q[3];
sx q[3];
rz(-0.86627427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3806216) q[0];
sx q[0];
rz(-0.80046099) q[0];
sx q[0];
rz(0.77734429) q[0];
rz(2.6457973) q[1];
sx q[1];
rz(-2.7340041) q[1];
sx q[1];
rz(0.19283238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58932006) q[0];
sx q[0];
rz(-2.5077132) q[0];
sx q[0];
rz(-0.53323563) q[0];
x q[1];
rz(0.76754153) q[2];
sx q[2];
rz(-1.7450252) q[2];
sx q[2];
rz(1.7310639) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65102386) q[1];
sx q[1];
rz(-1.2747163) q[1];
sx q[1];
rz(-0.091551642) q[1];
rz(-pi) q[2];
rz(1.1689831) q[3];
sx q[3];
rz(-1.7208368) q[3];
sx q[3];
rz(-0.91317859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5567646) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(-0.83621109) q[2];
rz(-2.1861475) q[3];
sx q[3];
rz(-1.4232057) q[3];
sx q[3];
rz(-2.6098765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29109508) q[0];
sx q[0];
rz(-1.2245155) q[0];
sx q[0];
rz(-3.1078872) q[0];
rz(2.2233502) q[1];
sx q[1];
rz(-0.70437175) q[1];
sx q[1];
rz(-0.69923002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13051899) q[0];
sx q[0];
rz(-1.0193745) q[0];
sx q[0];
rz(0.60451492) q[0];
x q[1];
rz(2.9543058) q[2];
sx q[2];
rz(-1.5178871) q[2];
sx q[2];
rz(-2.1250181) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.81046178) q[1];
sx q[1];
rz(-1.394636) q[1];
sx q[1];
rz(-0.44682746) q[1];
rz(-pi) q[2];
rz(-0.27856234) q[3];
sx q[3];
rz(-1.8980366) q[3];
sx q[3];
rz(2.183941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.08192) q[2];
sx q[2];
rz(-0.6531738) q[2];
sx q[2];
rz(3.0461404) q[2];
rz(1.0517612) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(2.2473647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28602257) q[0];
sx q[0];
rz(-2.6595071) q[0];
sx q[0];
rz(1.1676189) q[0];
rz(2.7883912) q[1];
sx q[1];
rz(-0.51274061) q[1];
sx q[1];
rz(2.8599427) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324556) q[0];
sx q[0];
rz(-1.0971709) q[0];
sx q[0];
rz(-0.098129674) q[0];
x q[1];
rz(0.64729624) q[2];
sx q[2];
rz(-1.4855097) q[2];
sx q[2];
rz(-0.46009395) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.901501) q[1];
sx q[1];
rz(-0.20040837) q[1];
sx q[1];
rz(-0.83205454) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8696712) q[3];
sx q[3];
rz(-2.4574349) q[3];
sx q[3];
rz(1.1604094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99792751) q[2];
sx q[2];
rz(-1.2219656) q[2];
sx q[2];
rz(1.8188875) q[2];
rz(-1.5023242) q[3];
sx q[3];
rz(-2.0570677) q[3];
sx q[3];
rz(-0.098202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.99763501) q[0];
sx q[0];
rz(-1.8956381) q[0];
sx q[0];
rz(-0.27398807) q[0];
rz(0.59448376) q[1];
sx q[1];
rz(-1.4130054) q[1];
sx q[1];
rz(0.53057539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45121058) q[0];
sx q[0];
rz(-1.57461) q[0];
sx q[0];
rz(-0.03937748) q[0];
x q[1];
rz(-1.6470026) q[2];
sx q[2];
rz(-2.1222024) q[2];
sx q[2];
rz(1.4879256) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1889607) q[1];
sx q[1];
rz(-0.41302339) q[1];
sx q[1];
rz(2.7160591) q[1];
x q[2];
rz(2.2835571) q[3];
sx q[3];
rz(-1.4467561) q[3];
sx q[3];
rz(0.8197561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24255594) q[2];
sx q[2];
rz(-1.2071004) q[2];
sx q[2];
rz(-2.8846018) q[2];
rz(1.1993923) q[3];
sx q[3];
rz(-1.896984) q[3];
sx q[3];
rz(1.5132343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5892107) q[0];
sx q[0];
rz(-0.95443812) q[0];
sx q[0];
rz(2.6115665) q[0];
rz(1.5026622) q[1];
sx q[1];
rz(-1.1549779) q[1];
sx q[1];
rz(0.93856215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5230704) q[0];
sx q[0];
rz(-1.1022727) q[0];
sx q[0];
rz(2.4836088) q[0];
x q[1];
rz(2.4468719) q[2];
sx q[2];
rz(-0.18604569) q[2];
sx q[2];
rz(-3.0228066) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0238513) q[1];
sx q[1];
rz(-1.3582025) q[1];
sx q[1];
rz(-2.0076058) q[1];
rz(-1.7750562) q[3];
sx q[3];
rz(-0.59511772) q[3];
sx q[3];
rz(-0.68347733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75371257) q[2];
sx q[2];
rz(-0.91999274) q[2];
sx q[2];
rz(-2.06125) q[2];
rz(2.3267817) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(1.4174392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1431047) q[0];
sx q[0];
rz(-1.4838706) q[0];
sx q[0];
rz(1.0888354) q[0];
rz(-3.0740652) q[1];
sx q[1];
rz(-1.8636401) q[1];
sx q[1];
rz(-2.3823104) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35847202) q[0];
sx q[0];
rz(-1.6896588) q[0];
sx q[0];
rz(-1.3141743) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6443321) q[2];
sx q[2];
rz(-0.25575519) q[2];
sx q[2];
rz(1.7921247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49803653) q[1];
sx q[1];
rz(-0.78460556) q[1];
sx q[1];
rz(-0.30679171) q[1];
rz(1.3410232) q[3];
sx q[3];
rz(-1.9073351) q[3];
sx q[3];
rz(2.0391109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8074983) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(-2.5661772) q[2];
rz(-0.56257644) q[3];
sx q[3];
rz(-2.176599) q[3];
sx q[3];
rz(1.3728728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060487735) q[0];
sx q[0];
rz(-1.2860379) q[0];
sx q[0];
rz(2.2338569) q[0];
rz(-2.0545215) q[1];
sx q[1];
rz(-1.1320976) q[1];
sx q[1];
rz(-0.017398106) q[1];
rz(-1.7123863) q[2];
sx q[2];
rz(-1.7238486) q[2];
sx q[2];
rz(0.42949745) q[2];
rz(1.1997052) q[3];
sx q[3];
rz(-2.1368847) q[3];
sx q[3];
rz(0.18659244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
