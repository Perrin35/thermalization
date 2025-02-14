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
rz(1.0951618) q[0];
sx q[0];
rz(-2.996063) q[0];
sx q[0];
rz(-0.50330436) q[0];
rz(-1.9637928) q[1];
sx q[1];
rz(4.7363321) q[1];
sx q[1];
rz(7.9793032) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58848042) q[0];
sx q[0];
rz(-2.7403826) q[0];
sx q[0];
rz(0.11102872) q[0];
rz(-pi) q[1];
rz(-2.3906227) q[2];
sx q[2];
rz(-1.079796) q[2];
sx q[2];
rz(-3.1042539) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7912476) q[1];
sx q[1];
rz(-2.4295632) q[1];
sx q[1];
rz(-2.6256034) q[1];
rz(-1.5863933) q[3];
sx q[3];
rz(-1.6199779) q[3];
sx q[3];
rz(1.552207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6426927) q[2];
sx q[2];
rz(-1.91012) q[2];
sx q[2];
rz(1.6373681) q[2];
rz(-2.4046992) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(-0.98114291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.40642834) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(-0.51097393) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(-0.32646349) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24169479) q[0];
sx q[0];
rz(-0.37942552) q[0];
sx q[0];
rz(2.3510758) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88739354) q[2];
sx q[2];
rz(-1.7345725) q[2];
sx q[2];
rz(2.0205409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0357246) q[1];
sx q[1];
rz(-0.99109036) q[1];
sx q[1];
rz(0.3649664) q[1];
x q[2];
rz(-2.191698) q[3];
sx q[3];
rz(-2.6767459) q[3];
sx q[3];
rz(0.45528938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9713126) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(1.8453321) q[2];
rz(-2.7235459) q[3];
sx q[3];
rz(-2.1635735) q[3];
sx q[3];
rz(-0.70438284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.068785) q[0];
sx q[0];
rz(-2.7588221) q[0];
sx q[0];
rz(-0.54247722) q[0];
rz(1.3756649) q[1];
sx q[1];
rz(-0.41789564) q[1];
sx q[1];
rz(3.1105522) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.870793) q[0];
sx q[0];
rz(-2.7794771) q[0];
sx q[0];
rz(0.55678456) q[0];
rz(-pi) q[1];
rz(0.91340567) q[2];
sx q[2];
rz(-1.7511466) q[2];
sx q[2];
rz(3.0424398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87796383) q[1];
sx q[1];
rz(-2.1641556) q[1];
sx q[1];
rz(0.29747648) q[1];
x q[2];
rz(-1.5543082) q[3];
sx q[3];
rz(-2.1257536) q[3];
sx q[3];
rz(2.6944427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5151908) q[2];
sx q[2];
rz(-2.4050737) q[2];
sx q[2];
rz(0.82702965) q[2];
rz(-2.6229897) q[3];
sx q[3];
rz(-1.1122455) q[3];
sx q[3];
rz(-1.5145068) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9943635) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(-1.9299141) q[0];
rz(-1.0481102) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(1.3406219) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3163534) q[0];
sx q[0];
rz(-1.8779426) q[0];
sx q[0];
rz(-2.4841043) q[0];
x q[1];
rz(0.2366613) q[2];
sx q[2];
rz(-1.4259031) q[2];
sx q[2];
rz(-1.3587111) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27599469) q[1];
sx q[1];
rz(-1.0038687) q[1];
sx q[1];
rz(1.2204942) q[1];
rz(-pi) q[2];
rz(0.31130917) q[3];
sx q[3];
rz(-1.9640067) q[3];
sx q[3];
rz(-1.9462727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2010605) q[2];
sx q[2];
rz(-2.9643855) q[2];
sx q[2];
rz(-2.000467) q[2];
rz(-1.6107669) q[3];
sx q[3];
rz(-2.174236) q[3];
sx q[3];
rz(-2.358986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56213266) q[0];
sx q[0];
rz(-1.1701522) q[0];
sx q[0];
rz(-0.60254565) q[0];
rz(0.58019477) q[1];
sx q[1];
rz(-1.5126901) q[1];
sx q[1];
rz(0.7811195) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9255722) q[0];
sx q[0];
rz(-1.4013774) q[0];
sx q[0];
rz(1.2792619) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76483043) q[2];
sx q[2];
rz(-1.0862319) q[2];
sx q[2];
rz(-2.7861905) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.14706068) q[1];
sx q[1];
rz(-2.414738) q[1];
sx q[1];
rz(-1.7009853) q[1];
rz(-pi) q[2];
rz(1.7479728) q[3];
sx q[3];
rz(-0.85292888) q[3];
sx q[3];
rz(-2.126463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.96131229) q[2];
sx q[2];
rz(-2.0013516) q[2];
sx q[2];
rz(-0.22892924) q[2];
rz(-1.8217575) q[3];
sx q[3];
rz(-1.4819375) q[3];
sx q[3];
rz(0.84407097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22885403) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(-1.9629021) q[0];
rz(-1.9121869) q[1];
sx q[1];
rz(-2.7472159) q[1];
sx q[1];
rz(-1.2961402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73979171) q[0];
sx q[0];
rz(-1.8515046) q[0];
sx q[0];
rz(-0.14174353) q[0];
rz(-0.8528233) q[2];
sx q[2];
rz(-2.1666424) q[2];
sx q[2];
rz(-2.9862474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.072159616) q[1];
sx q[1];
rz(-2.0516703) q[1];
sx q[1];
rz(2.062501) q[1];
rz(-pi) q[2];
rz(-0.18978203) q[3];
sx q[3];
rz(-0.49473195) q[3];
sx q[3];
rz(-0.035898048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65427762) q[2];
sx q[2];
rz(-0.9026022) q[2];
sx q[2];
rz(-0.26228341) q[2];
rz(-0.30019635) q[3];
sx q[3];
rz(-1.9191977) q[3];
sx q[3];
rz(2.9330971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7009785) q[0];
sx q[0];
rz(-2.746026) q[0];
sx q[0];
rz(-1.3025008) q[0];
rz(2.473096) q[1];
sx q[1];
rz(-1.3553456) q[1];
sx q[1];
rz(-1.0350593) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67643205) q[0];
sx q[0];
rz(-2.6817828) q[0];
sx q[0];
rz(-2.0681429) q[0];
rz(1.1167489) q[2];
sx q[2];
rz(-0.40698689) q[2];
sx q[2];
rz(1.5429516) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8866402) q[1];
sx q[1];
rz(-1.4952085) q[1];
sx q[1];
rz(-1.2264538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97888246) q[3];
sx q[3];
rz(-0.91875263) q[3];
sx q[3];
rz(-1.0197717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0237026) q[2];
sx q[2];
rz(-1.002545) q[2];
sx q[2];
rz(-1.6560076) q[2];
rz(-0.28247908) q[3];
sx q[3];
rz(-1.7829203) q[3];
sx q[3];
rz(-2.7070467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3561803) q[0];
sx q[0];
rz(-2.0996576) q[0];
sx q[0];
rz(0.27454141) q[0];
rz(-1.2046332) q[1];
sx q[1];
rz(-0.58012539) q[1];
sx q[1];
rz(-2.5801632) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0208541) q[0];
sx q[0];
rz(-2.0971813) q[0];
sx q[0];
rz(2.3096183) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7700069) q[2];
sx q[2];
rz(-0.71096651) q[2];
sx q[2];
rz(1.154605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2632713) q[1];
sx q[1];
rz(-2.6215963) q[1];
sx q[1];
rz(-2.2161198) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8772914) q[3];
sx q[3];
rz(-1.9076548) q[3];
sx q[3];
rz(-1.565153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56208912) q[2];
sx q[2];
rz(-2.0359437) q[2];
sx q[2];
rz(-1.7378463) q[2];
rz(1.7956519) q[3];
sx q[3];
rz(-2.3965049) q[3];
sx q[3];
rz(0.48817202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76029921) q[0];
sx q[0];
rz(-0.60307044) q[0];
sx q[0];
rz(2.2316933) q[0];
rz(1.9445885) q[1];
sx q[1];
rz(-2.0884114) q[1];
sx q[1];
rz(-0.88814703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5407046) q[0];
sx q[0];
rz(-1.6645169) q[0];
sx q[0];
rz(1.4613495) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4817564) q[2];
sx q[2];
rz(-1.6183231) q[2];
sx q[2];
rz(-0.30872503) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55938086) q[1];
sx q[1];
rz(-2.191698) q[1];
sx q[1];
rz(1.0889174) q[1];
rz(1.2751332) q[3];
sx q[3];
rz(-1.5348892) q[3];
sx q[3];
rz(1.1334238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1278648) q[2];
sx q[2];
rz(-1.700054) q[2];
sx q[2];
rz(-1.0825276) q[2];
rz(0.23076375) q[3];
sx q[3];
rz(-1.328238) q[3];
sx q[3];
rz(-0.48399353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312663) q[0];
sx q[0];
rz(-2.3895097) q[0];
sx q[0];
rz(-0.58151522) q[0];
rz(-0.61559081) q[1];
sx q[1];
rz(-0.94771996) q[1];
sx q[1];
rz(-1.8448578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6433577) q[0];
sx q[0];
rz(-1.5603934) q[0];
sx q[0];
rz(1.5889945) q[0];
rz(-2.7974772) q[2];
sx q[2];
rz(-1.5617579) q[2];
sx q[2];
rz(-1.8923541) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3677926) q[1];
sx q[1];
rz(-1.3145216) q[1];
sx q[1];
rz(-2.0092177) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9711791) q[3];
sx q[3];
rz(-2.2030911) q[3];
sx q[3];
rz(2.1454772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7221308) q[2];
sx q[2];
rz(-3.0007134) q[2];
sx q[2];
rz(-0.44931832) q[2];
rz(2.8214473) q[3];
sx q[3];
rz(-1.3195427) q[3];
sx q[3];
rz(1.2464574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.37353361) q[0];
sx q[0];
rz(-2.4644485) q[0];
sx q[0];
rz(1.2039536) q[0];
rz(-1.3846579) q[1];
sx q[1];
rz(-2.90381) q[1];
sx q[1];
rz(-0.79868383) q[1];
rz(2.8821386) q[2];
sx q[2];
rz(-0.96895192) q[2];
sx q[2];
rz(-0.12018991) q[2];
rz(-2.2438335) q[3];
sx q[3];
rz(-2.5644292) q[3];
sx q[3];
rz(0.226365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
