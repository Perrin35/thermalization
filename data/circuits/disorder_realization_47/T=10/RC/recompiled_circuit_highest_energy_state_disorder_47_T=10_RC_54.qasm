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
rz(-2.0464309) q[0];
sx q[0];
rz(-0.14552966) q[0];
sx q[0];
rz(-2.6382883) q[0];
rz(-1.9637928) q[1];
sx q[1];
rz(-1.5468532) q[1];
sx q[1];
rz(-1.4454747) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58848042) q[0];
sx q[0];
rz(-2.7403826) q[0];
sx q[0];
rz(-0.11102872) q[0];
x q[1];
rz(0.75096994) q[2];
sx q[2];
rz(-1.079796) q[2];
sx q[2];
rz(0.037338749) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7912476) q[1];
sx q[1];
rz(-2.4295632) q[1];
sx q[1];
rz(-0.51598926) q[1];
rz(1.5551994) q[3];
sx q[3];
rz(-1.6199779) q[3];
sx q[3];
rz(1.552207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6426927) q[2];
sx q[2];
rz(-1.91012) q[2];
sx q[2];
rz(-1.5042245) q[2];
rz(0.73689342) q[3];
sx q[3];
rz(-1.5259909) q[3];
sx q[3];
rz(-2.1604497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7351643) q[0];
sx q[0];
rz(-2.6728215) q[0];
sx q[0];
rz(2.6306187) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(-0.32646349) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58566203) q[0];
sx q[0];
rz(-1.8343975) q[0];
sx q[0];
rz(-1.8469514) q[0];
rz(-0.88739354) q[2];
sx q[2];
rz(-1.7345725) q[2];
sx q[2];
rz(-1.1210517) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.105868) q[1];
sx q[1];
rz(-2.1505023) q[1];
sx q[1];
rz(0.3649664) q[1];
rz(-pi) q[2];
rz(-2.8577153) q[3];
sx q[3];
rz(-1.9440158) q[3];
sx q[3];
rz(2.9220327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9713126) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(-1.8453321) q[2];
rz(0.41804677) q[3];
sx q[3];
rz(-2.1635735) q[3];
sx q[3];
rz(2.4372098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.068785) q[0];
sx q[0];
rz(-0.3827706) q[0];
sx q[0];
rz(-0.54247722) q[0];
rz(-1.3756649) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(-0.03104041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.870793) q[0];
sx q[0];
rz(-0.36211553) q[0];
sx q[0];
rz(2.5848081) q[0];
rz(0.91340567) q[2];
sx q[2];
rz(-1.3904461) q[2];
sx q[2];
rz(-3.0424398) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.52306226) q[1];
sx q[1];
rz(-1.8162604) q[1];
sx q[1];
rz(2.1851468) q[1];
rz(-pi) q[2];
rz(-3.1150041) q[3];
sx q[3];
rz(-2.5864161) q[3];
sx q[3];
rz(0.41586499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.62640181) q[2];
sx q[2];
rz(-0.73651892) q[2];
sx q[2];
rz(0.82702965) q[2];
rz(2.6229897) q[3];
sx q[3];
rz(-2.0293472) q[3];
sx q[3];
rz(-1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9943635) q[0];
sx q[0];
rz(-1.3059068) q[0];
sx q[0];
rz(1.9299141) q[0];
rz(1.0481102) q[1];
sx q[1];
rz(-0.71980372) q[1];
sx q[1];
rz(1.3406219) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82523924) q[0];
sx q[0];
rz(-1.2636501) q[0];
sx q[0];
rz(-2.4841043) q[0];
rz(-0.55669341) q[2];
sx q[2];
rz(-0.27678267) q[2];
sx q[2];
rz(2.3900696) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27599469) q[1];
sx q[1];
rz(-1.0038687) q[1];
sx q[1];
rz(-1.9210985) q[1];
rz(2.8302835) q[3];
sx q[3];
rz(-1.177586) q[3];
sx q[3];
rz(1.1953199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94053215) q[2];
sx q[2];
rz(-0.17720711) q[2];
sx q[2];
rz(-2.000467) q[2];
rz(-1.6107669) q[3];
sx q[3];
rz(-0.96735668) q[3];
sx q[3];
rz(-0.78260666) q[3];
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
rz(0.56213266) q[0];
sx q[0];
rz(-1.1701522) q[0];
sx q[0];
rz(0.60254565) q[0];
rz(-2.5613979) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(-0.7811195) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8373708) q[0];
sx q[0];
rz(-1.2835555) q[0];
sx q[0];
rz(0.17673136) q[0];
rz(-pi) q[1];
rz(2.2011312) q[2];
sx q[2];
rz(-2.2302719) q[2];
sx q[2];
rz(-2.3465921) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32051099) q[1];
sx q[1];
rz(-0.85144224) q[1];
sx q[1];
rz(-3.0266552) q[1];
x q[2];
rz(-2.9424465) q[3];
sx q[3];
rz(-0.73561397) q[3];
sx q[3];
rz(-0.74935952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1802804) q[2];
sx q[2];
rz(-1.140241) q[2];
sx q[2];
rz(-0.22892924) q[2];
rz(-1.3198352) q[3];
sx q[3];
rz(-1.4819375) q[3];
sx q[3];
rz(-0.84407097) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22885403) q[0];
sx q[0];
rz(-1.5541394) q[0];
sx q[0];
rz(1.9629021) q[0];
rz(-1.2294058) q[1];
sx q[1];
rz(-2.7472159) q[1];
sx q[1];
rz(-1.8454525) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3501007) q[0];
sx q[0];
rz(-1.7069567) q[0];
sx q[0];
rz(-1.2873935) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3713114) q[2];
sx q[2];
rz(-0.89777032) q[2];
sx q[2];
rz(1.1550624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78627045) q[1];
sx q[1];
rz(-2.4679524) q[1];
sx q[1];
rz(2.4060529) q[1];
rz(2.9518106) q[3];
sx q[3];
rz(-0.49473195) q[3];
sx q[3];
rz(-0.035898048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.487315) q[2];
sx q[2];
rz(-0.9026022) q[2];
sx q[2];
rz(-0.26228341) q[2];
rz(-0.30019635) q[3];
sx q[3];
rz(-1.2223949) q[3];
sx q[3];
rz(-2.9330971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.4406141) q[0];
sx q[0];
rz(-0.3955667) q[0];
sx q[0];
rz(1.3025008) q[0];
rz(2.473096) q[1];
sx q[1];
rz(-1.786247) q[1];
sx q[1];
rz(-2.1065333) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44159206) q[0];
sx q[0];
rz(-1.3574575) q[0];
sx q[0];
rz(1.9812897) q[0];
x q[1];
rz(1.1167489) q[2];
sx q[2];
rz(-0.40698689) q[2];
sx q[2];
rz(1.5429516) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0332843) q[1];
sx q[1];
rz(-0.35221812) q[1];
sx q[1];
rz(-1.3501106) q[1];
rz(-pi) q[2];
rz(2.5104654) q[3];
sx q[3];
rz(-0.85048095) q[3];
sx q[3];
rz(-2.9583954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0237026) q[2];
sx q[2];
rz(-2.1390476) q[2];
sx q[2];
rz(1.6560076) q[2];
rz(0.28247908) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(0.43454596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3561803) q[0];
sx q[0];
rz(-2.0996576) q[0];
sx q[0];
rz(-2.8670512) q[0];
rz(1.2046332) q[1];
sx q[1];
rz(-0.58012539) q[1];
sx q[1];
rz(2.5801632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12073853) q[0];
sx q[0];
rz(-1.0444114) q[0];
sx q[0];
rz(2.3096183) q[0];
rz(-pi) q[1];
rz(-2.7700069) q[2];
sx q[2];
rz(-0.71096651) q[2];
sx q[2];
rz(1.154605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8703979) q[1];
sx q[1];
rz(-1.8742838) q[1];
sx q[1];
rz(1.1417901) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92989489) q[3];
sx q[3];
rz(-0.42503438) q[3];
sx q[3];
rz(-0.89034058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5795035) q[2];
sx q[2];
rz(-2.0359437) q[2];
sx q[2];
rz(1.7378463) q[2];
rz(1.3459407) q[3];
sx q[3];
rz(-2.3965049) q[3];
sx q[3];
rz(2.6534206) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3812934) q[0];
sx q[0];
rz(-2.5385222) q[0];
sx q[0];
rz(-0.90989939) q[0];
rz(-1.1970041) q[1];
sx q[1];
rz(-1.0531813) q[1];
sx q[1];
rz(-2.2534456) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040374856) q[0];
sx q[0];
rz(-1.4618317) q[0];
sx q[0];
rz(0.094281406) q[0];
rz(-1.5106701) q[2];
sx q[2];
rz(-0.91183582) q[2];
sx q[2];
rz(1.9163641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3067813) q[1];
sx q[1];
rz(-1.1843006) q[1];
sx q[1];
rz(-2.4624834) q[1];
rz(-0.03753438) q[3];
sx q[3];
rz(-1.8662631) q[3];
sx q[3];
rz(2.7151544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0137279) q[2];
sx q[2];
rz(-1.700054) q[2];
sx q[2];
rz(-2.0590651) q[2];
rz(0.23076375) q[3];
sx q[3];
rz(-1.8133546) q[3];
sx q[3];
rz(-2.6575991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312663) q[0];
sx q[0];
rz(-2.3895097) q[0];
sx q[0];
rz(-0.58151522) q[0];
rz(0.61559081) q[1];
sx q[1];
rz(-0.94771996) q[1];
sx q[1];
rz(-1.2967348) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.588284) q[0];
sx q[0];
rz(-0.020961449) q[0];
sx q[0];
rz(1.0514489) q[0];
x q[1];
rz(3.1148071) q[2];
sx q[2];
rz(-2.7973632) q[2];
sx q[2];
rz(-0.29634288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8402214) q[1];
sx q[1];
rz(-2.6379963) q[1];
sx q[1];
rz(-1.0177769) q[1];
rz(-0.6720242) q[3];
sx q[3];
rz(-1.250953) q[3];
sx q[3];
rz(-2.8120188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7221308) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(-0.44931832) q[2];
rz(-0.32014534) q[3];
sx q[3];
rz(-1.82205) q[3];
sx q[3];
rz(-1.2464574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37353361) q[0];
sx q[0];
rz(-2.4644485) q[0];
sx q[0];
rz(1.2039536) q[0];
rz(1.3846579) q[1];
sx q[1];
rz(-0.23778267) q[1];
sx q[1];
rz(2.3429088) q[1];
rz(-1.9282707) q[2];
sx q[2];
rz(-2.4926179) q[2];
sx q[2];
rz(0.31821584) q[2];
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
