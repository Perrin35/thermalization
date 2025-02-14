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
rz(-0.43871969) q[0];
sx q[0];
rz(-2.5320142) q[0];
sx q[0];
rz(0.69019812) q[0];
rz(-1.8742427) q[1];
sx q[1];
rz(-2.6522418) q[1];
sx q[1];
rz(-2.7971921) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716838) q[0];
sx q[0];
rz(-2.8835886) q[0];
sx q[0];
rz(2.1841316) q[0];
rz(-pi) q[1];
rz(0.51235389) q[2];
sx q[2];
rz(-0.47253451) q[2];
sx q[2];
rz(-0.81250459) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2799543) q[1];
sx q[1];
rz(-0.32740232) q[1];
sx q[1];
rz(-0.91307098) q[1];
rz(-pi) q[2];
rz(2.9696483) q[3];
sx q[3];
rz(-2.056582) q[3];
sx q[3];
rz(0.43648374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.043618) q[2];
sx q[2];
rz(-2.2100885) q[2];
sx q[2];
rz(2.0841058) q[2];
rz(2.1438694) q[3];
sx q[3];
rz(-1.7505373) q[3];
sx q[3];
rz(-1.9690431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69720307) q[0];
sx q[0];
rz(-2.3219705) q[0];
sx q[0];
rz(-2.3890553) q[0];
rz(1.8684111) q[1];
sx q[1];
rz(-1.1538785) q[1];
sx q[1];
rz(0.85535991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46525644) q[0];
sx q[0];
rz(-1.4727061) q[0];
sx q[0];
rz(-1.3560182) q[0];
rz(-pi) q[1];
rz(2.8933002) q[2];
sx q[2];
rz(-0.96370164) q[2];
sx q[2];
rz(0.33252972) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8458057) q[1];
sx q[1];
rz(-1.9545022) q[1];
sx q[1];
rz(-1.3944666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5961038) q[3];
sx q[3];
rz(-1.3719333) q[3];
sx q[3];
rz(-2.1391099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7684043) q[2];
sx q[2];
rz(-1.6104001) q[2];
sx q[2];
rz(1.9739523) q[2];
rz(-0.097955616) q[3];
sx q[3];
rz(-2.8601213) q[3];
sx q[3];
rz(1.9907985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34404594) q[0];
sx q[0];
rz(-2.5757289) q[0];
sx q[0];
rz(1.4669482) q[0];
rz(3.103718) q[1];
sx q[1];
rz(-1.9297618) q[1];
sx q[1];
rz(-1.1143335) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1087462) q[0];
sx q[0];
rz(-1.4583197) q[0];
sx q[0];
rz(1.9839909) q[0];
x q[1];
rz(-2.3856106) q[2];
sx q[2];
rz(-2.4376483) q[2];
sx q[2];
rz(-0.24105016) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73015651) q[1];
sx q[1];
rz(-1.1056285) q[1];
sx q[1];
rz(-0.48480715) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30348482) q[3];
sx q[3];
rz(-1.0110098) q[3];
sx q[3];
rz(-1.8523077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.36180878) q[2];
sx q[2];
rz(-0.59541687) q[2];
sx q[2];
rz(-0.53100604) q[2];
rz(-0.94046721) q[3];
sx q[3];
rz(-2.2898424) q[3];
sx q[3];
rz(-1.4550335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4130212) q[0];
sx q[0];
rz(-2.72609) q[0];
sx q[0];
rz(0.79237932) q[0];
rz(-3.0089695) q[1];
sx q[1];
rz(-2.4515929) q[1];
sx q[1];
rz(-0.92540583) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6974119) q[0];
sx q[0];
rz(-1.7562997) q[0];
sx q[0];
rz(2.9836487) q[0];
rz(-pi) q[1];
rz(1.6907352) q[2];
sx q[2];
rz(-0.45443568) q[2];
sx q[2];
rz(1.1392405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5645743) q[1];
sx q[1];
rz(-2.1547103) q[1];
sx q[1];
rz(0.24078943) q[1];
x q[2];
rz(-0.53823353) q[3];
sx q[3];
rz(-0.58659121) q[3];
sx q[3];
rz(-2.339102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35599071) q[2];
sx q[2];
rz(-1.8428558) q[2];
sx q[2];
rz(-1.7388657) q[2];
rz(-1.2518903) q[3];
sx q[3];
rz(-1.3024878) q[3];
sx q[3];
rz(-0.29786626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.70016015) q[0];
sx q[0];
rz(-0.97235632) q[0];
sx q[0];
rz(-0.98883072) q[0];
rz(1.5160457) q[1];
sx q[1];
rz(-1.6700309) q[1];
sx q[1];
rz(0.54738799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.669466) q[0];
sx q[0];
rz(-1.5385813) q[0];
sx q[0];
rz(-2.6526582) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7774974) q[2];
sx q[2];
rz(-1.4949706) q[2];
sx q[2];
rz(-0.66616601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.061364) q[1];
sx q[1];
rz(-1.974611) q[1];
sx q[1];
rz(0.81527995) q[1];
x q[2];
rz(-2.8071466) q[3];
sx q[3];
rz(-2.8502101) q[3];
sx q[3];
rz(-2.4117985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0872515) q[2];
sx q[2];
rz(-0.50830278) q[2];
sx q[2];
rz(-0.15092078) q[2];
rz(1.5058676) q[3];
sx q[3];
rz(-1.8412291) q[3];
sx q[3];
rz(-1.6747564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48581377) q[0];
sx q[0];
rz(-1.9447615) q[0];
sx q[0];
rz(2.6659513) q[0];
rz(0.42426839) q[1];
sx q[1];
rz(-1.905966) q[1];
sx q[1];
rz(1.3729399) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17775336) q[0];
sx q[0];
rz(-0.31597695) q[0];
sx q[0];
rz(2.55992) q[0];
rz(0.67694725) q[2];
sx q[2];
rz(-0.71093762) q[2];
sx q[2];
rz(-0.18260156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.331363) q[1];
sx q[1];
rz(-0.64953564) q[1];
sx q[1];
rz(-1.002921) q[1];
rz(-pi) q[2];
rz(2.0128355) q[3];
sx q[3];
rz(-1.6762432) q[3];
sx q[3];
rz(-0.46322185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0515392) q[2];
sx q[2];
rz(-2.1681483) q[2];
sx q[2];
rz(-0.092183979) q[2];
rz(1.1719545) q[3];
sx q[3];
rz(-0.53297526) q[3];
sx q[3];
rz(-0.91641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6072657) q[0];
sx q[0];
rz(-2.3054275) q[0];
sx q[0];
rz(-1.0323866) q[0];
rz(0.1296002) q[1];
sx q[1];
rz(-1.383129) q[1];
sx q[1];
rz(1.3547156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98785214) q[0];
sx q[0];
rz(-0.82621413) q[0];
sx q[0];
rz(2.4620442) q[0];
rz(-2.7849136) q[2];
sx q[2];
rz(-1.2497777) q[2];
sx q[2];
rz(-0.73250801) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9224841) q[1];
sx q[1];
rz(-1.0661844) q[1];
sx q[1];
rz(-0.64448661) q[1];
rz(-1.0951429) q[3];
sx q[3];
rz(-1.5009592) q[3];
sx q[3];
rz(1.2418408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20595343) q[2];
sx q[2];
rz(-0.26198584) q[2];
sx q[2];
rz(1.6550753) q[2];
rz(0.67241159) q[3];
sx q[3];
rz(-2.0629864) q[3];
sx q[3];
rz(-0.068517223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9392149) q[0];
sx q[0];
rz(-0.12513932) q[0];
sx q[0];
rz(2.8676721) q[0];
rz(2.9778453) q[1];
sx q[1];
rz(-1.3122908) q[1];
sx q[1];
rz(0.40072498) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3253097) q[0];
sx q[0];
rz(-2.0643215) q[0];
sx q[0];
rz(-1.2801583) q[0];
rz(1.2873307) q[2];
sx q[2];
rz(-2.7230883) q[2];
sx q[2];
rz(-1.1039656) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.54998794) q[1];
sx q[1];
rz(-2.2428136) q[1];
sx q[1];
rz(0.023663533) q[1];
rz(1.9951882) q[3];
sx q[3];
rz(-1.8901069) q[3];
sx q[3];
rz(0.13892787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.87216941) q[2];
sx q[2];
rz(-2.9672406) q[2];
sx q[2];
rz(1.4738458) q[2];
rz(-3.040124) q[3];
sx q[3];
rz(-1.2227819) q[3];
sx q[3];
rz(-1.3946704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2673016) q[0];
sx q[0];
rz(-0.75031459) q[0];
sx q[0];
rz(-0.0016203298) q[0];
rz(0.56456176) q[1];
sx q[1];
rz(-1.8058585) q[1];
sx q[1];
rz(2.6394305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82715568) q[0];
sx q[0];
rz(-2.1232044) q[0];
sx q[0];
rz(-2.3236426) q[0];
x q[1];
rz(-0.61567581) q[2];
sx q[2];
rz(-2.4936495) q[2];
sx q[2];
rz(-2.9790762) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6356023) q[1];
sx q[1];
rz(-2.5924304) q[1];
sx q[1];
rz(-0.378774) q[1];
rz(-1.2099482) q[3];
sx q[3];
rz(-0.54665414) q[3];
sx q[3];
rz(0.15693675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0994215) q[2];
sx q[2];
rz(-1.1974988) q[2];
sx q[2];
rz(-1.8966804) q[2];
rz(-0.20243195) q[3];
sx q[3];
rz(-1.3916241) q[3];
sx q[3];
rz(0.99948731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56364432) q[0];
sx q[0];
rz(-0.36643323) q[0];
sx q[0];
rz(1.7250489) q[0];
rz(1.1997403) q[1];
sx q[1];
rz(-1.9681294) q[1];
sx q[1];
rz(-2.8660668) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5274297) q[0];
sx q[0];
rz(-2.6350281) q[0];
sx q[0];
rz(0.81679438) q[0];
rz(-pi) q[1];
rz(-2.7904835) q[2];
sx q[2];
rz(-1.4867696) q[2];
sx q[2];
rz(0.42421266) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4135189) q[1];
sx q[1];
rz(-2.6797543) q[1];
sx q[1];
rz(1.6175265) q[1];
x q[2];
rz(1.4472792) q[3];
sx q[3];
rz(-0.97573167) q[3];
sx q[3];
rz(1.7473011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.766091) q[2];
sx q[2];
rz(-1.4644863) q[2];
sx q[2];
rz(2.7601833) q[2];
rz(2.8219847) q[3];
sx q[3];
rz(-2.3242798) q[3];
sx q[3];
rz(-0.6680502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850591) q[0];
sx q[0];
rz(-1.1548797) q[0];
sx q[0];
rz(-0.67697939) q[0];
rz(1.001724) q[1];
sx q[1];
rz(-1.4717419) q[1];
sx q[1];
rz(-0.91632661) q[1];
rz(1.6940633) q[2];
sx q[2];
rz(-0.90926778) q[2];
sx q[2];
rz(-2.8127083) q[2];
rz(0.65883815) q[3];
sx q[3];
rz(-1.659395) q[3];
sx q[3];
rz(-0.79959662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
