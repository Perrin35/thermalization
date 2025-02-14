OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(1.9429053) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(-0.79678798) q[1];
sx q[1];
rz(0.30153433) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806592) q[0];
sx q[0];
rz(-1.7552688) q[0];
sx q[0];
rz(-0.16452275) q[0];
x q[1];
rz(-1.0204591) q[2];
sx q[2];
rz(-0.34780234) q[2];
sx q[2];
rz(-2.2344799) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.49820013) q[1];
sx q[1];
rz(-0.82005703) q[1];
sx q[1];
rz(3.0489075) q[1];
rz(-pi) q[2];
rz(1.1490378) q[3];
sx q[3];
rz(-0.43614498) q[3];
sx q[3];
rz(-2.8835322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4528759) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(1.1348881) q[2];
rz(-0.12198837) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(-0.056099135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-1.6720471) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(3.1392414) q[0];
rz(2.7408842) q[1];
sx q[1];
rz(-1.7727163) q[1];
sx q[1];
rz(-0.8561264) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6295638) q[0];
sx q[0];
rz(-2.1855547) q[0];
sx q[0];
rz(0.31350664) q[0];
x q[1];
rz(-1.0925757) q[2];
sx q[2];
rz(-0.56861955) q[2];
sx q[2];
rz(-0.61517119) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5634671) q[1];
sx q[1];
rz(-1.620786) q[1];
sx q[1];
rz(1.0432788) q[1];
x q[2];
rz(-0.56613381) q[3];
sx q[3];
rz(-1.2581035) q[3];
sx q[3];
rz(-1.4110402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.037228435) q[2];
sx q[2];
rz(-1.2040441) q[2];
sx q[2];
rz(-0.47002235) q[2];
rz(-0.59703279) q[3];
sx q[3];
rz(-0.11490122) q[3];
sx q[3];
rz(-1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700579) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(-0.22802995) q[0];
rz(2.6920964) q[1];
sx q[1];
rz(-2.1411965) q[1];
sx q[1];
rz(-2.1953348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3458918) q[0];
sx q[0];
rz(-1.8270703) q[0];
sx q[0];
rz(1.9836224) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65785711) q[2];
sx q[2];
rz(-2.0924753) q[2];
sx q[2];
rz(1.9470095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6201278) q[1];
sx q[1];
rz(-1.2343654) q[1];
sx q[1];
rz(1.7640339) q[1];
rz(-pi) q[2];
rz(1.2872996) q[3];
sx q[3];
rz(-1.4186267) q[3];
sx q[3];
rz(-0.37934549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9049282) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(-1.6620592) q[2];
rz(-0.18105257) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(2.2553867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36412305) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(0.89853483) q[0];
rz(0.50615519) q[1];
sx q[1];
rz(-1.2944784) q[1];
sx q[1];
rz(-1.4599962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.579297) q[0];
sx q[0];
rz(-1.0387207) q[0];
sx q[0];
rz(-1.4946372) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6674472) q[2];
sx q[2];
rz(-2.3398826) q[2];
sx q[2];
rz(0.8566044) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0467028) q[1];
sx q[1];
rz(-3.0071435) q[1];
sx q[1];
rz(-1.0195169) q[1];
rz(-pi) q[2];
rz(-1.4921032) q[3];
sx q[3];
rz(-0.65641145) q[3];
sx q[3];
rz(0.66167964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8755181) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(-1.6391222) q[2];
rz(1.8858887) q[3];
sx q[3];
rz(-1.4016822) q[3];
sx q[3];
rz(-0.046253117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2332377) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(0.18049845) q[0];
rz(1.7620697) q[1];
sx q[1];
rz(-2.5738398) q[1];
sx q[1];
rz(2.0464121) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9126606) q[0];
sx q[0];
rz(-1.7062418) q[0];
sx q[0];
rz(1.5082466) q[0];
rz(-pi) q[1];
rz(0.081549598) q[2];
sx q[2];
rz(-2.3429541) q[2];
sx q[2];
rz(0.35432409) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.74873456) q[1];
sx q[1];
rz(-1.7104516) q[1];
sx q[1];
rz(2.3631494) q[1];
rz(-0.88222701) q[3];
sx q[3];
rz(-1.467657) q[3];
sx q[3];
rz(2.8571199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77202648) q[2];
sx q[2];
rz(-2.2942746) q[2];
sx q[2];
rz(-1.0687211) q[2];
rz(-0.42701834) q[3];
sx q[3];
rz(-0.87978274) q[3];
sx q[3];
rz(-1.2515742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0291979) q[0];
sx q[0];
rz(-0.91570941) q[0];
sx q[0];
rz(0.29119626) q[0];
rz(-1.4578106) q[1];
sx q[1];
rz(-1.0271415) q[1];
sx q[1];
rz(2.2529032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4584258) q[0];
sx q[0];
rz(-2.1447276) q[0];
sx q[0];
rz(2.5794585) q[0];
x q[1];
rz(-2.3843147) q[2];
sx q[2];
rz(-1.5520397) q[2];
sx q[2];
rz(2.4919486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3515489) q[1];
sx q[1];
rz(-2.1218532) q[1];
sx q[1];
rz(-2.181777) q[1];
rz(-2.3360152) q[3];
sx q[3];
rz(-2.71851) q[3];
sx q[3];
rz(-0.42832908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39151829) q[2];
sx q[2];
rz(-1.3041648) q[2];
sx q[2];
rz(2.1938426) q[2];
rz(2.9446972) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(0.31015629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8295558) q[0];
sx q[0];
rz(-1.1126248) q[0];
sx q[0];
rz(0.46992508) q[0];
rz(1.4620818) q[1];
sx q[1];
rz(-1.3009678) q[1];
sx q[1];
rz(1.7376815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3199661) q[0];
sx q[0];
rz(-1.9187154) q[0];
sx q[0];
rz(-0.97396429) q[0];
x q[1];
rz(-0.20848863) q[2];
sx q[2];
rz(-2.0402153) q[2];
sx q[2];
rz(2.4698348) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1199855) q[1];
sx q[1];
rz(-2.48263) q[1];
sx q[1];
rz(-2.1849385) q[1];
rz(2.6440355) q[3];
sx q[3];
rz(-0.64421875) q[3];
sx q[3];
rz(-1.8807172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64138428) q[2];
sx q[2];
rz(-2.310014) q[2];
sx q[2];
rz(0.10759648) q[2];
rz(0.61947668) q[3];
sx q[3];
rz(-1.2958801) q[3];
sx q[3];
rz(1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7463995) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(0.25303823) q[0];
rz(0.088317618) q[1];
sx q[1];
rz(-2.2679236) q[1];
sx q[1];
rz(-2.1926682) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2491566) q[0];
sx q[0];
rz(-1.8674441) q[0];
sx q[0];
rz(-2.9846624) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19472943) q[2];
sx q[2];
rz(-1.1695332) q[2];
sx q[2];
rz(1.8068318) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9946549) q[1];
sx q[1];
rz(-0.48284621) q[1];
sx q[1];
rz(-1.4551244) q[1];
rz(-pi) q[2];
rz(2.3036495) q[3];
sx q[3];
rz(-0.70819127) q[3];
sx q[3];
rz(1.4670682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(2.8308947) q[2];
rz(0.24142309) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(-2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017515) q[0];
sx q[0];
rz(-0.40626353) q[0];
sx q[0];
rz(1.9061506) q[0];
rz(2.6605117) q[1];
sx q[1];
rz(-0.95293871) q[1];
sx q[1];
rz(1.7135886) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3612551) q[0];
sx q[0];
rz(-1.5819823) q[0];
sx q[0];
rz(-0.7594603) q[0];
x q[1];
rz(-2.253654) q[2];
sx q[2];
rz(-0.35338923) q[2];
sx q[2];
rz(1.8977752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2216827) q[1];
sx q[1];
rz(-2.7836728) q[1];
sx q[1];
rz(-0.81166761) q[1];
rz(-pi) q[2];
rz(0.72254027) q[3];
sx q[3];
rz(-1.3219537) q[3];
sx q[3];
rz(-1.794508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1741751) q[2];
sx q[2];
rz(-2.316541) q[2];
sx q[2];
rz(3.1136801) q[2];
rz(-0.86999718) q[3];
sx q[3];
rz(-2.3614707) q[3];
sx q[3];
rz(-1.1869259) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025539909) q[0];
sx q[0];
rz(-1.5486131) q[0];
sx q[0];
rz(-2.0822339) q[0];
rz(1.7268044) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(0.47763225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519417) q[0];
sx q[0];
rz(-2.383259) q[0];
sx q[0];
rz(0.21679057) q[0];
rz(-pi) q[1];
rz(-1.9066493) q[2];
sx q[2];
rz(-1.5689465) q[2];
sx q[2];
rz(0.84650485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79523008) q[1];
sx q[1];
rz(-1.5168961) q[1];
sx q[1];
rz(1.6832215) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60250256) q[3];
sx q[3];
rz(-0.91711603) q[3];
sx q[3];
rz(2.9127171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9719505) q[2];
sx q[2];
rz(-0.98196882) q[2];
sx q[2];
rz(0.36894813) q[2];
rz(1.2667027) q[3];
sx q[3];
rz(-0.72487512) q[3];
sx q[3];
rz(-2.1784311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030768) q[0];
sx q[0];
rz(-1.2110447) q[0];
sx q[0];
rz(-1.3883653) q[0];
rz(0.44454642) q[1];
sx q[1];
rz(-2.3311756) q[1];
sx q[1];
rz(-0.12745007) q[1];
rz(1.2426022) q[2];
sx q[2];
rz(-1.9933619) q[2];
sx q[2];
rz(0.91771916) q[2];
rz(-3.1034341) q[3];
sx q[3];
rz(-0.46783075) q[3];
sx q[3];
rz(-0.086188407) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
