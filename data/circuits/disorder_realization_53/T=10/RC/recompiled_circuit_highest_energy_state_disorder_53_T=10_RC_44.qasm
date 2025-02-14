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
rz(-0.80225575) q[0];
sx q[0];
rz(-1.7576317) q[0];
sx q[0];
rz(1.2686165) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(-3.0825244) q[1];
sx q[1];
rz(-1.7360092) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168306) q[0];
sx q[0];
rz(-0.66436902) q[0];
sx q[0];
rz(-0.82358255) q[0];
rz(-pi) q[1];
rz(0.67063801) q[2];
sx q[2];
rz(-1.7522289) q[2];
sx q[2];
rz(1.8701613) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6224712) q[1];
sx q[1];
rz(-0.85825413) q[1];
sx q[1];
rz(1.0153517) q[1];
rz(0.89449331) q[3];
sx q[3];
rz(-1.749012) q[3];
sx q[3];
rz(1.7591998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1543701) q[2];
sx q[2];
rz(-1.1673678) q[2];
sx q[2];
rz(0.78987375) q[2];
rz(-0.023898276) q[3];
sx q[3];
rz(-1.3393211) q[3];
sx q[3];
rz(2.6051615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8928878) q[0];
sx q[0];
rz(-1.6510115) q[0];
sx q[0];
rz(1.254068) q[0];
rz(-1.6528543) q[1];
sx q[1];
rz(-2.2658927) q[1];
sx q[1];
rz(2.658433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.457204) q[0];
sx q[0];
rz(-1.7596608) q[0];
sx q[0];
rz(0.052653814) q[0];
x q[1];
rz(-1.05422) q[2];
sx q[2];
rz(-2.2369875) q[2];
sx q[2];
rz(-1.4201128) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3846783) q[1];
sx q[1];
rz(-1.9325629) q[1];
sx q[1];
rz(-1.2902618) q[1];
rz(2.8511413) q[3];
sx q[3];
rz(-1.402352) q[3];
sx q[3];
rz(-3.0149533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.88863215) q[2];
sx q[2];
rz(-1.6855719) q[2];
sx q[2];
rz(1.5144833) q[2];
rz(3.0857981) q[3];
sx q[3];
rz(-1.5972219) q[3];
sx q[3];
rz(1.4896711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2273939) q[0];
sx q[0];
rz(-0.44454235) q[0];
sx q[0];
rz(-0.12810531) q[0];
rz(-0.57614342) q[1];
sx q[1];
rz(-2.4100401) q[1];
sx q[1];
rz(-0.76549706) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.573581) q[0];
sx q[0];
rz(-2.4666602) q[0];
sx q[0];
rz(-2.5715067) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9472935) q[2];
sx q[2];
rz(-2.2222179) q[2];
sx q[2];
rz(-2.5848856) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.89807) q[1];
sx q[1];
rz(-2.0546034) q[1];
sx q[1];
rz(1.0462199) q[1];
rz(1.5963022) q[3];
sx q[3];
rz(-1.4273941) q[3];
sx q[3];
rz(2.0633782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3196044) q[2];
sx q[2];
rz(-2.2687948) q[2];
sx q[2];
rz(0.86307159) q[2];
rz(-0.34267628) q[3];
sx q[3];
rz(-0.42049146) q[3];
sx q[3];
rz(-1.6531403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2494025) q[0];
sx q[0];
rz(-2.1718195) q[0];
sx q[0];
rz(-2.6595907) q[0];
rz(2.8328698) q[1];
sx q[1];
rz(-0.32544193) q[1];
sx q[1];
rz(-2.5293005) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9252602) q[0];
sx q[0];
rz(-2.0811006) q[0];
sx q[0];
rz(-2.5863674) q[0];
x q[1];
rz(-0.014043645) q[2];
sx q[2];
rz(-1.4870475) q[2];
sx q[2];
rz(1.75911) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5911124) q[1];
sx q[1];
rz(-0.96516536) q[1];
sx q[1];
rz(3.0870599) q[1];
rz(-pi) q[2];
rz(2.9923066) q[3];
sx q[3];
rz(-1.5773492) q[3];
sx q[3];
rz(0.1803151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65773949) q[2];
sx q[2];
rz(-2.190399) q[2];
sx q[2];
rz(-1.308002) q[2];
rz(-0.43404239) q[3];
sx q[3];
rz(-1.7900034) q[3];
sx q[3];
rz(-1.2691809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9798715) q[0];
sx q[0];
rz(-1.5890108) q[0];
sx q[0];
rz(2.8671434) q[0];
rz(1.6203923) q[1];
sx q[1];
rz(-1.043964) q[1];
sx q[1];
rz(1.257198) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1311943) q[0];
sx q[0];
rz(-1.3290414) q[0];
sx q[0];
rz(-0.25126033) q[0];
rz(-2.7726452) q[2];
sx q[2];
rz(-1.0668179) q[2];
sx q[2];
rz(1.9874017) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.32361832) q[1];
sx q[1];
rz(-0.58992851) q[1];
sx q[1];
rz(-1.9650311) q[1];
x q[2];
rz(1.9433504) q[3];
sx q[3];
rz(-0.74001678) q[3];
sx q[3];
rz(-1.3510902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7460798) q[2];
sx q[2];
rz(-1.505625) q[2];
sx q[2];
rz(0.5557605) q[2];
rz(0.79103509) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(-1.6382431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898734) q[0];
sx q[0];
rz(-1.4396311) q[0];
sx q[0];
rz(-2.4295501) q[0];
rz(-0.60246077) q[1];
sx q[1];
rz(-0.41350499) q[1];
sx q[1];
rz(0.079040225) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1940546) q[0];
sx q[0];
rz(-1.506516) q[0];
sx q[0];
rz(-0.044919515) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84509648) q[2];
sx q[2];
rz(-3.0281986) q[2];
sx q[2];
rz(1.5886943) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5309294) q[1];
sx q[1];
rz(-2.3641415) q[1];
sx q[1];
rz(-0.31458293) q[1];
rz(0.2980663) q[3];
sx q[3];
rz(-1.4722927) q[3];
sx q[3];
rz(0.57386604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.041212335) q[2];
sx q[2];
rz(-2.2104287) q[2];
sx q[2];
rz(-0.97990123) q[2];
rz(-1.529871) q[3];
sx q[3];
rz(-1.2414705) q[3];
sx q[3];
rz(2.8821168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3945382) q[0];
sx q[0];
rz(-1.2815481) q[0];
sx q[0];
rz(0.10093149) q[0];
rz(-2.586567) q[1];
sx q[1];
rz(-2.2075768) q[1];
sx q[1];
rz(-1.2947882) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7930113) q[0];
sx q[0];
rz(-1.3378578) q[0];
sx q[0];
rz(-0.88653867) q[0];
rz(-1.8651269) q[2];
sx q[2];
rz(-1.5698676) q[2];
sx q[2];
rz(1.9850736) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2699822) q[1];
sx q[1];
rz(-1.1448323) q[1];
sx q[1];
rz(-1.4130089) q[1];
rz(-1.0247308) q[3];
sx q[3];
rz(-0.52123681) q[3];
sx q[3];
rz(1.8841528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9635705) q[2];
sx q[2];
rz(-2.2691085) q[2];
sx q[2];
rz(0.24173173) q[2];
rz(-0.47197765) q[3];
sx q[3];
rz(-2.9265672) q[3];
sx q[3];
rz(1.5836466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3374775) q[0];
sx q[0];
rz(-2.1601456) q[0];
sx q[0];
rz(2.5313983) q[0];
rz(2.923851) q[1];
sx q[1];
rz(-2.1861031) q[1];
sx q[1];
rz(2.311923) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0141896) q[0];
sx q[0];
rz(-0.9547736) q[0];
sx q[0];
rz(0.69973972) q[0];
x q[1];
rz(-2.81189) q[2];
sx q[2];
rz(-1.2296943) q[2];
sx q[2];
rz(1.4280048) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45308847) q[1];
sx q[1];
rz(-1.9769393) q[1];
sx q[1];
rz(-1.7726834) q[1];
x q[2];
rz(-3.025981) q[3];
sx q[3];
rz(-2.278961) q[3];
sx q[3];
rz(-2.1747605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3024451) q[2];
sx q[2];
rz(-1.382901) q[2];
sx q[2];
rz(-2.006532) q[2];
rz(-2.6531687) q[3];
sx q[3];
rz(-1.4819744) q[3];
sx q[3];
rz(1.9058913) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14588533) q[0];
sx q[0];
rz(-1.1111525) q[0];
sx q[0];
rz(2.1671894) q[0];
rz(-0.79565945) q[1];
sx q[1];
rz(-0.37745044) q[1];
sx q[1];
rz(0.61765751) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5159588) q[0];
sx q[0];
rz(-0.083481073) q[0];
sx q[0];
rz(0.030103695) q[0];
rz(-pi) q[1];
rz(-0.52613132) q[2];
sx q[2];
rz(-1.9205893) q[2];
sx q[2];
rz(-2.443231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61231919) q[1];
sx q[1];
rz(-2.5041083) q[1];
sx q[1];
rz(0.46302621) q[1];
rz(0.97891221) q[3];
sx q[3];
rz(-2.5911461) q[3];
sx q[3];
rz(1.5456256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1254897) q[2];
sx q[2];
rz(-1.2697271) q[2];
sx q[2];
rz(1.7529091) q[2];
rz(-1.8359418) q[3];
sx q[3];
rz(-0.7235705) q[3];
sx q[3];
rz(1.8587662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1284803) q[0];
sx q[0];
rz(-2.4025669) q[0];
sx q[0];
rz(-2.541743) q[0];
rz(2.5685617) q[1];
sx q[1];
rz(-2.532798) q[1];
sx q[1];
rz(1.0240239) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0759648) q[0];
sx q[0];
rz(-2.5602523) q[0];
sx q[0];
rz(0.79728787) q[0];
rz(0.055209514) q[2];
sx q[2];
rz(-0.6533567) q[2];
sx q[2];
rz(0.094815985) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9309195) q[1];
sx q[1];
rz(-1.1751047) q[1];
sx q[1];
rz(-1.4798505) q[1];
x q[2];
rz(1.7217595) q[3];
sx q[3];
rz(-1.9486041) q[3];
sx q[3];
rz(-1.0491766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92274252) q[2];
sx q[2];
rz(-1.7538193) q[2];
sx q[2];
rz(2.3663523) q[2];
rz(1.5357337) q[3];
sx q[3];
rz(-1.9550867) q[3];
sx q[3];
rz(1.2576013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6105462) q[0];
sx q[0];
rz(-0.95120593) q[0];
sx q[0];
rz(-1.9134941) q[0];
rz(1.3950521) q[1];
sx q[1];
rz(-1.6680622) q[1];
sx q[1];
rz(-0.39573085) q[1];
rz(-0.57714085) q[2];
sx q[2];
rz(-2.6147523) q[2];
sx q[2];
rz(-0.025512248) q[2];
rz(1.6041605) q[3];
sx q[3];
rz(-2.3725474) q[3];
sx q[3];
rz(-0.65548246) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
