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
rz(2.3393369) q[0];
sx q[0];
rz(-1.383961) q[0];
sx q[0];
rz(-1.2686165) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(6.3422536) q[1];
sx q[1];
rz(8.0191945) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85044139) q[0];
sx q[0];
rz(-2.0401401) q[0];
sx q[0];
rz(-2.6525081) q[0];
rz(0.67063801) q[2];
sx q[2];
rz(-1.3893638) q[2];
sx q[2];
rz(1.2714314) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66623233) q[1];
sx q[1];
rz(-1.1603429) q[1];
sx q[1];
rz(-2.3478541) q[1];
rz(-pi) q[2];
rz(0.89449331) q[3];
sx q[3];
rz(-1.3925806) q[3];
sx q[3];
rz(-1.7591998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9872226) q[2];
sx q[2];
rz(-1.9742249) q[2];
sx q[2];
rz(2.3517189) q[2];
rz(0.023898276) q[3];
sx q[3];
rz(-1.8022715) q[3];
sx q[3];
rz(-0.53643119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24870482) q[0];
sx q[0];
rz(-1.4905812) q[0];
sx q[0];
rz(1.254068) q[0];
rz(1.4887384) q[1];
sx q[1];
rz(-0.87569991) q[1];
sx q[1];
rz(-2.658433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2452908) q[0];
sx q[0];
rz(-1.622513) q[0];
sx q[0];
rz(-1.381676) q[0];
rz(-pi) q[1];
x q[1];
rz(1.05422) q[2];
sx q[2];
rz(-0.90460515) q[2];
sx q[2];
rz(-1.4201128) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4401958) q[1];
sx q[1];
rz(-0.45397511) q[1];
sx q[1];
rz(2.5099436) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7464569) q[3];
sx q[3];
rz(-1.2845728) q[3];
sx q[3];
rz(1.3940879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88863215) q[2];
sx q[2];
rz(-1.6855719) q[2];
sx q[2];
rz(-1.6271094) q[2];
rz(-3.0857981) q[3];
sx q[3];
rz(-1.5972219) q[3];
sx q[3];
rz(1.6519215) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91419879) q[0];
sx q[0];
rz(-0.44454235) q[0];
sx q[0];
rz(3.0134873) q[0];
rz(-0.57614342) q[1];
sx q[1];
rz(-0.73155254) q[1];
sx q[1];
rz(-2.3760956) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.573581) q[0];
sx q[0];
rz(-2.4666602) q[0];
sx q[0];
rz(2.5715067) q[0];
x q[1];
rz(1.9472935) q[2];
sx q[2];
rz(-2.2222179) q[2];
sx q[2];
rz(-2.5848856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24352267) q[1];
sx q[1];
rz(-1.0869893) q[1];
sx q[1];
rz(-2.0953728) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14344826) q[3];
sx q[3];
rz(-1.5455523) q[3];
sx q[3];
rz(2.6526566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3196044) q[2];
sx q[2];
rz(-2.2687948) q[2];
sx q[2];
rz(0.86307159) q[2];
rz(0.34267628) q[3];
sx q[3];
rz(-2.7211012) q[3];
sx q[3];
rz(1.4884523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2494025) q[0];
sx q[0];
rz(-0.96977314) q[0];
sx q[0];
rz(-0.48200193) q[0];
rz(0.30872289) q[1];
sx q[1];
rz(-0.32544193) q[1];
sx q[1];
rz(-0.6122922) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68744874) q[0];
sx q[0];
rz(-0.73546919) q[0];
sx q[0];
rz(-0.81540458) q[0];
x q[1];
rz(1.7365513) q[2];
sx q[2];
rz(-3.0566772) q[2];
sx q[2];
rz(-1.9254534) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.051381342) q[1];
sx q[1];
rz(-1.5259698) q[1];
sx q[1];
rz(2.1771237) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5774229) q[3];
sx q[3];
rz(-1.4215135) q[3];
sx q[3];
rz(1.7501259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.65773949) q[2];
sx q[2];
rz(-2.190399) q[2];
sx q[2];
rz(-1.308002) q[2];
rz(2.7075503) q[3];
sx q[3];
rz(-1.3515892) q[3];
sx q[3];
rz(1.2691809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.16172116) q[0];
sx q[0];
rz(-1.5890108) q[0];
sx q[0];
rz(-0.27444926) q[0];
rz(1.5212003) q[1];
sx q[1];
rz(-1.043964) q[1];
sx q[1];
rz(1.8843947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3109821) q[0];
sx q[0];
rz(-2.794696) q[0];
sx q[0];
rz(-0.78126379) q[0];
rz(-pi) q[1];
rz(2.1047701) q[2];
sx q[2];
rz(-1.2495012) q[2];
sx q[2];
rz(-0.23201135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3538227) q[1];
sx q[1];
rz(-1.0313883) q[1];
sx q[1];
rz(0.25168519) q[1];
rz(-pi) q[2];
rz(-2.2755695) q[3];
sx q[3];
rz(-1.8187722) q[3];
sx q[3];
rz(-0.06125227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7460798) q[2];
sx q[2];
rz(-1.505625) q[2];
sx q[2];
rz(-2.5858322) q[2];
rz(0.79103509) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(-1.6382431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24285862) q[0];
sx q[0];
rz(-1.7019615) q[0];
sx q[0];
rz(-0.71204251) q[0];
rz(0.60246077) q[1];
sx q[1];
rz(-0.41350499) q[1];
sx q[1];
rz(-0.079040225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9475381) q[0];
sx q[0];
rz(-1.506516) q[0];
sx q[0];
rz(3.0966731) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4858133) q[2];
sx q[2];
rz(-1.6459609) q[2];
sx q[2];
rz(2.4369881) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.95922071) q[1];
sx q[1];
rz(-2.3010265) q[1];
sx q[1];
rz(1.275179) q[1];
rz(-1.4677804) q[3];
sx q[3];
rz(-1.8673737) q[3];
sx q[3];
rz(2.1144583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.041212335) q[2];
sx q[2];
rz(-0.93116394) q[2];
sx q[2];
rz(0.97990123) q[2];
rz(-1.529871) q[3];
sx q[3];
rz(-1.9001222) q[3];
sx q[3];
rz(-2.8821168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74705446) q[0];
sx q[0];
rz(-1.8600445) q[0];
sx q[0];
rz(-3.0406612) q[0];
rz(-0.55502564) q[1];
sx q[1];
rz(-2.2075768) q[1];
sx q[1];
rz(1.2947882) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7930113) q[0];
sx q[0];
rz(-1.8037348) q[0];
sx q[0];
rz(2.255054) q[0];
rz(-1.2764658) q[2];
sx q[2];
rz(-1.571725) q[2];
sx q[2];
rz(1.9850736) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87161049) q[1];
sx q[1];
rz(-1.9967604) q[1];
sx q[1];
rz(-1.4130089) q[1];
rz(-2.1168618) q[3];
sx q[3];
rz(-2.6203558) q[3];
sx q[3];
rz(-1.2574399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1780221) q[2];
sx q[2];
rz(-0.87248412) q[2];
sx q[2];
rz(2.8998609) q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80411512) q[0];
sx q[0];
rz(-2.1601456) q[0];
sx q[0];
rz(2.5313983) q[0];
rz(-0.21774165) q[1];
sx q[1];
rz(-2.1861031) q[1];
sx q[1];
rz(2.311923) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9834845) q[0];
sx q[0];
rz(-2.2451441) q[0];
sx q[0];
rz(2.308962) q[0];
x q[1];
rz(0.83135624) q[2];
sx q[2];
rz(-2.6717253) q[2];
sx q[2];
rz(2.510315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.45308847) q[1];
sx q[1];
rz(-1.1646534) q[1];
sx q[1];
rz(-1.7726834) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4368966) q[3];
sx q[3];
rz(-2.4256676) q[3];
sx q[3];
rz(-1.9980824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3024451) q[2];
sx q[2];
rz(-1.382901) q[2];
sx q[2];
rz(-1.1350606) q[2];
rz(-2.6531687) q[3];
sx q[3];
rz(-1.4819744) q[3];
sx q[3];
rz(-1.2357014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14588533) q[0];
sx q[0];
rz(-1.1111525) q[0];
sx q[0];
rz(-2.1671894) q[0];
rz(-0.79565945) q[1];
sx q[1];
rz(-0.37745044) q[1];
sx q[1];
rz(0.61765751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6256339) q[0];
sx q[0];
rz(-3.0581116) q[0];
sx q[0];
rz(0.030103695) q[0];
x q[1];
rz(-2.5133649) q[2];
sx q[2];
rz(-2.5190353) q[2];
sx q[2];
rz(0.33930919) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1682081) q[1];
sx q[1];
rz(-1.0092371) q[1];
sx q[1];
rz(1.8902768) q[1];
rz(-0.32989721) q[3];
sx q[3];
rz(-1.1217692) q[3];
sx q[3];
rz(0.87775081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.016102942) q[2];
sx q[2];
rz(-1.8718655) q[2];
sx q[2];
rz(-1.3886836) q[2];
rz(-1.8359418) q[3];
sx q[3];
rz(-2.4180222) q[3];
sx q[3];
rz(-1.8587662) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0131123) q[0];
sx q[0];
rz(-0.73902577) q[0];
sx q[0];
rz(2.541743) q[0];
rz(-0.57303095) q[1];
sx q[1];
rz(-2.532798) q[1];
sx q[1];
rz(-2.1175687) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0759648) q[0];
sx q[0];
rz(-0.5813404) q[0];
sx q[0];
rz(0.79728787) q[0];
rz(-pi) q[1];
rz(0.65262087) q[2];
sx q[2];
rz(-1.5372477) q[2];
sx q[2];
rz(-1.5198358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8166067) q[1];
sx q[1];
rz(-1.6546975) q[1];
sx q[1];
rz(0.3971667) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7598466) q[3];
sx q[3];
rz(-1.7110398) q[3];
sx q[3];
rz(0.5776757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92274252) q[2];
sx q[2];
rz(-1.3877733) q[2];
sx q[2];
rz(0.77524033) q[2];
rz(1.5357337) q[3];
sx q[3];
rz(-1.1865059) q[3];
sx q[3];
rz(1.8839914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53104644) q[0];
sx q[0];
rz(-2.1903867) q[0];
sx q[0];
rz(1.2280986) q[0];
rz(1.7465406) q[1];
sx q[1];
rz(-1.4735305) q[1];
sx q[1];
rz(2.7458618) q[1];
rz(2.6880245) q[2];
sx q[2];
rz(-1.8487052) q[2];
sx q[2];
rz(-2.1089274) q[2];
rz(-2.3395634) q[3];
sx q[3];
rz(-1.5475954) q[3];
sx q[3];
rz(-2.2022998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
