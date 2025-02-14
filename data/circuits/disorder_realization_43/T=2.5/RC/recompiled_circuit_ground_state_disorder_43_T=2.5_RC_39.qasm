OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3620152) q[0];
sx q[0];
rz(0.22688046) q[0];
sx q[0];
rz(14.332834) q[0];
rz(-0.094376266) q[1];
sx q[1];
rz(-2.2278995) q[1];
sx q[1];
rz(2.8606666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.323712) q[0];
sx q[0];
rz(-1.7294939) q[0];
sx q[0];
rz(-2.236332) q[0];
rz(1.2035349) q[2];
sx q[2];
rz(-0.31031552) q[2];
sx q[2];
rz(2.4825437) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.62796616) q[1];
sx q[1];
rz(-1.1556323) q[1];
sx q[1];
rz(-2.8529851) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5689718) q[3];
sx q[3];
rz(-1.2702281) q[3];
sx q[3];
rz(0.68651783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0278339) q[2];
sx q[2];
rz(-1.4049302) q[2];
sx q[2];
rz(-3.0244381) q[2];
rz(2.8095918) q[3];
sx q[3];
rz(-0.74688512) q[3];
sx q[3];
rz(-2.3565256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.4431045) q[0];
sx q[0];
rz(-0.76258689) q[0];
sx q[0];
rz(2.7681328) q[0];
rz(0.39730486) q[1];
sx q[1];
rz(-1.9970857) q[1];
sx q[1];
rz(0.73748803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0126368) q[0];
sx q[0];
rz(-2.1018493) q[0];
sx q[0];
rz(-0.80967243) q[0];
x q[1];
rz(0.55464427) q[2];
sx q[2];
rz(-2.0337542) q[2];
sx q[2];
rz(-2.2210768) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.26230815) q[1];
sx q[1];
rz(-1.4170987) q[1];
sx q[1];
rz(0.55262312) q[1];
rz(-pi) q[2];
rz(2.3839398) q[3];
sx q[3];
rz(-1.831358) q[3];
sx q[3];
rz(0.16989947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.071659) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(-2.7640479) q[2];
rz(-2.8318882) q[3];
sx q[3];
rz(-1.0891424) q[3];
sx q[3];
rz(-1.6449876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6250703) q[0];
sx q[0];
rz(-2.3023038) q[0];
sx q[0];
rz(1.2155493) q[0];
rz(-1.0141605) q[1];
sx q[1];
rz(-1.4233669) q[1];
sx q[1];
rz(-2.1713712) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.559186) q[0];
sx q[0];
rz(-2.0699887) q[0];
sx q[0];
rz(3.0954719) q[0];
rz(-2.2367291) q[2];
sx q[2];
rz(-2.3653125) q[2];
sx q[2];
rz(-1.5011476) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41170317) q[1];
sx q[1];
rz(-1.7474993) q[1];
sx q[1];
rz(-1.2248301) q[1];
rz(-1.6861451) q[3];
sx q[3];
rz(-1.2209425) q[3];
sx q[3];
rz(0.11125166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0936475) q[2];
sx q[2];
rz(-2.2990172) q[2];
sx q[2];
rz(0.45316163) q[2];
rz(1.9122745) q[3];
sx q[3];
rz(-1.8639576) q[3];
sx q[3];
rz(2.9696828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9446843) q[0];
sx q[0];
rz(-0.6854282) q[0];
sx q[0];
rz(0.36566439) q[0];
rz(-2.2293495) q[1];
sx q[1];
rz(-1.279) q[1];
sx q[1];
rz(-2.7361187) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49783865) q[0];
sx q[0];
rz(-2.9037243) q[0];
sx q[0];
rz(-1.0435304) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22511668) q[2];
sx q[2];
rz(-2.3729513) q[2];
sx q[2];
rz(0.74736881) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1846611) q[1];
sx q[1];
rz(-1.4029364) q[1];
sx q[1];
rz(-2.6015758) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.502127) q[3];
sx q[3];
rz(-1.829095) q[3];
sx q[3];
rz(2.4909508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.071986467) q[2];
sx q[2];
rz(-1.7013841) q[2];
sx q[2];
rz(1.5427422) q[2];
rz(0.37374464) q[3];
sx q[3];
rz(-1.475324) q[3];
sx q[3];
rz(1.6811949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54295802) q[0];
sx q[0];
rz(-2.3668508) q[0];
sx q[0];
rz(2.4454818) q[0];
rz(1.2593345) q[1];
sx q[1];
rz(-1.101661) q[1];
sx q[1];
rz(-0.36516821) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2717948) q[0];
sx q[0];
rz(-2.4569017) q[0];
sx q[0];
rz(2.4224046) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2155207) q[2];
sx q[2];
rz(-0.50130166) q[2];
sx q[2];
rz(-2.9013366) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2272353) q[1];
sx q[1];
rz(-1.0252756) q[1];
sx q[1];
rz(-1.5724214) q[1];
rz(-2.5138084) q[3];
sx q[3];
rz(-1.2300228) q[3];
sx q[3];
rz(2.9140811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.460707) q[2];
sx q[2];
rz(-1.978771) q[2];
sx q[2];
rz(-0.17670259) q[2];
rz(-1.7548615) q[3];
sx q[3];
rz(-1.0597798) q[3];
sx q[3];
rz(2.7669014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6143167) q[0];
sx q[0];
rz(-0.85955954) q[0];
sx q[0];
rz(1.543462) q[0];
rz(-1.3044926) q[1];
sx q[1];
rz(-1.0544798) q[1];
sx q[1];
rz(0.80610448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.971834) q[0];
sx q[0];
rz(-0.8358347) q[0];
sx q[0];
rz(1.2494837) q[0];
rz(0.97053846) q[2];
sx q[2];
rz(-1.1379062) q[2];
sx q[2];
rz(2.8899756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66530217) q[1];
sx q[1];
rz(-0.21033439) q[1];
sx q[1];
rz(-1.3306985) q[1];
rz(-2.0982355) q[3];
sx q[3];
rz(-2.0813848) q[3];
sx q[3];
rz(3.1289738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1690037) q[2];
sx q[2];
rz(-2.6469595) q[2];
sx q[2];
rz(3.0779823) q[2];
rz(-0.58498597) q[3];
sx q[3];
rz(-1.7413185) q[3];
sx q[3];
rz(-1.6271648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99120283) q[0];
sx q[0];
rz(-1.7504033) q[0];
sx q[0];
rz(-0.72917953) q[0];
rz(-1.179262) q[1];
sx q[1];
rz(-0.74960342) q[1];
sx q[1];
rz(-2.8470305) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33074327) q[0];
sx q[0];
rz(-1.0746351) q[0];
sx q[0];
rz(-1.2788354) q[0];
rz(1.8053889) q[2];
sx q[2];
rz(-0.44371191) q[2];
sx q[2];
rz(3.0598874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8940436) q[1];
sx q[1];
rz(-1.4534543) q[1];
sx q[1];
rz(-2.3168922) q[1];
rz(-2.2410349) q[3];
sx q[3];
rz(-1.2646741) q[3];
sx q[3];
rz(-2.8365119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6509167) q[2];
sx q[2];
rz(-1.7016405) q[2];
sx q[2];
rz(0.73224625) q[2];
rz(-0.8693153) q[3];
sx q[3];
rz(-1.9711875) q[3];
sx q[3];
rz(1.4066345) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.231584) q[0];
sx q[0];
rz(-1.5554447) q[0];
sx q[0];
rz(-0.79767942) q[0];
rz(-0.32294598) q[1];
sx q[1];
rz(-2.1452417) q[1];
sx q[1];
rz(1.7332044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7777611) q[0];
sx q[0];
rz(-2.1182502) q[0];
sx q[0];
rz(1.8108032) q[0];
rz(-pi) q[1];
rz(1.1118719) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(-1.9314507) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47695404) q[1];
sx q[1];
rz(-2.3424397) q[1];
sx q[1];
rz(-0.33070143) q[1];
rz(-pi) q[2];
rz(0.50813913) q[3];
sx q[3];
rz(-2.7781825) q[3];
sx q[3];
rz(-1.5809218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8352167) q[2];
sx q[2];
rz(-1.110346) q[2];
sx q[2];
rz(-0.19732538) q[2];
rz(-2.3399682) q[3];
sx q[3];
rz(-2.1620731) q[3];
sx q[3];
rz(0.42158034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1590969) q[0];
sx q[0];
rz(-1.4747341) q[0];
sx q[0];
rz(2.2267447) q[0];
rz(-2.9313056) q[1];
sx q[1];
rz(-0.57840127) q[1];
sx q[1];
rz(-0.68967825) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8473916) q[0];
sx q[0];
rz(-0.24450824) q[0];
sx q[0];
rz(-0.39347009) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.44697) q[2];
sx q[2];
rz(-3.0231907) q[2];
sx q[2];
rz(1.8121383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0090078) q[1];
sx q[1];
rz(-2.0783536) q[1];
sx q[1];
rz(-1.0222438) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34387572) q[3];
sx q[3];
rz(-1.0707885) q[3];
sx q[3];
rz(2.9150034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0476394) q[2];
sx q[2];
rz(-2.2932105) q[2];
sx q[2];
rz(0.32996714) q[2];
rz(2.0983569) q[3];
sx q[3];
rz(-1.7236575) q[3];
sx q[3];
rz(0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1107776) q[0];
sx q[0];
rz(-1.8106221) q[0];
sx q[0];
rz(-1.0571085) q[0];
rz(-2.8874176) q[1];
sx q[1];
rz(-1.9994241) q[1];
sx q[1];
rz(-2.4370297) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1833558) q[0];
sx q[0];
rz(-0.57794774) q[0];
sx q[0];
rz(2.1653752) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.061526394) q[2];
sx q[2];
rz(-1.3869993) q[2];
sx q[2];
rz(-1.4169324) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5479991) q[1];
sx q[1];
rz(-1.3915466) q[1];
sx q[1];
rz(2.0229193) q[1];
x q[2];
rz(2.4289565) q[3];
sx q[3];
rz(-0.27240005) q[3];
sx q[3];
rz(-2.2505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5083984) q[2];
sx q[2];
rz(-1.9693547) q[2];
sx q[2];
rz(-2.634826) q[2];
rz(2.9554328) q[3];
sx q[3];
rz(-2.9045744) q[3];
sx q[3];
rz(-0.53562927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3770461) q[0];
sx q[0];
rz(-1.6896387) q[0];
sx q[0];
rz(1.1402546) q[0];
rz(-0.90691943) q[1];
sx q[1];
rz(-0.74105558) q[1];
sx q[1];
rz(1.8190609) q[1];
rz(2.604031) q[2];
sx q[2];
rz(-0.96502177) q[2];
sx q[2];
rz(1.6966664) q[2];
rz(1.838253) q[3];
sx q[3];
rz(-1.8831913) q[3];
sx q[3];
rz(2.5358653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
