OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.59453073) q[0];
sx q[0];
rz(-1.1214331) q[0];
sx q[0];
rz(-2.9601331) q[0];
rz(2.060086) q[1];
sx q[1];
rz(-0.67343155) q[1];
sx q[1];
rz(-1.0531309) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5232789) q[0];
sx q[0];
rz(-1.3594472) q[0];
sx q[0];
rz(1.1002512) q[0];
rz(-0.62280957) q[2];
sx q[2];
rz(-1.8138759) q[2];
sx q[2];
rz(-2.958975) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6932106) q[1];
sx q[1];
rz(-0.795185) q[1];
sx q[1];
rz(-2.8661212) q[1];
rz(2.9795322) q[3];
sx q[3];
rz(-1.2352304) q[3];
sx q[3];
rz(-3.0367362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9782605) q[2];
sx q[2];
rz(-2.0330727) q[2];
sx q[2];
rz(-1.367761) q[2];
rz(-1.0129499) q[3];
sx q[3];
rz(-0.84665853) q[3];
sx q[3];
rz(-3.0701385) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.098009) q[0];
sx q[0];
rz(-1.4494891) q[0];
sx q[0];
rz(2.5464771) q[0];
rz(2.0960506) q[1];
sx q[1];
rz(-1.4141934) q[1];
sx q[1];
rz(1.6275303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94905108) q[0];
sx q[0];
rz(-2.4727614) q[0];
sx q[0];
rz(2.2632954) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3081231) q[2];
sx q[2];
rz(-1.5916628) q[2];
sx q[2];
rz(1.2534864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0552057) q[1];
sx q[1];
rz(-1.7654164) q[1];
sx q[1];
rz(-2.9662532) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6731095) q[3];
sx q[3];
rz(-1.521109) q[3];
sx q[3];
rz(-1.697584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65511584) q[2];
sx q[2];
rz(-2.1003508) q[2];
sx q[2];
rz(2.3454323) q[2];
rz(-2.1697309) q[3];
sx q[3];
rz(-2.4367417) q[3];
sx q[3];
rz(2.9698353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3765091) q[0];
sx q[0];
rz(-0.77873814) q[0];
sx q[0];
rz(0.064095108) q[0];
rz(-0.31072101) q[1];
sx q[1];
rz(-1.4704082) q[1];
sx q[1];
rz(-1.4583189) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653939) q[0];
sx q[0];
rz(-1.6178693) q[0];
sx q[0];
rz(-1.7182299) q[0];
rz(-0.65285039) q[2];
sx q[2];
rz(-2.0712426) q[2];
sx q[2];
rz(-2.1476538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8628143) q[1];
sx q[1];
rz(-1.1618975) q[1];
sx q[1];
rz(-1.6294953) q[1];
x q[2];
rz(1.5536669) q[3];
sx q[3];
rz(-1.8855842) q[3];
sx q[3];
rz(1.1312315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1594499) q[2];
sx q[2];
rz(-0.84656707) q[2];
sx q[2];
rz(2.4625051) q[2];
rz(1.3726161) q[3];
sx q[3];
rz(-1.8434098) q[3];
sx q[3];
rz(2.1931271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452334) q[0];
sx q[0];
rz(-0.56931749) q[0];
sx q[0];
rz(1.1244208) q[0];
rz(-1.9212978) q[1];
sx q[1];
rz(-1.1923469) q[1];
sx q[1];
rz(-1.6569998) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254018) q[0];
sx q[0];
rz(-1.6214217) q[0];
sx q[0];
rz(-3.1111858) q[0];
rz(0.19214432) q[2];
sx q[2];
rz(-2.0803183) q[2];
sx q[2];
rz(-1.1332878) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0255819) q[1];
sx q[1];
rz(-0.75041795) q[1];
sx q[1];
rz(-2.2376899) q[1];
rz(-pi) q[2];
rz(-1.2403537) q[3];
sx q[3];
rz(-2.2664865) q[3];
sx q[3];
rz(0.89360039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0174039) q[2];
sx q[2];
rz(-1.2839395) q[2];
sx q[2];
rz(2.1253288) q[2];
rz(-1.4034363) q[3];
sx q[3];
rz(-1.6442464) q[3];
sx q[3];
rz(-2.0531634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6338585) q[0];
sx q[0];
rz(-1.8554747) q[0];
sx q[0];
rz(2.741709) q[0];
rz(1.1625066) q[1];
sx q[1];
rz(-1.3299273) q[1];
sx q[1];
rz(-0.17366017) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0608873) q[0];
sx q[0];
rz(-1.6433435) q[0];
sx q[0];
rz(-1.5425496) q[0];
x q[1];
rz(2.6661751) q[2];
sx q[2];
rz(-0.43736514) q[2];
sx q[2];
rz(1.6944483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9160737) q[1];
sx q[1];
rz(-1.5991296) q[1];
sx q[1];
rz(1.6760992) q[1];
rz(-2.8097314) q[3];
sx q[3];
rz(-1.3734986) q[3];
sx q[3];
rz(-0.33533898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.48012039) q[2];
sx q[2];
rz(-1.9813333) q[2];
sx q[2];
rz(2.373467) q[2];
rz(2.2875732) q[3];
sx q[3];
rz(-1.7212399) q[3];
sx q[3];
rz(2.1693726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1059234) q[0];
sx q[0];
rz(-2.8978455) q[0];
sx q[0];
rz(1.6712028) q[0];
rz(0.51180965) q[1];
sx q[1];
rz(-2.6302331) q[1];
sx q[1];
rz(-1.9981729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91771942) q[0];
sx q[0];
rz(-1.6544764) q[0];
sx q[0];
rz(2.9907945) q[0];
rz(-pi) q[1];
rz(-2.5383699) q[2];
sx q[2];
rz(-1.8458741) q[2];
sx q[2];
rz(3.0657257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9547564) q[1];
sx q[1];
rz(-1.7446767) q[1];
sx q[1];
rz(-1.1980921) q[1];
rz(-pi) q[2];
rz(-1.0877803) q[3];
sx q[3];
rz(-0.59637585) q[3];
sx q[3];
rz(-3.0697825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4986971) q[2];
sx q[2];
rz(-0.84474793) q[2];
sx q[2];
rz(1.6112304) q[2];
rz(-1.4536084) q[3];
sx q[3];
rz(-1.0691103) q[3];
sx q[3];
rz(2.8924275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221508) q[0];
sx q[0];
rz(-2.3972153) q[0];
sx q[0];
rz(-1.7013593) q[0];
rz(-0.72921905) q[1];
sx q[1];
rz(-1.1518642) q[1];
sx q[1];
rz(-1.1332606) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5232552) q[0];
sx q[0];
rz(-1.528307) q[0];
sx q[0];
rz(2.8190814) q[0];
x q[1];
rz(-3.1256345) q[2];
sx q[2];
rz(-0.71231132) q[2];
sx q[2];
rz(2.6238837) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28753528) q[1];
sx q[1];
rz(-2.4234728) q[1];
sx q[1];
rz(1.992804) q[1];
rz(-pi) q[2];
rz(1.4134737) q[3];
sx q[3];
rz(-1.3120578) q[3];
sx q[3];
rz(1.3820005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40522727) q[2];
sx q[2];
rz(-2.0183125) q[2];
sx q[2];
rz(-2.6531632) q[2];
rz(1.3119665) q[3];
sx q[3];
rz(-0.043881504) q[3];
sx q[3];
rz(2.5002938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064780386) q[0];
sx q[0];
rz(-1.2925873) q[0];
sx q[0];
rz(3.0704165) q[0];
rz(-0.03216234) q[1];
sx q[1];
rz(-1.8036489) q[1];
sx q[1];
rz(-1.9326928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4874408) q[0];
sx q[0];
rz(-1.3432661) q[0];
sx q[0];
rz(2.3339416) q[0];
rz(-pi) q[1];
rz(-2.9572763) q[2];
sx q[2];
rz(-1.8039928) q[2];
sx q[2];
rz(1.2554982) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.68122411) q[1];
sx q[1];
rz(-1.7011233) q[1];
sx q[1];
rz(3.1254915) q[1];
rz(-pi) q[2];
rz(-2.5197221) q[3];
sx q[3];
rz(-1.8727881) q[3];
sx q[3];
rz(-2.5627476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9341087) q[2];
sx q[2];
rz(-0.19342962) q[2];
sx q[2];
rz(-1.5709546) q[2];
rz(-0.87336826) q[3];
sx q[3];
rz(-1.4040754) q[3];
sx q[3];
rz(-0.35203716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677925) q[0];
sx q[0];
rz(-1.6163102) q[0];
sx q[0];
rz(2.8299676) q[0];
rz(2.3198126) q[1];
sx q[1];
rz(-2.5506134) q[1];
sx q[1];
rz(1.6315546) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85743839) q[0];
sx q[0];
rz(-2.7630002) q[0];
sx q[0];
rz(-2.5749102) q[0];
rz(-pi) q[1];
rz(2.9456003) q[2];
sx q[2];
rz(-1.2336944) q[2];
sx q[2];
rz(0.54346426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.77126399) q[1];
sx q[1];
rz(-2.1762098) q[1];
sx q[1];
rz(1.3243115) q[1];
rz(-pi) q[2];
rz(1.3872836) q[3];
sx q[3];
rz(-1.7997051) q[3];
sx q[3];
rz(-0.51975091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2733549) q[2];
sx q[2];
rz(-0.53871012) q[2];
sx q[2];
rz(-0.81595016) q[2];
rz(2.6319035) q[3];
sx q[3];
rz(-2.3563801) q[3];
sx q[3];
rz(-1.2379439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8907392) q[0];
sx q[0];
rz(-1.4287404) q[0];
sx q[0];
rz(0.2302641) q[0];
rz(0.62581217) q[1];
sx q[1];
rz(-2.1964549) q[1];
sx q[1];
rz(2.4831916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9068245) q[0];
sx q[0];
rz(-1.0848197) q[0];
sx q[0];
rz(0.66396873) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41215956) q[2];
sx q[2];
rz(-2.4240652) q[2];
sx q[2];
rz(-1.621643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2645996) q[1];
sx q[1];
rz(-1.0942642) q[1];
sx q[1];
rz(-0.18696733) q[1];
x q[2];
rz(0.53604605) q[3];
sx q[3];
rz(-1.2671766) q[3];
sx q[3];
rz(0.21751465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.6752424) q[2];
sx q[2];
rz(-2.3190976) q[2];
sx q[2];
rz(0.33774439) q[2];
rz(2.1181469) q[3];
sx q[3];
rz(-1.3600391) q[3];
sx q[3];
rz(1.2158998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52453775) q[0];
sx q[0];
rz(-2.239997) q[0];
sx q[0];
rz(-0.064185113) q[0];
rz(1.2383923) q[1];
sx q[1];
rz(-1.7108142) q[1];
sx q[1];
rz(-1.6731813) q[1];
rz(-0.17910437) q[2];
sx q[2];
rz(-2.1894107) q[2];
sx q[2];
rz(-1.4084963) q[2];
rz(-1.6347377) q[3];
sx q[3];
rz(-0.99935525) q[3];
sx q[3];
rz(2.9730848) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
