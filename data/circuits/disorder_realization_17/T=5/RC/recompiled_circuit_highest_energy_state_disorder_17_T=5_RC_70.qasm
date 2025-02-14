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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(-2.4315779) q[0];
rz(2.8275936) q[1];
sx q[1];
rz(-2.2065838) q[1];
sx q[1];
rz(1.3318055) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074866991) q[0];
sx q[0];
rz(-2.7829889) q[0];
sx q[0];
rz(2.4551366) q[0];
x q[1];
rz(0.34320404) q[2];
sx q[2];
rz(-1.8095952) q[2];
sx q[2];
rz(-0.78087872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1430681) q[1];
sx q[1];
rz(-1.8170763) q[1];
sx q[1];
rz(-2.1889436) q[1];
x q[2];
rz(-0.93721849) q[3];
sx q[3];
rz(-1.211418) q[3];
sx q[3];
rz(2.0546953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72267246) q[2];
sx q[2];
rz(-1.924943) q[2];
sx q[2];
rz(0.81532064) q[2];
rz(-0.80960649) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0631436) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(0.94509205) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5623215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54540578) q[0];
sx q[0];
rz(-0.89841398) q[0];
sx q[0];
rz(2.5491712) q[0];
rz(-pi) q[1];
rz(0.42065545) q[2];
sx q[2];
rz(-1.1569945) q[2];
sx q[2];
rz(2.5977573) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89797276) q[1];
sx q[1];
rz(-2.0531516) q[1];
sx q[1];
rz(-1.4528758) q[1];
x q[2];
rz(-0.48247997) q[3];
sx q[3];
rz(-2.2594995) q[3];
sx q[3];
rz(0.83079332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.14757806) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(2.6914524) q[2];
rz(-1.133793) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(-0.97761893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56871539) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(-1.7549134) q[1];
sx q[1];
rz(-1.7117056) q[1];
sx q[1];
rz(-0.6001572) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48334405) q[0];
sx q[0];
rz(-1.5243825) q[0];
sx q[0];
rz(-0.29057403) q[0];
rz(-pi) q[1];
rz(0.05608989) q[2];
sx q[2];
rz(-2.1768835) q[2];
sx q[2];
rz(-2.1457399) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8474947) q[1];
sx q[1];
rz(-2.1408242) q[1];
sx q[1];
rz(-1.9600201) q[1];
rz(-pi) q[2];
rz(-1.3998109) q[3];
sx q[3];
rz(-1.0043084) q[3];
sx q[3];
rz(-1.8734219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0967789) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(-2.0951994) q[2];
rz(-0.3977631) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(-3.0090295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9699049) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(-2.0890253) q[0];
rz(1.196208) q[1];
sx q[1];
rz(-1.9320678) q[1];
sx q[1];
rz(-0.41950163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2118111) q[0];
sx q[0];
rz(-2.503184) q[0];
sx q[0];
rz(2.6891461) q[0];
x q[1];
rz(0.93650903) q[2];
sx q[2];
rz(-1.3872996) q[2];
sx q[2];
rz(-1.2681792) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4070523) q[1];
sx q[1];
rz(-2.190935) q[1];
sx q[1];
rz(2.6728515) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.64873) q[3];
sx q[3];
rz(-2.5774084) q[3];
sx q[3];
rz(-1.9436556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18998751) q[2];
sx q[2];
rz(-1.1185458) q[2];
sx q[2];
rz(1.0901964) q[2];
rz(0.53127855) q[3];
sx q[3];
rz(-0.71610206) q[3];
sx q[3];
rz(-1.5268911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7200318) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(1.74362) q[0];
rz(0.13294237) q[1];
sx q[1];
rz(-1.839919) q[1];
sx q[1];
rz(-0.001210777) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069790445) q[0];
sx q[0];
rz(-1.9513357) q[0];
sx q[0];
rz(-1.5881722) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0693477) q[2];
sx q[2];
rz(-2.8184888) q[2];
sx q[2];
rz(-2.1954775) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7442786) q[1];
sx q[1];
rz(-2.4220253) q[1];
sx q[1];
rz(-0.47081635) q[1];
x q[2];
rz(1.4210971) q[3];
sx q[3];
rz(-1.5970917) q[3];
sx q[3];
rz(2.5291075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2289537) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(-0.21052989) q[2];
rz(-1.0362961) q[3];
sx q[3];
rz(-2.3573037) q[3];
sx q[3];
rz(-1.2150631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9482816) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(-0.40570983) q[0];
rz(1.7871208) q[1];
sx q[1];
rz(-2.3190119) q[1];
sx q[1];
rz(2.2192661) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852823) q[0];
sx q[0];
rz(-1.4643702) q[0];
sx q[0];
rz(-1.2805416) q[0];
x q[1];
rz(-0.10293691) q[2];
sx q[2];
rz(-2.2349226) q[2];
sx q[2];
rz(-0.31598202) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2353224) q[1];
sx q[1];
rz(-0.962469) q[1];
sx q[1];
rz(0.13923213) q[1];
rz(-pi) q[2];
rz(-1.9590553) q[3];
sx q[3];
rz(-1.045533) q[3];
sx q[3];
rz(-1.0356667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(-2.7396438) q[2];
rz(-1.8534144) q[3];
sx q[3];
rz(-0.86789075) q[3];
sx q[3];
rz(-1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61889082) q[0];
sx q[0];
rz(-1.1033449) q[0];
sx q[0];
rz(-0.064362137) q[0];
rz(-2.3035658) q[1];
sx q[1];
rz(-2.3152654) q[1];
sx q[1];
rz(-1.3305957) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42444077) q[0];
sx q[0];
rz(-1.0006051) q[0];
sx q[0];
rz(1.7198404) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0688017) q[2];
sx q[2];
rz(-1.0235707) q[2];
sx q[2];
rz(2.4233873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8992675) q[1];
sx q[1];
rz(-2.0914255) q[1];
sx q[1];
rz(-1.8249056) q[1];
rz(-pi) q[2];
rz(1.3279541) q[3];
sx q[3];
rz(-0.62627072) q[3];
sx q[3];
rz(1.0427208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7360721) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(-2.4580477) q[2];
rz(-0.73379597) q[3];
sx q[3];
rz(-1.7283231) q[3];
sx q[3];
rz(0.84764135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7965294) q[0];
sx q[0];
rz(-2.0579484) q[0];
sx q[0];
rz(1.9684568) q[0];
rz(0.37551156) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(1.0367905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6818559) q[0];
sx q[0];
rz(-1.5899183) q[0];
sx q[0];
rz(-0.023650344) q[0];
x q[1];
rz(1.0388184) q[2];
sx q[2];
rz(-0.3178645) q[2];
sx q[2];
rz(1.7828072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6822676) q[1];
sx q[1];
rz(-2.3258491) q[1];
sx q[1];
rz(1.7641279) q[1];
x q[2];
rz(-2.310318) q[3];
sx q[3];
rz(-2.4711802) q[3];
sx q[3];
rz(1.1408653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.57572395) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(-2.4110528) q[2];
rz(-2.9324487) q[3];
sx q[3];
rz(-0.52339619) q[3];
sx q[3];
rz(1.2701579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.65649477) q[0];
sx q[0];
rz(-2.7182343) q[0];
sx q[0];
rz(-0.04743162) q[0];
rz(-2.9196396) q[1];
sx q[1];
rz(-1.0739948) q[1];
sx q[1];
rz(2.2304631) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.52235) q[0];
sx q[0];
rz(-1.5372835) q[0];
sx q[0];
rz(1.5977809) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5510727) q[2];
sx q[2];
rz(-0.64177401) q[2];
sx q[2];
rz(2.3373147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4621689) q[1];
sx q[1];
rz(-2.7738681) q[1];
sx q[1];
rz(-3.1064139) q[1];
x q[2];
rz(2.1371331) q[3];
sx q[3];
rz(-1.3977504) q[3];
sx q[3];
rz(-2.6030948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2719416) q[2];
sx q[2];
rz(-0.38001529) q[2];
sx q[2];
rz(-0.16858777) q[2];
rz(0.21909675) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(-1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1472226) q[0];
sx q[0];
rz(-2.8012025) q[0];
sx q[0];
rz(-0.45743531) q[0];
rz(1.9153197) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(-1.7785243) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5209893) q[0];
sx q[0];
rz(-2.531157) q[0];
sx q[0];
rz(3.0362398) q[0];
rz(2.5117433) q[2];
sx q[2];
rz(-1.7694663) q[2];
sx q[2];
rz(-2.3185286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5879092) q[1];
sx q[1];
rz(-1.739555) q[1];
sx q[1];
rz(-2.7588899) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4665589) q[3];
sx q[3];
rz(-0.6851894) q[3];
sx q[3];
rz(-2.4327421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3739796) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(-2.6895831) q[2];
rz(2.7376145) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(1.1261136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6954738) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(2.6192464) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(-0.28235565) q[2];
sx q[2];
rz(-1.8179802) q[2];
sx q[2];
rz(2.5674934) q[2];
rz(2.6330144) q[3];
sx q[3];
rz(-2.8420035) q[3];
sx q[3];
rz(-3.0723078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
