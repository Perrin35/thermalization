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
rz(-2.383411) q[0];
sx q[0];
rz(-0.36769205) q[0];
sx q[0];
rz(2.0453069) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(-2.0356324) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5263514) q[0];
sx q[0];
rz(-2.5509074) q[0];
sx q[0];
rz(-2.3870941) q[0];
x q[1];
rz(-2.8739086) q[2];
sx q[2];
rz(-1.4776728) q[2];
sx q[2];
rz(-0.47506079) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3333862) q[1];
sx q[1];
rz(-1.5667586) q[1];
sx q[1];
rz(-1.6081078) q[1];
rz(-pi) q[2];
x q[2];
rz(0.023882377) q[3];
sx q[3];
rz(-2.8571354) q[3];
sx q[3];
rz(-2.6989355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7652863) q[2];
sx q[2];
rz(-0.47069612) q[2];
sx q[2];
rz(0.37962309) q[2];
rz(-1.6342573) q[3];
sx q[3];
rz(-2.0416656) q[3];
sx q[3];
rz(2.1716993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955502) q[0];
sx q[0];
rz(-0.33691418) q[0];
sx q[0];
rz(0.90091339) q[0];
rz(0.13149978) q[1];
sx q[1];
rz(-1.9410746) q[1];
sx q[1];
rz(-0.44201717) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95831224) q[0];
sx q[0];
rz(-1.4889476) q[0];
sx q[0];
rz(-0.32572066) q[0];
x q[1];
rz(-1.5698264) q[2];
sx q[2];
rz(-1.5754149) q[2];
sx q[2];
rz(2.6498891) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.17168853) q[1];
sx q[1];
rz(-2.0398047) q[1];
sx q[1];
rz(2.5086839) q[1];
rz(2.7412299) q[3];
sx q[3];
rz(-1.6737079) q[3];
sx q[3];
rz(-2.3617552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3758731) q[2];
sx q[2];
rz(-1.4292382) q[2];
sx q[2];
rz(2.5326552) q[2];
rz(-0.40306148) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(-0.5016996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32560638) q[0];
sx q[0];
rz(-2.7222962) q[0];
sx q[0];
rz(2.0035279) q[0];
rz(-1.5886547) q[1];
sx q[1];
rz(-0.76858968) q[1];
sx q[1];
rz(2.3345711) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3024705) q[0];
sx q[0];
rz(-3.0387276) q[0];
sx q[0];
rz(-0.26786049) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4078232) q[2];
sx q[2];
rz(-2.0103243) q[2];
sx q[2];
rz(1.7874315) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0453607) q[1];
sx q[1];
rz(-1.3825455) q[1];
sx q[1];
rz(1.4493311) q[1];
rz(1.9495973) q[3];
sx q[3];
rz(-1.8378432) q[3];
sx q[3];
rz(-0.35940816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1300065) q[2];
sx q[2];
rz(-2.5277972) q[2];
sx q[2];
rz(1.2137132) q[2];
rz(-0.064149292) q[3];
sx q[3];
rz(-1.4870653) q[3];
sx q[3];
rz(-3.1411662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59378687) q[0];
sx q[0];
rz(-1.2603899) q[0];
sx q[0];
rz(-2.2547145) q[0];
rz(-0.15874323) q[1];
sx q[1];
rz(-0.49247772) q[1];
sx q[1];
rz(1.278272) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072479362) q[0];
sx q[0];
rz(-2.2769391) q[0];
sx q[0];
rz(-0.98346498) q[0];
x q[1];
rz(1.099894) q[2];
sx q[2];
rz(-1.0702025) q[2];
sx q[2];
rz(-0.72906993) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1310744) q[1];
sx q[1];
rz(-1.5107039) q[1];
sx q[1];
rz(0.87803234) q[1];
rz(-pi) q[2];
rz(-1.3298755) q[3];
sx q[3];
rz(-2.5735825) q[3];
sx q[3];
rz(0.78838511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13901237) q[2];
sx q[2];
rz(-0.85550344) q[2];
sx q[2];
rz(0.95376897) q[2];
rz(1.9468797) q[3];
sx q[3];
rz(-2.298893) q[3];
sx q[3];
rz(0.4755303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139528) q[0];
sx q[0];
rz(-0.71641818) q[0];
sx q[0];
rz(0.59999505) q[0];
rz(-0.68296877) q[1];
sx q[1];
rz(-2.3536847) q[1];
sx q[1];
rz(-2.8111828) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70008695) q[0];
sx q[0];
rz(-1.7528025) q[0];
sx q[0];
rz(-0.84199961) q[0];
x q[1];
rz(-0.16322132) q[2];
sx q[2];
rz(-1.7931869) q[2];
sx q[2];
rz(1.4537653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7327795) q[1];
sx q[1];
rz(-0.52460743) q[1];
sx q[1];
rz(-1.9636159) q[1];
x q[2];
rz(-2.7269269) q[3];
sx q[3];
rz(-2.2288481) q[3];
sx q[3];
rz(2.0882704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7052475) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(2.516563) q[2];
rz(-0.69508067) q[3];
sx q[3];
rz(-1.2806226) q[3];
sx q[3];
rz(-3.0888016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2136252) q[0];
sx q[0];
rz(-1.1518421) q[0];
sx q[0];
rz(1.7684162) q[0];
rz(1.7534509) q[1];
sx q[1];
rz(-0.87166798) q[1];
sx q[1];
rz(1.6067827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4563303) q[0];
sx q[0];
rz(-1.7182516) q[0];
sx q[0];
rz(-0.035650226) q[0];
x q[1];
rz(2.0892145) q[2];
sx q[2];
rz(-0.3568584) q[2];
sx q[2];
rz(-1.269358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.116356) q[1];
sx q[1];
rz(-2.3400809) q[1];
sx q[1];
rz(-0.58057433) q[1];
x q[2];
rz(-3.0032773) q[3];
sx q[3];
rz(-1.173291) q[3];
sx q[3];
rz(-2.60633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9390949) q[2];
sx q[2];
rz(-1.9715344) q[2];
sx q[2];
rz(-1.5411752) q[2];
rz(2.1206858) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(-0.30103621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0114667) q[0];
sx q[0];
rz(-3.0045894) q[0];
sx q[0];
rz(2.2504508) q[0];
rz(2.4961684) q[1];
sx q[1];
rz(-1.3963457) q[1];
sx q[1];
rz(-2.3088764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4629991) q[0];
sx q[0];
rz(-2.6581253) q[0];
sx q[0];
rz(-2.8449758) q[0];
rz(1.4786167) q[2];
sx q[2];
rz(-0.62609172) q[2];
sx q[2];
rz(0.41958671) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11800217) q[1];
sx q[1];
rz(-2.0609792) q[1];
sx q[1];
rz(1.791819) q[1];
rz(-pi) q[2];
rz(-1.9644968) q[3];
sx q[3];
rz(-2.0767005) q[3];
sx q[3];
rz(-0.33123744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.58925313) q[2];
sx q[2];
rz(-0.66698843) q[2];
sx q[2];
rz(-1.5935295) q[2];
rz(-0.42406905) q[3];
sx q[3];
rz(-1.721902) q[3];
sx q[3];
rz(-0.29485318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9645204) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(0.82343423) q[0];
rz(1.9442762) q[1];
sx q[1];
rz(-1.6540534) q[1];
sx q[1];
rz(-1.6808602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4005918) q[0];
sx q[0];
rz(-1.680458) q[0];
sx q[0];
rz(1.8909258) q[0];
rz(-pi) q[1];
rz(-0.16847289) q[2];
sx q[2];
rz(-2.1572067) q[2];
sx q[2];
rz(-2.7975067) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.15732652) q[1];
sx q[1];
rz(-1.5283902) q[1];
sx q[1];
rz(0.66161109) q[1];
x q[2];
rz(2.7815205) q[3];
sx q[3];
rz(-1.8723893) q[3];
sx q[3];
rz(2.2693279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0494277) q[2];
sx q[2];
rz(-0.66114134) q[2];
sx q[2];
rz(-0.1235505) q[2];
rz(-2.3030247) q[3];
sx q[3];
rz(-2.0288012) q[3];
sx q[3];
rz(-0.53028321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0274444) q[0];
sx q[0];
rz(-2.1338978) q[0];
sx q[0];
rz(1.7281519) q[0];
rz(0.89883262) q[1];
sx q[1];
rz(-0.90679589) q[1];
sx q[1];
rz(-2.5487505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4384189) q[0];
sx q[0];
rz(-0.98383622) q[0];
sx q[0];
rz(-0.86721768) q[0];
x q[1];
rz(-1.4845092) q[2];
sx q[2];
rz(-1.3893428) q[2];
sx q[2];
rz(-1.4149416) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.457691) q[1];
sx q[1];
rz(-1.8967034) q[1];
sx q[1];
rz(-1.1363455) q[1];
x q[2];
rz(2.1778706) q[3];
sx q[3];
rz(-0.98120171) q[3];
sx q[3];
rz(-2.251925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8574519) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(2.9202374) q[2];
rz(-1.0434693) q[3];
sx q[3];
rz(-0.9404434) q[3];
sx q[3];
rz(0.38844696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8429883) q[0];
sx q[0];
rz(-0.28268155) q[0];
sx q[0];
rz(0.84841949) q[0];
rz(0.09659718) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(2.3883147) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8814113) q[0];
sx q[0];
rz(-0.24883606) q[0];
sx q[0];
rz(1.1538065) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3123872) q[2];
sx q[2];
rz(-1.5463788) q[2];
sx q[2];
rz(-1.3361734) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92434525) q[1];
sx q[1];
rz(-2.0567523) q[1];
sx q[1];
rz(0.43083453) q[1];
rz(-pi) q[2];
rz(-0.71452272) q[3];
sx q[3];
rz(-1.2163289) q[3];
sx q[3];
rz(2.9936341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0640556) q[2];
sx q[2];
rz(-2.3040743) q[2];
sx q[2];
rz(2.688664) q[2];
rz(0.60106599) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(-2.1452904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30408981) q[0];
sx q[0];
rz(-1.6927728) q[0];
sx q[0];
rz(0.46020831) q[0];
rz(2.9651463) q[1];
sx q[1];
rz(-2.9063168) q[1];
sx q[1];
rz(1.6520687) q[1];
rz(3.1398793) q[2];
sx q[2];
rz(-0.34863498) q[2];
sx q[2];
rz(0.44193535) q[2];
rz(2.413977) q[3];
sx q[3];
rz(-1.1066827) q[3];
sx q[3];
rz(-0.91409693) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
