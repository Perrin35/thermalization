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
rz(-2.5341852) q[0];
sx q[0];
rz(-2.6976801) q[0];
sx q[0];
rz(2.8191415) q[0];
rz(-1.9880265) q[1];
sx q[1];
rz(-1.1918951) q[1];
sx q[1];
rz(0.21790394) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76989749) q[0];
sx q[0];
rz(-1.2697233) q[0];
sx q[0];
rz(-0.45905827) q[0];
rz(-2.430021) q[2];
sx q[2];
rz(-0.24453881) q[2];
sx q[2];
rz(-2.0681579) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7806381) q[1];
sx q[1];
rz(-1.8285969) q[1];
sx q[1];
rz(-1.6901438) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2530768) q[3];
sx q[3];
rz(-1.2938283) q[3];
sx q[3];
rz(-0.4539629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0107062) q[2];
sx q[2];
rz(-0.43121269) q[2];
sx q[2];
rz(2.9475589) q[2];
rz(-0.38328299) q[3];
sx q[3];
rz(-0.88519874) q[3];
sx q[3];
rz(-1.6398199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5782769) q[0];
sx q[0];
rz(-1.0551772) q[0];
sx q[0];
rz(1.1042327) q[0];
rz(2.4597994) q[1];
sx q[1];
rz(-1.3673404) q[1];
sx q[1];
rz(1.7368447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8786312) q[0];
sx q[0];
rz(-1.6978953) q[0];
sx q[0];
rz(-0.2574347) q[0];
rz(-pi) q[1];
rz(2.4787729) q[2];
sx q[2];
rz(-2.347555) q[2];
sx q[2];
rz(2.8341849) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42336418) q[1];
sx q[1];
rz(-1.2585139) q[1];
sx q[1];
rz(-0.52257089) q[1];
x q[2];
rz(1.4561171) q[3];
sx q[3];
rz(-1.8240098) q[3];
sx q[3];
rz(-0.51451433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76689395) q[2];
sx q[2];
rz(-2.1465492) q[2];
sx q[2];
rz(1.9083171) q[2];
rz(0.79902664) q[3];
sx q[3];
rz(-0.34501758) q[3];
sx q[3];
rz(0.12921216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1610573) q[0];
sx q[0];
rz(-1.072847) q[0];
sx q[0];
rz(2.5824353) q[0];
rz(2.0237538) q[1];
sx q[1];
rz(-1.585377) q[1];
sx q[1];
rz(0.1112172) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98830279) q[0];
sx q[0];
rz(-0.39228253) q[0];
sx q[0];
rz(-2.8610703) q[0];
rz(0.36685852) q[2];
sx q[2];
rz(-2.48403) q[2];
sx q[2];
rz(2.2012432) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7349941) q[1];
sx q[1];
rz(-1.3822228) q[1];
sx q[1];
rz(-2.2833634) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7484557) q[3];
sx q[3];
rz(-1.434978) q[3];
sx q[3];
rz(-0.47752646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3126276) q[2];
sx q[2];
rz(-1.6964922) q[2];
sx q[2];
rz(1.9435389) q[2];
rz(-1.4231921) q[3];
sx q[3];
rz(-1.6953902) q[3];
sx q[3];
rz(-0.47553441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6296185) q[0];
sx q[0];
rz(-2.8422575) q[0];
sx q[0];
rz(2.3053115) q[0];
rz(-1.5760999) q[1];
sx q[1];
rz(-1.0944159) q[1];
sx q[1];
rz(-0.85087585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5880575) q[0];
sx q[0];
rz(-1.6041821) q[0];
sx q[0];
rz(-2.1436611) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0298898) q[2];
sx q[2];
rz(-2.0238681) q[2];
sx q[2];
rz(-2.0889995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8703813) q[1];
sx q[1];
rz(-2.3496501) q[1];
sx q[1];
rz(0.19510896) q[1];
rz(-1.1034043) q[3];
sx q[3];
rz(-2.7041743) q[3];
sx q[3];
rz(1.9269772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0389082) q[2];
sx q[2];
rz(-0.85997283) q[2];
sx q[2];
rz(2.6343708) q[2];
rz(-0.66633362) q[3];
sx q[3];
rz(-2.3907876) q[3];
sx q[3];
rz(0.91608086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12793334) q[0];
sx q[0];
rz(-1.6935231) q[0];
sx q[0];
rz(-1.5284982) q[0];
rz(1.5616034) q[1];
sx q[1];
rz(-2.0785619) q[1];
sx q[1];
rz(-2.9772421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30475475) q[0];
sx q[0];
rz(-1.4753818) q[0];
sx q[0];
rz(1.4897904) q[0];
x q[1];
rz(-2.3037203) q[2];
sx q[2];
rz(-2.3114486) q[2];
sx q[2];
rz(1.6760507) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4867465) q[1];
sx q[1];
rz(-2.2426412) q[1];
sx q[1];
rz(-0.44507546) q[1];
rz(-pi) q[2];
rz(1.8643531) q[3];
sx q[3];
rz(-1.3428215) q[3];
sx q[3];
rz(-0.10693947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5932172) q[2];
sx q[2];
rz(-2.2943353) q[2];
sx q[2];
rz(0.66620052) q[2];
rz(-0.14144746) q[3];
sx q[3];
rz(-1.5513159) q[3];
sx q[3];
rz(-0.9945873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2675466) q[0];
sx q[0];
rz(-2.5135437) q[0];
sx q[0];
rz(0.48126599) q[0];
rz(0.73356837) q[1];
sx q[1];
rz(-1.8736519) q[1];
sx q[1];
rz(3.0487294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5294612) q[0];
sx q[0];
rz(-2.4680637) q[0];
sx q[0];
rz(-2.3271534) q[0];
rz(3.1029557) q[2];
sx q[2];
rz(-1.5687571) q[2];
sx q[2];
rz(0.80261999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9290849) q[1];
sx q[1];
rz(-1.9783535) q[1];
sx q[1];
rz(0.69504884) q[1];
rz(-pi) q[2];
rz(-2.8792297) q[3];
sx q[3];
rz(-1.922058) q[3];
sx q[3];
rz(-0.92811229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7218472) q[2];
sx q[2];
rz(-1.6215308) q[2];
sx q[2];
rz(1.2255555) q[2];
rz(0.078977481) q[3];
sx q[3];
rz(-2.0073399) q[3];
sx q[3];
rz(-1.8592161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52800286) q[0];
sx q[0];
rz(-0.71857518) q[0];
sx q[0];
rz(-0.80793107) q[0];
rz(-1.6174053) q[1];
sx q[1];
rz(-2.9831191) q[1];
sx q[1];
rz(-2.4375516) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28939279) q[0];
sx q[0];
rz(-1.5796229) q[0];
sx q[0];
rz(1.5509964) q[0];
x q[1];
rz(0.0041517082) q[2];
sx q[2];
rz(-1.5791681) q[2];
sx q[2];
rz(-2.3821751) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6345811) q[1];
sx q[1];
rz(-1.5935362) q[1];
sx q[1];
rz(-2.1918767) q[1];
x q[2];
rz(0.95583393) q[3];
sx q[3];
rz(-1.7059688) q[3];
sx q[3];
rz(-1.9036024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.070270553) q[2];
sx q[2];
rz(-2.4532048) q[2];
sx q[2];
rz(-2.2570611) q[2];
rz(-2.6279348) q[3];
sx q[3];
rz(-1.1762041) q[3];
sx q[3];
rz(2.6562712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1270444) q[0];
sx q[0];
rz(-2.3764648) q[0];
sx q[0];
rz(2.2807518) q[0];
rz(-0.035482081) q[1];
sx q[1];
rz(-1.9396962) q[1];
sx q[1];
rz(-1.6532345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4299773) q[0];
sx q[0];
rz(-1.1722783) q[0];
sx q[0];
rz(-2.7048955) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8143009) q[2];
sx q[2];
rz(-1.3483682) q[2];
sx q[2];
rz(0.13130638) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4379641) q[1];
sx q[1];
rz(-1.1878106) q[1];
sx q[1];
rz(-0.18937892) q[1];
rz(-pi) q[2];
x q[2];
rz(2.866167) q[3];
sx q[3];
rz(-1.4793494) q[3];
sx q[3];
rz(-0.25642727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52582467) q[2];
sx q[2];
rz(-2.0066228) q[2];
sx q[2];
rz(-1.6305249) q[2];
rz(0.75119558) q[3];
sx q[3];
rz(-0.54371756) q[3];
sx q[3];
rz(-2.3403919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57665956) q[0];
sx q[0];
rz(-1.3929921) q[0];
sx q[0];
rz(-2.9826953) q[0];
rz(1.502602) q[1];
sx q[1];
rz(-1.8105806) q[1];
sx q[1];
rz(1.166689) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.027307) q[0];
sx q[0];
rz(-2.5059359) q[0];
sx q[0];
rz(0.42814769) q[0];
x q[1];
rz(2.8620637) q[2];
sx q[2];
rz(-1.88481) q[2];
sx q[2];
rz(-1.90998) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.78913375) q[1];
sx q[1];
rz(-0.041445565) q[1];
sx q[1];
rz(1.9876473) q[1];
x q[2];
rz(1.6004815) q[3];
sx q[3];
rz(-1.5497713) q[3];
sx q[3];
rz(2.0898745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73929536) q[2];
sx q[2];
rz(-0.6604971) q[2];
sx q[2];
rz(-0.26563773) q[2];
rz(1.9499251) q[3];
sx q[3];
rz(-1.8040413) q[3];
sx q[3];
rz(-2.5875097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3091076) q[0];
sx q[0];
rz(-2.5256248) q[0];
sx q[0];
rz(-0.68514222) q[0];
rz(-2.5849672) q[1];
sx q[1];
rz(-1.4645422) q[1];
sx q[1];
rz(-2.305078) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3292824) q[0];
sx q[0];
rz(-0.54566452) q[0];
sx q[0];
rz(-0.87581234) q[0];
rz(2.0903953) q[2];
sx q[2];
rz(-1.8406036) q[2];
sx q[2];
rz(-2.4003911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3565607) q[1];
sx q[1];
rz(-1.1379269) q[1];
sx q[1];
rz(0.69286107) q[1];
x q[2];
rz(0.1565069) q[3];
sx q[3];
rz(-1.8409074) q[3];
sx q[3];
rz(-2.6685985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69540018) q[2];
sx q[2];
rz(-1.7170898) q[2];
sx q[2];
rz(0.85644537) q[2];
rz(-0.19927464) q[3];
sx q[3];
rz(-2.1592185) q[3];
sx q[3];
rz(-0.91163951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0675426) q[0];
sx q[0];
rz(-2.1265125) q[0];
sx q[0];
rz(1.4763005) q[0];
rz(-0.53786565) q[1];
sx q[1];
rz(-1.8712416) q[1];
sx q[1];
rz(1.3457294) q[1];
rz(-2.3026148) q[2];
sx q[2];
rz(-2.2175023) q[2];
sx q[2];
rz(-1.6432696) q[2];
rz(-0.27502455) q[3];
sx q[3];
rz(-2.2755819) q[3];
sx q[3];
rz(-1.6439846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
