OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0215065) q[0];
sx q[0];
rz(-1.9061631) q[0];
sx q[0];
rz(-1.0943227) q[0];
rz(-1.2878296) q[2];
sx q[2];
rz(-1.7093061) q[2];
sx q[2];
rz(1.1793009) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2687298) q[1];
sx q[1];
rz(-1.1464719) q[1];
sx q[1];
rz(-0.55117589) q[1];
x q[2];
rz(2.5039623) q[3];
sx q[3];
rz(-2.5507567) q[3];
sx q[3];
rz(1.5096111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(-0.93531936) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(-2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.4555567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6360977) q[0];
sx q[0];
rz(-2.539145) q[0];
sx q[0];
rz(-1.2732182) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7669719) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(0.74707109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.93745366) q[1];
sx q[1];
rz(-1.0480282) q[1];
sx q[1];
rz(0.015603113) q[1];
x q[2];
rz(0.060828408) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(-2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(-0.02877409) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(1.172539) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4370255) q[0];
sx q[0];
rz(-2.4798658) q[0];
sx q[0];
rz(2.6252803) q[0];
x q[1];
rz(0.34611361) q[2];
sx q[2];
rz(-0.78595224) q[2];
sx q[2];
rz(-1.9217938) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8923924) q[1];
sx q[1];
rz(-1.3274267) q[1];
sx q[1];
rz(-2.1057486) q[1];
rz(0.058733744) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(2.5854923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(0.89200154) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38055414) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(-0.52350837) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2098171) q[0];
sx q[0];
rz(-2.1486001) q[0];
sx q[0];
rz(-2.0156167) q[0];
rz(-pi) q[1];
rz(-1.4857616) q[2];
sx q[2];
rz(-1.2418613) q[2];
sx q[2];
rz(2.0563682) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47141155) q[1];
sx q[1];
rz(-1.3172611) q[1];
sx q[1];
rz(-0.86004911) q[1];
rz(0.97949667) q[3];
sx q[3];
rz(-1.6846091) q[3];
sx q[3];
rz(2.1554961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(0.28856746) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.6500641) q[0];
rz(-2.2619757) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0846227) q[0];
sx q[0];
rz(-1.3966494) q[0];
sx q[0];
rz(-1.6790381) q[0];
x q[1];
rz(-1.5993824) q[2];
sx q[2];
rz(-2.8386142) q[2];
sx q[2];
rz(2.6945393) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3757513) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(2.8326616) q[1];
rz(1.4335853) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(-1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0483635) q[0];
sx q[0];
rz(-1.7943802) q[0];
sx q[0];
rz(2.6327052) q[0];
rz(-0.09128696) q[2];
sx q[2];
rz(-1.316615) q[2];
sx q[2];
rz(-2.9448178) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9565935) q[1];
sx q[1];
rz(-1.5041659) q[1];
sx q[1];
rz(1.1207629) q[1];
x q[2];
rz(2.1720042) q[3];
sx q[3];
rz(-1.117327) q[3];
sx q[3];
rz(2.3525402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5086223) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-2.3197876) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069698378) q[0];
sx q[0];
rz(-1.0487723) q[0];
sx q[0];
rz(-2.4126023) q[0];
x q[1];
rz(-2.1173382) q[2];
sx q[2];
rz(-2.3356236) q[2];
sx q[2];
rz(-0.96217996) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6268839) q[1];
sx q[1];
rz(-0.85299546) q[1];
sx q[1];
rz(-1.8741329) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.20007) q[3];
sx q[3];
rz(-1.0154187) q[3];
sx q[3];
rz(-2.6851418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41641763) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(2.3186671) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.8364505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7628521) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(2.5255894) q[0];
x q[1];
rz(2.7140076) q[2];
sx q[2];
rz(-2.3901849) q[2];
sx q[2];
rz(-2.4497355) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9040363) q[1];
sx q[1];
rz(-0.95103969) q[1];
sx q[1];
rz(1.1035641) q[1];
rz(-pi) q[2];
rz(0.023530258) q[3];
sx q[3];
rz(-1.9100034) q[3];
sx q[3];
rz(2.6587405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1398853) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.4902327) q[2];
rz(1.0772609) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(-0.3219147) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-0.70294356) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87899938) q[0];
sx q[0];
rz(-1.2125373) q[0];
sx q[0];
rz(-2.7260289) q[0];
rz(-pi) q[1];
rz(1.0614971) q[2];
sx q[2];
rz(-1.6054389) q[2];
sx q[2];
rz(-2.408037) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.34300464) q[1];
sx q[1];
rz(-1.7354021) q[1];
sx q[1];
rz(2.1144457) q[1];
rz(-pi) q[2];
rz(2.1861107) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(0.51481065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90074173) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(1.8019603) q[2];
rz(2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(-0.46863619) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7077431) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(-1.3462523) q[0];
x q[1];
rz(-2.2655728) q[2];
sx q[2];
rz(-2.6555736) q[2];
sx q[2];
rz(-1.3181869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.48060265) q[1];
sx q[1];
rz(-0.27462474) q[1];
sx q[1];
rz(1.8250699) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81810276) q[3];
sx q[3];
rz(-2.5460498) q[3];
sx q[3];
rz(0.18225741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(-0.11463595) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4951915) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-2.2767699) q[2];
sx q[2];
rz(-1.0100126) q[2];
sx q[2];
rz(1.0457912) q[2];
rz(-0.078483742) q[3];
sx q[3];
rz(-0.92072903) q[3];
sx q[3];
rz(-2.846684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];