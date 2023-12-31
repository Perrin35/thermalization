OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3044843) q[0];
sx q[0];
rz(-1.6882856) q[0];
sx q[0];
rz(-0.31153554) q[0];
rz(-0.43752924) q[1];
sx q[1];
rz(-1.8234) q[1];
sx q[1];
rz(-2.5826366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564724) q[0];
sx q[0];
rz(-2.7010242) q[0];
sx q[0];
rz(1.23929) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5582325) q[2];
sx q[2];
rz(-2.1280834) q[2];
sx q[2];
rz(-1.81665) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.641687) q[1];
sx q[1];
rz(-0.99398621) q[1];
sx q[1];
rz(-0.66024248) q[1];
rz(-pi) q[2];
rz(-0.55922237) q[3];
sx q[3];
rz(-0.50563522) q[3];
sx q[3];
rz(-1.3666183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7493593) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(-0.63670811) q[2];
rz(0.84896815) q[3];
sx q[3];
rz(-0.62148062) q[3];
sx q[3];
rz(0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37671509) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(2.9887181) q[0];
rz(-2.3846467) q[1];
sx q[1];
rz(-1.5870973) q[1];
sx q[1];
rz(2.1551932) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9372285) q[0];
sx q[0];
rz(-1.5110656) q[0];
sx q[0];
rz(1.6110957) q[0];
rz(-1.0436922) q[2];
sx q[2];
rz(-2.2615848) q[2];
sx q[2];
rz(1.8156798) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6978554) q[1];
sx q[1];
rz(-1.9471696) q[1];
sx q[1];
rz(-1.8259551) q[1];
rz(-pi) q[2];
rz(-1.8996703) q[3];
sx q[3];
rz(-1.6992913) q[3];
sx q[3];
rz(-2.0957029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.219517) q[2];
sx q[2];
rz(-2.3584649) q[2];
rz(3.1230208) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(-2.1799178) q[0];
rz(-2.7812474) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(-3.0128984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0823682) q[0];
sx q[0];
rz(-1.6186065) q[0];
sx q[0];
rz(-1.6533018) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80758904) q[2];
sx q[2];
rz(-1.979504) q[2];
sx q[2];
rz(-1.7847716) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8752746) q[1];
sx q[1];
rz(-2.6058063) q[1];
sx q[1];
rz(-2.5712625) q[1];
x q[2];
rz(1.4314753) q[3];
sx q[3];
rz(-2.3686667) q[3];
sx q[3];
rz(-1.3611925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0814357) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(-1.241768) q[2];
rz(-2.5545819) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(-0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.509165) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(-1.084491) q[0];
rz(1.4831316) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(-0.09253563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70150347) q[0];
sx q[0];
rz(-1.2883696) q[0];
sx q[0];
rz(-0.28042067) q[0];
rz(-pi) q[1];
rz(-2.7961568) q[2];
sx q[2];
rz(-1.1159117) q[2];
sx q[2];
rz(-2.5330184) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42379984) q[1];
sx q[1];
rz(-1.6289627) q[1];
sx q[1];
rz(0.34323005) q[1];
rz(1.362364) q[3];
sx q[3];
rz(-2.4760893) q[3];
sx q[3];
rz(-0.044737577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(2.8095424) q[2];
rz(1.0559233) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32245359) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(0.21155587) q[0];
rz(-1.3062723) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(0.64770118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0928597) q[0];
sx q[0];
rz(-1.6824241) q[0];
sx q[0];
rz(2.9812921) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7933502) q[2];
sx q[2];
rz(-1.4255382) q[2];
sx q[2];
rz(-0.20656221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0764717) q[1];
sx q[1];
rz(-1.7674801) q[1];
sx q[1];
rz(-1.0573439) q[1];
x q[2];
rz(-1.5557489) q[3];
sx q[3];
rz(-1.5516557) q[3];
sx q[3];
rz(-1.0552989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56090474) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(2.4482751) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.7047434) q[3];
sx q[3];
rz(0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(1.90907) q[0];
rz(-2.0690074) q[1];
sx q[1];
rz(-2.0697846) q[1];
sx q[1];
rz(-0.17428621) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3586853) q[0];
sx q[0];
rz(-1.8795965) q[0];
sx q[0];
rz(-0.82682825) q[0];
x q[1];
rz(-0.66138791) q[2];
sx q[2];
rz(-0.97860133) q[2];
sx q[2];
rz(-2.4370898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6433405) q[1];
sx q[1];
rz(-2.2689515) q[1];
sx q[1];
rz(-0.50486418) q[1];
rz(-pi) q[2];
rz(-1.2234736) q[3];
sx q[3];
rz(-2.8109549) q[3];
sx q[3];
rz(-0.91352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8217414) q[2];
sx q[2];
rz(-0.80703002) q[2];
sx q[2];
rz(-2.9439587) q[2];
rz(-0.28891426) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(-1.4060085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0777247) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(-0.20859627) q[0];
rz(-0.96616191) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(-1.6360412) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9150328) q[0];
sx q[0];
rz(-2.0953) q[0];
sx q[0];
rz(-2.4302308) q[0];
rz(-pi) q[1];
rz(1.534329) q[2];
sx q[2];
rz(-2.0009396) q[2];
sx q[2];
rz(-1.9311116) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.437285) q[1];
sx q[1];
rz(-1.1204801) q[1];
sx q[1];
rz(2.1561161) q[1];
rz(-0.286245) q[3];
sx q[3];
rz(-1.4713333) q[3];
sx q[3];
rz(-3.0486097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(3.0440142) q[2];
rz(1.7476667) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(-2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7512648) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(-2.6275997) q[0];
rz(-0.12318525) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(0.93200144) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829464) q[0];
sx q[0];
rz(-0.46447771) q[0];
sx q[0];
rz(1.8770201) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9281689) q[2];
sx q[2];
rz(-2.4744611) q[2];
sx q[2];
rz(0.92598976) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2145558) q[1];
sx q[1];
rz(-0.80074691) q[1];
sx q[1];
rz(-0.75896778) q[1];
rz(-pi) q[2];
rz(1.7979513) q[3];
sx q[3];
rz(-0.72941581) q[3];
sx q[3];
rz(-1.9279355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0294068) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(-2.712148) q[2];
rz(1.2094234) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(2.4485574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747303) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(-1.8027579) q[0];
rz(-2.4354637) q[1];
sx q[1];
rz(-1.8920205) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50981748) q[0];
sx q[0];
rz(-1.063856) q[0];
sx q[0];
rz(-2.1419873) q[0];
rz(-1.4016897) q[2];
sx q[2];
rz(-2.7204614) q[2];
sx q[2];
rz(-0.087547628) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95987684) q[1];
sx q[1];
rz(-2.1023589) q[1];
sx q[1];
rz(0.2172443) q[1];
x q[2];
rz(1.0334942) q[3];
sx q[3];
rz(-2.6712382) q[3];
sx q[3];
rz(-2.5126517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(-0.48428145) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(-0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1442239) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(2.9123059) q[0];
rz(-2.7067822) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(-0.71892175) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839448) q[0];
sx q[0];
rz(-2.5074208) q[0];
sx q[0];
rz(1.7303403) q[0];
x q[1];
rz(-0.11833338) q[2];
sx q[2];
rz(-0.98630691) q[2];
sx q[2];
rz(-3.0362533) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1246008) q[1];
sx q[1];
rz(-1.7182087) q[1];
sx q[1];
rz(-2.7684545) q[1];
rz(2.6370254) q[3];
sx q[3];
rz(-2.1310398) q[3];
sx q[3];
rz(-1.954078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3832613) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(2.3948005) q[2];
rz(2.2693999) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(-2.1993568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.223021) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(2.7643798) q[1];
sx q[1];
rz(-1.4706392) q[1];
sx q[1];
rz(-2.6249862) q[1];
rz(1.2582851) q[2];
sx q[2];
rz(-2.0682813) q[2];
sx q[2];
rz(-2.993519) q[2];
rz(2.4242998) q[3];
sx q[3];
rz(-1.6834696) q[3];
sx q[3];
rz(-3.0739741) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
