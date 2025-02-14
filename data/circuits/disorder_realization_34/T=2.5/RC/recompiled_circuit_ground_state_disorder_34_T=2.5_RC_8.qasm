OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3736149) q[0];
sx q[0];
rz(-0.64426214) q[0];
sx q[0];
rz(2.2777519) q[0];
rz(-1.5683501) q[1];
sx q[1];
rz(-1.22998) q[1];
sx q[1];
rz(-0.65043515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.602055) q[0];
sx q[0];
rz(-1.3663174) q[0];
sx q[0];
rz(0.50907593) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8312248) q[2];
sx q[2];
rz(-1.5987607) q[2];
sx q[2];
rz(1.8486277) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1623154) q[1];
sx q[1];
rz(-0.82274306) q[1];
sx q[1];
rz(-0.85461275) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6780544) q[3];
sx q[3];
rz(-1.4441731) q[3];
sx q[3];
rz(2.7789314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2245366) q[2];
sx q[2];
rz(-1.7660331) q[2];
sx q[2];
rz(1.7898233) q[2];
rz(0.81870643) q[3];
sx q[3];
rz(-1.4592146) q[3];
sx q[3];
rz(1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390035) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(-0.57574058) q[0];
rz(2.2585244) q[1];
sx q[1];
rz(-1.6022316) q[1];
sx q[1];
rz(1.8227122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3991845) q[0];
sx q[0];
rz(-0.78982805) q[0];
sx q[0];
rz(2.7862048) q[0];
x q[1];
rz(-0.024818161) q[2];
sx q[2];
rz(-0.78486004) q[2];
sx q[2];
rz(2.047124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12815753) q[1];
sx q[1];
rz(-2.2655106) q[1];
sx q[1];
rz(0.46626587) q[1];
rz(2.6448375) q[3];
sx q[3];
rz(-2.7803034) q[3];
sx q[3];
rz(0.75888097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5004702) q[2];
sx q[2];
rz(-2.1239855) q[2];
sx q[2];
rz(-2.9024331) q[2];
rz(0.33254361) q[3];
sx q[3];
rz(-1.6859237) q[3];
sx q[3];
rz(2.0502162) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7964433) q[0];
sx q[0];
rz(-2.394016) q[0];
sx q[0];
rz(0.24459608) q[0];
rz(-0.39464513) q[1];
sx q[1];
rz(-0.60401812) q[1];
sx q[1];
rz(1.9371202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8926579) q[0];
sx q[0];
rz(-1.5818074) q[0];
sx q[0];
rz(-1.5599773) q[0];
x q[1];
rz(-1.081624) q[2];
sx q[2];
rz(-1.2252639) q[2];
sx q[2];
rz(-0.45128497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2358346) q[1];
sx q[1];
rz(-1.2819965) q[1];
sx q[1];
rz(-2.1101084) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2368771) q[3];
sx q[3];
rz(-2.534233) q[3];
sx q[3];
rz(-1.5752618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34387732) q[2];
sx q[2];
rz(-2.158973) q[2];
sx q[2];
rz(-0.24125153) q[2];
rz(1.7734211) q[3];
sx q[3];
rz(-1.3894812) q[3];
sx q[3];
rz(-2.1696137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57206804) q[0];
sx q[0];
rz(-1.854874) q[0];
sx q[0];
rz(1.1014112) q[0];
rz(-1.4934348) q[1];
sx q[1];
rz(-2.8547574) q[1];
sx q[1];
rz(-2.1410087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1959609) q[0];
sx q[0];
rz(-1.1544265) q[0];
sx q[0];
rz(-3.065057) q[0];
x q[1];
rz(1.5865636) q[2];
sx q[2];
rz(-2.5852499) q[2];
sx q[2];
rz(-0.76130262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2615255) q[1];
sx q[1];
rz(-1.6717922) q[1];
sx q[1];
rz(1.5965881) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8832753) q[3];
sx q[3];
rz(-0.44033209) q[3];
sx q[3];
rz(-0.60763121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3558041) q[2];
sx q[2];
rz(-1.641909) q[2];
sx q[2];
rz(-2.918952) q[2];
rz(1.6629793) q[3];
sx q[3];
rz(-0.40545774) q[3];
sx q[3];
rz(1.2610029) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8292238) q[0];
sx q[0];
rz(-2.0033422) q[0];
sx q[0];
rz(2.9421222) q[0];
rz(1.624674) q[1];
sx q[1];
rz(-1.805504) q[1];
sx q[1];
rz(-0.46474251) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366531) q[0];
sx q[0];
rz(-0.31855983) q[0];
sx q[0];
rz(-2.0037988) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6301441) q[2];
sx q[2];
rz(-1.7918824) q[2];
sx q[2];
rz(0.57002744) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.81632751) q[1];
sx q[1];
rz(-1.8598598) q[1];
sx q[1];
rz(1.2802034) q[1];
rz(-pi) q[2];
rz(0.0031849234) q[3];
sx q[3];
rz(-0.69744686) q[3];
sx q[3];
rz(1.7980597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0988079) q[2];
sx q[2];
rz(-1.837919) q[2];
sx q[2];
rz(-0.25351563) q[2];
rz(2.2809196) q[3];
sx q[3];
rz(-0.52152514) q[3];
sx q[3];
rz(-1.1389987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.072170243) q[0];
sx q[0];
rz(-1.5266029) q[0];
sx q[0];
rz(1.3602863) q[0];
rz(0.6036497) q[1];
sx q[1];
rz(-0.79704469) q[1];
sx q[1];
rz(2.6545677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384772) q[0];
sx q[0];
rz(-0.91775187) q[0];
sx q[0];
rz(-1.4794635) q[0];
x q[1];
rz(-1.4757882) q[2];
sx q[2];
rz(-0.83432799) q[2];
sx q[2];
rz(2.468994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7051866) q[1];
sx q[1];
rz(-1.7366221) q[1];
sx q[1];
rz(1.9351744) q[1];
x q[2];
rz(-0.2305687) q[3];
sx q[3];
rz(-1.9631223) q[3];
sx q[3];
rz(2.2986064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46502486) q[2];
sx q[2];
rz(-2.1239231) q[2];
sx q[2];
rz(-2.6505995) q[2];
rz(0.53389126) q[3];
sx q[3];
rz(-0.54755727) q[3];
sx q[3];
rz(3.0965613) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8462867) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(1.8701766) q[0];
rz(2.2170587) q[1];
sx q[1];
rz(-1.9969767) q[1];
sx q[1];
rz(-1.0251934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7780925) q[0];
sx q[0];
rz(-0.94787593) q[0];
sx q[0];
rz(1.8353375) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94550262) q[2];
sx q[2];
rz(-1.7472113) q[2];
sx q[2];
rz(-2.8636065) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89910903) q[1];
sx q[1];
rz(-0.61990737) q[1];
sx q[1];
rz(-1.0044881) q[1];
rz(-2.4348172) q[3];
sx q[3];
rz(-1.3825582) q[3];
sx q[3];
rz(2.1911603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2123432) q[2];
sx q[2];
rz(-0.29444567) q[2];
sx q[2];
rz(0.92246169) q[2];
rz(2.755002) q[3];
sx q[3];
rz(-1.2873193) q[3];
sx q[3];
rz(-1.6283901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59231049) q[0];
sx q[0];
rz(-0.03454241) q[0];
sx q[0];
rz(0.70796815) q[0];
rz(2.6225846) q[1];
sx q[1];
rz(-0.52878562) q[1];
sx q[1];
rz(-0.66600287) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7836664) q[0];
sx q[0];
rz(-0.54364288) q[0];
sx q[0];
rz(-1.3915791) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8911036) q[2];
sx q[2];
rz(-1.4605117) q[2];
sx q[2];
rz(-3.0406893) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.085364522) q[1];
sx q[1];
rz(-0.66909583) q[1];
sx q[1];
rz(3.1345033) q[1];
rz(1.610393) q[3];
sx q[3];
rz(-2.5259113) q[3];
sx q[3];
rz(-0.99056584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4697326) q[2];
sx q[2];
rz(-1.9817151) q[2];
sx q[2];
rz(-1.0603909) q[2];
rz(-1.7139858) q[3];
sx q[3];
rz(-1.5644045) q[3];
sx q[3];
rz(-0.76712999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0107182) q[0];
sx q[0];
rz(-2.3227782) q[0];
sx q[0];
rz(2.4406216) q[0];
rz(2.2531807) q[1];
sx q[1];
rz(-2.7263434) q[1];
sx q[1];
rz(1.7078687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2727391) q[0];
sx q[0];
rz(-2.0827385) q[0];
sx q[0];
rz(-1.4466132) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5436243) q[2];
sx q[2];
rz(-1.9765696) q[2];
sx q[2];
rz(-2.9494065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6810602) q[1];
sx q[1];
rz(-1.5892481) q[1];
sx q[1];
rz(-3.1204719) q[1];
rz(-pi) q[2];
rz(0.69207003) q[3];
sx q[3];
rz(-1.0502397) q[3];
sx q[3];
rz(-0.51560054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.99011699) q[2];
sx q[2];
rz(-1.9141804) q[2];
sx q[2];
rz(1.6647476) q[2];
rz(-2.9448523) q[3];
sx q[3];
rz(-1.1367831) q[3];
sx q[3];
rz(2.2018532) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6075341) q[0];
sx q[0];
rz(-1.5787831) q[0];
sx q[0];
rz(-2.0205355) q[0];
rz(2.6634482) q[1];
sx q[1];
rz(-0.68111626) q[1];
sx q[1];
rz(0.71472439) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8281855) q[0];
sx q[0];
rz(-1.3696693) q[0];
sx q[0];
rz(0.49230663) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0214731) q[2];
sx q[2];
rz(-1.2397814) q[2];
sx q[2];
rz(2.2084055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1063544) q[1];
sx q[1];
rz(-2.0788621) q[1];
sx q[1];
rz(0.50915995) q[1];
rz(-pi) q[2];
x q[2];
rz(2.274005) q[3];
sx q[3];
rz(-1.0699268) q[3];
sx q[3];
rz(-1.3102368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62632442) q[2];
sx q[2];
rz(-1.3268665) q[2];
sx q[2];
rz(-1.932762) q[2];
rz(-1.6995466) q[3];
sx q[3];
rz(-0.88651005) q[3];
sx q[3];
rz(-3.076021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6616853) q[0];
sx q[0];
rz(-1.8330782) q[0];
sx q[0];
rz(-0.080396419) q[0];
rz(-2.5936364) q[1];
sx q[1];
rz(-2.4520271) q[1];
sx q[1];
rz(1.3507623) q[1];
rz(-1.1153658) q[2];
sx q[2];
rz(-2.6276988) q[2];
sx q[2];
rz(-2.217271) q[2];
rz(0.73765909) q[3];
sx q[3];
rz(-0.59288965) q[3];
sx q[3];
rz(2.2137143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
