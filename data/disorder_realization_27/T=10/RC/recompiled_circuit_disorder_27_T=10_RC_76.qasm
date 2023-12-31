OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(-1.385153) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0692682) q[0];
sx q[0];
rz(-0.82075483) q[0];
sx q[0];
rz(0.44278352) q[0];
rz(-pi) q[1];
rz(2.7117549) q[2];
sx q[2];
rz(-0.59525437) q[2];
sx q[2];
rz(1.0560448) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79917158) q[1];
sx q[1];
rz(-0.28563269) q[1];
sx q[1];
rz(-0.61113961) q[1];
rz(-2.4888943) q[3];
sx q[3];
rz(-1.1678809) q[3];
sx q[3];
rz(-2.9626915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2549071) q[2];
sx q[2];
rz(-0.79780769) q[2];
sx q[2];
rz(-0.20516667) q[2];
rz(-0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(1.1024968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7339864) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(-0.45390391) q[0];
rz(-1.0247963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(-1.227238) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1821546) q[0];
sx q[0];
rz(-1.4417138) q[0];
sx q[0];
rz(3.0653619) q[0];
rz(2.6229157) q[2];
sx q[2];
rz(-1.9762602) q[2];
sx q[2];
rz(2.5395218) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2482359) q[1];
sx q[1];
rz(-1.8289939) q[1];
sx q[1];
rz(2.8398819) q[1];
x q[2];
rz(0.20083986) q[3];
sx q[3];
rz(-1.7240932) q[3];
sx q[3];
rz(-0.64985819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(-0.56742898) q[2];
rz(-2.7764017) q[3];
sx q[3];
rz(-1.4130211) q[3];
sx q[3];
rz(-0.96810961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.658618) q[0];
sx q[0];
rz(-2.5768319) q[0];
sx q[0];
rz(-0.89865249) q[0];
rz(2.1458416) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(0.333289) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62045287) q[0];
sx q[0];
rz(-2.1846909) q[0];
sx q[0];
rz(0.99434538) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.097924175) q[2];
sx q[2];
rz(-1.3946748) q[2];
sx q[2];
rz(0.91913659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.052735141) q[1];
sx q[1];
rz(-0.6328859) q[1];
sx q[1];
rz(-1.3120552) q[1];
rz(-2.1786147) q[3];
sx q[3];
rz(-0.64219785) q[3];
sx q[3];
rz(-1.1841139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68625346) q[2];
sx q[2];
rz(-1.3505961) q[2];
sx q[2];
rz(-1.1509482) q[2];
rz(2.3006556) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(1.1545198) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6999917) q[0];
sx q[0];
rz(-1.561152) q[0];
sx q[0];
rz(2.4568795) q[0];
rz(2.1060064) q[1];
sx q[1];
rz(-2.6338449) q[1];
sx q[1];
rz(-1.9365786) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689777) q[0];
sx q[0];
rz(-0.87991558) q[0];
sx q[0];
rz(-0.0432424) q[0];
rz(-0.92653805) q[2];
sx q[2];
rz(-1.7346003) q[2];
sx q[2];
rz(1.0545727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64441427) q[1];
sx q[1];
rz(-2.2481611) q[1];
sx q[1];
rz(0.35269423) q[1];
x q[2];
rz(3.0699176) q[3];
sx q[3];
rz(-0.86719162) q[3];
sx q[3];
rz(3.049831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24923199) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(2.7704346) q[2];
rz(1.4012198) q[3];
sx q[3];
rz(-0.6597844) q[3];
sx q[3];
rz(-1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-0.78614569) q[0];
sx q[0];
rz(-0.13312419) q[0];
rz(-2.1482824) q[1];
sx q[1];
rz(-1.7555833) q[1];
sx q[1];
rz(0.55508074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.426794) q[0];
sx q[0];
rz(-1.5580651) q[0];
sx q[0];
rz(1.8861594) q[0];
rz(-pi) q[1];
rz(1.7246036) q[2];
sx q[2];
rz(-0.4193192) q[2];
sx q[2];
rz(2.9749982) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3913369) q[1];
sx q[1];
rz(-1.5515944) q[1];
sx q[1];
rz(-2.0160497) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57226945) q[3];
sx q[3];
rz(-3.045445) q[3];
sx q[3];
rz(0.26086807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.83539) q[2];
sx q[2];
rz(-2.1357048) q[2];
sx q[2];
rz(0.13892697) q[2];
rz(-0.94240087) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(-2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5979364) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(0.57975769) q[0];
rz(-3.014091) q[1];
sx q[1];
rz(-1.9516877) q[1];
sx q[1];
rz(-1.5396083) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0149536) q[0];
sx q[0];
rz(-2.3987028) q[0];
sx q[0];
rz(1.2680306) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0110353) q[2];
sx q[2];
rz(-1.6160384) q[2];
sx q[2];
rz(0.10938489) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74486596) q[1];
sx q[1];
rz(-0.52280451) q[1];
sx q[1];
rz(-2.1641939) q[1];
x q[2];
rz(2.8508938) q[3];
sx q[3];
rz(-0.906956) q[3];
sx q[3];
rz(1.78474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5876028) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(-2.8721151) q[2];
rz(-0.23412165) q[3];
sx q[3];
rz(-2.6168489) q[3];
sx q[3];
rz(3.0814734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44678974) q[0];
sx q[0];
rz(-0.61820522) q[0];
sx q[0];
rz(0.68429464) q[0];
rz(3.0220095) q[1];
sx q[1];
rz(-1.8493098) q[1];
sx q[1];
rz(0.51876846) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8914688) q[0];
sx q[0];
rz(-0.84202535) q[0];
sx q[0];
rz(-2.9617589) q[0];
rz(0.22612818) q[2];
sx q[2];
rz(-0.78352189) q[2];
sx q[2];
rz(2.4353611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59473945) q[1];
sx q[1];
rz(-2.2097128) q[1];
sx q[1];
rz(-0.98313318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.026168907) q[3];
sx q[3];
rz(-1.9713638) q[3];
sx q[3];
rz(-2.8699584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2542904) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(-0.56345144) q[2];
rz(-3.0900132) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96034399) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(-1.3990336) q[0];
rz(-2.3545806) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(0.74434892) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3314914) q[0];
sx q[0];
rz(-1.5163251) q[0];
sx q[0];
rz(-2.9936552) q[0];
x q[1];
rz(-1.3049576) q[2];
sx q[2];
rz(-0.53005866) q[2];
sx q[2];
rz(2.838138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3202312) q[1];
sx q[1];
rz(-1.3470955) q[1];
sx q[1];
rz(0.51364586) q[1];
rz(-0.45458557) q[3];
sx q[3];
rz(-1.5919519) q[3];
sx q[3];
rz(-1.8225614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7156334) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(-0.60950935) q[2];
rz(2.4842747) q[3];
sx q[3];
rz(-0.64353639) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-3.0472026) q[0];
sx q[0];
rz(-1.6375861) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(-2.3666568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79561728) q[0];
sx q[0];
rz(-0.91751639) q[0];
sx q[0];
rz(0.4972636) q[0];
rz(-pi) q[1];
rz(-2.0140892) q[2];
sx q[2];
rz(-1.7121676) q[2];
sx q[2];
rz(3.1183426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7271125) q[1];
sx q[1];
rz(-0.4948805) q[1];
sx q[1];
rz(0.40463573) q[1];
rz(-pi) q[2];
rz(1.3516515) q[3];
sx q[3];
rz(-1.7100167) q[3];
sx q[3];
rz(2.5185891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.187414) q[2];
sx q[2];
rz(-2.9118907) q[2];
sx q[2];
rz(0.15787086) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-2.082943) q[3];
sx q[3];
rz(-1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(2.9272595) q[0];
rz(-0.65746039) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(2.0956031) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56453088) q[0];
sx q[0];
rz(-2.7663167) q[0];
sx q[0];
rz(2.8481086) q[0];
rz(-0.59098737) q[2];
sx q[2];
rz(-1.5948442) q[2];
sx q[2];
rz(1.9901333) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54088456) q[1];
sx q[1];
rz(-2.4661459) q[1];
sx q[1];
rz(2.1727968) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5973741) q[3];
sx q[3];
rz(-1.9760895) q[3];
sx q[3];
rz(-1.6534896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7252698) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(1.0894758) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.243306) q[3];
sx q[3];
rz(-0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(-1.5325585) q[1];
sx q[1];
rz(-1.4668203) q[1];
sx q[1];
rz(-1.1046881) q[1];
rz(-0.40907787) q[2];
sx q[2];
rz(-1.8553875) q[2];
sx q[2];
rz(2.6864048) q[2];
rz(-3.0388721) q[3];
sx q[3];
rz(-0.552388) q[3];
sx q[3];
rz(-1.296476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
