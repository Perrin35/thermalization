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
rz(0.25843698) q[0];
sx q[0];
rz(2.8539477) q[0];
sx q[0];
rz(10.532425) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(0.62825656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1028672) q[0];
sx q[0];
rz(-0.21351457) q[0];
sx q[0];
rz(-0.98189129) q[0];
rz(-pi) q[1];
rz(1.1789178) q[2];
sx q[2];
rz(-0.88681839) q[2];
sx q[2];
rz(2.1563403) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3648277) q[1];
sx q[1];
rz(-1.7994295) q[1];
sx q[1];
rz(1.5336114) q[1];
rz(-pi) q[2];
rz(2.7006458) q[3];
sx q[3];
rz(-2.0618871) q[3];
sx q[3];
rz(-2.2381312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85311741) q[2];
sx q[2];
rz(-1.709781) q[2];
sx q[2];
rz(-2.8004069) q[2];
rz(-2.2513385) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(0.53513479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5219236) q[0];
sx q[0];
rz(-0.6837908) q[0];
sx q[0];
rz(-0.71440119) q[0];
rz(-1.5419386) q[1];
sx q[1];
rz(-0.49595141) q[1];
sx q[1];
rz(3.0468859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3886914) q[0];
sx q[0];
rz(-1.4764278) q[0];
sx q[0];
rz(1.2497942) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8049967) q[2];
sx q[2];
rz(-2.1456652) q[2];
sx q[2];
rz(-0.47913715) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73037877) q[1];
sx q[1];
rz(-1.8068562) q[1];
sx q[1];
rz(2.7884931) q[1];
rz(-pi) q[2];
rz(-2.6323039) q[3];
sx q[3];
rz(-1.4319766) q[3];
sx q[3];
rz(-1.2535439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9793205) q[2];
sx q[2];
rz(-1.2459735) q[2];
sx q[2];
rz(-2.6675513) q[2];
rz(-2.3794409) q[3];
sx q[3];
rz(-0.63665974) q[3];
sx q[3];
rz(-0.4711841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2786461) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(0.85292029) q[0];
rz(-0.22311738) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(-2.6257637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2775622) q[0];
sx q[0];
rz(-1.5169965) q[0];
sx q[0];
rz(3.0911616) q[0];
rz(-pi) q[1];
rz(-2.2081096) q[2];
sx q[2];
rz(-1.3636944) q[2];
sx q[2];
rz(-1.1400956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22166907) q[1];
sx q[1];
rz(-0.95429568) q[1];
sx q[1];
rz(-2.4024582) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87275858) q[3];
sx q[3];
rz(-1.8379194) q[3];
sx q[3];
rz(1.9655379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3333266) q[2];
sx q[2];
rz(-0.9767248) q[2];
sx q[2];
rz(0.094956368) q[2];
rz(2.700108) q[3];
sx q[3];
rz(-1.7850826) q[3];
sx q[3];
rz(-1.1536185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3917291) q[0];
sx q[0];
rz(-2.5666105) q[0];
sx q[0];
rz(2.0854501) q[0];
rz(-0.9437584) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(-1.6897078) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75987923) q[0];
sx q[0];
rz(-0.96764123) q[0];
sx q[0];
rz(0.85079129) q[0];
x q[1];
rz(-1.1805492) q[2];
sx q[2];
rz(-2.0023291) q[2];
sx q[2];
rz(0.081588946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46717146) q[1];
sx q[1];
rz(-1.3302263) q[1];
sx q[1];
rz(1.5340541) q[1];
rz(-pi) q[2];
x q[2];
rz(0.068030997) q[3];
sx q[3];
rz(-0.11091954) q[3];
sx q[3];
rz(-2.7601931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(-0.70017868) q[2];
rz(-2.2705196) q[3];
sx q[3];
rz(-1.0659404) q[3];
sx q[3];
rz(-1.7843436) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(2.6569271) q[0];
rz(-2.2635745) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(2.7427618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1171378) q[0];
sx q[0];
rz(-1.0255735) q[0];
sx q[0];
rz(-1.412789) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9100175) q[2];
sx q[2];
rz(-0.89618635) q[2];
sx q[2];
rz(-2.2774027) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.63020413) q[1];
sx q[1];
rz(-0.92899073) q[1];
sx q[1];
rz(-0.46496572) q[1];
rz(-pi) q[2];
rz(2.2696617) q[3];
sx q[3];
rz(-2.0060349) q[3];
sx q[3];
rz(2.1520681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47348076) q[2];
sx q[2];
rz(-1.6010189) q[2];
sx q[2];
rz(-2.2637746) q[2];
rz(2.7569438) q[3];
sx q[3];
rz(-0.84482241) q[3];
sx q[3];
rz(1.9627242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28410742) q[0];
sx q[0];
rz(-2.6326038) q[0];
sx q[0];
rz(2.8711163) q[0];
rz(-0.30666223) q[1];
sx q[1];
rz(-0.51391927) q[1];
sx q[1];
rz(-2.564863) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79918843) q[0];
sx q[0];
rz(-1.4831838) q[0];
sx q[0];
rz(-1.3848928) q[0];
x q[1];
rz(-2.7790766) q[2];
sx q[2];
rz(-0.91372817) q[2];
sx q[2];
rz(0.18489472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57177793) q[1];
sx q[1];
rz(-1.0228436) q[1];
sx q[1];
rz(-2.673719) q[1];
rz(1.5848198) q[3];
sx q[3];
rz(-1.2360667) q[3];
sx q[3];
rz(1.0098135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3660672) q[2];
sx q[2];
rz(-2.0719353) q[2];
sx q[2];
rz(0.97037399) q[2];
rz(-2.7905285) q[3];
sx q[3];
rz(-0.23215663) q[3];
sx q[3];
rz(2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2358667) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(-2.6406777) q[0];
rz(-0.72055912) q[1];
sx q[1];
rz(-2.2961398) q[1];
sx q[1];
rz(2.2038584) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7739627) q[0];
sx q[0];
rz(-1.367462) q[0];
sx q[0];
rz(-3.0995661) q[0];
rz(-pi) q[1];
rz(-1.9483614) q[2];
sx q[2];
rz(-2.015471) q[2];
sx q[2];
rz(-1.7583762) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5590458) q[1];
sx q[1];
rz(-1.9314879) q[1];
sx q[1];
rz(-0.8949032) q[1];
rz(-pi) q[2];
rz(1.5608801) q[3];
sx q[3];
rz(-0.88396931) q[3];
sx q[3];
rz(-1.6520985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0096036) q[2];
sx q[2];
rz(-2.9961573) q[2];
sx q[2];
rz(2.2930938) q[2];
rz(2.2961473) q[3];
sx q[3];
rz(-0.65631056) q[3];
sx q[3];
rz(1.2257303) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41035143) q[0];
sx q[0];
rz(-0.9587962) q[0];
sx q[0];
rz(-0.89286667) q[0];
rz(-0.061575312) q[1];
sx q[1];
rz(-2.4696746) q[1];
sx q[1];
rz(0.16199131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0688871) q[0];
sx q[0];
rz(-1.6729682) q[0];
sx q[0];
rz(-0.01553681) q[0];
rz(-pi) q[1];
rz(-3.0062469) q[2];
sx q[2];
rz(-0.70643808) q[2];
sx q[2];
rz(1.634623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8422441) q[1];
sx q[1];
rz(-0.33978251) q[1];
sx q[1];
rz(2.7187125) q[1];
x q[2];
rz(0.83655436) q[3];
sx q[3];
rz(-3.0020107) q[3];
sx q[3];
rz(2.9849646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.050345) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(0.72489911) q[2];
rz(2.8643518) q[3];
sx q[3];
rz(-2.1767949) q[3];
sx q[3];
rz(-0.086517081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8898833) q[0];
sx q[0];
rz(-2.4885663) q[0];
sx q[0];
rz(3.1157893) q[0];
rz(1.39894) q[1];
sx q[1];
rz(-1.6394697) q[1];
sx q[1];
rz(0.10765156) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4336595) q[0];
sx q[0];
rz(-1.7042394) q[0];
sx q[0];
rz(0.99876499) q[0];
x q[1];
rz(-2.3978762) q[2];
sx q[2];
rz(-0.60600835) q[2];
sx q[2];
rz(-1.5569937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3947666) q[1];
sx q[1];
rz(-2.0438384) q[1];
sx q[1];
rz(-1.6113144) q[1];
rz(-1.6844382) q[3];
sx q[3];
rz(-1.5885308) q[3];
sx q[3];
rz(2.4271698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.671635) q[2];
sx q[2];
rz(-1.2592955) q[2];
sx q[2];
rz(-2.0476445) q[2];
rz(0.57540244) q[3];
sx q[3];
rz(-0.83991528) q[3];
sx q[3];
rz(2.0247139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6082918) q[0];
sx q[0];
rz(-0.58654439) q[0];
sx q[0];
rz(-2.4285512) q[0];
rz(-2.9167922) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(1.0932659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21205344) q[0];
sx q[0];
rz(-2.0418379) q[0];
sx q[0];
rz(-0.77917288) q[0];
x q[1];
rz(-0.20168882) q[2];
sx q[2];
rz(-1.8400157) q[2];
sx q[2];
rz(-1.9112916) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.097979989) q[1];
sx q[1];
rz(-0.47046767) q[1];
sx q[1];
rz(-0.2474252) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3087464) q[3];
sx q[3];
rz(-0.42632494) q[3];
sx q[3];
rz(-0.48708068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0010058086) q[2];
sx q[2];
rz(-0.69585496) q[2];
sx q[2];
rz(-0.33983964) q[2];
rz(3.0991683) q[3];
sx q[3];
rz(-1.2845311) q[3];
sx q[3];
rz(2.5137918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9545659) q[0];
sx q[0];
rz(-1.621959) q[0];
sx q[0];
rz(1.8931615) q[0];
rz(2.3595702) q[1];
sx q[1];
rz(-0.93688688) q[1];
sx q[1];
rz(-0.72601906) q[1];
rz(2.7693979) q[2];
sx q[2];
rz(-1.2584465) q[2];
sx q[2];
rz(-1.0312205) q[2];
rz(-1.9966765) q[3];
sx q[3];
rz(-1.799634) q[3];
sx q[3];
rz(-1.4992731) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
