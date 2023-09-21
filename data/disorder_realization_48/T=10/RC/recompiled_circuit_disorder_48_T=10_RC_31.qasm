OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(-2.350783) q[0];
sx q[0];
rz(2.8074582) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96555644) q[0];
sx q[0];
rz(-2.3210225) q[0];
sx q[0];
rz(1.6198938) q[0];
rz(-pi) q[1];
rz(-0.46484868) q[2];
sx q[2];
rz(-1.5343622) q[2];
sx q[2];
rz(-0.51793098) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76518607) q[1];
sx q[1];
rz(-1.1131439) q[1];
sx q[1];
rz(-0.54995579) q[1];
x q[2];
rz(-2.9524607) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(-1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1537271) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(0.38744774) q[0];
rz(-0.93915835) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(1.4025677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70799202) q[0];
sx q[0];
rz(-0.029768243) q[0];
sx q[0];
rz(-3.004651) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8005373) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(-1.7797433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1837511) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(0.5387696) q[1];
rz(-pi) q[2];
rz(-2.5094633) q[3];
sx q[3];
rz(-1.6475793) q[3];
sx q[3];
rz(1.1002822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(2.823901) q[2];
rz(-2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3690255) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-0.40107045) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553592) q[0];
sx q[0];
rz(-1.2463352) q[0];
sx q[0];
rz(-3.0326891) q[0];
x q[1];
rz(-2.4972649) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(-2.6459141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.229278) q[1];
sx q[1];
rz(-0.73017987) q[1];
sx q[1];
rz(0.0095403949) q[1];
rz(-pi) q[2];
rz(1.1776351) q[3];
sx q[3];
rz(-1.3893681) q[3];
sx q[3];
rz(-1.3822615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5473189) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926086) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.697631) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(-0.25340432) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93638203) q[0];
sx q[0];
rz(-2.10996) q[0];
sx q[0];
rz(-0.11986952) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6035945) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6771486) q[1];
sx q[1];
rz(-0.943999) q[1];
sx q[1];
rz(-1.6897175) q[1];
rz(1.6729309) q[3];
sx q[3];
rz(-0.92389744) q[3];
sx q[3];
rz(1.981786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0306586) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(2.3542662) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(0.45809349) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375181) q[0];
sx q[0];
rz(-2.6989557) q[0];
sx q[0];
rz(-3.1361561) q[0];
rz(-pi) q[1];
rz(2.6890254) q[2];
sx q[2];
rz(-1.0107702) q[2];
sx q[2];
rz(0.82816154) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.749436) q[1];
sx q[1];
rz(-1.6877618) q[1];
sx q[1];
rz(-0.82526308) q[1];
rz(0.094893806) q[3];
sx q[3];
rz(-1.0720383) q[3];
sx q[3];
rz(2.9961078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(2.9120973) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(-2.2221185) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.2671635) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8454682) q[0];
sx q[0];
rz(-0.86588973) q[0];
sx q[0];
rz(-2.092093) q[0];
rz(1.7153347) q[2];
sx q[2];
rz(-0.83206165) q[2];
sx q[2];
rz(1.9369672) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25176469) q[1];
sx q[1];
rz(-1.3740731) q[1];
sx q[1];
rz(-1.2377435) q[1];
rz(-pi) q[2];
rz(2.8980428) q[3];
sx q[3];
rz(-2.0381513) q[3];
sx q[3];
rz(-1.5817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1798114) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(-3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(2.069058) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(2.0297208) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5892964) q[0];
sx q[0];
rz(-2.6483905) q[0];
sx q[0];
rz(-1.769355) q[0];
x q[1];
rz(2.6107437) q[2];
sx q[2];
rz(-1.8804272) q[2];
sx q[2];
rz(-1.2209148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8060311) q[1];
sx q[1];
rz(-0.30004382) q[1];
sx q[1];
rz(0.26144822) q[1];
x q[2];
rz(2.0734378) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(0.84514602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325571) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6801493) q[0];
sx q[0];
rz(-1.4247243) q[0];
sx q[0];
rz(-2.2938674) q[0];
rz(2.5078012) q[2];
sx q[2];
rz(-2.3049424) q[2];
sx q[2];
rz(-1.6939236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8440486) q[1];
sx q[1];
rz(-0.58931671) q[1];
sx q[1];
rz(-2.938016) q[1];
rz(2.6492277) q[3];
sx q[3];
rz(-2.1015321) q[3];
sx q[3];
rz(-0.20519557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(2.0027347) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(-0.019502217) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6396128) q[0];
sx q[0];
rz(-1.6507971) q[0];
sx q[0];
rz(1.3638858) q[0];
rz(-2.9901864) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(-2.1422447) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2724185) q[1];
sx q[1];
rz(-0.65009102) q[1];
sx q[1];
rz(-1.8723349) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0476258) q[3];
sx q[3];
rz(-2.7928824) q[3];
sx q[3];
rz(-2.204012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-0.48450255) q[0];
rz(1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.9932995) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5507817) q[0];
sx q[0];
rz(-1.592144) q[0];
sx q[0];
rz(3.000893) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8414367) q[2];
sx q[2];
rz(-2.488392) q[2];
sx q[2];
rz(0.73210994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9705829) q[1];
sx q[1];
rz(-2.2323771) q[1];
sx q[1];
rz(0.2172825) q[1];
x q[2];
rz(2.9989472) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(2.8838317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(3.0129516) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(0.49114901) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
rz(-2.1309489) q[3];
sx q[3];
rz(-1.5394566) q[3];
sx q[3];
rz(1.4395366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
