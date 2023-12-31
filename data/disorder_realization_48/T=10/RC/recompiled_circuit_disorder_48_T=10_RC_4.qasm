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
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5028635) q[0];
sx q[0];
rz(-1.5348866) q[0];
sx q[0];
rz(-2.3907651) q[0];
rz(-pi) q[1];
rz(1.5300418) q[2];
sx q[2];
rz(-2.0353122) q[2];
sx q[2];
rz(-2.0704616) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4301181) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-2.3858566) q[1];
rz(0.18913194) q[3];
sx q[3];
rz(-1.3526219) q[3];
sx q[3];
rz(1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-2.5640326) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(0.66453385) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(0.93915835) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(1.739025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70799202) q[0];
sx q[0];
rz(-0.029768243) q[0];
sx q[0];
rz(-0.13694163) q[0];
rz(-pi) q[1];
rz(-2.5820929) q[2];
sx q[2];
rz(-1.7569949) q[2];
sx q[2];
rz(-2.6452243) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1075588) q[1];
sx q[1];
rz(-2.4817913) q[1];
sx q[1];
rz(0.69114139) q[1];
x q[2];
rz(-1.6658695) q[3];
sx q[3];
rz(-0.94082309) q[3];
sx q[3];
rz(0.41439393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(2.6638022) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(-2.7405222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0862335) q[0];
sx q[0];
rz(-1.2463352) q[0];
sx q[0];
rz(-3.0326891) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7602073) q[2];
sx q[2];
rz(-2.2064798) q[2];
sx q[2];
rz(0.96178255) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2420826) q[1];
sx q[1];
rz(-0.84065719) q[1];
sx q[1];
rz(1.5622557) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1242261) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(-0.22180804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.697631) q[0];
rz(1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(0.25340432) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70595104) q[0];
sx q[0];
rz(-0.55103978) q[0];
sx q[0];
rz(-1.3735229) q[0];
rz(-pi) q[1];
rz(3.1242712) q[2];
sx q[2];
rz(-1.0850731) q[2];
sx q[2];
rz(0.7893562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4761915) q[1];
sx q[1];
rz(-0.63648495) q[1];
sx q[1];
rz(2.979216) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0074189) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(-1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-1.0095989) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71516365) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(2.6834992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040745) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(-3.1361561) q[0];
x q[1];
rz(0.96197084) q[2];
sx q[2];
rz(-1.1912727) q[2];
sx q[2];
rz(0.48987197) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.304368) q[1];
sx q[1];
rz(-0.75290426) q[1];
sx q[1];
rz(1.7423082) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3985653) q[3];
sx q[3];
rz(-2.6346364) q[3];
sx q[3];
rz(-2.7996922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1725585) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(-0.0028006639) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(-1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7125268) q[0];
sx q[0];
rz(-2.292284) q[0];
sx q[0];
rz(-0.52961403) q[0];
rz(2.3976372) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(-0.26847408) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.25176469) q[1];
sx q[1];
rz(-1.3740731) q[1];
sx q[1];
rz(1.9038492) q[1];
x q[2];
rz(1.091286) q[3];
sx q[3];
rz(-1.3538085) q[3];
sx q[3];
rz(3.019141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1798114) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(2.0297208) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8060914) q[0];
sx q[0];
rz(-1.6643235) q[0];
sx q[0];
rz(-1.085824) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56356168) q[2];
sx q[2];
rz(-0.60699082) q[2];
sx q[2];
rz(-0.8286455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3355616) q[1];
sx q[1];
rz(-2.8415488) q[1];
sx q[1];
rz(-2.8801444) q[1];
x q[2];
rz(1.0681549) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(-0.84514602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325571) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(0.73079601) q[0];
rz(0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(2.3866167) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2728111) q[0];
sx q[0];
rz(-0.73505721) q[0];
sx q[0];
rz(-1.3520157) q[0];
rz(-pi) q[1];
rz(-0.72889363) q[2];
sx q[2];
rz(-2.0260099) q[2];
sx q[2];
rz(2.560937) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4432115) q[1];
sx q[1];
rz(-1.4581919) q[1];
sx q[1];
rz(0.57971445) q[1];
x q[2];
rz(-0.8927535) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(2.1219818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-0.019502217) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0855904) q[0];
sx q[0];
rz(-1.3645571) q[0];
sx q[0];
rz(0.081736728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15140622) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(-0.99934794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6448977) q[1];
sx q[1];
rz(-0.95458191) q[1];
sx q[1];
rz(-2.9195021) q[1];
x q[2];
rz(-1.2583624) q[3];
sx q[3];
rz(-1.4133246) q[3];
sx q[3];
rz(2.0563994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(-2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(-1.5415812) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(0.48450255) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(-1.9932995) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0120221) q[0];
sx q[0];
rz(-2.9992933) q[0];
sx q[0];
rz(-2.9905031) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8414367) q[2];
sx q[2];
rz(-2.488392) q[2];
sx q[2];
rz(2.4094827) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82603589) q[1];
sx q[1];
rz(-0.6912187) q[1];
sx q[1];
rz(1.300699) q[1];
rz(-pi) q[2];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(2.7764376) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0317595) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(-2.1784492) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(-1.1744432) q[2];
sx q[2];
rz(-0.94559961) q[2];
sx q[2];
rz(-0.0018975817) q[2];
rz(-0.036988463) q[3];
sx q[3];
rz(-2.130641) q[3];
sx q[3];
rz(-0.11161042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
