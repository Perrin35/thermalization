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
rz(-1.7582769) q[0];
sx q[0];
rz(-1.7051237) q[0];
sx q[0];
rz(-0.96487784) q[0];
rz(2.3896253) q[1];
sx q[1];
rz(-2.7120092) q[1];
sx q[1];
rz(1.0493976) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2478188) q[0];
sx q[0];
rz(-2.4501094) q[0];
sx q[0];
rz(2.1201503) q[0];
rz(2.9994316) q[2];
sx q[2];
rz(-1.249525) q[2];
sx q[2];
rz(2.4427593) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5738028) q[1];
sx q[1];
rz(-1.9490697) q[1];
sx q[1];
rz(-0.90444629) q[1];
rz(-pi) q[2];
rz(1.4780294) q[3];
sx q[3];
rz(-0.99348196) q[3];
sx q[3];
rz(2.4924459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.45399484) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(3.122984) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(-2.4936567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(1.7740771) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(-0.53502214) q[0];
rz(2.6016443) q[1];
sx q[1];
rz(-0.63553634) q[1];
sx q[1];
rz(-1.7417057) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2311418) q[0];
sx q[0];
rz(-2.5325091) q[0];
sx q[0];
rz(0.85394359) q[0];
x q[1];
rz(-0.56324701) q[2];
sx q[2];
rz(-2.0958825) q[2];
sx q[2];
rz(-1.4525177) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3185841) q[1];
sx q[1];
rz(-0.46097791) q[1];
sx q[1];
rz(0.58580841) q[1];
x q[2];
rz(-1.0008282) q[3];
sx q[3];
rz(-1.0335575) q[3];
sx q[3];
rz(-1.9443823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0672368) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(0.92607099) q[2];
rz(-0.98008424) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(-0.66942352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7922908) q[0];
sx q[0];
rz(-1.8632357) q[0];
sx q[0];
rz(-1.5128304) q[0];
rz(2.8577562) q[1];
sx q[1];
rz(-2.2161039) q[1];
sx q[1];
rz(2.0294752) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4130198) q[0];
sx q[0];
rz(-1.5778825) q[0];
sx q[0];
rz(-1.5779693) q[0];
rz(-0.9588963) q[2];
sx q[2];
rz(-3.009353) q[2];
sx q[2];
rz(0.43741495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49431153) q[1];
sx q[1];
rz(-1.0930645) q[1];
sx q[1];
rz(1.3603489) q[1];
rz(-pi) q[2];
rz(2.1792382) q[3];
sx q[3];
rz(-2.7076004) q[3];
sx q[3];
rz(0.8784465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5230368) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(2.0396566) q[2];
rz(2.2526422) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(-0.86047188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717644) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(-2.4523822) q[0];
rz(2.8269732) q[1];
sx q[1];
rz(-2.2210821) q[1];
sx q[1];
rz(1.4124195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1293761) q[0];
sx q[0];
rz(-2.3123992) q[0];
sx q[0];
rz(-2.0022805) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41823776) q[2];
sx q[2];
rz(-2.9091638) q[2];
sx q[2];
rz(-1.4168036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38217332) q[1];
sx q[1];
rz(-1.5133817) q[1];
sx q[1];
rz(2.7613954) q[1];
rz(-1.0869157) q[3];
sx q[3];
rz(-1.3984507) q[3];
sx q[3];
rz(0.92007557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41669258) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(-0.55595428) q[2];
rz(1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(2.8506193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3274662) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(-0.59447527) q[0];
rz(0.56198436) q[1];
sx q[1];
rz(-2.2832506) q[1];
sx q[1];
rz(-2.1655703) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9851889) q[0];
sx q[0];
rz(-1.1431343) q[0];
sx q[0];
rz(0.39910631) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1107509) q[2];
sx q[2];
rz(-0.75428666) q[2];
sx q[2];
rz(-0.13067836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.70006338) q[1];
sx q[1];
rz(-2.2231327) q[1];
sx q[1];
rz(1.1928012) q[1];
rz(-2.8401883) q[3];
sx q[3];
rz(-1.3632953) q[3];
sx q[3];
rz(-1.0338155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48656616) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(2.5308934) q[2];
rz(-0.30580172) q[3];
sx q[3];
rz(-2.1761201) q[3];
sx q[3];
rz(1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.82764757) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(1.6754643) q[0];
rz(-0.77955359) q[1];
sx q[1];
rz(-1.164271) q[1];
sx q[1];
rz(2.4028042) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9290498) q[0];
sx q[0];
rz(-1.3669786) q[0];
sx q[0];
rz(-1.8575559) q[0];
x q[1];
rz(-3.1081057) q[2];
sx q[2];
rz(-2.8346363) q[2];
sx q[2];
rz(0.54881664) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14911095) q[1];
sx q[1];
rz(-0.68435366) q[1];
sx q[1];
rz(2.1696287) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8035369) q[3];
sx q[3];
rz(-2.6430686) q[3];
sx q[3];
rz(2.5132708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61591992) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(2.2789148) q[2];
rz(-0.45977965) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(-1.1427243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.90726844) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(1.020485) q[0];
rz(0.057295784) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(2.3540156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6660992) q[0];
sx q[0];
rz(-1.3331873) q[0];
sx q[0];
rz(2.4852024) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25813132) q[2];
sx q[2];
rz(-1.6399334) q[2];
sx q[2];
rz(3.0237232) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5912093) q[1];
sx q[1];
rz(-2.0616643) q[1];
sx q[1];
rz(2.9390013) q[1];
rz(-0.98290409) q[3];
sx q[3];
rz(-1.5421151) q[3];
sx q[3];
rz(0.3405638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82137498) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(-2.5569432) q[2];
rz(-0.65230495) q[3];
sx q[3];
rz(-1.9125331) q[3];
sx q[3];
rz(1.8023796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24340165) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(-0.32824326) q[0];
rz(0.39168656) q[1];
sx q[1];
rz(-1.1513386) q[1];
sx q[1];
rz(1.4788871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2770689) q[0];
sx q[0];
rz(-1.796699) q[0];
sx q[0];
rz(-0.27033131) q[0];
x q[1];
rz(2.4166862) q[2];
sx q[2];
rz(-0.6577684) q[2];
sx q[2];
rz(-2.7810514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6819671) q[1];
sx q[1];
rz(-2.8468968) q[1];
sx q[1];
rz(-1.2380283) q[1];
rz(-pi) q[2];
rz(-1.201466) q[3];
sx q[3];
rz(-1.7672667) q[3];
sx q[3];
rz(0.35614355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0216003) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(0.78869406) q[2];
rz(-2.9005519) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929844) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(-2.1797144) q[0];
rz(0.23421639) q[1];
sx q[1];
rz(-1.1385671) q[1];
sx q[1];
rz(2.403517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2099521) q[0];
sx q[0];
rz(-2.6745213) q[0];
sx q[0];
rz(-0.54801236) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5066824) q[2];
sx q[2];
rz(-1.2037841) q[2];
sx q[2];
rz(-0.5762595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.187584) q[1];
sx q[1];
rz(-2.6212647) q[1];
sx q[1];
rz(-2.1854758) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8515737) q[3];
sx q[3];
rz(-0.92377907) q[3];
sx q[3];
rz(-3.0386277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6568079) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(-1.2602932) q[2];
rz(1.8079181) q[3];
sx q[3];
rz(-2.1402054) q[3];
sx q[3];
rz(-2.8792152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(2.2743478) q[0];
rz(1.0137089) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(-1.0702466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4152756) q[0];
sx q[0];
rz(-1.7219825) q[0];
sx q[0];
rz(3.0180406) q[0];
x q[1];
rz(-0.81999166) q[2];
sx q[2];
rz(-2.242616) q[2];
sx q[2];
rz(-1.2206248) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6079191) q[1];
sx q[1];
rz(-1.1864788) q[1];
sx q[1];
rz(0.55748765) q[1];
rz(3.002043) q[3];
sx q[3];
rz(-1.680239) q[3];
sx q[3];
rz(2.8982996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16537198) q[2];
sx q[2];
rz(-0.82087159) q[2];
sx q[2];
rz(-1.9099859) q[2];
rz(-1.6407137) q[3];
sx q[3];
rz(-0.21017635) q[3];
sx q[3];
rz(0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(2.7267743) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(2.6957569) q[2];
sx q[2];
rz(-1.5210963) q[2];
sx q[2];
rz(-1.7901909) q[2];
rz(0.75178643) q[3];
sx q[3];
rz(-1.1204168) q[3];
sx q[3];
rz(0.71747019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
