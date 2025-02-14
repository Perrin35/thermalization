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
rz(2.4647291) q[0];
sx q[0];
rz(1.6114177) q[0];
sx q[0];
rz(9.7220698) q[0];
rz(-2.3625506) q[1];
sx q[1];
rz(-2.9809597) q[1];
sx q[1];
rz(-1.4359441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69383088) q[0];
sx q[0];
rz(-2.075747) q[0];
sx q[0];
rz(2.2044066) q[0];
x q[1];
rz(1.3117242) q[2];
sx q[2];
rz(-1.6951778) q[2];
sx q[2];
rz(2.4476515) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42987016) q[1];
sx q[1];
rz(-0.87583576) q[1];
sx q[1];
rz(1.9942787) q[1];
rz(1.7706031) q[3];
sx q[3];
rz(-1.8733896) q[3];
sx q[3];
rz(-1.2452919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5982738) q[2];
sx q[2];
rz(-2.7136549) q[2];
sx q[2];
rz(2.9738026) q[2];
rz(-1.1281891) q[3];
sx q[3];
rz(-1.5690683) q[3];
sx q[3];
rz(1.9523841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0447277) q[0];
sx q[0];
rz(-0.6701349) q[0];
sx q[0];
rz(1.0750394) q[0];
rz(-1.7754405) q[1];
sx q[1];
rz(-1.7551883) q[1];
sx q[1];
rz(-1.8278106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9843861) q[0];
sx q[0];
rz(-1.3493378) q[0];
sx q[0];
rz(-0.33672543) q[0];
rz(-1.7402252) q[2];
sx q[2];
rz(-1.7544978) q[2];
sx q[2];
rz(1.6970762) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3071482) q[1];
sx q[1];
rz(-1.6475369) q[1];
sx q[1];
rz(2.2952324) q[1];
x q[2];
rz(1.9017392) q[3];
sx q[3];
rz(-1.6972491) q[3];
sx q[3];
rz(-3.0883873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4904867) q[2];
sx q[2];
rz(-1.7650812) q[2];
sx q[2];
rz(0.42438486) q[2];
rz(2.4091447) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(-1.6781767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80377793) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(2.6388229) q[0];
rz(0.94364199) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(2.3263993) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73683263) q[0];
sx q[0];
rz(-1.6175999) q[0];
sx q[0];
rz(0.18808774) q[0];
x q[1];
rz(-1.4645459) q[2];
sx q[2];
rz(-1.4521404) q[2];
sx q[2];
rz(0.066421631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6289754) q[1];
sx q[1];
rz(-1.830258) q[1];
sx q[1];
rz(1.2995059) q[1];
rz(-0.56166537) q[3];
sx q[3];
rz(-0.9540073) q[3];
sx q[3];
rz(1.3514618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0649071) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(1.1265075) q[2];
rz(-1.3569776) q[3];
sx q[3];
rz(-2.6419736) q[3];
sx q[3];
rz(-3.0016541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55946881) q[0];
sx q[0];
rz(-2.421565) q[0];
sx q[0];
rz(-2.7816787) q[0];
rz(-0.24078807) q[1];
sx q[1];
rz(-1.1477995) q[1];
sx q[1];
rz(-0.37318939) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8561962) q[0];
sx q[0];
rz(-1.6699978) q[0];
sx q[0];
rz(-1.7147934) q[0];
rz(2.7794851) q[2];
sx q[2];
rz(-2.547764) q[2];
sx q[2];
rz(-2.1613742) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9993101) q[1];
sx q[1];
rz(-1.9110084) q[1];
sx q[1];
rz(-0.55026023) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3803102) q[3];
sx q[3];
rz(-1.130874) q[3];
sx q[3];
rz(-0.89118496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.084426247) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(0.90844321) q[2];
rz(-2.0716136) q[3];
sx q[3];
rz(-1.9570276) q[3];
sx q[3];
rz(2.158304) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1835566) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(0.68963447) q[1];
sx q[1];
rz(-0.90877405) q[1];
sx q[1];
rz(-0.77401179) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9871983) q[0];
sx q[0];
rz(-0.78332213) q[0];
sx q[0];
rz(-0.67449595) q[0];
rz(-pi) q[1];
rz(-1.6025324) q[2];
sx q[2];
rz(-2.3661315) q[2];
sx q[2];
rz(-2.8808978) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7237295) q[1];
sx q[1];
rz(-0.47342074) q[1];
sx q[1];
rz(-2.0000877) q[1];
rz(-pi) q[2];
rz(1.6004531) q[3];
sx q[3];
rz(-2.1184485) q[3];
sx q[3];
rz(-0.95049196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8959877) q[2];
sx q[2];
rz(-1.7481123) q[2];
sx q[2];
rz(2.4020933) q[2];
rz(-1.6351522) q[3];
sx q[3];
rz(-2.2061901) q[3];
sx q[3];
rz(0.80772775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6226115) q[0];
sx q[0];
rz(-2.3464572) q[0];
sx q[0];
rz(1.125289) q[0];
rz(0.13825026) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(-1.5743871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042717) q[0];
sx q[0];
rz(-1.9372276) q[0];
sx q[0];
rz(1.2211114) q[0];
x q[1];
rz(1.2786516) q[2];
sx q[2];
rz(-1.9685479) q[2];
sx q[2];
rz(1.5270555) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9003657) q[1];
sx q[1];
rz(-1.8292184) q[1];
sx q[1];
rz(1.6994121) q[1];
x q[2];
rz(-2.8254824) q[3];
sx q[3];
rz(-2.5694642) q[3];
sx q[3];
rz(-2.8993949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77341998) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(-1.0991905) q[2];
rz(0.98569551) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(-2.709008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2987591) q[0];
sx q[0];
rz(-0.99501139) q[0];
sx q[0];
rz(0.80879912) q[0];
rz(-0.66566268) q[1];
sx q[1];
rz(-0.97949615) q[1];
sx q[1];
rz(2.1601802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2546651) q[0];
sx q[0];
rz(-1.5136295) q[0];
sx q[0];
rz(-2.1421823) q[0];
rz(-pi) q[1];
rz(-2.7000325) q[2];
sx q[2];
rz(-2.5924263) q[2];
sx q[2];
rz(1.5541981) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7715523) q[1];
sx q[1];
rz(-1.1070894) q[1];
sx q[1];
rz(-0.15180219) q[1];
rz(-pi) q[2];
rz(1.3990226) q[3];
sx q[3];
rz(-2.4618759) q[3];
sx q[3];
rz(1.8006067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0770646) q[2];
sx q[2];
rz(-1.6669824) q[2];
sx q[2];
rz(2.7725753) q[2];
rz(-1.7391694) q[3];
sx q[3];
rz(-2.2544421) q[3];
sx q[3];
rz(1.8203452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8333261) q[0];
sx q[0];
rz(-0.82887355) q[0];
sx q[0];
rz(0.33018026) q[0];
rz(0.45686832) q[1];
sx q[1];
rz(-2.7062682) q[1];
sx q[1];
rz(-0.59828573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.79222) q[0];
sx q[0];
rz(-2.0061919) q[0];
sx q[0];
rz(-2.4382082) q[0];
x q[1];
rz(3.0706579) q[2];
sx q[2];
rz(-1.1026376) q[2];
sx q[2];
rz(2.2092961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1674211) q[1];
sx q[1];
rz(-1.9037013) q[1];
sx q[1];
rz(1.3281338) q[1];
rz(0.71514327) q[3];
sx q[3];
rz(-2.3480881) q[3];
sx q[3];
rz(-1.8737829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8699708) q[2];
sx q[2];
rz(-0.87646708) q[2];
sx q[2];
rz(-3.0769707) q[2];
rz(-1.4501075) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(1.5456642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2794063) q[0];
sx q[0];
rz(-2.8074844) q[0];
sx q[0];
rz(2.620328) q[0];
rz(-2.0817256) q[1];
sx q[1];
rz(-1.9797378) q[1];
sx q[1];
rz(1.0313787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8143896) q[0];
sx q[0];
rz(-1.4966244) q[0];
sx q[0];
rz(-2.402284) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8487186) q[2];
sx q[2];
rz(-1.4559064) q[2];
sx q[2];
rz(3.1174768) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65341144) q[1];
sx q[1];
rz(-2.0926884) q[1];
sx q[1];
rz(-2.8074991) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6645722) q[3];
sx q[3];
rz(-1.4833772) q[3];
sx q[3];
rz(-1.9893579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3970268) q[2];
sx q[2];
rz(-2.4137745) q[2];
sx q[2];
rz(0.087372027) q[2];
rz(1.9258063) q[3];
sx q[3];
rz(-1.4601424) q[3];
sx q[3];
rz(2.9714835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755197) q[0];
sx q[0];
rz(-0.89101321) q[0];
sx q[0];
rz(1.1391621) q[0];
rz(0.49925223) q[1];
sx q[1];
rz(-1.6923994) q[1];
sx q[1];
rz(1.3334691) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77908726) q[0];
sx q[0];
rz(-1.6769857) q[0];
sx q[0];
rz(1.4800735) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0917909) q[2];
sx q[2];
rz(-1.3161471) q[2];
sx q[2];
rz(1.7083502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.002961) q[1];
sx q[1];
rz(-1.3182782) q[1];
sx q[1];
rz(1.7934402) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3995497) q[3];
sx q[3];
rz(-2.2931406) q[3];
sx q[3];
rz(2.6076041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.501005) q[2];
sx q[2];
rz(-1.86684) q[2];
sx q[2];
rz(-2.8583543) q[2];
rz(-1.9404274) q[3];
sx q[3];
rz(-0.675942) q[3];
sx q[3];
rz(-1.0454319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82861154) q[0];
sx q[0];
rz(-0.77684488) q[0];
sx q[0];
rz(-2.0364398) q[0];
rz(-0.10440566) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(-2.7945065) q[2];
sx q[2];
rz(-0.53205678) q[2];
sx q[2];
rz(2.2126641) q[2];
rz(1.1558044) q[3];
sx q[3];
rz(-0.94528755) q[3];
sx q[3];
rz(-0.62936924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
