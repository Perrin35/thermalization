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
rz(-1.530175) q[0];
sx q[0];
rz(-0.29729182) q[0];
rz(0.77904207) q[1];
sx q[1];
rz(-0.160633) q[1];
sx q[1];
rz(1.4359441) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061062) q[0];
sx q[0];
rz(-2.1155042) q[0];
sx q[0];
rz(2.5404055) q[0];
x q[1];
rz(-1.8298684) q[2];
sx q[2];
rz(-1.4464149) q[2];
sx q[2];
rz(0.69394116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2816726) q[1];
sx q[1];
rz(-1.249673) q[1];
sx q[1];
rz(-0.74076498) q[1];
rz(-pi) q[2];
rz(0.30835704) q[3];
sx q[3];
rz(-1.7614109) q[3];
sx q[3];
rz(2.7558143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5433189) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(-2.9738026) q[2];
rz(-2.0134036) q[3];
sx q[3];
rz(-1.5690683) q[3];
sx q[3];
rz(1.1892085) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096864916) q[0];
sx q[0];
rz(-0.6701349) q[0];
sx q[0];
rz(1.0750394) q[0];
rz(1.7754405) q[1];
sx q[1];
rz(-1.3864044) q[1];
sx q[1];
rz(-1.8278106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9843861) q[0];
sx q[0];
rz(-1.3493378) q[0];
sx q[0];
rz(-0.33672543) q[0];
rz(1.4013675) q[2];
sx q[2];
rz(-1.3870948) q[2];
sx q[2];
rz(1.4445164) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7915276) q[1];
sx q[1];
rz(-0.72775562) q[1];
sx q[1];
rz(1.6863053) q[1];
x q[2];
rz(-1.2398534) q[3];
sx q[3];
rz(-1.4443436) q[3];
sx q[3];
rz(-0.053205333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65110597) q[2];
sx q[2];
rz(-1.3765114) q[2];
sx q[2];
rz(-0.42438486) q[2];
rz(-2.4091447) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(-1.463416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80377793) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(2.6388229) q[0];
rz(-2.1979507) q[1];
sx q[1];
rz(-2.4991401) q[1];
sx q[1];
rz(0.81519333) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73683263) q[0];
sx q[0];
rz(-1.5239927) q[0];
sx q[0];
rz(-2.9535049) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6770467) q[2];
sx q[2];
rz(-1.6894522) q[2];
sx q[2];
rz(-0.066421631) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3383662) q[1];
sx q[1];
rz(-2.7684281) q[1];
sx q[1];
rz(0.79014059) q[1];
x q[2];
rz(2.5799273) q[3];
sx q[3];
rz(-2.1875854) q[3];
sx q[3];
rz(-1.3514618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0649071) q[2];
sx q[2];
rz(-1.6904597) q[2];
sx q[2];
rz(1.1265075) q[2];
rz(1.3569776) q[3];
sx q[3];
rz(-0.49961909) q[3];
sx q[3];
rz(0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821238) q[0];
sx q[0];
rz(-2.421565) q[0];
sx q[0];
rz(2.7816787) q[0];
rz(2.9008046) q[1];
sx q[1];
rz(-1.1477995) q[1];
sx q[1];
rz(-0.37318939) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2853965) q[0];
sx q[0];
rz(-1.4715949) q[0];
sx q[0];
rz(-1.7147934) q[0];
rz(2.7794851) q[2];
sx q[2];
rz(-2.547764) q[2];
sx q[2];
rz(-2.1613742) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2265985) q[1];
sx q[1];
rz(-1.0553331) q[1];
sx q[1];
rz(1.1771918) q[1];
rz(-pi) q[2];
rz(-1.3803102) q[3];
sx q[3];
rz(-2.0107186) q[3];
sx q[3];
rz(2.2504077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0571664) q[2];
sx q[2];
rz(-1.4320222) q[2];
sx q[2];
rz(0.90844321) q[2];
rz(-2.0716136) q[3];
sx q[3];
rz(-1.1845651) q[3];
sx q[3];
rz(0.98328868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95803607) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(-0.91304427) q[0];
rz(-0.68963447) q[1];
sx q[1];
rz(-2.2328186) q[1];
sx q[1];
rz(-0.77401179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15439437) q[0];
sx q[0];
rz(-2.3582705) q[0];
sx q[0];
rz(2.4670967) q[0];
rz(-pi) q[1];
rz(0.795587) q[2];
sx q[2];
rz(-1.5485816) q[2];
sx q[2];
rz(-1.8088248) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7237295) q[1];
sx q[1];
rz(-2.6681719) q[1];
sx q[1];
rz(1.141505) q[1];
x q[2];
rz(1.6004531) q[3];
sx q[3];
rz(-2.1184485) q[3];
sx q[3];
rz(-0.95049196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.245605) q[2];
sx q[2];
rz(-1.7481123) q[2];
sx q[2];
rz(-0.73949933) q[2];
rz(1.5064404) q[3];
sx q[3];
rz(-2.2061901) q[3];
sx q[3];
rz(-2.3338649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5189811) q[0];
sx q[0];
rz(-0.7951355) q[0];
sx q[0];
rz(-2.0163037) q[0];
rz(-3.0033424) q[1];
sx q[1];
rz(-1.6743276) q[1];
sx q[1];
rz(1.5743871) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7049887) q[0];
sx q[0];
rz(-1.2452176) q[0];
sx q[0];
rz(-2.7537936) q[0];
x q[1];
rz(-0.60092314) q[2];
sx q[2];
rz(-2.6527361) q[2];
sx q[2];
rz(2.2747441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7789845) q[1];
sx q[1];
rz(-1.6951188) q[1];
sx q[1];
rz(-0.26047867) q[1];
x q[2];
rz(-1.7683783) q[3];
sx q[3];
rz(-2.1113331) q[3];
sx q[3];
rz(0.61321248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.77341998) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(-2.0424021) q[2];
rz(0.98569551) q[3];
sx q[3];
rz(-1.516195) q[3];
sx q[3];
rz(2.709008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8428335) q[0];
sx q[0];
rz(-0.99501139) q[0];
sx q[0];
rz(0.80879912) q[0];
rz(-2.47593) q[1];
sx q[1];
rz(-2.1620965) q[1];
sx q[1];
rz(2.1601802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8869276) q[0];
sx q[0];
rz(-1.6279632) q[0];
sx q[0];
rz(-0.99941038) q[0];
rz(0.44156011) q[2];
sx q[2];
rz(-0.54916635) q[2];
sx q[2];
rz(1.5873945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26907193) q[1];
sx q[1];
rz(-1.4351294) q[1];
sx q[1];
rz(1.1024464) q[1];
rz(-pi) q[2];
rz(2.2432765) q[3];
sx q[3];
rz(-1.6784462) q[3];
sx q[3];
rz(3.0458991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0770646) q[2];
sx q[2];
rz(-1.4746102) q[2];
sx q[2];
rz(0.36901739) q[2];
rz(-1.4024233) q[3];
sx q[3];
rz(-0.88715059) q[3];
sx q[3];
rz(1.8203452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8333261) q[0];
sx q[0];
rz(-2.3127191) q[0];
sx q[0];
rz(2.8114124) q[0];
rz(0.45686832) q[1];
sx q[1];
rz(-2.7062682) q[1];
sx q[1];
rz(-0.59828573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5766524) q[0];
sx q[0];
rz(-2.1974753) q[0];
sx q[0];
rz(2.1184854) q[0];
x q[1];
rz(3.0706579) q[2];
sx q[2];
rz(-2.0389551) q[2];
sx q[2];
rz(-2.2092961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9741716) q[1];
sx q[1];
rz(-1.2378913) q[1];
sx q[1];
rz(1.3281338) q[1];
rz(-pi) q[2];
x q[2];
rz(0.65450683) q[3];
sx q[3];
rz(-1.084436) q[3];
sx q[3];
rz(2.8974722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.8699708) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(-3.0769707) q[2];
rz(1.4501075) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(1.5959285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2794063) q[0];
sx q[0];
rz(-0.33410826) q[0];
sx q[0];
rz(0.52126467) q[0];
rz(2.0817256) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(-2.1102139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17612621) q[0];
sx q[0];
rz(-0.83399189) q[0];
sx q[0];
rz(1.6710207) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3802558) q[2];
sx q[2];
rz(-0.31399841) q[2];
sx q[2];
rz(1.2316201) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.65341144) q[1];
sx q[1];
rz(-2.0926884) q[1];
sx q[1];
rz(2.8074991) q[1];
rz(-2.9529753) q[3];
sx q[3];
rz(-0.48435703) q[3];
sx q[3];
rz(-0.25121197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3970268) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(3.0542206) q[2];
rz(1.9258063) q[3];
sx q[3];
rz(-1.6814503) q[3];
sx q[3];
rz(0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(2.4755197) q[0];
sx q[0];
rz(-2.2505794) q[0];
sx q[0];
rz(1.1391621) q[0];
rz(-2.6423404) q[1];
sx q[1];
rz(-1.4491932) q[1];
sx q[1];
rz(1.8081236) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3402417) q[0];
sx q[0];
rz(-1.4805859) q[0];
sx q[0];
rz(0.10662459) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.75976) q[2];
sx q[2];
rz(-0.25936959) q[2];
sx q[2];
rz(1.237902) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.002961) q[1];
sx q[1];
rz(-1.3182782) q[1];
sx q[1];
rz(-1.3481524) q[1];
rz(-0.69656792) q[3];
sx q[3];
rz(-1.0390716) q[3];
sx q[3];
rz(1.5817489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6405876) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(-0.28323832) q[2];
rz(-1.9404274) q[3];
sx q[3];
rz(-0.675942) q[3];
sx q[3];
rz(-1.0454319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82861154) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(-0.10440566) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(-0.5055867) q[2];
sx q[2];
rz(-1.7442295) q[2];
sx q[2];
rz(-2.8019047) q[2];
rz(1.9857882) q[3];
sx q[3];
rz(-2.1963051) q[3];
sx q[3];
rz(2.5122234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
