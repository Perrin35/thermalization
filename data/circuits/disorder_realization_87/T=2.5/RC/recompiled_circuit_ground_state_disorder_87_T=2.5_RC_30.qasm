OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(-0.41271451) q[0];
sx q[0];
rz(0.8362008) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(-3.0795842) q[1];
sx q[1];
rz(-2.1332512) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8275534) q[0];
sx q[0];
rz(-1.1974821) q[0];
sx q[0];
rz(0.084777431) q[0];
rz(2.600619) q[2];
sx q[2];
rz(-0.43593513) q[2];
sx q[2];
rz(1.774314) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.086585933) q[1];
sx q[1];
rz(-1.0957967) q[1];
sx q[1];
rz(-1.7473298) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.020691607) q[3];
sx q[3];
rz(-1.3776099) q[3];
sx q[3];
rz(0.26259804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3081554) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(2.4885139) q[2];
rz(3.0270882) q[3];
sx q[3];
rz(-1.3936035) q[3];
sx q[3];
rz(1.0623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36650518) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.3436226) q[0];
rz(-1.9812745) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(0.62058273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4575093) q[0];
sx q[0];
rz(-2.7433222) q[0];
sx q[0];
rz(2.0837559) q[0];
x q[1];
rz(2.0257904) q[2];
sx q[2];
rz(-0.7209076) q[2];
sx q[2];
rz(-1.7146595) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1329886) q[1];
sx q[1];
rz(-0.78215608) q[1];
sx q[1];
rz(0.52495606) q[1];
rz(-pi) q[2];
rz(0.95497953) q[3];
sx q[3];
rz(-1.338025) q[3];
sx q[3];
rz(-1.3835088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6985942) q[2];
sx q[2];
rz(-1.5195165) q[2];
sx q[2];
rz(0.2300187) q[2];
rz(-1.5864141) q[3];
sx q[3];
rz(-1.9929726) q[3];
sx q[3];
rz(0.27568278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8860633) q[0];
sx q[0];
rz(-2.7524188) q[0];
sx q[0];
rz(2.0607167) q[0];
rz(-2.4983662) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(-0.08531514) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9179419) q[0];
sx q[0];
rz(-2.1408014) q[0];
sx q[0];
rz(0.78582363) q[0];
x q[1];
rz(-0.32344748) q[2];
sx q[2];
rz(-2.5451042) q[2];
sx q[2];
rz(-1.310854) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9602365) q[1];
sx q[1];
rz(-1.7367558) q[1];
sx q[1];
rz(1.4376498) q[1];
rz(-pi) q[2];
rz(-3.078104) q[3];
sx q[3];
rz(-0.62219071) q[3];
sx q[3];
rz(2.6486462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2049415) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(0.19634518) q[2];
rz(-0.28228545) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(-1.0452247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.008217) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(-2.7561482) q[0];
rz(-1.4133981) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(-2.8916496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37873822) q[0];
sx q[0];
rz(-2.1702499) q[0];
sx q[0];
rz(2.6361639) q[0];
rz(-pi) q[1];
rz(-1.6310299) q[2];
sx q[2];
rz(-1.4830198) q[2];
sx q[2];
rz(-2.6812832) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.085332) q[1];
sx q[1];
rz(-1.0823432) q[1];
sx q[1];
rz(-1.7621098) q[1];
x q[2];
rz(-0.79367359) q[3];
sx q[3];
rz(-1.4303615) q[3];
sx q[3];
rz(-2.2190203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46127737) q[2];
sx q[2];
rz(-1.8467434) q[2];
sx q[2];
rz(-0.2363905) q[2];
rz(1.8810898) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(-0.67729706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054955) q[0];
sx q[0];
rz(-1.6039055) q[0];
sx q[0];
rz(-2.6564823) q[0];
rz(-1.1803892) q[1];
sx q[1];
rz(-1.5472658) q[1];
sx q[1];
rz(2.3050883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0083975365) q[0];
sx q[0];
rz(-0.48497619) q[0];
sx q[0];
rz(0.72686813) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2274687) q[2];
sx q[2];
rz(-2.6466922) q[2];
sx q[2];
rz(-2.6656541) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0750351) q[1];
sx q[1];
rz(-1.6187985) q[1];
sx q[1];
rz(-0.90086277) q[1];
rz(-pi) q[2];
rz(2.740871) q[3];
sx q[3];
rz(-1.6288789) q[3];
sx q[3];
rz(-0.42089128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.89852077) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(0.76954976) q[2];
rz(-0.92929333) q[3];
sx q[3];
rz(-2.010689) q[3];
sx q[3];
rz(0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418075) q[0];
sx q[0];
rz(-2.6564044) q[0];
sx q[0];
rz(-0.76882452) q[0];
rz(1.7898412) q[1];
sx q[1];
rz(-2.6685346) q[1];
sx q[1];
rz(0.88968712) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3442197) q[0];
sx q[0];
rz(-1.780172) q[0];
sx q[0];
rz(-1.9599509) q[0];
rz(-pi) q[1];
rz(-1.6151516) q[2];
sx q[2];
rz(-1.5721448) q[2];
sx q[2];
rz(-1.1464455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1353768) q[1];
sx q[1];
rz(-1.6172503) q[1];
sx q[1];
rz(1.7170639) q[1];
rz(1.6857288) q[3];
sx q[3];
rz(-1.1224261) q[3];
sx q[3];
rz(1.861972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6512904) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(1.5550522) q[2];
rz(-1.2706612) q[3];
sx q[3];
rz(-0.41897604) q[3];
sx q[3];
rz(-3.0095625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9689869) q[0];
sx q[0];
rz(-1.1382599) q[0];
sx q[0];
rz(1.1440811) q[0];
rz(-2.3106958) q[1];
sx q[1];
rz(-1.3583207) q[1];
sx q[1];
rz(2.3544618) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09296552) q[0];
sx q[0];
rz(-2.9529609) q[0];
sx q[0];
rz(0.79690551) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5412124) q[2];
sx q[2];
rz(-2.1809309) q[2];
sx q[2];
rz(-1.0224563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28879189) q[1];
sx q[1];
rz(-2.0388985) q[1];
sx q[1];
rz(-2.1062984) q[1];
rz(-pi) q[2];
rz(2.5821286) q[3];
sx q[3];
rz(-1.9064404) q[3];
sx q[3];
rz(1.0192724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1181011) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-0.033163158) q[2];
rz(-1.5432594) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24932662) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(1.7484885) q[0];
rz(-2.6847367) q[1];
sx q[1];
rz(-1.9970048) q[1];
sx q[1];
rz(-2.0410062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7129214) q[0];
sx q[0];
rz(-1.5690737) q[0];
sx q[0];
rz(3.0655906) q[0];
x q[1];
rz(0.69533657) q[2];
sx q[2];
rz(-1.0461263) q[2];
sx q[2];
rz(2.6161043) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.17328429) q[1];
sx q[1];
rz(-1.5368764) q[1];
sx q[1];
rz(-2.2380301) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32842095) q[3];
sx q[3];
rz(-2.7304318) q[3];
sx q[3];
rz(-0.31842768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0988203) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(1.5607321) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.5127425) q[3];
sx q[3];
rz(2.1055351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7324657) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(0.29148802) q[0];
rz(1.2358707) q[1];
sx q[1];
rz(-1.196685) q[1];
sx q[1];
rz(3.018697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0451806) q[0];
sx q[0];
rz(-1.9171606) q[0];
sx q[0];
rz(-1.6254495) q[0];
x q[1];
rz(-2.9192186) q[2];
sx q[2];
rz(-1.8463328) q[2];
sx q[2];
rz(2.0982519) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0258515) q[1];
sx q[1];
rz(-2.9065955) q[1];
sx q[1];
rz(-2.0629289) q[1];
x q[2];
rz(0.0060283383) q[3];
sx q[3];
rz(-1.42618) q[3];
sx q[3];
rz(-2.5580689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5128532) q[2];
sx q[2];
rz(-1.7491919) q[2];
sx q[2];
rz(1.0635618) q[2];
rz(-2.9641446) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(-1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7040831) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(-1.4310687) q[0];
rz(0.60072947) q[1];
sx q[1];
rz(-0.44980106) q[1];
sx q[1];
rz(-2.5240555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3356614) q[0];
sx q[0];
rz(-1.3552203) q[0];
sx q[0];
rz(-0.12713253) q[0];
x q[1];
rz(2.6170116) q[2];
sx q[2];
rz(-2.2426105) q[2];
sx q[2];
rz(1.8730522) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75467089) q[1];
sx q[1];
rz(-0.53334177) q[1];
sx q[1];
rz(1.6120595) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5687739) q[3];
sx q[3];
rz(-1.8316934) q[3];
sx q[3];
rz(0.12810055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90126976) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(-2.1642302) q[2];
rz(1.0824925) q[3];
sx q[3];
rz(-0.87477028) q[3];
sx q[3];
rz(-0.055518363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8025773) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(-0.045724178) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(0.81188079) q[2];
sx q[2];
rz(-2.4627081) q[2];
sx q[2];
rz(-2.2344786) q[2];
rz(-2.4847538) q[3];
sx q[3];
rz(-1.7056864) q[3];
sx q[3];
rz(1.0581072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
