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
rz(-2.3053919) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(6.3451938) q[1];
sx q[1];
rz(11.558029) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22577408) q[0];
sx q[0];
rz(-1.4918707) q[0];
sx q[0];
rz(1.945334) q[0];
x q[1];
rz(-0.37990976) q[2];
sx q[2];
rz(-1.3515944) q[2];
sx q[2];
rz(0.70218147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5656149) q[1];
sx q[1];
rz(-1.7276141) q[1];
sx q[1];
rz(-0.48139907) q[1];
rz(-pi) q[2];
x q[2];
rz(3.120901) q[3];
sx q[3];
rz(-1.3776099) q[3];
sx q[3];
rz(0.26259804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36650518) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.7979701) q[0];
rz(1.1603181) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(0.62058273) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40784971) q[0];
sx q[0];
rz(-1.3792999) q[0];
sx q[0];
rz(1.922185) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2389675) q[2];
sx q[2];
rz(-1.2764954) q[2];
sx q[2];
rz(-2.6455392) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0479792) q[1];
sx q[1];
rz(-1.9318214) q[1];
sx q[1];
rz(-2.4314636) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1813404) q[3];
sx q[3];
rz(-0.6529633) q[3];
sx q[3];
rz(3.0137526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4429984) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(2.911574) q[2];
rz(1.5551785) q[3];
sx q[3];
rz(-1.14862) q[3];
sx q[3];
rz(-0.27568278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8860633) q[0];
sx q[0];
rz(-0.38917381) q[0];
sx q[0];
rz(2.0607167) q[0];
rz(-0.64322645) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(-3.0562775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9935308) q[0];
sx q[0];
rz(-0.9328273) q[0];
sx q[0];
rz(-2.30739) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32344748) q[2];
sx q[2];
rz(-2.5451042) q[2];
sx q[2];
rz(-1.310854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6416723) q[1];
sx q[1];
rz(-2.9292078) q[1];
sx q[1];
rz(2.4714064) q[1];
x q[2];
rz(1.5253228) q[3];
sx q[3];
rz(-2.1915428) q[3];
sx q[3];
rz(-2.5705702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9366511) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(0.19634518) q[2];
rz(-2.8593072) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(1.0452247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.008217) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(-0.38544449) q[0];
rz(-1.4133981) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(0.24994303) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2521533) q[0];
sx q[0];
rz(-1.1595386) q[0];
sx q[0];
rz(2.2338339) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6310299) q[2];
sx q[2];
rz(-1.4830198) q[2];
sx q[2];
rz(0.46030948) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.476772) q[1];
sx q[1];
rz(-2.6198434) q[1];
sx q[1];
rz(0.34362766) q[1];
x q[2];
rz(-1.3718666) q[3];
sx q[3];
rz(-0.78708157) q[3];
sx q[3];
rz(2.6347292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46127737) q[2];
sx q[2];
rz(-1.8467434) q[2];
sx q[2];
rz(-0.2363905) q[2];
rz(1.8810898) q[3];
sx q[3];
rz(-0.25006306) q[3];
sx q[3];
rz(-2.4642956) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054955) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(2.6564823) q[0];
rz(1.9612034) q[1];
sx q[1];
rz(-1.5943269) q[1];
sx q[1];
rz(-2.3050883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.361995) q[0];
sx q[0];
rz(-1.2149724) q[0];
sx q[0];
rz(-1.9076288) q[0];
x q[1];
rz(2.8233158) q[2];
sx q[2];
rz(-1.1851382) q[2];
sx q[2];
rz(1.9463271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4662469) q[1];
sx q[1];
rz(-0.9017749) q[1];
sx q[1];
rz(0.061208486) q[1];
rz(-pi) q[2];
rz(-1.6338634) q[3];
sx q[3];
rz(-1.1707889) q[3];
sx q[3];
rz(-2.0162752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2430719) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(-0.76954976) q[2];
rz(2.2122993) q[3];
sx q[3];
rz(-2.010689) q[3];
sx q[3];
rz(0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9418075) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(-0.76882452) q[0];
rz(1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(-2.2519055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.453196) q[0];
sx q[0];
rz(-1.9510117) q[0];
sx q[0];
rz(-2.9158457) q[0];
x q[1];
rz(1.6011997) q[2];
sx q[2];
rz(-0.044375751) q[2];
sx q[2];
rz(0.45472431) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5714214) q[1];
sx q[1];
rz(-1.716905) q[1];
sx q[1];
rz(0.04695462) q[1];
x q[2];
rz(1.6857288) q[3];
sx q[3];
rz(-2.0191666) q[3];
sx q[3];
rz(1.2796206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49030226) q[2];
sx q[2];
rz(-1.7040665) q[2];
sx q[2];
rz(1.5550522) q[2];
rz(1.2706612) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(0.13203013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1726058) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(-1.9975115) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(-2.3544618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0486271) q[0];
sx q[0];
rz(-2.9529609) q[0];
sx q[0];
rz(-0.79690551) q[0];
x q[1];
rz(0.88655858) q[2];
sx q[2];
rz(-1.1349003) q[2];
sx q[2];
rz(-0.8800216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5980144) q[1];
sx q[1];
rz(-2.0435963) q[1];
sx q[1];
rz(0.531457) q[1];
x q[2];
rz(0.58148099) q[3];
sx q[3];
rz(-2.498501) q[3];
sx q[3];
rz(-3.0742925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0234915) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-3.1084295) q[2];
rz(1.5432594) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(-1.5020812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.892266) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(-1.3931042) q[0];
rz(-2.6847367) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(2.0410062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8580061) q[0];
sx q[0];
rz(-1.4947944) q[0];
sx q[0];
rz(1.5690687) q[0];
rz(-pi) q[1];
rz(2.4462561) q[2];
sx q[2];
rz(-1.0461263) q[2];
sx q[2];
rz(-2.6161043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7707899) q[1];
sx q[1];
rz(-2.2375771) q[1];
sx q[1];
rz(0.043170269) q[1];
x q[2];
rz(-1.4310775) q[3];
sx q[3];
rz(-1.9587687) q[3];
sx q[3];
rz(3.1041077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0988203) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(1.5607321) q[2];
rz(-2.3780195) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(2.1055351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40912691) q[0];
sx q[0];
rz(-0.75508535) q[0];
sx q[0];
rz(-0.29148802) q[0];
rz(1.905722) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(3.018697) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50704573) q[0];
sx q[0];
rz(-1.6222008) q[0];
sx q[0];
rz(2.794751) q[0];
rz(-2.2332623) q[2];
sx q[2];
rz(-0.35229063) q[2];
sx q[2];
rz(0.34991821) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5220241) q[1];
sx q[1];
rz(-1.777473) q[1];
sx q[1];
rz(0.11264888) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-1.5767617) q[3];
sx q[3];
rz(0.98814135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5128532) q[2];
sx q[2];
rz(-1.7491919) q[2];
sx q[2];
rz(-2.0780308) q[2];
rz(2.9641446) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.5408632) q[1];
sx q[1];
rz(-2.6917916) q[1];
sx q[1];
rz(-2.5240555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4040632) q[0];
sx q[0];
rz(-1.6949708) q[0];
sx q[0];
rz(-1.7880718) q[0];
x q[1];
rz(0.82771684) q[2];
sx q[2];
rz(-1.9736145) q[2];
sx q[2];
rz(0.043443505) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2899298) q[1];
sx q[1];
rz(-1.5917707) q[1];
sx q[1];
rz(2.1037654) q[1];
rz(-pi) q[2];
rz(-2.6838949) q[3];
sx q[3];
rz(-0.62333306) q[3];
sx q[3];
rz(1.8231572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2403229) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(2.1642302) q[2];
rz(-1.0824925) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(-0.055518363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8025773) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(3.0958685) q[1];
sx q[1];
rz(-2.3241691) q[1];
sx q[1];
rz(-1.587422) q[1];
rz(2.1004213) q[2];
sx q[2];
rz(-1.1239792) q[2];
sx q[2];
rz(-1.299528) q[2];
rz(2.9228899) q[3];
sx q[3];
rz(-0.66853157) q[3];
sx q[3];
rz(2.4562277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
