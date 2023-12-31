OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6264412) q[0];
sx q[0];
rz(-3.1080973) q[0];
sx q[0];
rz(-1.7749696) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(-1.4895952) q[1];
sx q[1];
rz(1.1319914) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223335) q[0];
sx q[0];
rz(-0.83689892) q[0];
sx q[0];
rz(2.0440621) q[0];
x q[1];
rz(-0.16146026) q[2];
sx q[2];
rz(-1.4070639) q[2];
sx q[2];
rz(1.8979567) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0238029) q[1];
sx q[1];
rz(-2.1315247) q[1];
sx q[1];
rz(2.6610713) q[1];
rz(-2.6996783) q[3];
sx q[3];
rz(-2.5406197) q[3];
sx q[3];
rz(2.0244983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91036096) q[2];
sx q[2];
rz(-1.8138764) q[2];
sx q[2];
rz(1.988391) q[2];
rz(-0.48405805) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(0.4593862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.83950481) q[0];
sx q[0];
rz(-0.33058259) q[0];
sx q[0];
rz(0.50305811) q[0];
rz(1.5867651) q[1];
sx q[1];
rz(-0.70650548) q[1];
sx q[1];
rz(0.15393004) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0921558) q[0];
sx q[0];
rz(-1.0679809) q[0];
sx q[0];
rz(-3.1252607) q[0];
x q[1];
rz(-0.92127992) q[2];
sx q[2];
rz(-1.3037455) q[2];
sx q[2];
rz(-0.18077476) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1068374) q[1];
sx q[1];
rz(-1.1498067) q[1];
sx q[1];
rz(-1.0616598) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5716343) q[3];
sx q[3];
rz(-1.0610233) q[3];
sx q[3];
rz(1.9364995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8872035) q[2];
sx q[2];
rz(-2.7837191) q[2];
sx q[2];
rz(-1.8187693) q[2];
rz(1.6555188) q[3];
sx q[3];
rz(-1.6008987) q[3];
sx q[3];
rz(-0.71050182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94674295) q[0];
sx q[0];
rz(-1.1179593) q[0];
sx q[0];
rz(0.90993607) q[0];
rz(-2.3643156) q[1];
sx q[1];
rz(-0.83559075) q[1];
sx q[1];
rz(-0.98532239) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9604208) q[0];
sx q[0];
rz(-2.4665678) q[0];
sx q[0];
rz(1.8280562) q[0];
x q[1];
rz(1.9820205) q[2];
sx q[2];
rz(-0.69934884) q[2];
sx q[2];
rz(1.5027836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31919033) q[1];
sx q[1];
rz(-1.5158347) q[1];
sx q[1];
rz(-2.8922006) q[1];
rz(0.76936929) q[3];
sx q[3];
rz(-2.3017075) q[3];
sx q[3];
rz(-1.4589256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.862792) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(1.2476236) q[2];
rz(2.6990081) q[3];
sx q[3];
rz(-1.7740039) q[3];
sx q[3];
rz(-0.32143337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9001532) q[0];
sx q[0];
rz(-1.0192008) q[0];
sx q[0];
rz(-2.0181657) q[0];
rz(-2.5627047) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(-1.7480063) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14560315) q[0];
sx q[0];
rz(-1.1568501) q[0];
sx q[0];
rz(-2.4664509) q[0];
x q[1];
rz(-0.43487866) q[2];
sx q[2];
rz(-1.7231427) q[2];
sx q[2];
rz(0.16532126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.77809282) q[1];
sx q[1];
rz(-2.2908195) q[1];
sx q[1];
rz(-2.5219445) q[1];
x q[2];
rz(0.23354236) q[3];
sx q[3];
rz(-0.84905784) q[3];
sx q[3];
rz(-0.28802179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0397772) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(0.0021136443) q[2];
rz(-0.56143108) q[3];
sx q[3];
rz(-1.0174454) q[3];
sx q[3];
rz(2.5533365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11557065) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(-1.9157238) q[0];
rz(1.6745802) q[1];
sx q[1];
rz(-1.8672698) q[1];
sx q[1];
rz(1.7747169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850692) q[0];
sx q[0];
rz(-1.9559304) q[0];
sx q[0];
rz(-1.7760081) q[0];
rz(-pi) q[1];
x q[1];
rz(0.032614313) q[2];
sx q[2];
rz(-2.3941052) q[2];
sx q[2];
rz(-0.05687296) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67147672) q[1];
sx q[1];
rz(-1.4911545) q[1];
sx q[1];
rz(-0.86298841) q[1];
rz(-pi) q[2];
rz(-0.59646888) q[3];
sx q[3];
rz(-2.5554552) q[3];
sx q[3];
rz(0.49367192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(-0.40536353) q[2];
rz(-0.46164414) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(-1.5464787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49155238) q[0];
sx q[0];
rz(-1.4251645) q[0];
sx q[0];
rz(2.6382085) q[0];
rz(-2.9227496) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(2.4898081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3394649) q[0];
sx q[0];
rz(-1.8101242) q[0];
sx q[0];
rz(1.375074) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82181828) q[2];
sx q[2];
rz(-1.6449252) q[2];
sx q[2];
rz(0.2766343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.3411322) q[1];
sx q[1];
rz(-2.9974077) q[1];
sx q[1];
rz(1.202281) q[1];
rz(-1.6739453) q[3];
sx q[3];
rz(-1.2892937) q[3];
sx q[3];
rz(-0.75479773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6289604) q[2];
sx q[2];
rz(-1.7444538) q[2];
sx q[2];
rz(-2.7499278) q[2];
rz(2.9351249) q[3];
sx q[3];
rz(-2.4075017) q[3];
sx q[3];
rz(-2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7713292) q[0];
sx q[0];
rz(-1.6917133) q[0];
sx q[0];
rz(2.1719334) q[0];
rz(-0.58352739) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-0.13024174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1657432) q[0];
sx q[0];
rz(-1.2184869) q[0];
sx q[0];
rz(1.5177112) q[0];
rz(-0.98353705) q[2];
sx q[2];
rz(-2.1526255) q[2];
sx q[2];
rz(-2.8772417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.817037) q[1];
sx q[1];
rz(-2.5845924) q[1];
sx q[1];
rz(0.42028285) q[1];
rz(-2.7534915) q[3];
sx q[3];
rz(-1.023205) q[3];
sx q[3];
rz(-0.91528085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(2.0557892) q[2];
rz(-3.0893677) q[3];
sx q[3];
rz(-1.6970535) q[3];
sx q[3];
rz(2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6770342) q[0];
sx q[0];
rz(-2.1573986) q[0];
sx q[0];
rz(-1.1664671) q[0];
rz(-0.72558609) q[1];
sx q[1];
rz(-1.8478994) q[1];
sx q[1];
rz(-2.3988147) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5398139) q[0];
sx q[0];
rz(-0.83866461) q[0];
sx q[0];
rz(1.3108805) q[0];
rz(-pi) q[1];
rz(-0.62959813) q[2];
sx q[2];
rz(-2.3592735) q[2];
sx q[2];
rz(-2.9462189) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2837977) q[1];
sx q[1];
rz(-2.2504914) q[1];
sx q[1];
rz(-2.6074175) q[1];
rz(-pi) q[2];
rz(2.5524213) q[3];
sx q[3];
rz(-2.1450451) q[3];
sx q[3];
rz(1.1816927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.6122931) q[2];
sx q[2];
rz(-0.26088866) q[2];
rz(1.1076814) q[3];
sx q[3];
rz(-1.2872144) q[3];
sx q[3];
rz(2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7850007) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(2.2055431) q[0];
rz(-0.014135663) q[1];
sx q[1];
rz(-1.8363258) q[1];
sx q[1];
rz(2.4900808) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12670853) q[0];
sx q[0];
rz(-1.5471022) q[0];
sx q[0];
rz(0.004304927) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4849309) q[2];
sx q[2];
rz(-1.0495532) q[2];
sx q[2];
rz(1.148664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7266255) q[1];
sx q[1];
rz(-1.8888998) q[1];
sx q[1];
rz(-0.12673641) q[1];
rz(-1.2048079) q[3];
sx q[3];
rz(-1.9848739) q[3];
sx q[3];
rz(-0.28596349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2953879) q[2];
sx q[2];
rz(-0.69869763) q[2];
sx q[2];
rz(2.6182168) q[2];
rz(-0.30803099) q[3];
sx q[3];
rz(-0.57294661) q[3];
sx q[3];
rz(1.3380922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4073407) q[0];
sx q[0];
rz(-0.14695209) q[0];
sx q[0];
rz(1.403632) q[0];
rz(-2.667528) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(1.7636991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5546075) q[0];
sx q[0];
rz(-2.4025612) q[0];
sx q[0];
rz(-2.2686195) q[0];
x q[1];
rz(-0.7147185) q[2];
sx q[2];
rz(-2.208459) q[2];
sx q[2];
rz(-1.6705461) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1174406) q[1];
sx q[1];
rz(-1.9630868) q[1];
sx q[1];
rz(-1.3082318) q[1];
rz(-pi) q[2];
rz(2.58425) q[3];
sx q[3];
rz(-1.7094269) q[3];
sx q[3];
rz(1.0728474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(1.4577929) q[3];
sx q[3];
rz(-1.0554353) q[3];
sx q[3];
rz(-0.84993258) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4703341) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(-2.509027) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(0.40614265) q[2];
sx q[2];
rz(-0.71229013) q[2];
sx q[2];
rz(-2.9980414) q[2];
rz(2.395527) q[3];
sx q[3];
rz(-0.88701556) q[3];
sx q[3];
rz(-1.5942667) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
