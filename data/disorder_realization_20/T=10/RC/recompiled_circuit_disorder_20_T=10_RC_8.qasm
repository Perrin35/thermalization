OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.5469172) q[0];
sx q[0];
rz(4.1424799) q[0];
sx q[0];
rz(9.6371798) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(1.8600872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40385383) q[0];
sx q[0];
rz(-1.7716584) q[0];
sx q[0];
rz(1.2234729) q[0];
x q[1];
rz(-1.5413061) q[2];
sx q[2];
rz(-0.87848308) q[2];
sx q[2];
rz(2.8147547) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.80860955) q[1];
sx q[1];
rz(-1.1472417) q[1];
sx q[1];
rz(0.18402789) q[1];
rz(-pi) q[2];
rz(-1.4420907) q[3];
sx q[3];
rz(-1.5081815) q[3];
sx q[3];
rz(-0.45941976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(2.4860399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(-0.15629388) q[0];
rz(-2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.6527536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3296559) q[0];
sx q[0];
rz(-0.57528472) q[0];
sx q[0];
rz(0.29559691) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82897562) q[2];
sx q[2];
rz(-2.1301221) q[2];
sx q[2];
rz(2.0916569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78755806) q[1];
sx q[1];
rz(-2.2984142) q[1];
sx q[1];
rz(2.3393199) q[1];
x q[2];
rz(0.1366027) q[3];
sx q[3];
rz(-1.540629) q[3];
sx q[3];
rz(-0.68340106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(-2.8448811) q[2];
rz(-2.8091649) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082224) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.189032) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2169164) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-2.6485373) q[0];
rz(2.9844935) q[2];
sx q[2];
rz(-0.25114775) q[2];
sx q[2];
rz(-0.37065333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58454734) q[1];
sx q[1];
rz(-2.097633) q[1];
sx q[1];
rz(0.79750632) q[1];
rz(-pi) q[2];
rz(1.1071113) q[3];
sx q[3];
rz(-2.2148805) q[3];
sx q[3];
rz(0.37582276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.2515602) q[2];
rz(2.6873612) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(-0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(2.5355133) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(-1.9569424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458762) q[0];
sx q[0];
rz(-1.3564975) q[0];
sx q[0];
rz(-1.3927668) q[0];
rz(-pi) q[1];
rz(-0.84817727) q[2];
sx q[2];
rz(-0.64372534) q[2];
sx q[2];
rz(-1.6456457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9222316) q[1];
sx q[1];
rz(-2.1582099) q[1];
sx q[1];
rz(-0.14031336) q[1];
rz(-pi) q[2];
rz(0.70373669) q[3];
sx q[3];
rz(-0.53547137) q[3];
sx q[3];
rz(-0.4243917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(2.502029) q[2];
rz(-2.284164) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(0.48979428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(0.52039352) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(0.72881126) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6729352) q[0];
sx q[0];
rz(-0.22738439) q[0];
sx q[0];
rz(2.6926281) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9786733) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(-1.7999072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5908969) q[1];
sx q[1];
rz(-1.7753074) q[1];
sx q[1];
rz(-1.900161) q[1];
x q[2];
rz(2.6008984) q[3];
sx q[3];
rz(-2.1139507) q[3];
sx q[3];
rz(2.3054996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0985428) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(-2.4318802) q[2];
rz(0.63306159) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(-3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0387886) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(0.88371712) q[0];
rz(0.9206413) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(2.9439435) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5699428) q[0];
sx q[0];
rz(-0.93402223) q[0];
sx q[0];
rz(-2.3417579) q[0];
x q[1];
rz(-1.6502041) q[2];
sx q[2];
rz(-1.6769771) q[2];
sx q[2];
rz(0.90027819) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9787394) q[1];
sx q[1];
rz(-2.0254685) q[1];
sx q[1];
rz(-1.5233558) q[1];
rz(-pi) q[2];
rz(-2.2570409) q[3];
sx q[3];
rz(-0.93730799) q[3];
sx q[3];
rz(-1.8201049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(0.075604288) q[2];
rz(1.4893701) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(-0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.41912115) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.1694318) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6299316) q[0];
sx q[0];
rz(-0.17782623) q[0];
sx q[0];
rz(-1.1849665) q[0];
x q[1];
rz(-0.39016907) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(-0.90261501) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.84050256) q[1];
sx q[1];
rz(-1.9060578) q[1];
sx q[1];
rz(-0.5718949) q[1];
rz(2.7536105) q[3];
sx q[3];
rz(-1.531732) q[3];
sx q[3];
rz(2.8475873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.823267) q[2];
sx q[2];
rz(-2.0264758) q[2];
sx q[2];
rz(-0.014483359) q[2];
rz(1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(-1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4938175) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(-2.0477707) q[0];
x q[1];
rz(2.2051815) q[2];
sx q[2];
rz(-0.80447703) q[2];
sx q[2];
rz(-0.22125439) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9131635) q[1];
sx q[1];
rz(-0.73691165) q[1];
sx q[1];
rz(-1.5494898) q[1];
rz(0.3433414) q[3];
sx q[3];
rz(-2.0715054) q[3];
sx q[3];
rz(-0.2688558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(-2.6628475) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95865059) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(1.4022934) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.2618077) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5926873) q[0];
sx q[0];
rz(-1.840001) q[0];
sx q[0];
rz(-2.6972428) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65323921) q[2];
sx q[2];
rz(-1.4047033) q[2];
sx q[2];
rz(1.3939898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7463747) q[1];
sx q[1];
rz(-2.4002541) q[1];
sx q[1];
rz(-2.7625601) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6120139) q[3];
sx q[3];
rz(-1.3929318) q[3];
sx q[3];
rz(-1.9999268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(-0.66649377) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(2.3898747) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(1.9816459) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(-2.0711526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72296952) q[0];
sx q[0];
rz(-1.7374991) q[0];
sx q[0];
rz(-2.3136501) q[0];
x q[1];
rz(-0.80600593) q[2];
sx q[2];
rz(-1.5984048) q[2];
sx q[2];
rz(2.1105786) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.479446) q[1];
sx q[1];
rz(-1.8188735) q[1];
sx q[1];
rz(-2.3723888) q[1];
x q[2];
rz(-1.7461734) q[3];
sx q[3];
rz(-1.8150107) q[3];
sx q[3];
rz(-2.0004686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3728309) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(-0.57957831) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51331818) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(1.0340446) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(-1.5288011) q[2];
sx q[2];
rz(-2.4091346) q[2];
sx q[2];
rz(-0.11382881) q[2];
rz(0.74450775) q[3];
sx q[3];
rz(-2.4662938) q[3];
sx q[3];
rz(0.288356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];