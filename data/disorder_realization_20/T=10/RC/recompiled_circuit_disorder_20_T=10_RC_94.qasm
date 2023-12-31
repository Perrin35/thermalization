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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0467487) q[0];
sx q[0];
rz(-1.9108512) q[0];
sx q[0];
rz(0.21324555) q[0];
rz(3.1060495) q[2];
sx q[2];
rz(-0.69283743) q[2];
sx q[2];
rz(2.7685744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9077397) q[1];
sx q[1];
rz(-2.6820161) q[1];
sx q[1];
rz(-1.9563667) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.063135677) q[3];
sx q[3];
rz(-1.4423443) q[3];
sx q[3];
rz(1.1032784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8193801) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(1.6248576) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(2.4860399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(0.15629388) q[0];
rz(-2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.6527536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8119368) q[0];
sx q[0];
rz(-2.5663079) q[0];
sx q[0];
rz(-2.8459957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3180914) q[2];
sx q[2];
rz(-2.245792) q[2];
sx q[2];
rz(-1.0456955) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21178791) q[1];
sx q[1];
rz(-1.025052) q[1];
sx q[1];
rz(-0.8916698) q[1];
rz(-pi) q[2];
rz(-3.00499) q[3];
sx q[3];
rz(-1.6009637) q[3];
sx q[3];
rz(0.68340106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2911825) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(-2.8448811) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333703) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(2.6112774) q[0];
rz(-2.5976394) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(1.189032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408802) q[0];
sx q[0];
rz(-1.2749201) q[0];
sx q[0];
rz(1.7353252) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9844935) q[2];
sx q[2];
rz(-2.8904449) q[2];
sx q[2];
rz(2.7709393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4424858) q[1];
sx q[1];
rz(-0.92256303) q[1];
sx q[1];
rz(2.4590765) q[1];
rz(-1.1071113) q[3];
sx q[3];
rz(-2.2148805) q[3];
sx q[3];
rz(2.7657699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(1.2515602) q[2];
rz(-2.6873612) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(-0.35513487) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3635451) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(-1.8006515) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(-1.1846503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23683324) q[0];
sx q[0];
rz(-1.7447115) q[0];
sx q[0];
rz(0.21763344) q[0];
rz(-2.6809533) q[2];
sx q[2];
rz(-2.0377635) q[2];
sx q[2];
rz(-2.4796782) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6726917) q[1];
sx q[1];
rz(-0.60201445) q[1];
sx q[1];
rz(-1.3637581) q[1];
x q[2];
rz(-1.9373478) q[3];
sx q[3];
rz(-1.1712211) q[3];
sx q[3];
rz(-1.2031581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9081395) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(2.502029) q[2];
rz(0.8574287) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5109167) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(-0.72881126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8048808) q[0];
sx q[0];
rz(-1.4727955) q[0];
sx q[0];
rz(2.9360807) q[0];
rz(-pi) q[1];
rz(-0.16291933) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(1.7999072) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55620533) q[1];
sx q[1];
rz(-0.3857179) q[1];
sx q[1];
rz(-2.1410039) q[1];
rz(2.2767931) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(-1.696123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0985428) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(0.70971242) q[2];
rz(-2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0387886) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(2.2578755) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(0.19764915) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55006856) q[0];
sx q[0];
rz(-0.95614377) q[0];
sx q[0];
rz(-2.385925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63981668) q[2];
sx q[2];
rz(-0.13249995) q[2];
sx q[2];
rz(2.8853531) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16285322) q[1];
sx q[1];
rz(-2.0254685) q[1];
sx q[1];
rz(1.6182369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88455172) q[3];
sx q[3];
rz(-0.93730799) q[3];
sx q[3];
rz(-1.3214878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(0.075604288) q[2];
rz(1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224715) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(-1.9627337) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(-1.9721608) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6299316) q[0];
sx q[0];
rz(-0.17782623) q[0];
sx q[0];
rz(-1.9566262) q[0];
x q[1];
rz(2.5240426) q[2];
sx q[2];
rz(-0.46717656) q[2];
sx q[2];
rz(-3.0385366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.84050256) q[1];
sx q[1];
rz(-1.2355348) q[1];
sx q[1];
rz(2.5696978) q[1];
rz(-pi) q[2];
rz(-1.5285989) q[3];
sx q[3];
rz(-1.9584667) q[3];
sx q[3];
rz(-1.2927511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-0.014483359) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(-2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977215) q[0];
sx q[0];
rz(-0.35035366) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(-1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6477752) q[0];
sx q[0];
rz(-1.4871948) q[0];
sx q[0];
rz(1.093822) q[0];
rz(0.93641113) q[2];
sx q[2];
rz(-2.3371156) q[2];
sx q[2];
rz(-0.22125439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22842912) q[1];
sx q[1];
rz(-2.404681) q[1];
sx q[1];
rz(-1.5494898) q[1];
x q[2];
rz(2.1222955) q[3];
sx q[3];
rz(-2.5428452) q[3];
sx q[3];
rz(2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(-0.95747581) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(0.47874513) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95865059) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(0.69933403) q[1];
sx q[1];
rz(-1.8383086) q[1];
sx q[1];
rz(1.8797849) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5489053) q[0];
sx q[0];
rz(-1.840001) q[0];
sx q[0];
rz(0.44434987) q[0];
rz(-2.4883534) q[2];
sx q[2];
rz(-1.4047033) q[2];
sx q[2];
rz(-1.7476029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.251235) q[1];
sx q[1];
rz(-2.2489378) q[1];
sx q[1];
rz(1.8974341) q[1];
x q[2];
rz(1.3654361) q[3];
sx q[3];
rz(-2.0911651) q[3];
sx q[3];
rz(-2.6092649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.70301473) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(0.66649377) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(1.1599468) q[0];
rz(-0.054140422) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(-2.0711526) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4186231) q[0];
sx q[0];
rz(-1.7374991) q[0];
sx q[0];
rz(2.3136501) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3355867) q[2];
sx q[2];
rz(-1.5431879) q[2];
sx q[2];
rz(-2.1105786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66026238) q[1];
sx q[1];
rz(-0.80033014) q[1];
sx q[1];
rz(0.34923133) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24786592) q[3];
sx q[3];
rz(-1.7409179) q[3];
sx q[3];
rz(2.6691013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3728309) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(1.6147511) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(2.4894864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51331818) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
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
rz(-0.53229971) q[3];
sx q[3];
rz(-1.1333864) q[3];
sx q[3];
rz(-0.65896853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
