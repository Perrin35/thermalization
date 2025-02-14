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
rz(2.4636318) q[0];
sx q[0];
rz(-2.8673661) q[0];
sx q[0];
rz(-1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(-1.6419799) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8052657) q[0];
sx q[0];
rz(-1.8417497) q[0];
sx q[0];
rz(-1.8150369) q[0];
rz(-1.2798645) q[2];
sx q[2];
rz(-2.3290344) q[2];
sx q[2];
rz(-1.9951374) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1823481) q[1];
sx q[1];
rz(-0.79234353) q[1];
sx q[1];
rz(-2.6807129) q[1];
rz(-pi) q[2];
rz(2.6628482) q[3];
sx q[3];
rz(-0.9968206) q[3];
sx q[3];
rz(-0.51306242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1715673) q[2];
sx q[2];
rz(-1.0415404) q[2];
sx q[2];
rz(-0.901326) q[2];
rz(2.6775635) q[3];
sx q[3];
rz(-1.2814859) q[3];
sx q[3];
rz(0.06981167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0483911) q[0];
sx q[0];
rz(-0.036660107) q[0];
sx q[0];
rz(-1.8638336) q[0];
rz(-0.094206421) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(-1.6069848) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8044455) q[0];
sx q[0];
rz(-2.5408816) q[0];
sx q[0];
rz(-1.5128193) q[0];
rz(1.8707451) q[2];
sx q[2];
rz(-0.19057628) q[2];
sx q[2];
rz(1.9924763) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1086191) q[1];
sx q[1];
rz(-2.0619644) q[1];
sx q[1];
rz(1.8673378) q[1];
rz(0.8022763) q[3];
sx q[3];
rz(-0.55219383) q[3];
sx q[3];
rz(-1.7084165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.719912) q[2];
sx q[2];
rz(-1.2016502) q[2];
sx q[2];
rz(-2.9812532) q[2];
rz(-2.4028589) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(-1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62683231) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(-0.50977388) q[0];
rz(-1.431142) q[1];
sx q[1];
rz(-1.9606083) q[1];
sx q[1];
rz(-2.3950155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9812766) q[0];
sx q[0];
rz(-1.0183987) q[0];
sx q[0];
rz(3.1180326) q[0];
rz(0.72991972) q[2];
sx q[2];
rz(-1.0906719) q[2];
sx q[2];
rz(1.4381222) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.67184) q[1];
sx q[1];
rz(-2.3418509) q[1];
sx q[1];
rz(1.0554764) q[1];
x q[2];
rz(0.88991092) q[3];
sx q[3];
rz(-1.2632252) q[3];
sx q[3];
rz(-2.6079082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15472445) q[2];
sx q[2];
rz(-0.26686033) q[2];
sx q[2];
rz(-1.9179087) q[2];
rz(1.0736505) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(0.8980208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2735485) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(2.7622188) q[0];
rz(-0.1768449) q[1];
sx q[1];
rz(-1.6601446) q[1];
sx q[1];
rz(-2.3462229) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0220227) q[0];
sx q[0];
rz(-1.2031462) q[0];
sx q[0];
rz(-3.0421542) q[0];
rz(1.2878809) q[2];
sx q[2];
rz(-0.77436111) q[2];
sx q[2];
rz(-2.0191569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2242608) q[1];
sx q[1];
rz(-2.5493109) q[1];
sx q[1];
rz(-2.0439953) q[1];
rz(-pi) q[2];
rz(3.0159608) q[3];
sx q[3];
rz(-1.6274618) q[3];
sx q[3];
rz(1.980933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7103601) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(1.7809407) q[2];
rz(-1.0294754) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(1.3220538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.65115702) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(-1.5486451) q[0];
rz(0.40052888) q[1];
sx q[1];
rz(-1.6533886) q[1];
sx q[1];
rz(-1.7318447) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.627305) q[0];
sx q[0];
rz(-0.76111932) q[0];
sx q[0];
rz(2.850015) q[0];
rz(0.39938853) q[2];
sx q[2];
rz(-0.3872954) q[2];
sx q[2];
rz(3.0687817) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6741989) q[1];
sx q[1];
rz(-2.0122897) q[1];
sx q[1];
rz(-0.91798325) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1928113) q[3];
sx q[3];
rz(-1.6262486) q[3];
sx q[3];
rz(0.90622073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8283525) q[2];
sx q[2];
rz(-1.4408377) q[2];
sx q[2];
rz(0.43133119) q[2];
rz(-3.0865772) q[3];
sx q[3];
rz(-2.8300245) q[3];
sx q[3];
rz(2.5855248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677652) q[0];
sx q[0];
rz(-0.25280935) q[0];
sx q[0];
rz(-2.1164236) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-2.5621474) q[1];
sx q[1];
rz(-0.90389171) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65176455) q[0];
sx q[0];
rz(-1.8790885) q[0];
sx q[0];
rz(-2.2405008) q[0];
rz(0.80861096) q[2];
sx q[2];
rz(-2.3249194) q[2];
sx q[2];
rz(-2.6220235) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0817944) q[1];
sx q[1];
rz(-2.3507893) q[1];
sx q[1];
rz(-0.38500144) q[1];
rz(0.76670209) q[3];
sx q[3];
rz(-0.23582102) q[3];
sx q[3];
rz(2.9028149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30770939) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(2.9252388) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-0.68550617) q[3];
sx q[3];
rz(-1.5007796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-2.0348771) q[0];
sx q[0];
rz(-2.4427781) q[0];
rz(-1.1837333) q[1];
sx q[1];
rz(-1.4212757) q[1];
sx q[1];
rz(2.2898477) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257024) q[0];
sx q[0];
rz(-0.53759241) q[0];
sx q[0];
rz(1.4966775) q[0];
rz(-0.57735195) q[2];
sx q[2];
rz(-1.5492348) q[2];
sx q[2];
rz(2.2039764) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4364383) q[1];
sx q[1];
rz(-0.90696224) q[1];
sx q[1];
rz(-1.1052119) q[1];
rz(-2.5301039) q[3];
sx q[3];
rz(-0.73966714) q[3];
sx q[3];
rz(-2.8304493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8375887) q[2];
sx q[2];
rz(-1.8303266) q[2];
sx q[2];
rz(0.50789976) q[2];
rz(2.1128283) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(2.540551) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4429338) q[0];
sx q[0];
rz(-1.3154987) q[0];
sx q[0];
rz(3.1319295) q[0];
rz(-1.6161605) q[1];
sx q[1];
rz(-1.3776255) q[1];
sx q[1];
rz(-0.38472167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5949109) q[0];
sx q[0];
rz(-1.4230886) q[0];
sx q[0];
rz(-1.8099996) q[0];
rz(-pi) q[1];
rz(-0.20236774) q[2];
sx q[2];
rz(-1.9500537) q[2];
sx q[2];
rz(-0.33563644) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8249121) q[1];
sx q[1];
rz(-2.7928565) q[1];
sx q[1];
rz(2.2168) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58458768) q[3];
sx q[3];
rz(-2.2798924) q[3];
sx q[3];
rz(0.60230909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5821417) q[2];
sx q[2];
rz(-0.95953512) q[2];
sx q[2];
rz(-0.01586308) q[2];
rz(1.9150241) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(-2.5198643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7568307) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(-1.0913947) q[0];
rz(-0.81820828) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(-0.6689201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92907897) q[0];
sx q[0];
rz(-2.2336279) q[0];
sx q[0];
rz(-0.79848358) q[0];
rz(-pi) q[1];
rz(-2.2833334) q[2];
sx q[2];
rz(-1.2405292) q[2];
sx q[2];
rz(-1.0077604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3260169) q[1];
sx q[1];
rz(-1.6417802) q[1];
sx q[1];
rz(2.1287588) q[1];
x q[2];
rz(-0.98812859) q[3];
sx q[3];
rz(-1.291953) q[3];
sx q[3];
rz(-3.0089965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6341256) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(0.96366209) q[2];
rz(-1.4108747) q[3];
sx q[3];
rz(-1.7129292) q[3];
sx q[3];
rz(-1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249579) q[0];
sx q[0];
rz(-2.5732915) q[0];
sx q[0];
rz(-1.4547263) q[0];
rz(1.1085054) q[1];
sx q[1];
rz(-1.6051555) q[1];
sx q[1];
rz(2.813521) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011400819) q[0];
sx q[0];
rz(-1.5131803) q[0];
sx q[0];
rz(-1.4559697) q[0];
rz(1.4544308) q[2];
sx q[2];
rz(-1.7504217) q[2];
sx q[2];
rz(-0.76919523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3076664) q[1];
sx q[1];
rz(-0.95912479) q[1];
sx q[1];
rz(-1.0990547) q[1];
rz(-pi) q[2];
x q[2];
rz(1.397473) q[3];
sx q[3];
rz(-1.4651235) q[3];
sx q[3];
rz(-0.25143501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46128094) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(-0.4168365) q[2];
rz(0.69878116) q[3];
sx q[3];
rz(-0.67658934) q[3];
sx q[3];
rz(-0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9417435) q[0];
sx q[0];
rz(-1.9244292) q[0];
sx q[0];
rz(1.695965) q[0];
rz(2.4327714) q[1];
sx q[1];
rz(-2.2292021) q[1];
sx q[1];
rz(3.0880047) q[1];
rz(-2.9403953) q[2];
sx q[2];
rz(-2.1977949) q[2];
sx q[2];
rz(2.0044873) q[2];
rz(-0.96434595) q[3];
sx q[3];
rz(-0.83899211) q[3];
sx q[3];
rz(-0.28750026) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
