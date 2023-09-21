OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5946755) q[0];
sx q[0];
rz(-1.0008873) q[0];
sx q[0];
rz(-0.21240182) q[0];
rz(0.71495932) q[1];
sx q[1];
rz(-2.3541048) q[1];
sx q[1];
rz(1.2815055) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4709269) q[0];
sx q[0];
rz(-2.7424194) q[0];
sx q[0];
rz(1.0317208) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69252695) q[2];
sx q[2];
rz(-1.5480969) q[2];
sx q[2];
rz(-1.2251309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80860955) q[1];
sx q[1];
rz(-1.994351) q[1];
sx q[1];
rz(-0.18402789) q[1];
x q[2];
rz(-2.0251861) q[3];
sx q[3];
rz(-2.9985399) q[3];
sx q[3];
rz(-1.5798626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(2.8939698) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(-2.9852988) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.4888391) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.70398753) q[2];
sx q[2];
rz(-0.96103243) q[2];
sx q[2];
rz(-3.0733382) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9298047) q[1];
sx q[1];
rz(-1.025052) q[1];
sx q[1];
rz(2.2499229) q[1];
rz(-1.6012472) q[3];
sx q[3];
rz(-1.4342562) q[3];
sx q[3];
rz(0.88324916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(-0.29671159) q[2];
rz(-2.8091649) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.6333703) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(-0.53031522) q[0];
rz(2.5976394) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(1.189032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30071242) q[0];
sx q[0];
rz(-1.2749201) q[0];
sx q[0];
rz(-1.4062675) q[0];
rz(1.5306773) q[2];
sx q[2];
rz(-1.8187858) q[2];
sx q[2];
rz(2.9330394) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51057286) q[1];
sx q[1];
rz(-2.2377308) q[1];
sx q[1];
rz(2.2651947) q[1];
rz(-pi) q[2];
rz(-0.69840188) q[3];
sx q[3];
rz(-1.2050556) q[3];
sx q[3];
rz(1.6549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(2.6873612) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.93901912) q[1];
sx q[1];
rz(-1.1846503) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458762) q[0];
sx q[0];
rz(-1.7850951) q[0];
sx q[0];
rz(1.7488259) q[0];
rz(-pi) q[1];
rz(-2.2934154) q[2];
sx q[2];
rz(-0.64372534) q[2];
sx q[2];
rz(-1.4959469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6726917) q[1];
sx q[1];
rz(-2.5395782) q[1];
sx q[1];
rz(1.7778346) q[1];
x q[2];
rz(-2.7167927) q[3];
sx q[3];
rz(-1.9072755) q[3];
sx q[3];
rz(-0.51589033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23345315) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(0.63956368) q[2];
rz(2.284164) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(2.6517984) q[3];
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
rz(-0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(0.72881126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8048808) q[0];
sx q[0];
rz(-1.4727955) q[0];
sx q[0];
rz(-2.9360807) q[0];
rz(0.16291933) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(1.3416854) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.049206991) q[1];
sx q[1];
rz(-1.8930463) q[1];
sx q[1];
rz(2.9258123) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1843188) q[3];
sx q[3];
rz(-1.1144708) q[3];
sx q[3];
rz(-0.43382713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0430498) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(0.70971242) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(0.13171296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0387886) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(2.2578755) q[0];
rz(-0.9206413) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(2.9439435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716498) q[0];
sx q[0];
rz(-2.2075704) q[0];
sx q[0];
rz(-2.3417579) q[0];
x q[1];
rz(-2.501776) q[2];
sx q[2];
rz(-3.0090927) q[2];
sx q[2];
rz(0.25623955) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9787394) q[1];
sx q[1];
rz(-2.0254685) q[1];
sx q[1];
rz(-1.5233558) q[1];
rz(-0.88455172) q[3];
sx q[3];
rz(-2.2042847) q[3];
sx q[3];
rz(1.3214878) q[3];
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
rz(-1.6522225) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(0.82908019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41912115) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(-3.0995195) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(-1.9721608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7021381) q[0];
sx q[0];
rz(-1.5041782) q[0];
sx q[0];
rz(-1.4057977) q[0];
rz(-pi) q[1];
rz(0.39016907) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(0.90261501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52160727) q[1];
sx q[1];
rz(-1.0343401) q[1];
sx q[1];
rz(1.9636088) q[1];
rz(-pi) q[2];
rz(1.6129937) q[3];
sx q[3];
rz(-1.183126) q[3];
sx q[3];
rz(1.2927511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(3.1271093) q[2];
rz(-1.621834) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(-0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0018205) q[0];
sx q[0];
rz(-0.35035366) q[0];
sx q[0];
rz(-2.5792504) q[0];
rz(-1.7100122) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(2.6838578) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2248174) q[0];
sx q[0];
rz(-2.6579034) q[0];
sx q[0];
rz(1.7513357) q[0];
x q[1];
rz(-0.87403239) q[2];
sx q[2];
rz(-2.0119785) q[2];
sx q[2];
rz(-0.87768427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3265877) q[1];
sx q[1];
rz(-1.5564789) q[1];
sx q[1];
rz(-0.83399764) q[1];
rz(2.1222955) q[3];
sx q[3];
rz(-2.5428452) q[3];
sx q[3];
rz(-0.37125722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(-2.8590554) q[2];
rz(2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(2.6628475) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1829421) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.7392993) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.8797849) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5489053) q[0];
sx q[0];
rz(-1.3015916) q[0];
sx q[0];
rz(-2.6972428) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3627522) q[2];
sx q[2];
rz(-2.2135452) q[2];
sx q[2];
rz(0.3026697) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0314434) q[1];
sx q[1];
rz(-1.823339) q[1];
sx q[1];
rz(2.4367711) q[1];
rz(-pi) q[2];
rz(1.7761566) q[3];
sx q[3];
rz(-2.0911651) q[3];
sx q[3];
rz(2.6092649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70301473) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(-0.66649377) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2465729) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(-1.0704401) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4186231) q[0];
sx q[0];
rz(-1.7374991) q[0];
sx q[0];
rz(-0.82794257) q[0];
rz(-pi) q[1];
rz(-1.6106597) q[2];
sx q[2];
rz(-2.3764052) q[2];
sx q[2];
rz(-0.51102343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66214661) q[1];
sx q[1];
rz(-1.8188735) q[1];
sx q[1];
rz(-0.76920385) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7461734) q[3];
sx q[3];
rz(-1.8150107) q[3];
sx q[3];
rz(-2.0004686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76876172) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(-1.5268415) q[2];
rz(-0.57957831) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6282745) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(1.0340446) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(-2.3028159) q[2];
sx q[2];
rz(-1.5427187) q[2];
sx q[2];
rz(-1.6533921) q[2];
rz(-2.0680239) q[3];
sx q[3];
rz(-2.0484925) q[3];
sx q[3];
rz(-1.9852553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];