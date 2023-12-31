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
rz(2.9291908) q[0];
rz(0.71495932) q[1];
sx q[1];
rz(3.9290805) q[1];
sx q[1];
rz(10.706283) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6706657) q[0];
sx q[0];
rz(-2.7424194) q[0];
sx q[0];
rz(-2.1098718) q[0];
x q[1];
rz(3.1060495) q[2];
sx q[2];
rz(-0.69283743) q[2];
sx q[2];
rz(-0.37301829) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4557588) q[1];
sx q[1];
rz(-1.7384006) q[1];
sx q[1];
rz(-2.0007677) q[1];
rz(-pi) q[2];
rz(0.063135677) q[3];
sx q[3];
rz(-1.6992484) q[3];
sx q[3];
rz(1.1032784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32221258) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(-1.5167351) q[2];
rz(2.8939698) q[3];
sx q[3];
rz(-1.473749) q[3];
sx q[3];
rz(0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(-2.9852988) q[0];
rz(-2.5022751) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.6527536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1600906) q[0];
sx q[0];
rz(-2.1182051) q[0];
sx q[0];
rz(1.3840958) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3180914) q[2];
sx q[2];
rz(-0.89580065) q[2];
sx q[2];
rz(-2.0958971) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21178791) q[1];
sx q[1];
rz(-1.025052) q[1];
sx q[1];
rz(0.8916698) q[1];
rz(2.923521) q[3];
sx q[3];
rz(-3.0017188) q[3];
sx q[3];
rz(-2.0381896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2911825) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(-0.29671159) q[2];
rz(0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.2338352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333703) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(0.5439533) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(-1.189032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408802) q[0];
sx q[0];
rz(-1.8666726) q[0];
sx q[0];
rz(-1.7353252) q[0];
x q[1];
rz(-1.5306773) q[2];
sx q[2];
rz(-1.3228068) q[2];
sx q[2];
rz(2.9330394) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51057286) q[1];
sx q[1];
rz(-0.90386183) q[1];
sx q[1];
rz(-2.2651947) q[1];
rz(-pi) q[2];
rz(1.1071113) q[3];
sx q[3];
rz(-0.92671219) q[3];
sx q[3];
rz(2.7657699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.236078) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.8900324) q[2];
rz(-2.6873612) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(1.1846503) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8458762) q[0];
sx q[0];
rz(-1.3564975) q[0];
sx q[0];
rz(-1.7488259) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84817727) q[2];
sx q[2];
rz(-2.4978673) q[2];
sx q[2];
rz(-1.4959469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9222316) q[1];
sx q[1];
rz(-0.9833828) q[1];
sx q[1];
rz(0.14031336) q[1];
rz(-pi) q[2];
rz(1.2042449) q[3];
sx q[3];
rz(-1.9703715) q[3];
sx q[3];
rz(1.2031581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(0.63956368) q[2];
rz(2.284164) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-0.48979428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5109167) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(-0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(0.72881126) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279014) q[0];
sx q[0];
rz(-1.7753082) q[0];
sx q[0];
rz(-1.4707028) q[0];
x q[1];
rz(-2.2949218) q[2];
sx q[2];
rz(-1.4482822) q[2];
sx q[2];
rz(0.12144897) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0923857) q[1];
sx q[1];
rz(-1.8930463) q[1];
sx q[1];
rz(2.9258123) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8647996) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(1.4454696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0430498) q[2];
sx q[2];
rz(-2.7816732) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0387886) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(-0.88371712) q[0];
rz(-2.2209514) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(0.19764915) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6186002) q[0];
sx q[0];
rz(-2.1654961) q[0];
sx q[0];
rz(-0.80070514) q[0];
rz(-pi) q[1];
rz(3.0350787) q[2];
sx q[2];
rz(-1.4918367) q[2];
sx q[2];
rz(-0.67895141) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27053988) q[1];
sx q[1];
rz(-0.45696837) q[1];
sx q[1];
rz(-3.0448826) q[1];
x q[2];
rz(-2.382155) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(-2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(3.0659884) q[2];
rz(1.4893701) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.41912115) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(0.042073123) q[0];
rz(1.178859) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(-1.9721608) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5116611) q[0];
sx q[0];
rz(-2.9637664) q[0];
sx q[0];
rz(1.1849665) q[0];
rz(-pi) q[1];
rz(-1.2866227) q[2];
sx q[2];
rz(-1.9467762) q[2];
sx q[2];
rz(-0.77501955) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6199854) q[1];
sx q[1];
rz(-1.0343401) q[1];
sx q[1];
rz(-1.1779838) q[1];
rz(2.7536105) q[3];
sx q[3];
rz(-1.6098607) q[3];
sx q[3];
rz(-2.8475873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3183257) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-3.1271093) q[2];
rz(1.621834) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-2.5792504) q[0];
rz(-1.4315804) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(-2.6838578) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.021488) q[0];
sx q[0];
rz(-1.095626) q[0];
sx q[0];
rz(3.0475463) q[0];
rz(-0.93641113) q[2];
sx q[2];
rz(-0.80447703) q[2];
sx q[2];
rz(2.9203383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25719682) q[1];
sx q[1];
rz(-2.307502) q[1];
sx q[1];
rz(-0.019330545) q[1];
rz(-1.0192972) q[3];
sx q[3];
rz(-0.5987474) q[3];
sx q[3];
rz(-2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(0.28253728) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(0.47874513) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95865059) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(1.4022934) q[0];
rz(-2.4422586) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.8797849) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5926873) q[0];
sx q[0];
rz(-1.840001) q[0];
sx q[0];
rz(-0.44434987) q[0];
rz(-pi) q[1];
rz(2.8724573) q[2];
sx q[2];
rz(-0.67101523) q[2];
sx q[2];
rz(3.105643) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1101493) q[1];
sx q[1];
rz(-1.823339) q[1];
sx q[1];
rz(2.4367711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34187596) q[3];
sx q[3];
rz(-2.5856527) q[3];
sx q[3];
rz(3.0059909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70301473) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(-1.6837439) q[2];
rz(-0.66649377) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(-0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(1.1599468) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(-2.0711526) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72296952) q[0];
sx q[0];
rz(-1.4040935) q[0];
sx q[0];
rz(-2.3136501) q[0];
x q[1];
rz(-1.530933) q[2];
sx q[2];
rz(-2.3764052) q[2];
sx q[2];
rz(0.51102343) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1420319) q[1];
sx q[1];
rz(-0.83082098) q[1];
sx q[1];
rz(1.9097411) q[1];
rz(-pi) q[2];
rz(-1.3954193) q[3];
sx q[3];
rz(-1.8150107) q[3];
sx q[3];
rz(2.0004686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3728309) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(-1.5268415) q[2];
rz(0.57957831) q[3];
sx q[3];
rz(-1.3468346) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51331818) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(2.1075481) q[1];
sx q[1];
rz(-0.68987344) q[1];
sx q[1];
rz(-0.72763163) q[1];
rz(1.6127916) q[2];
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
