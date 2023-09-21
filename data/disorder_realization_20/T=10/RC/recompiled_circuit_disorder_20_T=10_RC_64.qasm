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
rz(-2.0467487) q[0];
sx q[0];
rz(-1.9108512) q[0];
sx q[0];
rz(-2.9283471) q[0];
x q[1];
rz(-1.5413061) q[2];
sx q[2];
rz(-0.87848308) q[2];
sx q[2];
rz(-0.32683795) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80860955) q[1];
sx q[1];
rz(-1.1472417) q[1];
sx q[1];
rz(-2.9575648) q[1];
rz(-0.063135677) q[3];
sx q[3];
rz(-1.6992484) q[3];
sx q[3];
rz(-1.1032784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8193801) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.4780761) q[0];
sx q[0];
rz(2.9852988) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.4888391) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8119368) q[0];
sx q[0];
rz(-0.57528472) q[0];
sx q[0];
rz(-0.29559691) q[0];
x q[1];
rz(-2.4376051) q[2];
sx q[2];
rz(-0.96103243) q[2];
sx q[2];
rz(3.0733382) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.78755806) q[1];
sx q[1];
rz(-0.84317849) q[1];
sx q[1];
rz(0.80227279) q[1];
x q[2];
rz(2.923521) q[3];
sx q[3];
rz(-0.13987386) q[3];
sx q[3];
rz(2.0381896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(0.29671159) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(-1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082224) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(-0.5439533) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.189032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8231359) q[0];
sx q[0];
rz(-1.4134777) q[0];
sx q[0];
rz(-0.29968981) q[0];
x q[1];
rz(-0.15709917) q[2];
sx q[2];
rz(-2.8904449) q[2];
sx q[2];
rz(0.37065333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5570453) q[1];
sx q[1];
rz(-2.097633) q[1];
sx q[1];
rz(-2.3440863) q[1];
x q[2];
rz(1.1071113) q[3];
sx q[3];
rz(-0.92671219) q[3];
sx q[3];
rz(-0.37582276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.8900324) q[2];
rz(-2.6873612) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77804756) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(1.3409412) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-2.2025735) q[1];
sx q[1];
rz(1.1846503) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957164) q[0];
sx q[0];
rz(-1.7850951) q[0];
sx q[0];
rz(-1.3927668) q[0];
x q[1];
rz(2.0834288) q[2];
sx q[2];
rz(-1.1626273) q[2];
sx q[2];
rz(2.4525016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8682755) q[1];
sx q[1];
rz(-1.454121) q[1];
sx q[1];
rz(2.1627746) q[1];
rz(0.42479997) q[3];
sx q[3];
rz(-1.9072755) q[3];
sx q[3];
rz(2.6257023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.23345315) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(-2.502029) q[2];
rz(-0.8574287) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(-0.48979428) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(-0.52039352) q[0];
rz(-0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(0.72881126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6729352) q[0];
sx q[0];
rz(-0.22738439) q[0];
sx q[0];
rz(0.44896455) q[0];
rz(-2.9786733) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(1.3416854) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5908969) q[1];
sx q[1];
rz(-1.7753074) q[1];
sx q[1];
rz(1.900161) q[1];
rz(-pi) q[2];
rz(0.95727386) q[3];
sx q[3];
rz(-2.0271218) q[3];
sx q[3];
rz(-0.43382713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0430498) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(0.70971242) q[2];
rz(-0.63306159) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(0.13171296) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028041) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(-2.2578755) q[0];
rz(2.2209514) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716498) q[0];
sx q[0];
rz(-0.93402223) q[0];
sx q[0];
rz(-2.3417579) q[0];
rz(-1.6502041) q[2];
sx q[2];
rz(-1.6769771) q[2];
sx q[2];
rz(0.90027819) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8710528) q[1];
sx q[1];
rz(-0.45696837) q[1];
sx q[1];
rz(-0.096710042) q[1];
x q[2];
rz(2.382155) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.027823042) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(-1.4893701) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41912115) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.1584287) q[1];
sx q[1];
rz(-1.1694318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1202576) q[0];
sx q[0];
rz(-1.7354256) q[0];
sx q[0];
rz(-3.0740601) q[0];
x q[1];
rz(-1.85497) q[2];
sx q[2];
rz(-1.9467762) q[2];
sx q[2];
rz(0.77501955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6199854) q[1];
sx q[1];
rz(-2.1072525) q[1];
sx q[1];
rz(-1.1779838) q[1];
rz(-pi) q[2];
rz(-0.10294442) q[3];
sx q[3];
rz(-2.7517481) q[3];
sx q[3];
rz(-1.1815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-2.0264758) q[2];
sx q[2];
rz(0.014483359) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(0.55317944) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(2.5792504) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.3969914) q[1];
sx q[1];
rz(-0.45773488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91677529) q[0];
sx q[0];
rz(-2.6579034) q[0];
sx q[0];
rz(1.3902569) q[0];
rz(-pi) q[1];
rz(-2.5896795) q[2];
sx q[2];
rz(-2.1898824) q[2];
sx q[2];
rz(-1.0362831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.25719682) q[1];
sx q[1];
rz(-0.83409062) q[1];
sx q[1];
rz(3.1222621) q[1];
x q[2];
rz(0.3433414) q[3];
sx q[3];
rz(-2.0715054) q[3];
sx q[3];
rz(2.8727369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83595014) q[2];
sx q[2];
rz(-2.2017411) q[2];
sx q[2];
rz(-0.28253728) q[2];
rz(-2.1841168) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(-0.47874513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1829421) q[0];
sx q[0];
rz(-2.4562953) q[0];
sx q[0];
rz(-1.7392993) q[0];
rz(-0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.2618077) q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-0.26913531) q[2];
sx q[2];
rz(-2.4705774) q[2];
sx q[2];
rz(0.035949635) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.251235) q[1];
sx q[1];
rz(-0.8926549) q[1];
sx q[1];
rz(1.8974341) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7997167) q[3];
sx q[3];
rz(-2.5856527) q[3];
sx q[3];
rz(-3.0059909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4385779) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(1.6837439) q[2];
rz(-2.4750989) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(2.3898747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89501971) q[0];
sx q[0];
rz(-2.331215) q[0];
sx q[0];
rz(-1.9816459) q[0];
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
rz(-2.1150072) q[0];
sx q[0];
rz(-2.3837649) q[0];
sx q[0];
rz(-1.8146145) q[0];
rz(0.80600593) q[2];
sx q[2];
rz(-1.5984048) q[2];
sx q[2];
rz(1.031014) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66214661) q[1];
sx q[1];
rz(-1.8188735) q[1];
sx q[1];
rz(-0.76920385) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6108547) q[3];
sx q[3];
rz(-2.8419552) q[3];
sx q[3];
rz(0.50869298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76876172) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(2.5620143) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6282745) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(-1.6127916) q[2];
sx q[2];
rz(-0.73245807) q[2];
sx q[2];
rz(3.0277638) q[2];
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