OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(-2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(3.0537002) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1577507) q[0];
sx q[0];
rz(-1.112126) q[0];
sx q[0];
rz(1.0755324) q[0];
rz(-0.42638875) q[2];
sx q[2];
rz(-2.0585367) q[2];
sx q[2];
rz(-0.79977712) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76074782) q[1];
sx q[1];
rz(-1.5729841) q[1];
sx q[1];
rz(-0.50427498) q[1];
x q[2];
rz(1.8294931) q[3];
sx q[3];
rz(-1.5764578) q[3];
sx q[3];
rz(-2.2744846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0007881) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(-0.56935707) q[2];
rz(-1.6128444) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(2.5879481) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.9083317) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90399088) q[0];
sx q[0];
rz(-1.0407789) q[0];
sx q[0];
rz(-0.7941829) q[0];
rz(2.4670062) q[2];
sx q[2];
rz(-1.6709176) q[2];
sx q[2];
rz(1.1652201) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0765637) q[1];
sx q[1];
rz(-2.934888) q[1];
sx q[1];
rz(1.6480584) q[1];
rz(-0.29970701) q[3];
sx q[3];
rz(-2.4518659) q[3];
sx q[3];
rz(-1.0563861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0040940293) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(-0.0022350524) q[2];
rz(-0.8301174) q[3];
sx q[3];
rz(-0.79289645) q[3];
sx q[3];
rz(-1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8925791) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(-2.2900443) q[0];
rz(-0.7710723) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(-3.070389) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.203513) q[0];
sx q[0];
rz(-2.8305614) q[0];
sx q[0];
rz(-0.87840338) q[0];
rz(-1.115828) q[2];
sx q[2];
rz(-1.3340545) q[2];
sx q[2];
rz(0.67982212) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.066312) q[1];
sx q[1];
rz(-1.4367141) q[1];
sx q[1];
rz(2.4877497) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86665385) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(1.31124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8292134) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-2.6181347) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4098542) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(-1.4473787) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3576413) q[0];
sx q[0];
rz(-2.8581627) q[0];
sx q[0];
rz(0.95143239) q[0];
rz(1.5201735) q[2];
sx q[2];
rz(-0.17855893) q[2];
sx q[2];
rz(2.6176917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4059658) q[1];
sx q[1];
rz(-1.4117129) q[1];
sx q[1];
rz(-0.87079485) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36376603) q[3];
sx q[3];
rz(-1.3437265) q[3];
sx q[3];
rz(2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(2.9966667) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1458364) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(1.2129983) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4764403) q[0];
sx q[0];
rz(-1.9422918) q[0];
sx q[0];
rz(1.2907589) q[0];
x q[1];
rz(-0.80981363) q[2];
sx q[2];
rz(-2.3286616) q[2];
sx q[2];
rz(-0.44250689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1274384) q[1];
sx q[1];
rz(-1.452983) q[1];
sx q[1];
rz(-2.884306) q[1];
rz(-pi) q[2];
rz(1.8323684) q[3];
sx q[3];
rz(-1.3993169) q[3];
sx q[3];
rz(-1.24303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(-2.5115013) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(0.28636006) q[0];
rz(2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(0.58475959) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8557381) q[0];
sx q[0];
rz(-1.2585088) q[0];
sx q[0];
rz(2.5166442) q[0];
x q[1];
rz(0.61442394) q[2];
sx q[2];
rz(-1.5922058) q[2];
sx q[2];
rz(1.1504088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0271921) q[1];
sx q[1];
rz(-2.235734) q[1];
sx q[1];
rz(0.36892885) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29295178) q[3];
sx q[3];
rz(-0.17360273) q[3];
sx q[3];
rz(-1.1374744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16453234) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(-0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(0.54876304) q[0];
rz(-0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(-2.5351977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081189922) q[0];
sx q[0];
rz(-2.5572544) q[0];
sx q[0];
rz(-1.1330182) q[0];
x q[1];
rz(-0.90784448) q[2];
sx q[2];
rz(-1.2693451) q[2];
sx q[2];
rz(-0.19336685) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8555774) q[1];
sx q[1];
rz(-2.3708323) q[1];
sx q[1];
rz(-1.4090178) q[1];
rz(-pi) q[2];
rz(-1.1612732) q[3];
sx q[3];
rz(-1.2425353) q[3];
sx q[3];
rz(-1.3271774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(0.50841224) q[2];
rz(1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(-1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(2.3876277) q[0];
rz(-0.94999653) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(-0.22769134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75922155) q[0];
sx q[0];
rz(-1.641944) q[0];
sx q[0];
rz(1.1573777) q[0];
x q[1];
rz(0.44616206) q[2];
sx q[2];
rz(-1.9036306) q[2];
sx q[2];
rz(1.7762426) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88340064) q[1];
sx q[1];
rz(-2.0996143) q[1];
sx q[1];
rz(1.8316815) q[1];
rz(3.0759096) q[3];
sx q[3];
rz(-0.55995299) q[3];
sx q[3];
rz(-0.63510676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(-0.83855808) q[2];
rz(-0.36455425) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(2.3013039) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39695981) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(1.0546168) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.9931591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6156857) q[0];
sx q[0];
rz(-1.4988168) q[0];
sx q[0];
rz(-3.1176438) q[0];
rz(0.9903332) q[2];
sx q[2];
rz(-0.60539421) q[2];
sx q[2];
rz(-0.47547542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25363898) q[1];
sx q[1];
rz(-1.9983074) q[1];
sx q[1];
rz(2.0931787) q[1];
rz(1.5356482) q[3];
sx q[3];
rz(-1.9225549) q[3];
sx q[3];
rz(-2.7621321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-2.3633374) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(-2.4018438) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(2.443312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90726501) q[0];
sx q[0];
rz(-2.8303879) q[0];
sx q[0];
rz(-0.23647232) q[0];
rz(-1.6449498) q[2];
sx q[2];
rz(-1.5460795) q[2];
sx q[2];
rz(-0.93862426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2472898) q[1];
sx q[1];
rz(-1.3548046) q[1];
sx q[1];
rz(-2.0386019) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9000557) q[3];
sx q[3];
rz(-1.0165443) q[3];
sx q[3];
rz(0.17061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0239821) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(-2.3013766) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(-1.9742112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8828076) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(-0.029126833) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(-1.2226979) q[2];
sx q[2];
rz(-2.1445027) q[2];
sx q[2];
rz(2.4377844) q[2];
rz(-2.9914231) q[3];
sx q[3];
rz(-2.2554845) q[3];
sx q[3];
rz(3.1156202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
