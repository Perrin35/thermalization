OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6450206) q[0];
sx q[0];
rz(-1.8121769) q[0];
sx q[0];
rz(-0.42005959) q[0];
rz(-pi) q[1];
rz(2.1098233) q[2];
sx q[2];
rz(-1.4570149) q[2];
sx q[2];
rz(3.0381418) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9125036) q[1];
sx q[1];
rz(-0.87516811) q[1];
sx q[1];
rz(-1.2234664) q[1];
x q[2];
rz(0.31906268) q[3];
sx q[3];
rz(-2.3205119) q[3];
sx q[3];
rz(1.0533489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0212705) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.583064) q[2];
rz(-0.99672404) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(-2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(-3.0875207) q[0];
rz(-1.9460829) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(-2.6057459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5862522) q[0];
sx q[0];
rz(-0.99389168) q[0];
sx q[0];
rz(2.9931086) q[0];
rz(-pi) q[1];
rz(0.29157721) q[2];
sx q[2];
rz(-2.4107286) q[2];
sx q[2];
rz(0.72327327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6260813) q[1];
sx q[1];
rz(-1.9768081) q[1];
sx q[1];
rz(-1.8477693) q[1];
rz(-1.348098) q[3];
sx q[3];
rz(-2.9950812) q[3];
sx q[3];
rz(-3.1080064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(-2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(0.025807468) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2496101) q[0];
sx q[0];
rz(-1.550807) q[0];
sx q[0];
rz(1.32094) q[0];
rz(-pi) q[1];
rz(-1.8823207) q[2];
sx q[2];
rz(-1.6561964) q[2];
sx q[2];
rz(-0.66750079) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3734696) q[1];
sx q[1];
rz(-2.3489531) q[1];
sx q[1];
rz(2.5440689) q[1];
rz(-1.953655) q[3];
sx q[3];
rz(-1.448505) q[3];
sx q[3];
rz(1.0505291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(-2.2606842) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(2.6054629) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49318424) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(2.1302845) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23668004) q[2];
sx q[2];
rz(-1.5152021) q[2];
sx q[2];
rz(-2.6829164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.22524658) q[1];
sx q[1];
rz(-0.9496453) q[1];
sx q[1];
rz(-2.0398554) q[1];
x q[2];
rz(-0.67646497) q[3];
sx q[3];
rz(-1.3389265) q[3];
sx q[3];
rz(-1.3456618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6716016) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(2.0969351) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(-0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2351284) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-1.0513603) q[0];
rz(-1.4936739) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(-0.043118127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331677) q[0];
sx q[0];
rz(-0.63430099) q[0];
sx q[0];
rz(-1.2483291) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0897323) q[2];
sx q[2];
rz(-1.3216002) q[2];
sx q[2];
rz(3.0254226) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4700714) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(-0.77312153) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78955663) q[3];
sx q[3];
rz(-1.6766747) q[3];
sx q[3];
rz(-1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.9859001) q[0];
rz(2.39134) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(-1.0587143) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9583225) q[0];
sx q[0];
rz(-1.1407307) q[0];
sx q[0];
rz(1.3413315) q[0];
x q[1];
rz(-0.58745678) q[2];
sx q[2];
rz(-0.60411462) q[2];
sx q[2];
rz(0.47746745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7212972) q[1];
sx q[1];
rz(-1.2696206) q[1];
sx q[1];
rz(-2.9073614) q[1];
rz(-1.3771463) q[3];
sx q[3];
rz(-2.2691257) q[3];
sx q[3];
rz(-2.8805672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041615151) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(-1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(2.2999433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3232988) q[0];
sx q[0];
rz(-2.5677486) q[0];
sx q[0];
rz(-0.42563514) q[0];
rz(-pi) q[1];
rz(-0.26288962) q[2];
sx q[2];
rz(-2.4979696) q[2];
sx q[2];
rz(2.7123244) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0591653) q[1];
sx q[1];
rz(-1.4945684) q[1];
sx q[1];
rz(-3.0928897) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82960415) q[3];
sx q[3];
rz(-1.4083574) q[3];
sx q[3];
rz(0.12773578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(3.0916396) q[2];
rz(-2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(-0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(0.41762525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091052) q[0];
sx q[0];
rz(-2.0628477) q[0];
sx q[0];
rz(1.0743272) q[0];
x q[1];
rz(3.0268961) q[2];
sx q[2];
rz(-1.2071929) q[2];
sx q[2];
rz(0.39322688) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8466134) q[1];
sx q[1];
rz(-1.1049005) q[1];
sx q[1];
rz(-0.93211517) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8475259) q[3];
sx q[3];
rz(-1.9293647) q[3];
sx q[3];
rz(2.8420574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(2.2414482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(-2.5667403) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(0.94690698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6553584) q[0];
sx q[0];
rz(-1.8706733) q[0];
sx q[0];
rz(0.34098682) q[0];
rz(-3.0445381) q[2];
sx q[2];
rz(-1.6104873) q[2];
sx q[2];
rz(-1.5511712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8967594) q[1];
sx q[1];
rz(-1.4866801) q[1];
sx q[1];
rz(2.8068078) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19091786) q[3];
sx q[3];
rz(-2.0058161) q[3];
sx q[3];
rz(-2.2866979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-0.26556382) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(0.65504909) q[0];
rz(-0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8048332) q[0];
sx q[0];
rz(-1.5157489) q[0];
sx q[0];
rz(2.9677662) q[0];
rz(-1.8337433) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(-2.3761689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5035142) q[1];
sx q[1];
rz(-0.44736171) q[1];
sx q[1];
rz(-2.4380986) q[1];
rz(-pi) q[2];
x q[2];
rz(2.321407) q[3];
sx q[3];
rz(-0.60884848) q[3];
sx q[3];
rz(1.604515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(2.541686) q[2];
rz(-2.24263) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-0.22944336) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(-0.12702604) q[2];
sx q[2];
rz(-2.1913678) q[2];
sx q[2];
rz(0.36507228) q[2];
rz(-1.0944081) q[3];
sx q[3];
rz(-2.4088358) q[3];
sx q[3];
rz(0.81523304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];