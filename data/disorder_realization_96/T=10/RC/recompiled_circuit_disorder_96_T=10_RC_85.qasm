OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(3.0135305) q[0];
sx q[0];
rz(11.749) q[0];
rz(0.983239) q[1];
sx q[1];
rz(2.6020738) q[1];
sx q[1];
rz(10.625216) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1737328) q[0];
sx q[0];
rz(-1.9779441) q[0];
sx q[0];
rz(-1.3074387) q[0];
x q[1];
rz(1.7898516) q[2];
sx q[2];
rz(-0.54974216) q[2];
sx q[2];
rz(-1.4866536) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22908902) q[1];
sx q[1];
rz(-2.2664245) q[1];
sx q[1];
rz(-1.2234664) q[1];
rz(-pi) q[2];
rz(0.79521631) q[3];
sx q[3];
rz(-1.3391558) q[3];
sx q[3];
rz(-2.8455646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1203221) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(-2.1448686) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(-0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5834171) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(-1.1955098) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(-2.6057459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5862522) q[0];
sx q[0];
rz(-2.147701) q[0];
sx q[0];
rz(-2.9931086) q[0];
x q[1];
rz(0.70948647) q[2];
sx q[2];
rz(-1.7638793) q[2];
sx q[2];
rz(2.513934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0332196) q[1];
sx q[1];
rz(-0.48711005) q[1];
sx q[1];
rz(0.56652041) q[1];
rz(0.032580839) q[3];
sx q[3];
rz(-1.7136646) q[3];
sx q[3];
rz(-0.25861614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(-0.95834857) q[2];
rz(-3.0751394) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41985837) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(-0.45723215) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(0.025807468) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3993527) q[0];
sx q[0];
rz(-0.25063801) q[0];
sx q[0];
rz(1.4901194) q[0];
x q[1];
rz(-0.089695887) q[2];
sx q[2];
rz(-1.2604453) q[2];
sx q[2];
rz(-2.2657564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60359309) q[1];
sx q[1];
rz(-2.2003761) q[1];
sx q[1];
rz(1.0521207) q[1];
x q[2];
rz(-1.1879376) q[3];
sx q[3];
rz(-1.448505) q[3];
sx q[3];
rz(-1.0505291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1228483) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(-0.5853816) q[2];
rz(-2.9600926) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-0.62269354) q[0];
sx q[0];
rz(-0.17661072) q[0];
rz(-2.2606842) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(-2.6054629) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6484084) q[0];
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
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.199898) q[1];
sx q[1];
rz(-0.75921339) q[1];
sx q[1];
rz(-0.56337507) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7809308) q[3];
sx q[3];
rz(-2.4324527) q[3];
sx q[3];
rz(-3.0879471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46999103) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(-0.70703834) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(-2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(0.043118127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9263822) q[0];
sx q[0];
rz(-2.1676817) q[0];
sx q[0];
rz(-2.9125288) q[0];
rz(-pi) q[1];
rz(-1.0967069) q[2];
sx q[2];
rz(-2.5708963) q[2];
sx q[2];
rz(1.2793465) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4823618) q[1];
sx q[1];
rz(-1.8419957) q[1];
sx q[1];
rz(-1.2960474) q[1];
rz(-pi) q[2];
x q[2];
rz(2.99302) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(3.1032004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(-2.5514065) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5181638) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(-2.0828784) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48460618) q[0];
sx q[0];
rz(-1.3625506) q[0];
sx q[0];
rz(-2.7013742) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5541359) q[2];
sx q[2];
rz(-2.537478) q[2];
sx q[2];
rz(-2.6641252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25710479) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(2.2124955) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9162354) q[3];
sx q[3];
rz(-0.72031027) q[3];
sx q[3];
rz(2.5845137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(-1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041615151) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.51145) q[0];
rz(-1.7639683) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(-2.2999433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.323303) q[0];
sx q[0];
rz(-1.0536195) q[0];
sx q[0];
rz(1.3099758) q[0];
x q[1];
rz(2.878703) q[2];
sx q[2];
rz(-0.6436231) q[2];
sx q[2];
rz(0.4292683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51534286) q[1];
sx q[1];
rz(-1.6193577) q[1];
sx q[1];
rz(-1.6471144) q[1];
rz(-pi) q[2];
rz(1.8089201) q[3];
sx q[3];
rz(-2.3861285) q[3];
sx q[3];
rz(1.2681703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(-0.66155457) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(-1.7393973) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(0.41762525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2554889) q[0];
sx q[0];
rz(-0.68414738) q[0];
sx q[0];
rz(0.72649254) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11469658) q[2];
sx q[2];
rz(-1.2071929) q[2];
sx q[2];
rz(0.39322688) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.046125267) q[1];
sx q[1];
rz(-2.1324665) q[1];
sx q[1];
rz(0.559505) q[1];
x q[2];
rz(-1.1975708) q[3];
sx q[3];
rz(-1.2959359) q[3];
sx q[3];
rz(-1.1653792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79545704) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(2.0987089) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-2.2414482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8326571) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-0.62966627) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.31153) q[1];
sx q[1];
rz(2.1946857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7787951) q[0];
sx q[0];
rz(-0.4501833) q[0];
sx q[0];
rz(-0.74624004) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6106748) q[2];
sx q[2];
rz(-1.6677742) q[2];
sx q[2];
rz(-3.1181042) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7864101) q[1];
sx q[1];
rz(-1.9043515) q[1];
sx q[1];
rz(1.4817609) q[1];
x q[2];
rz(1.1287273) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(0.63463075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.8927195) q[2];
rz(-0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0749851) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(2.24522) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33675942) q[0];
sx q[0];
rz(-1.5157489) q[0];
sx q[0];
rz(-2.9677662) q[0];
x q[1];
rz(1.3078493) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(-2.3761689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2584553) q[1];
sx q[1];
rz(-1.2346134) q[1];
sx q[1];
rz(1.8717481) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82018567) q[3];
sx q[3];
rz(-0.60884848) q[3];
sx q[3];
rz(1.5370777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(2.541686) q[2];
rz(2.24263) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.8469289) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(2.1951998) q[2];
sx q[2];
rz(-1.6740435) q[2];
sx q[2];
rz(-1.2798535) q[2];
rz(-2.7502144) q[3];
sx q[3];
rz(-0.93422514) q[3];
sx q[3];
rz(1.422062) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
