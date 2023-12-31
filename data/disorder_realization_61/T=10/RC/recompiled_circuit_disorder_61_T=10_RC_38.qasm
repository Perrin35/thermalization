OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(-2.0079186) q[0];
sx q[0];
rz(-1.5925621) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(-0.54418286) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2092752) q[0];
sx q[0];
rz(-1.5341175) q[0];
sx q[0];
rz(0.067461405) q[0];
rz(3.0402096) q[2];
sx q[2];
rz(-1.6147436) q[2];
sx q[2];
rz(-1.6909388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6737353) q[1];
sx q[1];
rz(-0.89092548) q[1];
sx q[1];
rz(2.6152337) q[1];
rz(-pi) q[2];
rz(2.5892026) q[3];
sx q[3];
rz(-0.033621764) q[3];
sx q[3];
rz(0.53510964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0698174) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(1.3624181) q[2];
rz(0.028256265) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(-0.79022592) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64489275) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(1.9702966) q[0];
rz(-0.21121875) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(-1.7374977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728714) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(2.1268334) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8318769) q[2];
sx q[2];
rz(-2.6385912) q[2];
sx q[2];
rz(-0.75370698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6200871) q[1];
sx q[1];
rz(-1.5477991) q[1];
sx q[1];
rz(0.28316811) q[1];
x q[2];
rz(2.6614463) q[3];
sx q[3];
rz(-1.4501791) q[3];
sx q[3];
rz(-3.0062825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30963787) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(2.1976166) q[2];
rz(0.55654636) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(0.81958333) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(-2.5779285) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99479988) q[0];
sx q[0];
rz(-0.99780647) q[0];
sx q[0];
rz(-1.6550199) q[0];
rz(-pi) q[1];
rz(-0.32391874) q[2];
sx q[2];
rz(-2.516526) q[2];
sx q[2];
rz(-2.6697086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1395531) q[1];
sx q[1];
rz(-0.66322749) q[1];
sx q[1];
rz(-1.0242277) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1590957) q[3];
sx q[3];
rz(-2.1512096) q[3];
sx q[3];
rz(1.8713649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.50239572) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(2.0813265) q[2];
rz(3.0660196) q[3];
sx q[3];
rz(-2.2836756) q[3];
sx q[3];
rz(0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183427) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(-1.5441783) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(-2.4838122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13947978) q[0];
sx q[0];
rz(-1.6008899) q[0];
sx q[0];
rz(-0.17480236) q[0];
rz(-pi) q[1];
rz(2.8305956) q[2];
sx q[2];
rz(-1.1592835) q[2];
sx q[2];
rz(-2.3222773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4752136) q[1];
sx q[1];
rz(-0.5961601) q[1];
sx q[1];
rz(1.1071043) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8834881) q[3];
sx q[3];
rz(-1.8542284) q[3];
sx q[3];
rz(0.62844814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0458935) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(2.0783157) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(-2.1176594) q[0];
rz(0.78760415) q[1];
sx q[1];
rz(-1.5757631) q[1];
sx q[1];
rz(-0.62686282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89676566) q[0];
sx q[0];
rz(-1.5875495) q[0];
sx q[0];
rz(-1.5411975) q[0];
rz(2.3847694) q[2];
sx q[2];
rz(-0.42381091) q[2];
sx q[2];
rz(-0.28048453) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8265241) q[1];
sx q[1];
rz(-1.8693923) q[1];
sx q[1];
rz(1.8930757) q[1];
rz(3.0179126) q[3];
sx q[3];
rz(-1.1808174) q[3];
sx q[3];
rz(-2.9434413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91288599) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(2.1234925) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.3221778) q[3];
sx q[3];
rz(-0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6927032) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(-1.1992136) q[0];
rz(-1.1692283) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(-2.8170524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8400612) q[0];
sx q[0];
rz(-1.9948298) q[0];
sx q[0];
rz(1.016894) q[0];
rz(-pi) q[1];
rz(-0.89914497) q[2];
sx q[2];
rz(-2.2567344) q[2];
sx q[2];
rz(-1.9859973) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2734087) q[1];
sx q[1];
rz(-1.5955828) q[1];
sx q[1];
rz(2.1048057) q[1];
x q[2];
rz(1.2732029) q[3];
sx q[3];
rz(-1.1614359) q[3];
sx q[3];
rz(0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67363182) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(-1.8661873) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6773029) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(1.4087079) q[0];
rz(2.6858221) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(-1.9394402) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31583187) q[0];
sx q[0];
rz(-1.9681265) q[0];
sx q[0];
rz(-2.2560918) q[0];
x q[1];
rz(-1.1832778) q[2];
sx q[2];
rz(-2.2829208) q[2];
sx q[2];
rz(2.9930263) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4533838) q[1];
sx q[1];
rz(-1.6628478) q[1];
sx q[1];
rz(-1.6545014) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3863181) q[3];
sx q[3];
rz(-1.8309621) q[3];
sx q[3];
rz(-0.20862143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(-2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(1.5406066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071035944) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-1.113168) q[0];
rz(2.44599) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(1.6092469) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5743916) q[0];
sx q[0];
rz(-2.6041457) q[0];
sx q[0];
rz(2.6525932) q[0];
rz(-pi) q[1];
rz(-2.3959827) q[2];
sx q[2];
rz(-2.5157305) q[2];
sx q[2];
rz(-1.9288043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0678565) q[1];
sx q[1];
rz(-1.590775) q[1];
sx q[1];
rz(1.1205792) q[1];
rz(-pi) q[2];
rz(0.60204695) q[3];
sx q[3];
rz(-1.4730519) q[3];
sx q[3];
rz(-2.5412113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(0.85419401) q[2];
rz(-1.9130075) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(1.4860229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(2.4771931) q[0];
rz(-1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(-1.857035) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4391543) q[0];
sx q[0];
rz(-2.2674773) q[0];
sx q[0];
rz(2.7657763) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8753488) q[2];
sx q[2];
rz(-1.9387445) q[2];
sx q[2];
rz(-0.14012303) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7094972) q[1];
sx q[1];
rz(-2.2120683) q[1];
sx q[1];
rz(2.1367367) q[1];
rz(-2.7515718) q[3];
sx q[3];
rz(-1.5454925) q[3];
sx q[3];
rz(-0.43061531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61218843) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(-0.55076304) q[2];
rz(-2.4380056) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187014) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(-3.0974467) q[0];
rz(1.6126397) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(-2.3619161) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89833242) q[0];
sx q[0];
rz(-1.7120449) q[0];
sx q[0];
rz(2.6000573) q[0];
rz(-pi) q[1];
rz(-2.0585971) q[2];
sx q[2];
rz(-0.81760815) q[2];
sx q[2];
rz(-1.4873193) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70982198) q[1];
sx q[1];
rz(-2.5858472) q[1];
sx q[1];
rz(-0.35094786) q[1];
rz(-pi) q[2];
rz(2.3603504) q[3];
sx q[3];
rz(-1.2282073) q[3];
sx q[3];
rz(0.70171802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49446517) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(-0.70458448) q[3];
sx q[3];
rz(-1.5141809) q[3];
sx q[3];
rz(-0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4298532) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(-1.343887) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(2.1716739) q[2];
sx q[2];
rz(-2.6261332) q[2];
sx q[2];
rz(0.020608227) q[2];
rz(-0.6380973) q[3];
sx q[3];
rz(-1.6774072) q[3];
sx q[3];
rz(0.93117136) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
