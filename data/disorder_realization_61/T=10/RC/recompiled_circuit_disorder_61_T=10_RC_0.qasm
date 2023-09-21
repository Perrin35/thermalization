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
rz(2.5974098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.005851) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(-2.6430921) q[0];
x q[1];
rz(2.7317023) q[2];
sx q[2];
rz(-3.0311243) q[2];
sx q[2];
rz(-2.8540749) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6737353) q[1];
sx q[1];
rz(-0.89092548) q[1];
sx q[1];
rz(0.52635898) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5531494) q[3];
sx q[3];
rz(-1.5994161) q[3];
sx q[3];
rz(2.0538405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0698174) q[2];
sx q[2];
rz(-1.8775619) q[2];
sx q[2];
rz(-1.7791746) q[2];
rz(3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(-2.3513667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(-1.171296) q[0];
rz(-2.9303739) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(1.7374977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0728714) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(1.0147592) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48268433) q[2];
sx q[2];
rz(-1.4233372) q[2];
sx q[2];
rz(-2.5978616) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1120226) q[1];
sx q[1];
rz(-0.28407541) q[1];
sx q[1];
rz(3.0594538) q[1];
x q[2];
rz(2.6614463) q[3];
sx q[3];
rz(-1.6914136) q[3];
sx q[3];
rz(3.0062825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.6654623) q[2];
sx q[2];
rz(-2.1976166) q[2];
rz(-0.55654636) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(-1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.29207644) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(-1.7180432) q[0];
rz(-2.3220093) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(2.5779285) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99479988) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(-1.6550199) q[0];
x q[1];
rz(-1.3450422) q[2];
sx q[2];
rz(-0.98276897) q[2];
sx q[2];
rz(3.0622481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7970265) q[1];
sx q[1];
rz(-2.1246506) q[1];
sx q[1];
rz(-0.38573854) q[1];
x q[2];
rz(-1.9824969) q[3];
sx q[3];
rz(-0.99038306) q[3];
sx q[3];
rz(-1.8713649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6391969) q[2];
sx q[2];
rz(-2.0194619) q[2];
sx q[2];
rz(-1.0602661) q[2];
rz(-3.0660196) q[3];
sx q[3];
rz(-2.2836756) q[3];
sx q[3];
rz(3.0269572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.32325) q[0];
sx q[0];
rz(-2.8338354) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(1.5441783) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(-2.4838122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7155899) q[0];
sx q[0];
rz(-1.7455187) q[0];
sx q[0];
rz(-1.5402373) q[0];
rz(-pi) q[1];
rz(-1.1409608) q[2];
sx q[2];
rz(-1.2865215) q[2];
sx q[2];
rz(-0.62361275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2967589) q[1];
sx q[1];
rz(-1.3169603) q[1];
sx q[1];
rz(-2.1162507) q[1];
rz(-pi) q[2];
rz(-2.2903719) q[3];
sx q[3];
rz(-2.7606066) q[3];
sx q[3];
rz(1.3850497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0956991) q[2];
sx q[2];
rz(-2.7313488) q[2];
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
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(-2.1176594) q[0];
rz(-2.3539885) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(0.62686282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67353467) q[0];
sx q[0];
rz(-1.600391) q[0];
sx q[0];
rz(-3.1248321) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8711898) q[2];
sx q[2];
rz(-1.2671748) q[2];
sx q[2];
rz(-2.6189569) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6080731) q[1];
sx q[1];
rz(-2.7058209) q[1];
sx q[1];
rz(-0.79969745) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0179126) q[3];
sx q[3];
rz(-1.1808174) q[3];
sx q[3];
rz(-0.19815138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2287067) q[2];
sx q[2];
rz(-1.2330981) q[2];
sx q[2];
rz(-2.1234925) q[2];
rz(-2.5300238) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-0.95156804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4488895) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(1.1992136) q[0];
rz(-1.9723643) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(-0.32454023) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6231461) q[0];
sx q[0];
rz(-2.0707957) q[0];
sx q[0];
rz(-2.6536233) q[0];
x q[1];
rz(0.89914497) q[2];
sx q[2];
rz(-2.2567344) q[2];
sx q[2];
rz(-1.1555954) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4536344) q[1];
sx q[1];
rz(-1.0369685) q[1];
sx q[1];
rz(-0.028793528) q[1];
rz(-pi) q[2];
rz(1.8683897) q[3];
sx q[3];
rz(-1.1614359) q[3];
sx q[3];
rz(-0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4679608) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(1.0547137) q[3];
sx q[3];
rz(-0.79189363) q[3];
sx q[3];
rz(-1.595343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6773029) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(1.4087079) q[0];
rz(-2.6858221) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(1.9394402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5802655) q[0];
sx q[0];
rz(-2.1939477) q[0];
sx q[0];
rz(2.6448963) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9583148) q[2];
sx q[2];
rz(-0.85867184) q[2];
sx q[2];
rz(-0.14856635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68820885) q[1];
sx q[1];
rz(-1.6628478) q[1];
sx q[1];
rz(-1.4870912) q[1];
rz(-pi) q[2];
rz(1.3863181) q[3];
sx q[3];
rz(-1.8309621) q[3];
sx q[3];
rz(0.20862143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(2.2953575) q[2];
rz(-0.49514654) q[3];
sx q[3];
rz(-2.4954002) q[3];
sx q[3];
rz(-1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705567) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-1.113168) q[0];
rz(-0.69560266) q[1];
sx q[1];
rz(-2.7516987) q[1];
sx q[1];
rz(1.5323458) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217864) q[0];
sx q[0];
rz(-1.1018503) q[0];
sx q[0];
rz(-1.297834) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1148909) q[2];
sx q[2];
rz(-2.0156983) q[2];
sx q[2];
rz(-1.0783431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.53837) q[1];
sx q[1];
rz(-2.6909628) q[1];
sx q[1];
rz(-1.6166812) q[1];
rz(-pi) q[2];
rz(0.60204695) q[3];
sx q[3];
rz(-1.6685408) q[3];
sx q[3];
rz(-0.60038139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1476851) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(2.2873986) q[2];
rz(-1.2285852) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.4860229) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478741) q[0];
sx q[0];
rz(-2.6429208) q[0];
sx q[0];
rz(0.6643995) q[0];
rz(-1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(-1.857035) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4391543) q[0];
sx q[0];
rz(-2.2674773) q[0];
sx q[0];
rz(0.37581635) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9719073) q[2];
sx q[2];
rz(-0.45058695) q[2];
sx q[2];
rz(2.633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7582693) q[1];
sx q[1];
rz(-0.82793923) q[1];
sx q[1];
rz(2.5187056) q[1];
rz(0.066468863) q[3];
sx q[3];
rz(-2.7507938) q[3];
sx q[3];
rz(1.9399411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5294042) q[2];
sx q[2];
rz(-2.2470784) q[2];
sx q[2];
rz(-0.55076304) q[2];
rz(2.4380056) q[3];
sx q[3];
rz(-2.1313322) q[3];
sx q[3];
rz(1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7228912) q[0];
sx q[0];
rz(-0.35690618) q[0];
sx q[0];
rz(0.044145949) q[0];
rz(-1.6126397) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(2.3619161) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89833242) q[0];
sx q[0];
rz(-1.7120449) q[0];
sx q[0];
rz(-2.6000573) q[0];
rz(-pi) q[1];
rz(0.81515628) q[2];
sx q[2];
rz(-1.2218468) q[2];
sx q[2];
rz(-2.7100035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.838678) q[1];
sx q[1];
rz(-1.0524787) q[1];
sx q[1];
rz(-1.360449) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6727647) q[3];
sx q[3];
rz(-0.83823293) q[3];
sx q[3];
rz(-1.1956786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6471275) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(1.5987827) q[2];
rz(-0.70458448) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.71173944) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(-1.343887) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(1.1336397) q[2];
sx q[2];
rz(-1.2883678) q[2];
sx q[2];
rz(-1.012445) q[2];
rz(-2.963831) q[3];
sx q[3];
rz(-0.64571417) q[3];
sx q[3];
rz(-0.4971102) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];