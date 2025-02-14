OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(-1.1986873) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(-0.79678798) q[1];
sx q[1];
rz(0.30153433) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806592) q[0];
sx q[0];
rz(-1.3863239) q[0];
sx q[0];
rz(0.16452275) q[0];
rz(1.8704988) q[2];
sx q[2];
rz(-1.7500008) q[2];
sx q[2];
rz(1.186893) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49820013) q[1];
sx q[1];
rz(-2.3215356) q[1];
sx q[1];
rz(-3.0489075) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1490378) q[3];
sx q[3];
rz(-2.7054477) q[3];
sx q[3];
rz(-0.25806043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.6887168) q[2];
sx q[2];
rz(-2.0417002) q[2];
sx q[2];
rz(-1.1348881) q[2];
rz(0.12198837) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(0.056099135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6720471) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(3.1392414) q[0];
rz(2.7408842) q[1];
sx q[1];
rz(-1.3688764) q[1];
sx q[1];
rz(-2.2854663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6295638) q[0];
sx q[0];
rz(-0.95603795) q[0];
sx q[0];
rz(0.31350664) q[0];
rz(-1.0547423) q[2];
sx q[2];
rz(-1.8212089) q[2];
sx q[2];
rz(2.5977787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1198188) q[1];
sx q[1];
rz(-2.0975862) q[1];
sx q[1];
rz(3.0837545) q[1];
rz(-pi) q[2];
rz(1.9366185) q[3];
sx q[3];
rz(-2.1064039) q[3];
sx q[3];
rz(2.7887695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1043642) q[2];
sx q[2];
rz(-1.2040441) q[2];
sx q[2];
rz(-2.6715703) q[2];
rz(0.59703279) q[3];
sx q[3];
rz(-0.11490122) q[3];
sx q[3];
rz(-1.6980096) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7715348) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(0.22802995) q[0];
rz(2.6920964) q[1];
sx q[1];
rz(-2.1411965) q[1];
sx q[1];
rz(0.94625783) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3458918) q[0];
sx q[0];
rz(-1.8270703) q[0];
sx q[0];
rz(-1.9836224) q[0];
x q[1];
rz(-0.65785711) q[2];
sx q[2];
rz(-1.0491174) q[2];
sx q[2];
rz(-1.9470095) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1552004) q[1];
sx q[1];
rz(-2.7554535) q[1];
sx q[1];
rz(2.6394352) q[1];
rz(-pi) q[2];
rz(1.069293) q[3];
sx q[3];
rz(-0.32078241) q[3];
sx q[3];
rz(-2.4298885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9049282) q[2];
sx q[2];
rz(-1.6764574) q[2];
sx q[2];
rz(1.6620592) q[2];
rz(2.9605401) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(-0.886206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7774696) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(-2.2430578) q[0];
rz(-0.50615519) q[1];
sx q[1];
rz(-1.2944784) q[1];
sx q[1];
rz(-1.6815965) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7285883) q[0];
sx q[0];
rz(-2.6046136) q[0];
sx q[0];
rz(3.0130638) q[0];
rz(-pi) q[1];
rz(-2.3982749) q[2];
sx q[2];
rz(-1.2365336) q[2];
sx q[2];
rz(-1.0570132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0702563) q[1];
sx q[1];
rz(-1.6410636) q[1];
sx q[1];
rz(1.6855168) q[1];
x q[2];
rz(-0.060486501) q[3];
sx q[3];
rz(-2.2248259) q[3];
sx q[3];
rz(-0.56246434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8755181) q[2];
sx q[2];
rz(-2.6660599) q[2];
sx q[2];
rz(-1.5024705) q[2];
rz(-1.8858887) q[3];
sx q[3];
rz(-1.7399104) q[3];
sx q[3];
rz(3.0953395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9083549) q[0];
sx q[0];
rz(-2.4112356) q[0];
sx q[0];
rz(-2.9610942) q[0];
rz(1.3795229) q[1];
sx q[1];
rz(-0.5677529) q[1];
sx q[1];
rz(-1.0951805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4783673) q[0];
sx q[0];
rz(-2.992482) q[0];
sx q[0];
rz(0.43004934) q[0];
rz(1.6542475) q[2];
sx q[2];
rz(-0.7755643) q[2];
sx q[2];
rz(-2.6706539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95848786) q[1];
sx q[1];
rz(-2.339683) q[1];
sx q[1];
rz(-1.3758749) q[1];
rz(-pi) q[2];
rz(3.0083382) q[3];
sx q[3];
rz(-2.2550003) q[3];
sx q[3];
rz(-1.2017991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77202648) q[2];
sx q[2];
rz(-0.84731805) q[2];
sx q[2];
rz(2.0728716) q[2];
rz(0.42701834) q[3];
sx q[3];
rz(-2.2618099) q[3];
sx q[3];
rz(1.8900185) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0291979) q[0];
sx q[0];
rz(-2.2258832) q[0];
sx q[0];
rz(2.8503964) q[0];
rz(1.4578106) q[1];
sx q[1];
rz(-1.0271415) q[1];
sx q[1];
rz(0.88868946) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5425974) q[0];
sx q[0];
rz(-0.7805191) q[0];
sx q[0];
rz(-2.2602343) q[0];
rz(-0.75727792) q[2];
sx q[2];
rz(-1.5520397) q[2];
sx q[2];
rz(0.64964408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0093569) q[1];
sx q[1];
rz(-2.0814544) q[1];
sx q[1];
rz(2.4978994) q[1];
rz(-pi) q[2];
rz(-1.8848041) q[3];
sx q[3];
rz(-1.8591789) q[3];
sx q[3];
rz(0.42321229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7500744) q[2];
sx q[2];
rz(-1.3041648) q[2];
sx q[2];
rz(0.94775003) q[2];
rz(-0.1968955) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(0.31015629) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8295558) q[0];
sx q[0];
rz(-2.0289679) q[0];
sx q[0];
rz(-0.46992508) q[0];
rz(-1.6795109) q[1];
sx q[1];
rz(-1.8406248) q[1];
sx q[1];
rz(1.4039111) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023177308) q[0];
sx q[0];
rz(-1.0141393) q[0];
sx q[0];
rz(-2.7283637) q[0];
x q[1];
rz(1.1833625) q[2];
sx q[2];
rz(-2.63113) q[2];
sx q[2];
rz(-2.032377) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7498878) q[1];
sx q[1];
rz(-1.0467231) q[1];
sx q[1];
rz(-0.4197555) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2265731) q[3];
sx q[3];
rz(-2.1267517) q[3];
sx q[3];
rz(1.2840301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64138428) q[2];
sx q[2];
rz(-0.83157867) q[2];
sx q[2];
rz(0.10759648) q[2];
rz(2.522116) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(-1.7776325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.7463995) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-0.25303823) q[0];
rz(0.088317618) q[1];
sx q[1];
rz(-2.2679236) q[1];
sx q[1];
rz(-2.1926682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509453) q[0];
sx q[0];
rz(-1.4207736) q[0];
sx q[0];
rz(1.2706744) q[0];
x q[1];
rz(-1.9789655) q[2];
sx q[2];
rz(-1.7498831) q[2];
sx q[2];
rz(-2.8286752) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1469377) q[1];
sx q[1];
rz(-2.6587464) q[1];
sx q[1];
rz(-1.6864683) q[1];
x q[2];
rz(2.3036495) q[3];
sx q[3];
rz(-2.4334014) q[3];
sx q[3];
rz(1.6745245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7439612) q[2];
sx q[2];
rz(-0.88699114) q[2];
sx q[2];
rz(0.31069791) q[2];
rz(-0.24142309) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(2.0588622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13984118) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(1.235442) q[0];
rz(0.48108092) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(-1.4280041) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1989216) q[0];
sx q[0];
rz(-0.81139542) q[0];
sx q[0];
rz(-1.5553724) q[0];
rz(-pi) q[1];
rz(-1.8495103) q[2];
sx q[2];
rz(-1.3506417) q[2];
sx q[2];
rz(-2.1627592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2216827) q[1];
sx q[1];
rz(-0.35791985) q[1];
sx q[1];
rz(-0.81166761) q[1];
rz(2.7747267) q[3];
sx q[3];
rz(-2.3847724) q[3];
sx q[3];
rz(0.048710195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1741751) q[2];
sx q[2];
rz(-0.82505161) q[2];
sx q[2];
rz(0.027912557) q[2];
rz(2.2715955) q[3];
sx q[3];
rz(-0.78012192) q[3];
sx q[3];
rz(-1.9546668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025539909) q[0];
sx q[0];
rz(-1.5929796) q[0];
sx q[0];
rz(2.0822339) q[0];
rz(-1.7268044) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(2.6639604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1519417) q[0];
sx q[0];
rz(-2.383259) q[0];
sx q[0];
rz(-2.9248021) q[0];
rz(-pi) q[1];
rz(-1.9066493) q[2];
sx q[2];
rz(-1.5689465) q[2];
sx q[2];
rz(0.84650485) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3463626) q[1];
sx q[1];
rz(-1.5168961) q[1];
sx q[1];
rz(1.4583711) q[1];
x q[2];
rz(-2.2077363) q[3];
sx q[3];
rz(-0.85784405) q[3];
sx q[3];
rz(-1.0757004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.16964218) q[2];
sx q[2];
rz(-0.98196882) q[2];
sx q[2];
rz(0.36894813) q[2];
rz(-1.87489) q[3];
sx q[3];
rz(-0.72487512) q[3];
sx q[3];
rz(0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13851588) q[0];
sx q[0];
rz(-1.930548) q[0];
sx q[0];
rz(1.7532274) q[0];
rz(-0.44454642) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(-2.6981392) q[2];
sx q[2];
rz(-1.8691861) q[2];
sx q[2];
rz(-0.791823) q[2];
rz(-3.1034341) q[3];
sx q[3];
rz(-0.46783075) q[3];
sx q[3];
rz(-0.086188407) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
