OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.37261) q[0];
sx q[0];
rz(-0.045172673) q[0];
sx q[0];
rz(2.4700408) q[0];
rz(-0.99611941) q[1];
sx q[1];
rz(-2.3528407) q[1];
sx q[1];
rz(2.8532343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.45123) q[0];
sx q[0];
rz(-1.7685862) q[0];
sx q[0];
rz(3.1044699) q[0];
rz(-0.66418437) q[2];
sx q[2];
rz(-1.3553047) q[2];
sx q[2];
rz(-1.8813949) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8858998) q[1];
sx q[1];
rz(-1.608537) q[1];
sx q[1];
rz(-2.6934212) q[1];
x q[2];
rz(-2.6872271) q[3];
sx q[3];
rz(-1.4404669) q[3];
sx q[3];
rz(-2.5443175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4396189) q[2];
sx q[2];
rz(-2.7957323) q[2];
sx q[2];
rz(0.71887476) q[2];
rz(1.4465205) q[3];
sx q[3];
rz(-1.6547763) q[3];
sx q[3];
rz(1.0369302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0816536) q[0];
sx q[0];
rz(-0.43123284) q[0];
sx q[0];
rz(-0.28847873) q[0];
rz(-0.55229315) q[1];
sx q[1];
rz(-2.0497132) q[1];
sx q[1];
rz(1.9658032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3551153) q[0];
sx q[0];
rz(-1.2594495) q[0];
sx q[0];
rz(0.15736736) q[0];
rz(-pi) q[1];
rz(-1.3973049) q[2];
sx q[2];
rz(-1.6855006) q[2];
sx q[2];
rz(-2.3115273) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9811444) q[1];
sx q[1];
rz(-1.6949883) q[1];
sx q[1];
rz(-2.1500514) q[1];
x q[2];
rz(-1.6064241) q[3];
sx q[3];
rz(-3.1094915) q[3];
sx q[3];
rz(-0.42313448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4484078) q[2];
sx q[2];
rz(-1.8790481) q[2];
sx q[2];
rz(-0.022424879) q[2];
rz(1.4261931) q[3];
sx q[3];
rz(-1.2778927) q[3];
sx q[3];
rz(-2.2414331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65396032) q[0];
sx q[0];
rz(-1.910169) q[0];
sx q[0];
rz(0.76474977) q[0];
rz(1.359831) q[1];
sx q[1];
rz(-1.9197074) q[1];
sx q[1];
rz(1.8870032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5072637) q[0];
sx q[0];
rz(-1.6741236) q[0];
sx q[0];
rz(-1.4402585) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32439442) q[2];
sx q[2];
rz(-1.430871) q[2];
sx q[2];
rz(0.10332271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0618825) q[1];
sx q[1];
rz(-1.8593504) q[1];
sx q[1];
rz(2.6614019) q[1];
x q[2];
rz(-2.7990667) q[3];
sx q[3];
rz(-1.4493128) q[3];
sx q[3];
rz(3.1241724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4059056) q[2];
sx q[2];
rz(-0.76852208) q[2];
sx q[2];
rz(-3.1206257) q[2];
rz(1.5812801) q[3];
sx q[3];
rz(-2.0798648) q[3];
sx q[3];
rz(0.90644065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86879325) q[0];
sx q[0];
rz(-0.11196207) q[0];
sx q[0];
rz(1.3695166) q[0];
rz(0.70961332) q[1];
sx q[1];
rz(-1.8636924) q[1];
sx q[1];
rz(0.69724625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6402138) q[0];
sx q[0];
rz(-1.9778628) q[0];
sx q[0];
rz(1.2329007) q[0];
rz(-0.68402779) q[2];
sx q[2];
rz(-2.1993756) q[2];
sx q[2];
rz(-0.39234871) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5647506) q[1];
sx q[1];
rz(-0.71755845) q[1];
sx q[1];
rz(-2.0400042) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94631291) q[3];
sx q[3];
rz(-2.8528004) q[3];
sx q[3];
rz(-2.2594947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1856445) q[2];
sx q[2];
rz(-0.94719013) q[2];
sx q[2];
rz(-2.4626125) q[2];
rz(-2.0164356) q[3];
sx q[3];
rz(-1.1449287) q[3];
sx q[3];
rz(2.9185435) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0515902) q[0];
sx q[0];
rz(-1.2035878) q[0];
sx q[0];
rz(-0.077127174) q[0];
rz(1.9902825) q[1];
sx q[1];
rz(-1.5682033) q[1];
sx q[1];
rz(-0.77879771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4753195) q[0];
sx q[0];
rz(-0.059564807) q[0];
sx q[0];
rz(2.2044529) q[0];
rz(-2.0955032) q[2];
sx q[2];
rz(-0.32654844) q[2];
sx q[2];
rz(2.7401217) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8647789) q[1];
sx q[1];
rz(-1.2363696) q[1];
sx q[1];
rz(2.8825376) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5803189) q[3];
sx q[3];
rz(-1.408934) q[3];
sx q[3];
rz(-1.7603742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.516958) q[2];
sx q[2];
rz(-1.7004852) q[2];
sx q[2];
rz(1.1901633) q[2];
rz(0.55365753) q[3];
sx q[3];
rz(-2.5167969) q[3];
sx q[3];
rz(-2.4042118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5364285) q[0];
sx q[0];
rz(-1.9909415) q[0];
sx q[0];
rz(0.27106699) q[0];
rz(-1.0003264) q[1];
sx q[1];
rz(-1.2157636) q[1];
sx q[1];
rz(2.2727374) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0603179) q[0];
sx q[0];
rz(-2.3485614) q[0];
sx q[0];
rz(-0.6834553) q[0];
x q[1];
rz(-1.4914091) q[2];
sx q[2];
rz(-0.83680162) q[2];
sx q[2];
rz(0.47576093) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7160401) q[1];
sx q[1];
rz(-1.2231266) q[1];
sx q[1];
rz(-0.1802318) q[1];
rz(-pi) q[2];
rz(-0.47021659) q[3];
sx q[3];
rz(-2.3002671) q[3];
sx q[3];
rz(0.39455345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.102999) q[2];
sx q[2];
rz(-0.30240348) q[2];
sx q[2];
rz(1.137286) q[2];
rz(0.31050995) q[3];
sx q[3];
rz(-2.2313084) q[3];
sx q[3];
rz(-0.92946068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0248658) q[0];
sx q[0];
rz(-1.4288582) q[0];
sx q[0];
rz(2.0615935) q[0];
rz(0.96177167) q[1];
sx q[1];
rz(-0.56517833) q[1];
sx q[1];
rz(0.22294179) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7540383) q[0];
sx q[0];
rz(-1.6056653) q[0];
sx q[0];
rz(3.1341482) q[0];
rz(2.8330363) q[2];
sx q[2];
rz(-2.1962104) q[2];
sx q[2];
rz(2.5297414) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5364616) q[1];
sx q[1];
rz(-0.85857449) q[1];
sx q[1];
rz(-2.773316) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1268397) q[3];
sx q[3];
rz(-0.97501576) q[3];
sx q[3];
rz(-0.19150133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.96233931) q[2];
sx q[2];
rz(-0.27955678) q[2];
sx q[2];
rz(-1.1780098) q[2];
rz(2.8152605) q[3];
sx q[3];
rz(-0.81063619) q[3];
sx q[3];
rz(-2.0012205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2278263) q[0];
sx q[0];
rz(-0.95781177) q[0];
sx q[0];
rz(-0.80818278) q[0];
rz(-0.35762865) q[1];
sx q[1];
rz(-1.459815) q[1];
sx q[1];
rz(1.0825895) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9900695) q[0];
sx q[0];
rz(-0.7034348) q[0];
sx q[0];
rz(-1.8358873) q[0];
rz(-pi) q[1];
x q[1];
rz(2.888526) q[2];
sx q[2];
rz(-1.0326721) q[2];
sx q[2];
rz(2.5122364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1114166) q[1];
sx q[1];
rz(-1.7716265) q[1];
sx q[1];
rz(2.3371731) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3305835) q[3];
sx q[3];
rz(-2.5222416) q[3];
sx q[3];
rz(0.15796433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.064934405) q[2];
sx q[2];
rz(-0.4883464) q[2];
sx q[2];
rz(-1.9261544) q[2];
rz(-0.91228929) q[3];
sx q[3];
rz(-0.89022294) q[3];
sx q[3];
rz(0.54273763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2927581) q[0];
sx q[0];
rz(-2.500535) q[0];
sx q[0];
rz(-2.7178398) q[0];
rz(1.1774225) q[1];
sx q[1];
rz(-0.91554987) q[1];
sx q[1];
rz(2.7659069) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59058023) q[0];
sx q[0];
rz(-0.19119054) q[0];
sx q[0];
rz(-1.952233) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5440953) q[2];
sx q[2];
rz(-1.809568) q[2];
sx q[2];
rz(-1.2899866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2089952) q[1];
sx q[1];
rz(-1.0699341) q[1];
sx q[1];
rz(-0.074525699) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71640941) q[3];
sx q[3];
rz(-1.3103974) q[3];
sx q[3];
rz(-1.6668056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0287013) q[2];
sx q[2];
rz(-1.6501004) q[2];
sx q[2];
rz(0.69226199) q[2];
rz(2.4568457) q[3];
sx q[3];
rz(-2.1919577) q[3];
sx q[3];
rz(-0.14028604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.5570062) q[0];
sx q[0];
rz(-2.1742915) q[0];
sx q[0];
rz(-1.517357) q[0];
rz(-1.6035621) q[1];
sx q[1];
rz(-1.3471194) q[1];
sx q[1];
rz(-1.921152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081471215) q[0];
sx q[0];
rz(-0.94917008) q[0];
sx q[0];
rz(-3.0682878) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8062808) q[2];
sx q[2];
rz(-1.0808498) q[2];
sx q[2];
rz(1.3563434) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5511849) q[1];
sx q[1];
rz(-0.77279323) q[1];
sx q[1];
rz(0.40829746) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8655928) q[3];
sx q[3];
rz(-2.285897) q[3];
sx q[3];
rz(3.0218647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5662235) q[2];
sx q[2];
rz(-1.0074002) q[2];
sx q[2];
rz(2.2929906) q[2];
rz(2.1620915) q[3];
sx q[3];
rz(-1.5749911) q[3];
sx q[3];
rz(-3.1041253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033087) q[0];
sx q[0];
rz(-1.6642234) q[0];
sx q[0];
rz(-2.4432175) q[0];
rz(2.5820844) q[1];
sx q[1];
rz(-0.63897501) q[1];
sx q[1];
rz(-0.41313304) q[1];
rz(-0.32158659) q[2];
sx q[2];
rz(-1.814331) q[2];
sx q[2];
rz(-0.41650256) q[2];
rz(0.57942617) q[3];
sx q[3];
rz(-1.5204932) q[3];
sx q[3];
rz(1.6411171) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
