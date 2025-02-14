OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2990155) q[0];
sx q[0];
rz(-2.9691073) q[0];
sx q[0];
rz(-0.01699288) q[0];
rz(3.6530082) q[1];
sx q[1];
rz(3.6379171) q[1];
sx q[1];
rz(12.080893) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4776909) q[0];
sx q[0];
rz(-1.6833651) q[0];
sx q[0];
rz(-1.4590053) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83037776) q[2];
sx q[2];
rz(-0.86566209) q[2];
sx q[2];
rz(-0.0041088897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6852666) q[1];
sx q[1];
rz(-2.2783372) q[1];
sx q[1];
rz(-2.4601542) q[1];
x q[2];
rz(-1.2608246) q[3];
sx q[3];
rz(-1.4616218) q[3];
sx q[3];
rz(2.7742164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39712507) q[2];
sx q[2];
rz(-3.0391389) q[2];
sx q[2];
rz(-2.3490119) q[2];
rz(-0.98627311) q[3];
sx q[3];
rz(-1.4873742) q[3];
sx q[3];
rz(-3.0084394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2546286) q[0];
sx q[0];
rz(-1.6034842) q[0];
sx q[0];
rz(3.0358553) q[0];
rz(2.3836783) q[1];
sx q[1];
rz(-0.71280232) q[1];
sx q[1];
rz(-0.47592083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4315368) q[0];
sx q[0];
rz(-1.3366295) q[0];
sx q[0];
rz(-1.6706927) q[0];
x q[1];
rz(1.4272825) q[2];
sx q[2];
rz(-1.0821618) q[2];
sx q[2];
rz(1.8861063) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1913174) q[1];
sx q[1];
rz(-1.4498267) q[1];
sx q[1];
rz(-2.1410393) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.845221) q[3];
sx q[3];
rz(-2.1934273) q[3];
sx q[3];
rz(0.53106703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9953352) q[2];
sx q[2];
rz(-0.47440752) q[2];
sx q[2];
rz(1.8246626) q[2];
rz(-2.3138192) q[3];
sx q[3];
rz(-2.0483978) q[3];
sx q[3];
rz(-2.304346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4227609) q[0];
sx q[0];
rz(-1.7976924) q[0];
sx q[0];
rz(1.9047009) q[0];
rz(1.749136) q[1];
sx q[1];
rz(-1.4207276) q[1];
sx q[1];
rz(0.69033355) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25094098) q[0];
sx q[0];
rz(-1.1162045) q[0];
sx q[0];
rz(-0.62418749) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1074781) q[2];
sx q[2];
rz(-1.5743739) q[2];
sx q[2];
rz(-1.5114552) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.54385932) q[1];
sx q[1];
rz(-1.3784474) q[1];
sx q[1];
rz(1.7406157) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4776633) q[3];
sx q[3];
rz(-0.85278836) q[3];
sx q[3];
rz(-1.688886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51848015) q[2];
sx q[2];
rz(-1.6192351) q[2];
sx q[2];
rz(3.0677262) q[2];
rz(1.1582003) q[3];
sx q[3];
rz(-1.2169714) q[3];
sx q[3];
rz(-1.4069675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2288007) q[0];
sx q[0];
rz(-0.055963628) q[0];
sx q[0];
rz(2.9486616) q[0];
rz(1.8979161) q[1];
sx q[1];
rz(-1.7453777) q[1];
sx q[1];
rz(0.23304932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4503964) q[0];
sx q[0];
rz(-1.47627) q[0];
sx q[0];
rz(-0.14733845) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5265756) q[2];
sx q[2];
rz(-0.70927519) q[2];
sx q[2];
rz(1.0929293) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8687081) q[1];
sx q[1];
rz(-1.2603425) q[1];
sx q[1];
rz(1.5116879) q[1];
rz(-pi) q[2];
rz(-0.95121164) q[3];
sx q[3];
rz(-1.1780292) q[3];
sx q[3];
rz(0.58863607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.21726808) q[2];
sx q[2];
rz(-1.7622207) q[2];
sx q[2];
rz(-0.064083727) q[2];
rz(-2.890375) q[3];
sx q[3];
rz(-2.678674) q[3];
sx q[3];
rz(0.29138756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083704405) q[0];
sx q[0];
rz(-1.8479481) q[0];
sx q[0];
rz(-1.2842913) q[0];
rz(3.1341556) q[1];
sx q[1];
rz(-1.1859272) q[1];
sx q[1];
rz(2.1844905) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.590132) q[0];
sx q[0];
rz(-0.55319417) q[0];
sx q[0];
rz(0.93450089) q[0];
x q[1];
rz(0.96896521) q[2];
sx q[2];
rz(-1.0442248) q[2];
sx q[2];
rz(-0.32770448) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51077183) q[1];
sx q[1];
rz(-1.242852) q[1];
sx q[1];
rz(-1.6751218) q[1];
rz(-pi) q[2];
rz(1.383833) q[3];
sx q[3];
rz(-2.2174944) q[3];
sx q[3];
rz(-1.9892429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0127504) q[2];
sx q[2];
rz(-0.70657554) q[2];
sx q[2];
rz(-1.0303222) q[2];
rz(-2.4230867) q[3];
sx q[3];
rz(-1.0822783) q[3];
sx q[3];
rz(0.70657402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71501032) q[0];
sx q[0];
rz(-2.3937245) q[0];
sx q[0];
rz(-2.7744875) q[0];
rz(2.7857065) q[1];
sx q[1];
rz(-1.117319) q[1];
sx q[1];
rz(1.5918559) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132431) q[0];
sx q[0];
rz(-1.5322881) q[0];
sx q[0];
rz(0.74158494) q[0];
rz(2.0075304) q[2];
sx q[2];
rz(-1.6730089) q[2];
sx q[2];
rz(-1.5222766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15799668) q[1];
sx q[1];
rz(-1.3648197) q[1];
sx q[1];
rz(-2.7218444) q[1];
rz(0.14629062) q[3];
sx q[3];
rz(-2.6669569) q[3];
sx q[3];
rz(-0.65505799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.97673577) q[2];
sx q[2];
rz(-1.2931436) q[2];
sx q[2];
rz(-3.1381651) q[2];
rz(0.88611832) q[3];
sx q[3];
rz(-2.0485853) q[3];
sx q[3];
rz(1.6974983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7057328) q[0];
sx q[0];
rz(-1.4737361) q[0];
sx q[0];
rz(-2.4263897) q[0];
rz(-3.0192979) q[1];
sx q[1];
rz(-1.0194174) q[1];
sx q[1];
rz(-0.14762793) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0334761) q[0];
sx q[0];
rz(-1.8980935) q[0];
sx q[0];
rz(1.9294338) q[0];
x q[1];
rz(0.72509693) q[2];
sx q[2];
rz(-2.5864961) q[2];
sx q[2];
rz(1.8738418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1654037) q[1];
sx q[1];
rz(-2.2641085) q[1];
sx q[1];
rz(2.90392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.070402831) q[3];
sx q[3];
rz(-1.8569267) q[3];
sx q[3];
rz(-2.4918258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3132396) q[2];
sx q[2];
rz(-0.30892631) q[2];
sx q[2];
rz(0.1114791) q[2];
rz(-0.76505032) q[3];
sx q[3];
rz(-1.4234411) q[3];
sx q[3];
rz(-0.033871977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2400804) q[0];
sx q[0];
rz(-0.96918786) q[0];
sx q[0];
rz(-1.0028268) q[0];
rz(-0.78549939) q[1];
sx q[1];
rz(-1.5248884) q[1];
sx q[1];
rz(-1.921382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5400235) q[0];
sx q[0];
rz(-0.51443729) q[0];
sx q[0];
rz(-3.0995447) q[0];
x q[1];
rz(1.6613879) q[2];
sx q[2];
rz(-1.633051) q[2];
sx q[2];
rz(1.7287776) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.073382811) q[1];
sx q[1];
rz(-1.0676976) q[1];
sx q[1];
rz(-2.7791609) q[1];
x q[2];
rz(1.3429014) q[3];
sx q[3];
rz(-2.493495) q[3];
sx q[3];
rz(-0.15290393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9057374) q[2];
sx q[2];
rz(-1.3964272) q[2];
sx q[2];
rz(-3.1077969) q[2];
rz(-2.3679521) q[3];
sx q[3];
rz(-1.528911) q[3];
sx q[3];
rz(1.0627559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3646024) q[0];
sx q[0];
rz(-2.5211625) q[0];
sx q[0];
rz(-0.34238368) q[0];
rz(-1.4211593) q[1];
sx q[1];
rz(-0.8232638) q[1];
sx q[1];
rz(1.8470496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7086805) q[0];
sx q[0];
rz(-1.4231235) q[0];
sx q[0];
rz(2.7765034) q[0];
x q[1];
rz(2.2850288) q[2];
sx q[2];
rz(-2.1571549) q[2];
sx q[2];
rz(-2.7991653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0755989) q[1];
sx q[1];
rz(-2.3477049) q[1];
sx q[1];
rz(0.17236472) q[1];
rz(2.8043069) q[3];
sx q[3];
rz(-1.4082068) q[3];
sx q[3];
rz(-0.43460571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8483868) q[2];
sx q[2];
rz(-1.7561971) q[2];
sx q[2];
rz(2.6050341) q[2];
rz(0.41283354) q[3];
sx q[3];
rz(-2.6091913) q[3];
sx q[3];
rz(2.0280973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8263016) q[0];
sx q[0];
rz(-1.2696126) q[0];
sx q[0];
rz(1.6444561) q[0];
rz(-2.3313088) q[1];
sx q[1];
rz(-1.5565926) q[1];
sx q[1];
rz(0.84698814) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772781) q[0];
sx q[0];
rz(-1.7314096) q[0];
sx q[0];
rz(1.3993652) q[0];
rz(-2.2184847) q[2];
sx q[2];
rz(-2.0294242) q[2];
sx q[2];
rz(-3.0748526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21928043) q[1];
sx q[1];
rz(-2.4488291) q[1];
sx q[1];
rz(2.7119066) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8448006) q[3];
sx q[3];
rz(-1.4782259) q[3];
sx q[3];
rz(1.1082197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0392796) q[2];
sx q[2];
rz(-1.3156834) q[2];
sx q[2];
rz(0.63423356) q[2];
rz(-0.63747326) q[3];
sx q[3];
rz(-2.0135148) q[3];
sx q[3];
rz(-2.9243961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98916003) q[0];
sx q[0];
rz(-1.9518873) q[0];
sx q[0];
rz(1.963203) q[0];
rz(2.4329026) q[1];
sx q[1];
rz(-1.7955753) q[1];
sx q[1];
rz(1.2830455) q[1];
rz(1.218956) q[2];
sx q[2];
rz(-1.2165804) q[2];
sx q[2];
rz(2.658398) q[2];
rz(1.9509964) q[3];
sx q[3];
rz(-2.1722542) q[3];
sx q[3];
rz(0.62361591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
