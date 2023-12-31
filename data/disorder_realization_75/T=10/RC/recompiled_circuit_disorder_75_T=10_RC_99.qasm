OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50492935) q[0];
sx q[0];
rz(-1.3267656) q[0];
sx q[0];
rz(1.7457477) q[0];
rz(-pi) q[1];
rz(-2.6909157) q[2];
sx q[2];
rz(-0.55826954) q[2];
sx q[2];
rz(0.4980586) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.287392) q[1];
sx q[1];
rz(-2.5890447) q[1];
sx q[1];
rz(-1.470675) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8098104) q[3];
sx q[3];
rz(-2.3506769) q[3];
sx q[3];
rz(-3.0856109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(0.35787004) q[2];
rz(2.9499124) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0033922694) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(3.0733118) q[0];
rz(1.0066907) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(0.65111792) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40819528) q[0];
sx q[0];
rz(-1.3696284) q[0];
sx q[0];
rz(0.64543076) q[0];
x q[1];
rz(0.5257734) q[2];
sx q[2];
rz(-1.5915046) q[2];
sx q[2];
rz(-2.6244342) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.843833) q[1];
sx q[1];
rz(-0.69121541) q[1];
sx q[1];
rz(-0.83696951) q[1];
x q[2];
rz(2.6504374) q[3];
sx q[3];
rz(-1.3782116) q[3];
sx q[3];
rz(1.1112978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45289257) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(0.90332705) q[2];
rz(2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2925401) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(0.49355155) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(1.0292056) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619336) q[0];
sx q[0];
rz(-0.56549457) q[0];
sx q[0];
rz(-0.66753597) q[0];
x q[1];
rz(2.7846787) q[2];
sx q[2];
rz(-1.1733574) q[2];
sx q[2];
rz(0.39427653) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2521816) q[1];
sx q[1];
rz(-1.2671789) q[1];
sx q[1];
rz(-1.7905411) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0164099) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(-0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-2.2600007) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13215412) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(2.8892345) q[0];
x q[1];
rz(-0.84882952) q[2];
sx q[2];
rz(-1.0366882) q[2];
sx q[2];
rz(-1.3918849) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91729212) q[1];
sx q[1];
rz(-0.91887337) q[1];
sx q[1];
rz(-2.2071381) q[1];
rz(-pi) q[2];
rz(2.9922585) q[3];
sx q[3];
rz(-0.42156005) q[3];
sx q[3];
rz(-1.644852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(-2.0289452) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(1.5295193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(-3.024335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1299767) q[0];
sx q[0];
rz(-2.5419309) q[0];
sx q[0];
rz(0.2325124) q[0];
rz(-2.8357382) q[2];
sx q[2];
rz(-2.5161985) q[2];
sx q[2];
rz(1.2139699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.458141) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(-1.4139921) q[1];
rz(-1.6860784) q[3];
sx q[3];
rz(-2.5513463) q[3];
sx q[3];
rz(2.3595927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(-1.3797181) q[2];
rz(-1.951925) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(-0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(-0.70621079) q[0];
rz(-1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(-3.0336753) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7100922) q[0];
sx q[0];
rz(-1.8653231) q[0];
sx q[0];
rz(-0.72580238) q[0];
rz(-pi) q[1];
rz(2.8651587) q[2];
sx q[2];
rz(-0.041639608) q[2];
sx q[2];
rz(-2.4561938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95277112) q[1];
sx q[1];
rz(-1.2880039) q[1];
sx q[1];
rz(-2.4464843) q[1];
rz(-pi) q[2];
rz(-2.1351486) q[3];
sx q[3];
rz(-2.3688149) q[3];
sx q[3];
rz(0.3533065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32249054) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(1.8479944) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.368822) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(-3.1266881) q[0];
rz(-0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(2.3419535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470456) q[0];
sx q[0];
rz(-0.73091113) q[0];
sx q[0];
rz(-0.38696179) q[0];
rz(-pi) q[1];
rz(-2.8268379) q[2];
sx q[2];
rz(-1.2999279) q[2];
sx q[2];
rz(2.1115007) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.423973) q[1];
sx q[1];
rz(-1.4255925) q[1];
sx q[1];
rz(1.6270301) q[1];
rz(3.0040967) q[3];
sx q[3];
rz(-1.0097479) q[3];
sx q[3];
rz(1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(2.2195393) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9668982) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(2.877537) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(-2.82428) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70586328) q[0];
sx q[0];
rz(-1.1840491) q[0];
sx q[0];
rz(0.078588967) q[0];
rz(-pi) q[1];
rz(-1.9900471) q[2];
sx q[2];
rz(-2.7630685) q[2];
sx q[2];
rz(-0.35767698) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74097733) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(-0.11465794) q[1];
x q[2];
rz(2.4962037) q[3];
sx q[3];
rz(-1.7040952) q[3];
sx q[3];
rz(0.79750878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7060966) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(-2.8912985) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89649993) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(1.7171575) q[0];
rz(-2.6622488) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(3.0260578) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45142052) q[0];
sx q[0];
rz(-0.85309404) q[0];
sx q[0];
rz(0.79057981) q[0];
x q[1];
rz(-0.71474448) q[2];
sx q[2];
rz(-1.7822767) q[2];
sx q[2];
rz(1.5541058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2781196) q[1];
sx q[1];
rz(-1.1520471) q[1];
sx q[1];
rz(0.31660415) q[1];
rz(2.9980738) q[3];
sx q[3];
rz(-2.4402938) q[3];
sx q[3];
rz(-1.6539751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.7158562) q[2];
rz(-1.4005631) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(0.28276309) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(0.26915959) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.3795308) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36294286) q[0];
sx q[0];
rz(-1.6761707) q[0];
sx q[0];
rz(-1.7371348) q[0];
x q[1];
rz(-3.0616947) q[2];
sx q[2];
rz(-1.986984) q[2];
sx q[2];
rz(1.4518567) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3456612) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(-1.3962586) q[1];
rz(2.3926211) q[3];
sx q[3];
rz(-0.43991551) q[3];
sx q[3];
rz(-2.9797152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0925838) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.5157549) q[2];
rz(1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-2.1102171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9243069) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(-2.6869607) q[1];
sx q[1];
rz(-2.0352719) q[1];
sx q[1];
rz(-0.24771053) q[1];
rz(0.17299962) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(-1.7975939) q[3];
sx q[3];
rz(-0.26567017) q[3];
sx q[3];
rz(-2.9220823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
