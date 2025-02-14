OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5934264) q[0];
sx q[0];
rz(-0.46841535) q[0];
sx q[0];
rz(-3.0031437) q[0];
rz(-1.2568714) q[1];
sx q[1];
rz(-1.0705907) q[1];
sx q[1];
rz(-3.1201153) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81651743) q[0];
sx q[0];
rz(-1.7735574) q[0];
sx q[0];
rz(-1.2120423) q[0];
rz(0.98928605) q[2];
sx q[2];
rz(-2.044319) q[2];
sx q[2];
rz(-3.1030637) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8565377) q[1];
sx q[1];
rz(-1.7453339) q[1];
sx q[1];
rz(-0.15600295) q[1];
rz(1.2091694) q[3];
sx q[3];
rz(-2.1361793) q[3];
sx q[3];
rz(-0.20547444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5452177) q[2];
sx q[2];
rz(-0.45520982) q[2];
sx q[2];
rz(1.5111766) q[2];
rz(0.049985416) q[3];
sx q[3];
rz(-1.9129916) q[3];
sx q[3];
rz(-2.889192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.2216457) q[0];
sx q[0];
rz(-1.9939461) q[0];
sx q[0];
rz(0.037121437) q[0];
rz(3.1275753) q[1];
sx q[1];
rz(-2.5633096) q[1];
sx q[1];
rz(2.9055273) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0229313) q[0];
sx q[0];
rz(-0.076181024) q[0];
sx q[0];
rz(0.19019048) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0843956) q[2];
sx q[2];
rz(-0.91569967) q[2];
sx q[2];
rz(-0.72593216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29480793) q[1];
sx q[1];
rz(-2.4681512) q[1];
sx q[1];
rz(0.93017857) q[1];
x q[2];
rz(-0.97352435) q[3];
sx q[3];
rz(-2.7144538) q[3];
sx q[3];
rz(-1.1534302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1290805) q[2];
sx q[2];
rz(-1.7725638) q[2];
sx q[2];
rz(0.033626076) q[2];
rz(-0.29020894) q[3];
sx q[3];
rz(-2.2955194) q[3];
sx q[3];
rz(-2.7510551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.520312) q[0];
sx q[0];
rz(-0.97977591) q[0];
sx q[0];
rz(2.8424971) q[0];
rz(-1.6200804) q[1];
sx q[1];
rz(-0.027093096) q[1];
sx q[1];
rz(0.026195899) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1258419) q[0];
sx q[0];
rz(-1.5643018) q[0];
sx q[0];
rz(-1.5918713) q[0];
x q[1];
rz(-2.0355939) q[2];
sx q[2];
rz(-2.8299677) q[2];
sx q[2];
rz(2.0457884) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67945665) q[1];
sx q[1];
rz(-0.95210451) q[1];
sx q[1];
rz(-2.2652744) q[1];
rz(1.3407767) q[3];
sx q[3];
rz(-2.0072847) q[3];
sx q[3];
rz(-0.088584049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4574778) q[2];
sx q[2];
rz(-2.699615) q[2];
sx q[2];
rz(-2.5341865) q[2];
rz(-0.36951798) q[3];
sx q[3];
rz(-1.642546) q[3];
sx q[3];
rz(0.0079689715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58005106) q[0];
sx q[0];
rz(-2.9396368) q[0];
sx q[0];
rz(2.0107021) q[0];
rz(-2.789403) q[1];
sx q[1];
rz(-0.49030855) q[1];
sx q[1];
rz(0.21603781) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7894831) q[0];
sx q[0];
rz(-2.0228068) q[0];
sx q[0];
rz(-1.841196) q[0];
rz(-0.44682002) q[2];
sx q[2];
rz(-1.6761314) q[2];
sx q[2];
rz(0.0086431816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6747848) q[1];
sx q[1];
rz(-2.4744316) q[1];
sx q[1];
rz(2.2391367) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.400835) q[3];
sx q[3];
rz(-1.6584087) q[3];
sx q[3];
rz(2.5690998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.5242247) q[2];
sx q[2];
rz(-0.91297954) q[2];
sx q[2];
rz(-2.6101904) q[2];
rz(-1.4517387) q[3];
sx q[3];
rz(-0.51133358) q[3];
sx q[3];
rz(-1.0303729) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78469974) q[0];
sx q[0];
rz(-1.0837311) q[0];
sx q[0];
rz(2.8392131) q[0];
rz(1.1131635) q[1];
sx q[1];
rz(-1.0013564) q[1];
sx q[1];
rz(-2.0804292) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7776124) q[0];
sx q[0];
rz(-1.5060695) q[0];
sx q[0];
rz(1.6311247) q[0];
x q[1];
rz(2.7348379) q[2];
sx q[2];
rz(-2.8977721) q[2];
sx q[2];
rz(1.8610473) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6401599) q[1];
sx q[1];
rz(-1.6987015) q[1];
sx q[1];
rz(2.3588965) q[1];
rz(-pi) q[2];
rz(1.8228028) q[3];
sx q[3];
rz(-1.2815426) q[3];
sx q[3];
rz(0.94925971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82659668) q[2];
sx q[2];
rz(-1.7642517) q[2];
sx q[2];
rz(0.19485168) q[2];
rz(-0.74243122) q[3];
sx q[3];
rz(-0.73115474) q[3];
sx q[3];
rz(2.8661695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18405296) q[0];
sx q[0];
rz(-0.66127151) q[0];
sx q[0];
rz(-1.1141962) q[0];
rz(-0.30955744) q[1];
sx q[1];
rz(-0.93300262) q[1];
sx q[1];
rz(-1.385744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4969869) q[0];
sx q[0];
rz(-0.33503767) q[0];
sx q[0];
rz(1.0100288) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3676332) q[2];
sx q[2];
rz(-0.51649714) q[2];
sx q[2];
rz(1.9488283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4655814) q[1];
sx q[1];
rz(-1.3106924) q[1];
sx q[1];
rz(2.4838402) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4900804) q[3];
sx q[3];
rz(-0.88884547) q[3];
sx q[3];
rz(1.7891974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4625357) q[2];
sx q[2];
rz(-0.44687301) q[2];
sx q[2];
rz(-0.40979579) q[2];
rz(-3.0637686) q[3];
sx q[3];
rz(-1.2160559) q[3];
sx q[3];
rz(2.3791544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78793144) q[0];
sx q[0];
rz(-2.9910112) q[0];
sx q[0];
rz(-1.4713564) q[0];
rz(2.3279066) q[1];
sx q[1];
rz(-2.4206471) q[1];
sx q[1];
rz(2.241316) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1293869) q[0];
sx q[0];
rz(-3.0086522) q[0];
sx q[0];
rz(2.6690527) q[0];
x q[1];
rz(1.8531606) q[2];
sx q[2];
rz(-2.8430364) q[2];
sx q[2];
rz(-0.82080847) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1301816) q[1];
sx q[1];
rz(-0.26618845) q[1];
sx q[1];
rz(-1.5397416) q[1];
rz(-1.5057218) q[3];
sx q[3];
rz(-1.3378007) q[3];
sx q[3];
rz(-2.9386793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5208931) q[2];
sx q[2];
rz(-0.40105477) q[2];
sx q[2];
rz(-0.14191423) q[2];
rz(2.9698931) q[3];
sx q[3];
rz(-1.2407691) q[3];
sx q[3];
rz(-2.5920674) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63355821) q[0];
sx q[0];
rz(-1.9867851) q[0];
sx q[0];
rz(-1.5800193) q[0];
rz(-2.436807) q[1];
sx q[1];
rz(-2.2403658) q[1];
sx q[1];
rz(-2.8660692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47833035) q[0];
sx q[0];
rz(-1.5092634) q[0];
sx q[0];
rz(-1.6986153) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1528831) q[2];
sx q[2];
rz(-2.8934921) q[2];
sx q[2];
rz(0.54923344) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6478579) q[1];
sx q[1];
rz(-1.0714515) q[1];
sx q[1];
rz(-0.33555056) q[1];
x q[2];
rz(-0.66657127) q[3];
sx q[3];
rz(-1.2351994) q[3];
sx q[3];
rz(-0.94387142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0813109) q[2];
sx q[2];
rz(-1.5014481) q[2];
sx q[2];
rz(-0.29067972) q[2];
rz(-2.3447013) q[3];
sx q[3];
rz(-2.6698038) q[3];
sx q[3];
rz(0.98381388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31442916) q[0];
sx q[0];
rz(-1.4845347) q[0];
sx q[0];
rz(0.86781251) q[0];
rz(0.21226352) q[1];
sx q[1];
rz(-0.8808732) q[1];
sx q[1];
rz(0.27264047) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51202392) q[0];
sx q[0];
rz(-1.4958515) q[0];
sx q[0];
rz(-1.7864947) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81959447) q[2];
sx q[2];
rz(-1.4965881) q[2];
sx q[2];
rz(-1.5641664) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8729108) q[1];
sx q[1];
rz(-0.93645331) q[1];
sx q[1];
rz(-0.047206248) q[1];
rz(-pi) q[2];
rz(2.7271184) q[3];
sx q[3];
rz(-1.2948841) q[3];
sx q[3];
rz(-1.3434354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0923882) q[2];
sx q[2];
rz(-2.7907351) q[2];
sx q[2];
rz(1.1445047) q[2];
rz(-0.62659621) q[3];
sx q[3];
rz(-2.383039) q[3];
sx q[3];
rz(-2.8265317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-2.3939683) q[0];
sx q[0];
rz(-2.4691041) q[0];
sx q[0];
rz(0.59338635) q[0];
rz(2.4001832) q[1];
sx q[1];
rz(-0.64677042) q[1];
sx q[1];
rz(-1.5470362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5760959) q[0];
sx q[0];
rz(-1.4617697) q[0];
sx q[0];
rz(-1.7897357) q[0];
rz(-pi) q[1];
rz(-0.52026622) q[2];
sx q[2];
rz(-2.6718585) q[2];
sx q[2];
rz(3.1343012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64338291) q[1];
sx q[1];
rz(-2.3003182) q[1];
sx q[1];
rz(1.6642844) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1118499) q[3];
sx q[3];
rz(-1.3689201) q[3];
sx q[3];
rz(1.0591398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0934304) q[2];
sx q[2];
rz(-0.54485816) q[2];
sx q[2];
rz(-2.1989934) q[2];
rz(1.7132828) q[3];
sx q[3];
rz(-1.2497679) q[3];
sx q[3];
rz(1.9863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1930595) q[0];
sx q[0];
rz(-1.5567224) q[0];
sx q[0];
rz(-1.3514883) q[0];
rz(1.9652741) q[1];
sx q[1];
rz(-0.90570025) q[1];
sx q[1];
rz(-0.66014231) q[1];
rz(1.2930585) q[2];
sx q[2];
rz(-1.3996887) q[2];
sx q[2];
rz(-2.5005093) q[2];
rz(-1.8342212) q[3];
sx q[3];
rz(-1.0527123) q[3];
sx q[3];
rz(0.15409877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
