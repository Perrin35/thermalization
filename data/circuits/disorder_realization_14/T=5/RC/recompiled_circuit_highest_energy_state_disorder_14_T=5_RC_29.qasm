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
rz(1.8847213) q[1];
sx q[1];
rz(4.2121834) q[1];
sx q[1];
rz(9.4033006) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8801418) q[0];
sx q[0];
rz(-0.40991789) q[0];
sx q[0];
rz(2.1005125) q[0];
rz(-2.1523066) q[2];
sx q[2];
rz(-2.044319) q[2];
sx q[2];
rz(-3.1030637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8565377) q[1];
sx q[1];
rz(-1.3962588) q[1];
sx q[1];
rz(2.9855897) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59604597) q[3];
sx q[3];
rz(-1.2674244) q[3];
sx q[3];
rz(1.16538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5452177) q[2];
sx q[2];
rz(-0.45520982) q[2];
sx q[2];
rz(-1.5111766) q[2];
rz(-3.0916072) q[3];
sx q[3];
rz(-1.9129916) q[3];
sx q[3];
rz(0.25240067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2216457) q[0];
sx q[0];
rz(-1.1476465) q[0];
sx q[0];
rz(-3.1044712) q[0];
rz(-3.1275753) q[1];
sx q[1];
rz(-2.5633096) q[1];
sx q[1];
rz(-2.9055273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2136617) q[0];
sx q[0];
rz(-1.4959916) q[0];
sx q[0];
rz(1.585225) q[0];
x q[1];
rz(2.5949941) q[2];
sx q[2];
rz(-0.79397092) q[2];
sx q[2];
rz(-3.1304718) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8467847) q[1];
sx q[1];
rz(-0.67344147) q[1];
sx q[1];
rz(-0.93017857) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1680683) q[3];
sx q[3];
rz(-0.42713886) q[3];
sx q[3];
rz(-1.1534302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1290805) q[2];
sx q[2];
rz(-1.3690288) q[2];
sx q[2];
rz(0.033626076) q[2];
rz(-0.29020894) q[3];
sx q[3];
rz(-0.84607327) q[3];
sx q[3];
rz(-0.39053759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6212807) q[0];
sx q[0];
rz(-0.97977591) q[0];
sx q[0];
rz(-0.2990956) q[0];
rz(1.5215123) q[1];
sx q[1];
rz(-0.027093096) q[1];
sx q[1];
rz(0.026195899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25616249) q[0];
sx q[0];
rz(-3.1195398) q[0];
sx q[0];
rz(-1.2718448) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2904335) q[2];
sx q[2];
rz(-1.4329264) q[2];
sx q[2];
rz(-3.1119135) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67945665) q[1];
sx q[1];
rz(-0.95210451) q[1];
sx q[1];
rz(0.87631823) q[1];
rz(-pi) q[2];
rz(-2.6947909) q[3];
sx q[3];
rz(-1.7789156) q[3];
sx q[3];
rz(1.3835386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68411487) q[2];
sx q[2];
rz(-2.699615) q[2];
sx q[2];
rz(0.6074062) q[2];
rz(2.7720747) q[3];
sx q[3];
rz(-1.4990467) q[3];
sx q[3];
rz(3.1336237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5615416) q[0];
sx q[0];
rz(-2.9396368) q[0];
sx q[0];
rz(-1.1308905) q[0];
rz(-0.35218969) q[1];
sx q[1];
rz(-2.6512841) q[1];
sx q[1];
rz(-2.9255548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7894831) q[0];
sx q[0];
rz(-1.1187859) q[0];
sx q[0];
rz(1.841196) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6874996) q[2];
sx q[2];
rz(-1.1266303) q[2];
sx q[2];
rz(-1.5291052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32111383) q[1];
sx q[1];
rz(-2.0778837) q[1];
sx q[1];
rz(2.687518) q[1];
x q[2];
rz(1.0918255) q[3];
sx q[3];
rz(-0.19102016) q[3];
sx q[3];
rz(0.52680069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5242247) q[2];
sx q[2];
rz(-2.2286131) q[2];
sx q[2];
rz(-2.6101904) q[2];
rz(1.6898539) q[3];
sx q[3];
rz(-0.51133358) q[3];
sx q[3];
rz(-1.0303729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78469974) q[0];
sx q[0];
rz(-2.0578616) q[0];
sx q[0];
rz(2.8392131) q[0];
rz(-1.1131635) q[1];
sx q[1];
rz(-1.0013564) q[1];
sx q[1];
rz(2.0804292) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.115175) q[0];
sx q[0];
rz(-0.088453293) q[0];
sx q[0];
rz(-2.3923516) q[0];
x q[1];
rz(1.6689014) q[2];
sx q[2];
rz(-1.3472234) q[2];
sx q[2];
rz(1.6983216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1969152) q[1];
sx q[1];
rz(-0.79087559) q[1];
sx q[1];
rz(-2.9612035) q[1];
x q[2];
rz(0.29814675) q[3];
sx q[3];
rz(-1.8121208) q[3];
sx q[3];
rz(-0.69484792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82659668) q[2];
sx q[2];
rz(-1.3773409) q[2];
sx q[2];
rz(-0.19485168) q[2];
rz(2.3991614) q[3];
sx q[3];
rz(-2.4104379) q[3];
sx q[3];
rz(0.2754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18405296) q[0];
sx q[0];
rz(-2.4803211) q[0];
sx q[0];
rz(-2.0273965) q[0];
rz(2.8320352) q[1];
sx q[1];
rz(-2.20859) q[1];
sx q[1];
rz(1.385744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64460574) q[0];
sx q[0];
rz(-0.33503767) q[0];
sx q[0];
rz(-2.1315639) q[0];
rz(-pi) q[1];
rz(-1.7739595) q[2];
sx q[2];
rz(-2.6250955) q[2];
sx q[2];
rz(1.1927644) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6760113) q[1];
sx q[1];
rz(-1.8309002) q[1];
sx q[1];
rz(2.4838402) q[1];
rz(-pi) q[2];
rz(0.65151229) q[3];
sx q[3];
rz(-2.2527472) q[3];
sx q[3];
rz(-1.7891974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4625357) q[2];
sx q[2];
rz(-2.6947196) q[2];
sx q[2];
rz(2.7317969) q[2];
rz(3.0637686) q[3];
sx q[3];
rz(-1.9255368) q[3];
sx q[3];
rz(-0.76243824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78793144) q[0];
sx q[0];
rz(-2.9910112) q[0];
sx q[0];
rz(-1.6702363) q[0];
rz(-0.81368601) q[1];
sx q[1];
rz(-2.4206471) q[1];
sx q[1];
rz(2.241316) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532458) q[0];
sx q[0];
rz(-1.6890959) q[0];
sx q[0];
rz(-1.6315881) q[0];
rz(-pi) q[1];
rz(1.8581819) q[2];
sx q[2];
rz(-1.6528439) q[2];
sx q[2];
rz(1.0204741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.011411) q[1];
sx q[1];
rz(-2.8754042) q[1];
sx q[1];
rz(-1.5397416) q[1];
rz(-pi) q[2];
rz(-0.23347207) q[3];
sx q[3];
rz(-1.6341101) q[3];
sx q[3];
rz(-1.7586643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5208931) q[2];
sx q[2];
rz(-2.7405379) q[2];
sx q[2];
rz(-0.14191423) q[2];
rz(-2.9698931) q[3];
sx q[3];
rz(-1.9008235) q[3];
sx q[3];
rz(0.54952526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5080344) q[0];
sx q[0];
rz(-1.9867851) q[0];
sx q[0];
rz(-1.5800193) q[0];
rz(-0.70478565) q[1];
sx q[1];
rz(-2.2403658) q[1];
sx q[1];
rz(2.8660692) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0412236) q[0];
sx q[0];
rz(-1.4432206) q[0];
sx q[0];
rz(-0.062037717) q[0];
rz(-pi) q[1];
rz(1.5322379) q[2];
sx q[2];
rz(-1.3256494) q[2];
sx q[2];
rz(2.4347255) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91159484) q[1];
sx q[1];
rz(-1.2775289) q[1];
sx q[1];
rz(-2.0946343) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6279758) q[3];
sx q[3];
rz(-0.73459638) q[3];
sx q[3];
rz(0.23046042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0813109) q[2];
sx q[2];
rz(-1.6401446) q[2];
sx q[2];
rz(2.8509129) q[2];
rz(0.79689133) q[3];
sx q[3];
rz(-2.6698038) q[3];
sx q[3];
rz(-2.1577788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8271635) q[0];
sx q[0];
rz(-1.4845347) q[0];
sx q[0];
rz(-0.86781251) q[0];
rz(2.9293291) q[1];
sx q[1];
rz(-2.2607195) q[1];
sx q[1];
rz(-2.8689522) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6295687) q[0];
sx q[0];
rz(-1.6457412) q[0];
sx q[0];
rz(-1.7864947) q[0];
rz(-pi) q[1];
rz(-2.3219982) q[2];
sx q[2];
rz(-1.4965881) q[2];
sx q[2];
rz(-1.5774262) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8729108) q[1];
sx q[1];
rz(-2.2051393) q[1];
sx q[1];
rz(0.047206248) q[1];
x q[2];
rz(2.7271184) q[3];
sx q[3];
rz(-1.8467086) q[3];
sx q[3];
rz(1.3434354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.049204443) q[2];
sx q[2];
rz(-2.7907351) q[2];
sx q[2];
rz(-1.1445047) q[2];
rz(0.62659621) q[3];
sx q[3];
rz(-0.75855362) q[3];
sx q[3];
rz(0.315061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7476244) q[0];
sx q[0];
rz(-2.4691041) q[0];
sx q[0];
rz(-2.5482063) q[0];
rz(0.74140948) q[1];
sx q[1];
rz(-2.4948222) q[1];
sx q[1];
rz(-1.5470362) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5911212) q[0];
sx q[0];
rz(-0.24419366) q[0];
sx q[0];
rz(-2.0376192) q[0];
rz(1.3236078) q[2];
sx q[2];
rz(-1.9744248) q[2];
sx q[2];
rz(0.57838049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.64338291) q[1];
sx q[1];
rz(-2.3003182) q[1];
sx q[1];
rz(1.4773083) q[1];
rz(1.4265028) q[3];
sx q[3];
rz(-2.9375667) q[3];
sx q[3];
rz(-0.91183364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0934304) q[2];
sx q[2];
rz(-0.54485816) q[2];
sx q[2];
rz(2.1989934) q[2];
rz(1.4283098) q[3];
sx q[3];
rz(-1.2497679) q[3];
sx q[3];
rz(-1.9863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1930595) q[0];
sx q[0];
rz(-1.5848703) q[0];
sx q[0];
rz(1.7901044) q[0];
rz(1.1763186) q[1];
sx q[1];
rz(-2.2358924) q[1];
sx q[1];
rz(2.4814503) q[1];
rz(2.1331486) q[2];
sx q[2];
rz(-2.8165419) q[2];
sx q[2];
rz(-1.4680924) q[2];
rz(2.6082718) q[3];
sx q[3];
rz(-1.7989895) q[3];
sx q[3];
rz(-1.2839272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
