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
rz(1.5481663) q[0];
sx q[0];
rz(-2.6731773) q[0];
sx q[0];
rz(-0.13844891) q[0];
rz(-1.2568714) q[1];
sx q[1];
rz(-1.0705907) q[1];
sx q[1];
rz(0.021477403) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8801418) q[0];
sx q[0];
rz(-2.7316748) q[0];
sx q[0];
rz(-2.1005125) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8201222) q[2];
sx q[2];
rz(-0.73222762) q[2];
sx q[2];
rz(-0.9257462) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5925046) q[1];
sx q[1];
rz(-0.23356423) q[1];
sx q[1];
rz(-0.84850581) q[1];
x q[2];
rz(-1.9324233) q[3];
sx q[3];
rz(-1.0054133) q[3];
sx q[3];
rz(0.20547444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5452177) q[2];
sx q[2];
rz(-0.45520982) q[2];
sx q[2];
rz(1.5111766) q[2];
rz(3.0916072) q[3];
sx q[3];
rz(-1.9129916) q[3];
sx q[3];
rz(-0.25240067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91994691) q[0];
sx q[0];
rz(-1.9939461) q[0];
sx q[0];
rz(-0.037121437) q[0];
rz(-0.014017398) q[1];
sx q[1];
rz(-0.57828301) q[1];
sx q[1];
rz(-2.9055273) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35821298) q[0];
sx q[0];
rz(-1.556408) q[0];
sx q[0];
rz(-3.0667801) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5949941) q[2];
sx q[2];
rz(-2.3476217) q[2];
sx q[2];
rz(-3.1304718) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6746857) q[1];
sx q[1];
rz(-1.0471736) q[1];
sx q[1];
rz(-0.44498131) q[1];
rz(-pi) q[2];
rz(-2.8909952) q[3];
sx q[3];
rz(-1.9204233) q[3];
sx q[3];
rz(-2.6299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0125121) q[2];
sx q[2];
rz(-1.7725638) q[2];
sx q[2];
rz(-0.033626076) q[2];
rz(-2.8513837) q[3];
sx q[3];
rz(-0.84607327) q[3];
sx q[3];
rz(-2.7510551) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6212807) q[0];
sx q[0];
rz(-2.1618167) q[0];
sx q[0];
rz(-2.8424971) q[0];
rz(1.6200804) q[1];
sx q[1];
rz(-3.1144996) q[1];
sx q[1];
rz(-3.1153968) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0157508) q[0];
sx q[0];
rz(-1.5643018) q[0];
sx q[0];
rz(-1.5497213) q[0];
rz(-pi) q[1];
rz(-1.8511591) q[2];
sx q[2];
rz(-1.4329264) q[2];
sx q[2];
rz(-0.029679178) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.700251) q[1];
sx q[1];
rz(-2.1192351) q[1];
sx q[1];
rz(-0.74728181) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.800816) q[3];
sx q[3];
rz(-1.134308) q[3];
sx q[3];
rz(0.088584049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.68411487) q[2];
sx q[2];
rz(-2.699615) q[2];
sx q[2];
rz(-0.6074062) q[2];
rz(-2.7720747) q[3];
sx q[3];
rz(-1.642546) q[3];
sx q[3];
rz(-0.0079689715) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-0.49030855) q[1];
sx q[1];
rz(2.9255548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0433885) q[0];
sx q[0];
rz(-1.3281315) q[0];
sx q[0];
rz(-2.674874) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44682002) q[2];
sx q[2];
rz(-1.6761314) q[2];
sx q[2];
rz(0.0086431816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8204788) q[1];
sx q[1];
rz(-1.063709) q[1];
sx q[1];
rz(0.45407461) q[1];
x q[2];
rz(-1.7407577) q[3];
sx q[3];
rz(-1.483184) q[3];
sx q[3];
rz(2.5690998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.617368) q[2];
sx q[2];
rz(-0.91297954) q[2];
sx q[2];
rz(0.53140223) q[2];
rz(-1.4517387) q[3];
sx q[3];
rz(-0.51133358) q[3];
sx q[3];
rz(-1.0303729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78469974) q[0];
sx q[0];
rz(-2.0578616) q[0];
sx q[0];
rz(0.30237958) q[0];
rz(2.0284292) q[1];
sx q[1];
rz(-1.0013564) q[1];
sx q[1];
rz(-1.0611634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7776124) q[0];
sx q[0];
rz(-1.6355231) q[0];
sx q[0];
rz(-1.6311247) q[0];
x q[1];
rz(1.6689014) q[2];
sx q[2];
rz(-1.7943693) q[2];
sx q[2];
rz(-1.6983216) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94316846) q[1];
sx q[1];
rz(-2.3454002) q[1];
sx q[1];
rz(1.3913586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29814675) q[3];
sx q[3];
rz(-1.3294719) q[3];
sx q[3];
rz(0.69484792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.314996) q[2];
sx q[2];
rz(-1.7642517) q[2];
sx q[2];
rz(0.19485168) q[2];
rz(-2.3991614) q[3];
sx q[3];
rz(-0.73115474) q[3];
sx q[3];
rz(0.2754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.18405296) q[0];
sx q[0];
rz(-2.4803211) q[0];
sx q[0];
rz(-1.1141962) q[0];
rz(-0.30955744) q[1];
sx q[1];
rz(-2.20859) q[1];
sx q[1];
rz(-1.7558487) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615212) q[0];
sx q[0];
rz(-1.7465704) q[0];
sx q[0];
rz(1.2840791) q[0];
rz(-pi) q[1];
rz(1.7739595) q[2];
sx q[2];
rz(-0.51649714) q[2];
sx q[2];
rz(1.1927644) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6760113) q[1];
sx q[1];
rz(-1.3106924) q[1];
sx q[1];
rz(-2.4838402) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65151229) q[3];
sx q[3];
rz(-2.2527472) q[3];
sx q[3];
rz(1.7891974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.679057) q[2];
sx q[2];
rz(-0.44687301) q[2];
sx q[2];
rz(-2.7317969) q[2];
rz(-3.0637686) q[3];
sx q[3];
rz(-1.9255368) q[3];
sx q[3];
rz(0.76243824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78793144) q[0];
sx q[0];
rz(-2.9910112) q[0];
sx q[0];
rz(1.4713564) q[0];
rz(-0.81368601) q[1];
sx q[1];
rz(-2.4206471) q[1];
sx q[1];
rz(2.241316) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6532458) q[0];
sx q[0];
rz(-1.4524967) q[0];
sx q[0];
rz(-1.5100046) q[0];
rz(-pi) q[1];
x q[1];
rz(1.288432) q[2];
sx q[2];
rz(-2.8430364) q[2];
sx q[2];
rz(-2.3207842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1623692) q[1];
sx q[1];
rz(-1.3047393) q[1];
sx q[1];
rz(0.0084657808) q[1];
rz(-pi) q[2];
rz(-2.9081206) q[3];
sx q[3];
rz(-1.6341101) q[3];
sx q[3];
rz(1.7586643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6206996) q[2];
sx q[2];
rz(-2.7405379) q[2];
sx q[2];
rz(2.9996784) q[2];
rz(-0.17169954) q[3];
sx q[3];
rz(-1.9008235) q[3];
sx q[3];
rz(-0.54952526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5080344) q[0];
sx q[0];
rz(-1.1548076) q[0];
sx q[0];
rz(1.5800193) q[0];
rz(0.70478565) q[1];
sx q[1];
rz(-2.2403658) q[1];
sx q[1];
rz(0.27552342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4954105) q[0];
sx q[0];
rz(-2.9998064) q[0];
sx q[0];
rz(1.1205733) q[0];
x q[1];
rz(2.9887096) q[2];
sx q[2];
rz(-2.8934921) q[2];
sx q[2];
rz(2.5923592) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49373473) q[1];
sx q[1];
rz(-1.0714515) q[1];
sx q[1];
rz(-2.8060421) q[1];
x q[2];
rz(0.51361689) q[3];
sx q[3];
rz(-0.73459638) q[3];
sx q[3];
rz(-0.23046042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0813109) q[2];
sx q[2];
rz(-1.5014481) q[2];
sx q[2];
rz(2.8509129) q[2];
rz(2.3447013) q[3];
sx q[3];
rz(-0.47178888) q[3];
sx q[3];
rz(-2.1577788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8271635) q[0];
sx q[0];
rz(-1.6570579) q[0];
sx q[0];
rz(-2.2737801) q[0];
rz(-2.9293291) q[1];
sx q[1];
rz(-0.8808732) q[1];
sx q[1];
rz(-2.8689522) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51202392) q[0];
sx q[0];
rz(-1.4958515) q[0];
sx q[0];
rz(-1.7864947) q[0];
rz(-3.0402203) q[2];
sx q[2];
rz(-0.82216149) q[2];
sx q[2];
rz(-3.0790975) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9524549) q[1];
sx q[1];
rz(-2.5057372) q[1];
sx q[1];
rz(1.5067504) q[1];
rz(-pi) q[2];
rz(2.5288071) q[3];
sx q[3];
rz(-2.6481508) q[3];
sx q[3];
rz(0.32705826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0923882) q[2];
sx q[2];
rz(-0.35085756) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.64677042) q[1];
sx q[1];
rz(1.5470362) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5911212) q[0];
sx q[0];
rz(-2.897399) q[0];
sx q[0];
rz(1.1039735) q[0];
x q[1];
rz(-2.6213264) q[2];
sx q[2];
rz(-2.6718585) q[2];
sx q[2];
rz(0.007291468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4982097) q[1];
sx q[1];
rz(-0.84127448) q[1];
sx q[1];
rz(1.6642844) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.029742777) q[3];
sx q[3];
rz(-1.3689201) q[3];
sx q[3];
rz(-1.0591398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.048162248) q[2];
sx q[2];
rz(-0.54485816) q[2];
sx q[2];
rz(-2.1989934) q[2];
rz(-1.7132828) q[3];
sx q[3];
rz(-1.8918248) q[3];
sx q[3];
rz(-1.1552756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1930595) q[0];
sx q[0];
rz(-1.5848703) q[0];
sx q[0];
rz(1.7901044) q[0];
rz(1.9652741) q[1];
sx q[1];
rz(-0.90570025) q[1];
sx q[1];
rz(-0.66014231) q[1];
rz(-1.8485342) q[2];
sx q[2];
rz(-1.3996887) q[2];
sx q[2];
rz(-2.5005093) q[2];
rz(2.7130933) q[3];
sx q[3];
rz(-2.5658723) q[3];
sx q[3];
rz(0.65278237) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
