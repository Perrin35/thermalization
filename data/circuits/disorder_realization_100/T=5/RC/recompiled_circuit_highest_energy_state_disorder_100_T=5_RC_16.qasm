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
rz(-1.1779397) q[0];
sx q[0];
rz(-3.0484634) q[0];
sx q[0];
rz(-0.73214522) q[0];
rz(1.572345) q[1];
sx q[1];
rz(-2.2557926) q[1];
sx q[1];
rz(-0.43876171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9222835) q[0];
sx q[0];
rz(-1.6262591) q[0];
sx q[0];
rz(1.903141) q[0];
rz(-pi) q[1];
rz(1.8541502) q[2];
sx q[2];
rz(-1.8210951) q[2];
sx q[2];
rz(1.3586501) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6631515) q[1];
sx q[1];
rz(-1.3125129) q[1];
sx q[1];
rz(-1.5228935) q[1];
rz(-pi) q[2];
rz(0.097919959) q[3];
sx q[3];
rz(-2.681308) q[3];
sx q[3];
rz(-1.8079881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0977122) q[2];
sx q[2];
rz(-0.040412929) q[2];
sx q[2];
rz(-0.91862339) q[2];
rz(2.0888603) q[3];
sx q[3];
rz(-0.91679263) q[3];
sx q[3];
rz(2.2241433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87616462) q[0];
sx q[0];
rz(-0.81887236) q[0];
sx q[0];
rz(-0.87232605) q[0];
rz(2.1288952) q[1];
sx q[1];
rz(-1.5112977) q[1];
sx q[1];
rz(2.8025119) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10386183) q[0];
sx q[0];
rz(-1.6997689) q[0];
sx q[0];
rz(-0.43855389) q[0];
rz(-pi) q[1];
rz(0.85810272) q[2];
sx q[2];
rz(-2.4913308) q[2];
sx q[2];
rz(-1.5257143) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87552947) q[1];
sx q[1];
rz(-1.5782981) q[1];
sx q[1];
rz(-0.041635978) q[1];
rz(3.0172164) q[3];
sx q[3];
rz(-1.892883) q[3];
sx q[3];
rz(-2.6895804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8104711) q[2];
sx q[2];
rz(-1.9942185) q[2];
sx q[2];
rz(-2.5362711) q[2];
rz(-2.807054) q[3];
sx q[3];
rz(-0.3722705) q[3];
sx q[3];
rz(-3.065006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2452068) q[0];
sx q[0];
rz(-2.4672282) q[0];
sx q[0];
rz(1.6687923) q[0];
rz(0.32707602) q[1];
sx q[1];
rz(-1.3027124) q[1];
sx q[1];
rz(0.27007857) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9327346) q[0];
sx q[0];
rz(-1.3737824) q[0];
sx q[0];
rz(-0.52893649) q[0];
rz(-pi) q[1];
rz(3.0659778) q[2];
sx q[2];
rz(-1.9984198) q[2];
sx q[2];
rz(-0.84578925) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1535608) q[1];
sx q[1];
rz(-2.2940293) q[1];
sx q[1];
rz(2.5392697) q[1];
rz(1.5746501) q[3];
sx q[3];
rz(-1.4906297) q[3];
sx q[3];
rz(-1.7146669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0373056) q[2];
sx q[2];
rz(-1.9992) q[2];
sx q[2];
rz(1.2874862) q[2];
rz(-1.9262975) q[3];
sx q[3];
rz(-1.3853962) q[3];
sx q[3];
rz(-0.17770411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52008587) q[0];
sx q[0];
rz(-0.44317133) q[0];
sx q[0];
rz(2.8218414) q[0];
rz(2.7301835) q[1];
sx q[1];
rz(-1.5450059) q[1];
sx q[1];
rz(3.06126) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7112996) q[0];
sx q[0];
rz(-1.7691433) q[0];
sx q[0];
rz(-1.2567149) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3351203) q[2];
sx q[2];
rz(-1.4676759) q[2];
sx q[2];
rz(1.1191302) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87819609) q[1];
sx q[1];
rz(-2.2913763) q[1];
sx q[1];
rz(-2.7769068) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5307054) q[3];
sx q[3];
rz(-1.3502933) q[3];
sx q[3];
rz(1.0569416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2692928) q[2];
sx q[2];
rz(-0.23739693) q[2];
sx q[2];
rz(0.044895127) q[2];
rz(-0.53523713) q[3];
sx q[3];
rz(-1.5072631) q[3];
sx q[3];
rz(3.0285192) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.991268) q[0];
sx q[0];
rz(-2.4327705) q[0];
sx q[0];
rz(2.2659361) q[0];
rz(-0.21370299) q[1];
sx q[1];
rz(-2.6468266) q[1];
sx q[1];
rz(-1.5692086) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039883651) q[0];
sx q[0];
rz(-1.3613762) q[0];
sx q[0];
rz(1.8906192) q[0];
rz(-pi) q[1];
rz(-1.8000697) q[2];
sx q[2];
rz(-1.9976282) q[2];
sx q[2];
rz(-0.62139702) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9006648) q[1];
sx q[1];
rz(-0.87605175) q[1];
sx q[1];
rz(1.7068295) q[1];
x q[2];
rz(2.8205958) q[3];
sx q[3];
rz(-1.6602274) q[3];
sx q[3];
rz(3.1232219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.956942) q[2];
sx q[2];
rz(-2.9968379) q[2];
sx q[2];
rz(1.0682028) q[2];
rz(0.60643658) q[3];
sx q[3];
rz(-1.9304099) q[3];
sx q[3];
rz(0.010312168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40474263) q[0];
sx q[0];
rz(-0.24979845) q[0];
sx q[0];
rz(-1.3909719) q[0];
rz(-0.5237611) q[1];
sx q[1];
rz(-0.57558376) q[1];
sx q[1];
rz(1.9578804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8075631) q[0];
sx q[0];
rz(-1.5153756) q[0];
sx q[0];
rz(-1.8593273) q[0];
x q[1];
rz(-2.1864088) q[2];
sx q[2];
rz(-3.0490626) q[2];
sx q[2];
rz(-1.1902155) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5038915) q[1];
sx q[1];
rz(-0.63234416) q[1];
sx q[1];
rz(1.8296627) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2565485) q[3];
sx q[3];
rz(-2.9687299) q[3];
sx q[3];
rz(-2.1751753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4411053) q[2];
sx q[2];
rz(-2.5429071) q[2];
sx q[2];
rz(-1.0318476) q[2];
rz(1.4619689) q[3];
sx q[3];
rz(-2.4470191) q[3];
sx q[3];
rz(1.7803378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.66330355) q[0];
sx q[0];
rz(-2.7155868) q[0];
sx q[0];
rz(-1.4248832) q[0];
rz(2.5583963) q[1];
sx q[1];
rz(-1.2251264) q[1];
sx q[1];
rz(-0.057417631) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1844139) q[0];
sx q[0];
rz(-0.8033565) q[0];
sx q[0];
rz(0.0048385421) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0591828) q[2];
sx q[2];
rz(-1.0213189) q[2];
sx q[2];
rz(-1.8715931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.098328248) q[1];
sx q[1];
rz(-2.4803313) q[1];
sx q[1];
rz(0.28820451) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.613861) q[3];
sx q[3];
rz(-0.70042983) q[3];
sx q[3];
rz(-0.34175261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3164369) q[2];
sx q[2];
rz(-1.2978483) q[2];
sx q[2];
rz(-1.8334897) q[2];
rz(-1.8209275) q[3];
sx q[3];
rz(-0.84669176) q[3];
sx q[3];
rz(2.2804885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84701455) q[0];
sx q[0];
rz(-0.68462831) q[0];
sx q[0];
rz(1.8735877) q[0];
rz(-1.1216724) q[1];
sx q[1];
rz(-1.2753762) q[1];
sx q[1];
rz(0.65006382) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0065816) q[0];
sx q[0];
rz(-1.0502292) q[0];
sx q[0];
rz(2.6755946) q[0];
x q[1];
rz(1.0380657) q[2];
sx q[2];
rz(-2.1830171) q[2];
sx q[2];
rz(2.39447) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6942581) q[1];
sx q[1];
rz(-1.9617394) q[1];
sx q[1];
rz(2.5503134) q[1];
rz(-pi) q[2];
rz(-0.88553377) q[3];
sx q[3];
rz(-0.69708744) q[3];
sx q[3];
rz(-3.0892885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0012297) q[2];
sx q[2];
rz(-1.4823806) q[2];
sx q[2];
rz(0.39230997) q[2];
rz(0.20491925) q[3];
sx q[3];
rz(-2.6409918) q[3];
sx q[3];
rz(0.71995455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1392764) q[0];
sx q[0];
rz(-1.5520232) q[0];
sx q[0];
rz(-1.4578777) q[0];
rz(-0.32683364) q[1];
sx q[1];
rz(-1.9953597) q[1];
sx q[1];
rz(2.0976417) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21699587) q[0];
sx q[0];
rz(-1.1490273) q[0];
sx q[0];
rz(-0.20713465) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4307666) q[2];
sx q[2];
rz(-2.0256038) q[2];
sx q[2];
rz(2.8190913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9324016) q[1];
sx q[1];
rz(-1.4889917) q[1];
sx q[1];
rz(1.1751168) q[1];
rz(2.2208728) q[3];
sx q[3];
rz(-1.8333922) q[3];
sx q[3];
rz(0.85185862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8331208) q[2];
sx q[2];
rz(-2.1127508) q[2];
sx q[2];
rz(2.2157045) q[2];
rz(2.0738257) q[3];
sx q[3];
rz(-2.0177149) q[3];
sx q[3];
rz(1.3656176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2866216) q[0];
sx q[0];
rz(-1.4737031) q[0];
sx q[0];
rz(2.5545252) q[0];
rz(0.58285561) q[1];
sx q[1];
rz(-0.28935495) q[1];
sx q[1];
rz(-0.48097441) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0547377) q[0];
sx q[0];
rz(-1.361892) q[0];
sx q[0];
rz(0.41685327) q[0];
rz(-pi) q[1];
rz(-1.0985259) q[2];
sx q[2];
rz(-2.741217) q[2];
sx q[2];
rz(0.085299678) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1470589) q[1];
sx q[1];
rz(-1.9939341) q[1];
sx q[1];
rz(2.5446241) q[1];
x q[2];
rz(-1.770788) q[3];
sx q[3];
rz(-1.1781577) q[3];
sx q[3];
rz(2.5225366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9440072) q[2];
sx q[2];
rz(-0.53637594) q[2];
sx q[2];
rz(1.3295004) q[2];
rz(-2.3594989) q[3];
sx q[3];
rz(-0.32556459) q[3];
sx q[3];
rz(-1.1480924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6999577) q[0];
sx q[0];
rz(-1.4186207) q[0];
sx q[0];
rz(1.6765208) q[0];
rz(-1.5326473) q[1];
sx q[1];
rz(-1.1679222) q[1];
sx q[1];
rz(2.4293778) q[1];
rz(2.2198792) q[2];
sx q[2];
rz(-0.94759181) q[2];
sx q[2];
rz(0.0061719051) q[2];
rz(-2.8529975) q[3];
sx q[3];
rz(-2.6251818) q[3];
sx q[3];
rz(2.8040721) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
