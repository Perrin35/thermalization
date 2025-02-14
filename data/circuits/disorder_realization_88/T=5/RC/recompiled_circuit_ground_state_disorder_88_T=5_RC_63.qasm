OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.705536) q[0];
sx q[0];
rz(-2.2826865) q[0];
sx q[0];
rz(-0.98281759) q[0];
rz(4.3217826) q[1];
sx q[1];
rz(5.5569841) q[1];
sx q[1];
rz(10.32294) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2959864) q[0];
sx q[0];
rz(-1.6750511) q[0];
sx q[0];
rz(-2.8235675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55805331) q[2];
sx q[2];
rz(-0.61667569) q[2];
sx q[2];
rz(1.0215173) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8232728) q[1];
sx q[1];
rz(-1.987769) q[1];
sx q[1];
rz(0.61064536) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4709365) q[3];
sx q[3];
rz(-1.4263527) q[3];
sx q[3];
rz(0.58510548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1232221) q[2];
sx q[2];
rz(-1.3192588) q[2];
sx q[2];
rz(2.3940864) q[2];
rz(-1.8788762) q[3];
sx q[3];
rz(-1.6044173) q[3];
sx q[3];
rz(-3.1178764) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30775192) q[0];
sx q[0];
rz(-3.096014) q[0];
sx q[0];
rz(-1.0652834) q[0];
rz(-0.65159687) q[1];
sx q[1];
rz(-2.7862796) q[1];
sx q[1];
rz(1.6337055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9414174) q[0];
sx q[0];
rz(-1.2696517) q[0];
sx q[0];
rz(0.99267545) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9266195) q[2];
sx q[2];
rz(-1.1220555) q[2];
sx q[2];
rz(-2.2858278) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.91707857) q[1];
sx q[1];
rz(-1.1076756) q[1];
sx q[1];
rz(-0.45189894) q[1];
x q[2];
rz(-3.1129131) q[3];
sx q[3];
rz(-2.1818871) q[3];
sx q[3];
rz(3.1163586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6642586) q[2];
sx q[2];
rz(-1.4724255) q[2];
sx q[2];
rz(-0.092546917) q[2];
rz(-0.44068286) q[3];
sx q[3];
rz(-0.40658545) q[3];
sx q[3];
rz(-1.996076) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6630845) q[0];
sx q[0];
rz(-2.9496851) q[0];
sx q[0];
rz(-1.9051911) q[0];
rz(-0.23794404) q[1];
sx q[1];
rz(-1.6704208) q[1];
sx q[1];
rz(-2.1740289) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0210812) q[0];
sx q[0];
rz(-0.97645611) q[0];
sx q[0];
rz(-0.8404151) q[0];
rz(-1.9163048) q[2];
sx q[2];
rz(-1.0566525) q[2];
sx q[2];
rz(1.1622857) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0261126) q[1];
sx q[1];
rz(-2.7819595) q[1];
sx q[1];
rz(1.0138649) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0839663) q[3];
sx q[3];
rz(-1.9836243) q[3];
sx q[3];
rz(-2.0183764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6588916) q[2];
sx q[2];
rz(-2.7693558) q[2];
sx q[2];
rz(-2.2230395) q[2];
rz(-2.3679768) q[3];
sx q[3];
rz(-0.88671237) q[3];
sx q[3];
rz(-2.1913967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8755662) q[0];
sx q[0];
rz(-2.4022864) q[0];
sx q[0];
rz(0.29228041) q[0];
rz(-2.267011) q[1];
sx q[1];
rz(-0.62868172) q[1];
sx q[1];
rz(-2.1983217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2349512) q[0];
sx q[0];
rz(-0.72777339) q[0];
sx q[0];
rz(1.1207188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0393911) q[2];
sx q[2];
rz(-1.7649586) q[2];
sx q[2];
rz(-2.5388427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8448461) q[1];
sx q[1];
rz(-1.2592717) q[1];
sx q[1];
rz(-2.2492337) q[1];
x q[2];
rz(-0.99340474) q[3];
sx q[3];
rz(-1.3311989) q[3];
sx q[3];
rz(-1.8519608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5248519) q[2];
sx q[2];
rz(-1.1139694) q[2];
sx q[2];
rz(2.2576766) q[2];
rz(1.135745) q[3];
sx q[3];
rz(-0.92848778) q[3];
sx q[3];
rz(-2.4511724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8166872) q[0];
sx q[0];
rz(-3.0563834) q[0];
sx q[0];
rz(0.19164044) q[0];
rz(-1.2958255) q[1];
sx q[1];
rz(-1.9486267) q[1];
sx q[1];
rz(0.4506909) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5862059) q[0];
sx q[0];
rz(-1.662001) q[0];
sx q[0];
rz(0.93409286) q[0];
rz(-pi) q[1];
rz(1.6411317) q[2];
sx q[2];
rz(-1.8425111) q[2];
sx q[2];
rz(-0.64411417) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1289336) q[1];
sx q[1];
rz(-1.4988572) q[1];
sx q[1];
rz(-0.99314816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8487318) q[3];
sx q[3];
rz(-1.6335051) q[3];
sx q[3];
rz(0.69578275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4638764) q[2];
sx q[2];
rz(-1.2748984) q[2];
sx q[2];
rz(-2.1837168) q[2];
rz(3.0887582) q[3];
sx q[3];
rz(-2.5962679) q[3];
sx q[3];
rz(-2.7132645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4868454) q[0];
sx q[0];
rz(-1.0423648) q[0];
sx q[0];
rz(0.47252193) q[0];
rz(-0.76332244) q[1];
sx q[1];
rz(-0.7205874) q[1];
sx q[1];
rz(0.29019132) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5369072) q[0];
sx q[0];
rz(-1.3504656) q[0];
sx q[0];
rz(-2.4308886) q[0];
x q[1];
rz(0.31748469) q[2];
sx q[2];
rz(-2.0986663) q[2];
sx q[2];
rz(1.1880529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.063263254) q[1];
sx q[1];
rz(-1.6734048) q[1];
sx q[1];
rz(-2.4797477) q[1];
x q[2];
rz(1.3715586) q[3];
sx q[3];
rz(-1.3955978) q[3];
sx q[3];
rz(-0.74149473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7758238) q[2];
sx q[2];
rz(-2.2799728) q[2];
sx q[2];
rz(2.5422868) q[2];
rz(1.4186836) q[3];
sx q[3];
rz(-1.320188) q[3];
sx q[3];
rz(-2.1883709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80239427) q[0];
sx q[0];
rz(-1.9559487) q[0];
sx q[0];
rz(1.3108569) q[0];
rz(1.5015548) q[1];
sx q[1];
rz(-2.7412667) q[1];
sx q[1];
rz(0.60360533) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045403584) q[0];
sx q[0];
rz(-1.9902162) q[0];
sx q[0];
rz(2.9950401) q[0];
rz(-0.55159388) q[2];
sx q[2];
rz(-1.7579798) q[2];
sx q[2];
rz(-2.1647705) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43638602) q[1];
sx q[1];
rz(-1.6387741) q[1];
sx q[1];
rz(1.3297362) q[1];
rz(1.7619753) q[3];
sx q[3];
rz(-2.303225) q[3];
sx q[3];
rz(0.54938706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32460585) q[2];
sx q[2];
rz(-2.2548455) q[2];
sx q[2];
rz(-2.9465607) q[2];
rz(-2.4845691) q[3];
sx q[3];
rz(-1.7200836) q[3];
sx q[3];
rz(0.11369625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8103771) q[0];
sx q[0];
rz(-1.4137784) q[0];
sx q[0];
rz(0.072176607) q[0];
rz(-0.78308925) q[1];
sx q[1];
rz(-0.63242811) q[1];
sx q[1];
rz(-0.91792387) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95436256) q[0];
sx q[0];
rz(-1.7779121) q[0];
sx q[0];
rz(-0.48204473) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8391219) q[2];
sx q[2];
rz(-0.96728281) q[2];
sx q[2];
rz(-0.32121111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.57088045) q[1];
sx q[1];
rz(-1.8547131) q[1];
sx q[1];
rz(1.3234119) q[1];
x q[2];
rz(-1.0467256) q[3];
sx q[3];
rz(-0.48180605) q[3];
sx q[3];
rz(3.1022933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5607295) q[2];
sx q[2];
rz(-1.9274638) q[2];
sx q[2];
rz(1.6403991) q[2];
rz(2.255127) q[3];
sx q[3];
rz(-1.804988) q[3];
sx q[3];
rz(3.0361573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60370541) q[0];
sx q[0];
rz(-2.7925346) q[0];
sx q[0];
rz(-0.50759298) q[0];
rz(0.72788584) q[1];
sx q[1];
rz(-1.6517086) q[1];
sx q[1];
rz(-2.0354039) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17783252) q[0];
sx q[0];
rz(-0.42362693) q[0];
sx q[0];
rz(2.2209441) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7462677) q[2];
sx q[2];
rz(-1.2173436) q[2];
sx q[2];
rz(-2.8025016) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.008279) q[1];
sx q[1];
rz(-2.6402557) q[1];
sx q[1];
rz(-0.8322844) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2415798) q[3];
sx q[3];
rz(-2.4379895) q[3];
sx q[3];
rz(2.7629064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20405208) q[2];
sx q[2];
rz(-2.0378518) q[2];
sx q[2];
rz(0.44720116) q[2];
rz(-1.4372545) q[3];
sx q[3];
rz(-2.3279326) q[3];
sx q[3];
rz(1.6566431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5337885) q[0];
sx q[0];
rz(-2.0543126) q[0];
sx q[0];
rz(0.5150038) q[0];
rz(1.1732514) q[1];
sx q[1];
rz(-1.2898022) q[1];
sx q[1];
rz(2.692093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363578) q[0];
sx q[0];
rz(-0.32669386) q[0];
sx q[0];
rz(1.5004083) q[0];
x q[1];
rz(-0.95566383) q[2];
sx q[2];
rz(-2.0086096) q[2];
sx q[2];
rz(2.4085338) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.24893568) q[1];
sx q[1];
rz(-0.79467862) q[1];
sx q[1];
rz(-0.67130295) q[1];
x q[2];
rz(2.8963666) q[3];
sx q[3];
rz(-0.86638993) q[3];
sx q[3];
rz(-0.77936649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0939193) q[2];
sx q[2];
rz(-0.27457044) q[2];
sx q[2];
rz(0.74335113) q[2];
rz(-2.7798233) q[3];
sx q[3];
rz(-1.0310562) q[3];
sx q[3];
rz(-2.7638392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34306985) q[0];
sx q[0];
rz(-1.1953851) q[0];
sx q[0];
rz(0.60085798) q[0];
rz(2.9877904) q[1];
sx q[1];
rz(-1.3187131) q[1];
sx q[1];
rz(-2.9153894) q[1];
rz(-3.0085887) q[2];
sx q[2];
rz(-1.8756744) q[2];
sx q[2];
rz(-1.6185624) q[2];
rz(-1.947851) q[3];
sx q[3];
rz(-0.98529639) q[3];
sx q[3];
rz(-2.1724971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
