OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(3.1403132) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89844184) q[0];
sx q[0];
rz(-1.713322) q[0];
sx q[0];
rz(1.6210763) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7156419) q[2];
sx q[2];
rz(-2.2507239) q[2];
sx q[2];
rz(-1.3198864) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70667627) q[1];
sx q[1];
rz(-1.7129363) q[1];
sx q[1];
rz(0.040277004) q[1];
rz(2.6585456) q[3];
sx q[3];
rz(-0.31937283) q[3];
sx q[3];
rz(-1.7628302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(-1.1958896) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(-1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(-1.8700245) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.7659448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130106) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(-2.9257665) q[0];
x q[1];
rz(-1.8091082) q[2];
sx q[2];
rz(-1.6147699) q[2];
sx q[2];
rz(0.51648742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83982044) q[1];
sx q[1];
rz(-0.41026792) q[1];
sx q[1];
rz(0.46373414) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0166753) q[3];
sx q[3];
rz(-1.2614054) q[3];
sx q[3];
rz(-0.644899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5415444) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(1.906357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2152527) q[0];
sx q[0];
rz(-1.0010166) q[0];
sx q[0];
rz(1.1220054) q[0];
rz(1.3182993) q[2];
sx q[2];
rz(-1.9324979) q[2];
sx q[2];
rz(-2.8351438) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0864799) q[1];
sx q[1];
rz(-1.8243454) q[1];
sx q[1];
rz(-0.26992814) q[1];
rz(2.0577621) q[3];
sx q[3];
rz(-2.1420797) q[3];
sx q[3];
rz(-0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.76413313) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(-0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(2.3377989) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(-0.11985699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940995) q[0];
sx q[0];
rz(-0.96512981) q[0];
sx q[0];
rz(-0.043461965) q[0];
x q[1];
rz(0.36977936) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(-1.0880926) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4704628) q[1];
sx q[1];
rz(-1.5531335) q[1];
sx q[1];
rz(-1.8030241) q[1];
rz(0.5991163) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(-1.8653387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-3.0692549) q[2];
rz(2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(-1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4478093) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(-0.48686349) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.9015076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4362674) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(1.5572085) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8334332) q[2];
sx q[2];
rz(-1.8034168) q[2];
sx q[2];
rz(-2.5831985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33927321) q[1];
sx q[1];
rz(-2.5335651) q[1];
sx q[1];
rz(3.0250711) q[1];
rz(-pi) q[2];
rz(2.6579882) q[3];
sx q[3];
rz(-1.0217474) q[3];
sx q[3];
rz(-0.99398127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.038625) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(0.15771244) q[3];
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
rz(-pi/2) q[3];
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
rz(-2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(-0.99037209) q[0];
rz(0.05274996) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(-1.4809158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4617417) q[0];
sx q[0];
rz(-1.6438419) q[0];
sx q[0];
rz(1.1887656) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.00077498282) q[2];
sx q[2];
rz(-2.0159855) q[2];
sx q[2];
rz(-2.0489401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8112091) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(1.9622383) q[1];
rz(-pi) q[2];
x q[2];
rz(3.092993) q[3];
sx q[3];
rz(-2.2489293) q[3];
sx q[3];
rz(2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(-2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787191) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.4690171) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.8168824) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64667386) q[0];
sx q[0];
rz(-0.12932983) q[0];
sx q[0];
rz(2.9652251) q[0];
rz(-0.5092233) q[2];
sx q[2];
rz(-2.074713) q[2];
sx q[2];
rz(-0.06171209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8285671) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(1.9810956) q[1];
rz(-pi) q[2];
rz(0.91388254) q[3];
sx q[3];
rz(-1.2607288) q[3];
sx q[3];
rz(3.0403746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(0.46052128) q[0];
rz(-3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(-1.2896279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562718) q[0];
sx q[0];
rz(-0.7542146) q[0];
sx q[0];
rz(-0.69099364) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71596594) q[2];
sx q[2];
rz(-2.1858474) q[2];
sx q[2];
rz(-1.3032608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3370812) q[1];
sx q[1];
rz(-1.4599428) q[1];
sx q[1];
rz(-1.7588508) q[1];
rz(0.66289785) q[3];
sx q[3];
rz(-1.7015966) q[3];
sx q[3];
rz(0.42636426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9926247) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(3.0454214) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-2.8073231) q[0];
rz(1.9175247) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-2.8589378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196913) q[0];
sx q[0];
rz(-2.0265409) q[0];
sx q[0];
rz(2.4569608) q[0];
rz(-pi) q[1];
rz(-0.33090584) q[2];
sx q[2];
rz(-1.3661642) q[2];
sx q[2];
rz(-0.76588878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5987451) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(-2.7302242) q[1];
rz(1.4275527) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(-2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(0.7412509) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(0.50869554) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(0.68181109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012321) q[0];
sx q[0];
rz(-1.0245748) q[0];
sx q[0];
rz(-2.9451314) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3186458) q[2];
sx q[2];
rz(-2.1087077) q[2];
sx q[2];
rz(-0.049494628) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0360003) q[1];
sx q[1];
rz(-0.75165527) q[1];
sx q[1];
rz(-2.4113301) q[1];
rz(-0.52313157) q[3];
sx q[3];
rz(-1.0732068) q[3];
sx q[3];
rz(2.350972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(-2.8005023) q[2];
rz(1.0836481) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(-2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-0.63275679) q[2];
sx q[2];
rz(-2.3709595) q[2];
sx q[2];
rz(-0.79216935) q[2];
rz(0.46799216) q[3];
sx q[3];
rz(-0.42588009) q[3];
sx q[3];
rz(1.7298021) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
