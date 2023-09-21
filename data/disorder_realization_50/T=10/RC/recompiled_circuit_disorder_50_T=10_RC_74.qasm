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
rz(-0.0012794415) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476386) q[0];
sx q[0];
rz(-1.6205661) q[0];
sx q[0];
rz(0.14270356) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84471976) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(-2.6127882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1577658) q[1];
sx q[1];
rz(-0.14769927) q[1];
sx q[1];
rz(-1.8450792) q[1];
rz(0.28489057) q[3];
sx q[3];
rz(-1.4244411) q[3];
sx q[3];
rz(-0.26998664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.8700245) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(1.3756479) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401706) q[0];
sx q[0];
rz(-1.4835843) q[0];
sx q[0];
rz(2.7313822) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3865115) q[2];
sx q[2];
rz(-0.24225907) q[2];
sx q[2];
rz(2.266303) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.34054204) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(-1.3786475) q[1];
rz(-pi) q[2];
rz(-2.1167011) q[3];
sx q[3];
rz(-0.62666577) q[3];
sx q[3];
rz(0.46862593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-2.6518872) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.906357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2152527) q[0];
sx q[0];
rz(-1.0010166) q[0];
sx q[0];
rz(2.0195872) q[0];
x q[1];
rz(2.5580102) q[2];
sx q[2];
rz(-2.7036813) q[2];
sx q[2];
rz(-2.8180518) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5879844) q[1];
sx q[1];
rz(-1.3097035) q[1];
sx q[1];
rz(1.3081461) q[1];
rz(2.5127605) q[3];
sx q[3];
rz(-1.1662081) q[3];
sx q[3];
rz(1.7189327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.7049449) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-0.11985699) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0474931) q[0];
sx q[0];
rz(-0.96512981) q[0];
sx q[0];
rz(0.043461965) q[0];
x q[1];
rz(0.36977936) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(2.0535) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2377492) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(-0.018149908) q[1];
x q[2];
rz(0.71505349) q[3];
sx q[3];
rz(-2.4063769) q[3];
sx q[3];
rz(-0.27763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.9015076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70532521) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(1.5572085) q[0];
x q[1];
rz(1.3081595) q[2];
sx q[2];
rz(-1.3381759) q[2];
sx q[2];
rz(2.5831985) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.005849) q[1];
sx q[1];
rz(-1.6372576) q[1];
sx q[1];
rz(-0.60484109) q[1];
x q[2];
rz(-2.6579882) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(-0.99398127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.038625) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-2.4385578) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(-0.99037209) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.4809158) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2213343) q[0];
sx q[0];
rz(-1.9517559) q[0];
sx q[0];
rz(3.0628946) q[0];
rz(2.0159857) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(0.47847754) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33038352) q[1];
sx q[1];
rz(-2.5685709) q[1];
sx q[1];
rz(1.9622383) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89208608) q[3];
sx q[3];
rz(-1.6086372) q[3];
sx q[3];
rz(-0.91176646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(-0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(1.0143657) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(1.8168824) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0425456) q[0];
sx q[0];
rz(-1.5481661) q[0];
sx q[0];
rz(-0.12734539) q[0];
x q[1];
rz(-2.2947646) q[2];
sx q[2];
rz(-0.70038751) q[2];
sx q[2];
rz(2.2221153) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99609375) q[1];
sx q[1];
rz(-2.5719574) q[1];
sx q[1];
rz(0.74665229) q[1];
rz(0.38447325) q[3];
sx q[3];
rz(-2.1914346) q[3];
sx q[3];
rz(-1.7006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5801195) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(-0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.2896279) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136918) q[0];
sx q[0];
rz(-2.0223589) q[0];
sx q[0];
rz(-0.62664647) q[0];
rz(-pi) q[1];
rz(-0.71596594) q[2];
sx q[2];
rz(-2.1858474) q[2];
sx q[2];
rz(1.3032608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21266567) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(0.11282632) q[1];
rz(1.4054221) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(-2.0588493) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(1.224068) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7455505) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(1.0066443) q[0];
x q[1];
rz(-0.56844175) q[2];
sx q[2];
rz(-0.38707765) q[2];
sx q[2];
rz(-2.8708411) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54284755) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(-2.7302242) q[1];
rz(1.4275527) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(0.43977246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21215542) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(2.4003417) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(-3.0264061) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(-0.68181109) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17446974) q[0];
sx q[0];
rz(-0.57708626) q[0];
sx q[0];
rz(1.2601) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3186458) q[2];
sx q[2];
rz(-1.032885) q[2];
sx q[2];
rz(-0.049494628) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0444813) q[1];
sx q[1];
rz(-1.0978571) q[1];
sx q[1];
rz(2.5333235) q[1];
rz(-pi) q[2];
rz(-0.52313157) q[3];
sx q[3];
rz(-2.0683859) q[3];
sx q[3];
rz(-2.350972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(2.8005023) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(-0.17102374) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(2.534261) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(0.66424673) q[2];
sx q[2];
rz(-1.1462117) q[2];
sx q[2];
rz(0.29427634) q[2];
rz(1.7726462) q[3];
sx q[3];
rz(-1.9484083) q[3];
sx q[3];
rz(1.2231135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
