OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.39419898) q[0];
sx q[0];
rz(-0.49180254) q[0];
sx q[0];
rz(-2.9536182) q[0];
rz(-1.1176874) q[1];
sx q[1];
rz(-1.517065) q[1];
sx q[1];
rz(2.7741073) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7029019) q[0];
sx q[0];
rz(-2.5950948) q[0];
sx q[0];
rz(1.7706857) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3221402) q[2];
sx q[2];
rz(-0.50422943) q[2];
sx q[2];
rz(2.3766975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6403113) q[1];
sx q[1];
rz(-0.2959364) q[1];
sx q[1];
rz(2.179115) q[1];
rz(-pi) q[2];
rz(0.46842694) q[3];
sx q[3];
rz(-0.37852968) q[3];
sx q[3];
rz(0.98640501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1774896) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(-2.5906079) q[2];
rz(-1.3059113) q[3];
sx q[3];
rz(-1.4923613) q[3];
sx q[3];
rz(-1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(0.42981237) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(-2.205251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124356) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(-2.9108414) q[0];
rz(1.1510552) q[2];
sx q[2];
rz(-0.83101666) q[2];
sx q[2];
rz(1.8120399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97570005) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(2.9794934) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.02336054) q[3];
sx q[3];
rz(-1.0682032) q[3];
sx q[3];
rz(-2.392829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3669746) q[2];
sx q[2];
rz(-0.32704157) q[2];
sx q[2];
rz(2.7152087) q[2];
rz(-1.2373699) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.8957829) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(2.202503) q[0];
rz(-0.89871961) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(2.5476707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9868601) q[0];
sx q[0];
rz(-1.8911456) q[0];
sx q[0];
rz(2.8264168) q[0];
rz(-pi) q[1];
rz(-0.8823231) q[2];
sx q[2];
rz(-1.6154628) q[2];
sx q[2];
rz(0.67827144) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18332874) q[1];
sx q[1];
rz(-1.3576344) q[1];
sx q[1];
rz(-1.1413241) q[1];
rz(-pi) q[2];
rz(2.5251758) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(2.4544231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.64017355) q[2];
sx q[2];
rz(-2.4121425) q[2];
sx q[2];
rz(-1.4397941) q[2];
rz(-0.38763186) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3751635) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(2.373383) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(-0.75685135) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605646) q[0];
sx q[0];
rz(-1.2111944) q[0];
sx q[0];
rz(1.9268376) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0365305) q[2];
sx q[2];
rz(-1.5823936) q[2];
sx q[2];
rz(-0.84601814) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0087352) q[1];
sx q[1];
rz(-2.1454304) q[1];
sx q[1];
rz(-0.43032129) q[1];
rz(2.0714949) q[3];
sx q[3];
rz(-0.90161937) q[3];
sx q[3];
rz(-2.5239528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42671529) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(-1.654401) q[2];
rz(-2.5590844) q[3];
sx q[3];
rz(-1.094386) q[3];
sx q[3];
rz(0.55707651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-1.0852381) q[0];
sx q[0];
rz(-0.75772444) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-2.2133591) q[1];
sx q[1];
rz(2.0910738) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8331063) q[0];
sx q[0];
rz(-1.8252488) q[0];
sx q[0];
rz(1.893505) q[0];
x q[1];
rz(-0.96111416) q[2];
sx q[2];
rz(-0.69176199) q[2];
sx q[2];
rz(1.8209396) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1861021) q[1];
sx q[1];
rz(-1.0483861) q[1];
sx q[1];
rz(-0.91764692) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53253048) q[3];
sx q[3];
rz(-0.88056394) q[3];
sx q[3];
rz(-2.0257476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.918255) q[2];
sx q[2];
rz(-2.788322) q[2];
sx q[2];
rz(-2.5081432) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69960064) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-2.8884086) q[0];
rz(-1.6075915) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(-1.4621428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77308649) q[0];
sx q[0];
rz(-2.9162772) q[0];
sx q[0];
rz(0.36264514) q[0];
rz(-pi) q[1];
rz(-0.24121933) q[2];
sx q[2];
rz(-0.93572817) q[2];
sx q[2];
rz(-1.5347753) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0837005) q[1];
sx q[1];
rz(-1.18827) q[1];
sx q[1];
rz(1.3871357) q[1];
rz(-2.8263894) q[3];
sx q[3];
rz(-1.2425555) q[3];
sx q[3];
rz(0.70716508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9399461) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(-2.3366826) q[2];
rz(1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(-2.2139363) q[0];
rz(2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(-2.129508) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9575858) q[0];
sx q[0];
rz(-2.036096) q[0];
sx q[0];
rz(-0.9853978) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38883932) q[2];
sx q[2];
rz(-1.4423587) q[2];
sx q[2];
rz(1.625979) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5092897) q[1];
sx q[1];
rz(-2.1142694) q[1];
sx q[1];
rz(-1.6093045) q[1];
x q[2];
rz(-0.55862553) q[3];
sx q[3];
rz(-1.5010251) q[3];
sx q[3];
rz(-2.74673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(0.42993316) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(2.608192) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6124509) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(2.9274143) q[0];
rz(-2.0902436) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(2.8578551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9203628) q[0];
sx q[0];
rz(-1.8717248) q[0];
sx q[0];
rz(-1.6065341) q[0];
rz(0.19212888) q[2];
sx q[2];
rz(-1.2461975) q[2];
sx q[2];
rz(2.0786376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9323862) q[1];
sx q[1];
rz(-1.7665518) q[1];
sx q[1];
rz(-2.6372361) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.489336) q[3];
sx q[3];
rz(-2.1144923) q[3];
sx q[3];
rz(-0.039747681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(-1.4452176) q[2];
rz(1.5444267) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(-2.8022695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504869) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(-0.069256393) q[0];
rz(-1.6537369) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(-1.5725296) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5208961) q[0];
sx q[0];
rz(-2.3224152) q[0];
sx q[0];
rz(0.92185123) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32585085) q[2];
sx q[2];
rz(-1.648765) q[2];
sx q[2];
rz(-0.98380145) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.442765) q[1];
sx q[1];
rz(-2.8932533) q[1];
sx q[1];
rz(0.34766867) q[1];
rz(-0.2089573) q[3];
sx q[3];
rz(-2.2088802) q[3];
sx q[3];
rz(-2.0861422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1853603) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(3.0017079) q[2];
rz(0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-2.4998253) q[0];
rz(1.2311252) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(0.26783255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56775996) q[0];
sx q[0];
rz(-1.1288252) q[0];
sx q[0];
rz(2.3964336) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1688813) q[2];
sx q[2];
rz(-2.9209666) q[2];
sx q[2];
rz(2.0711183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2512868) q[1];
sx q[1];
rz(-2.7000513) q[1];
sx q[1];
rz(2.4114386) q[1];
rz(-pi) q[2];
rz(-2.0249428) q[3];
sx q[3];
rz(-0.59026679) q[3];
sx q[3];
rz(0.86436194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3315167) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(2.8752575) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(-0.77990445) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(2.4895888) q[2];
sx q[2];
rz(-2.3163788) q[2];
sx q[2];
rz(-0.083995081) q[2];
rz(2.8779463) q[3];
sx q[3];
rz(-2.4468138) q[3];
sx q[3];
rz(3.1082603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];