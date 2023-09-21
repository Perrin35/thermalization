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
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6715235) q[0];
sx q[0];
rz(-1.0363665) q[0];
sx q[0];
rz(3.021391) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1349749) q[2];
sx q[2];
rz(-1.0834603) q[2];
sx q[2];
rz(2.6589573) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0109947) q[1];
sx q[1];
rz(-1.8124609) q[1];
sx q[1];
rz(-0.17250891) q[1];
x q[2];
rz(-1.7484619) q[3];
sx q[3];
rz(-1.2347617) q[3];
sx q[3];
rz(0.48776585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1774896) q[2];
sx q[2];
rz(-0.51012817) q[2];
sx q[2];
rz(2.5906079) q[2];
rz(-1.3059113) q[3];
sx q[3];
rz(-1.4923613) q[3];
sx q[3];
rz(1.8252385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47857639) q[0];
sx q[0];
rz(-2.1415648) q[0];
sx q[0];
rz(-0.4719032) q[0];
rz(2.7117803) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(-0.93634161) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82114007) q[0];
sx q[0];
rz(-1.3409412) q[0];
sx q[0];
rz(1.4810522) q[0];
x q[1];
rz(-0.78511946) q[2];
sx q[2];
rz(-1.2650507) q[2];
sx q[2];
rz(-0.53346764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1658926) q[1];
sx q[1];
rz(-0.70043889) q[1];
sx q[1];
rz(0.16209929) q[1];
rz(-pi) q[2];
x q[2];
rz(1.068088) q[3];
sx q[3];
rz(-1.5503251) q[3];
sx q[3];
rz(-2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(-2.7152087) q[2];
rz(-1.9042227) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(-3.1085076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24580978) q[0];
sx q[0];
rz(-1.3170467) q[0];
sx q[0];
rz(2.202503) q[0];
rz(2.242873) q[1];
sx q[1];
rz(-0.4788613) q[1];
sx q[1];
rz(-0.59392196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.827841) q[0];
sx q[0];
rz(-1.8694287) q[0];
sx q[0];
rz(1.9065501) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8823231) q[2];
sx q[2];
rz(-1.5261298) q[2];
sx q[2];
rz(-2.4633212) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2908823) q[1];
sx q[1];
rz(-1.1516654) q[1];
sx q[1];
rz(-2.9078729) q[1];
x q[2];
rz(2.6731554) q[3];
sx q[3];
rz(-1.2611946) q[3];
sx q[3];
rz(2.7999511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64017355) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(-1.7017986) q[2];
rz(-2.7539608) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(-0.38813996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.3751635) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(0.50278062) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-0.50351024) q[1];
sx q[1];
rz(2.3847413) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2198974) q[0];
sx q[0];
rz(-1.9031525) q[0];
sx q[0];
rz(0.38145782) q[0];
rz(-2.0365305) q[2];
sx q[2];
rz(-1.5823936) q[2];
sx q[2];
rz(0.84601814) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13285747) q[1];
sx q[1];
rz(-0.99616226) q[1];
sx q[1];
rz(0.43032129) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73369153) q[3];
sx q[3];
rz(-1.1847704) q[3];
sx q[3];
rz(-1.8611849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7148774) q[2];
sx q[2];
rz(-1.2299512) q[2];
sx q[2];
rz(1.654401) q[2];
rz(-0.58250827) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(-2.5845161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61705631) q[0];
sx q[0];
rz(-2.7334088) q[0];
sx q[0];
rz(-0.88390669) q[0];
rz(0.44287037) q[2];
sx q[2];
rz(-2.1211229) q[2];
sx q[2];
rz(-2.0572822) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95549059) q[1];
sx q[1];
rz(-1.0483861) q[1];
sx q[1];
rz(-2.2239457) q[1];
rz(-0.53253048) q[3];
sx q[3];
rz(-2.2610287) q[3];
sx q[3];
rz(1.1158451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.918255) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(0.63344947) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(-0.25318405) q[0];
rz(-1.6075915) q[1];
sx q[1];
rz(-1.7065159) q[1];
sx q[1];
rz(1.4621428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77308649) q[0];
sx q[0];
rz(-2.9162772) q[0];
sx q[0];
rz(-2.7789475) q[0];
rz(-0.92163779) q[2];
sx q[2];
rz(-1.7643133) q[2];
sx q[2];
rz(2.9606539) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51965442) q[1];
sx q[1];
rz(-2.7192273) q[1];
sx q[1];
rz(-0.42610355) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2268279) q[3];
sx q[3];
rz(-1.2729537) q[3];
sx q[3];
rz(0.96836585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9399461) q[2];
sx q[2];
rz(-0.74735171) q[2];
sx q[2];
rz(2.3366826) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-1.0369119) q[3];
sx q[3];
rz(1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(-0.92765635) q[0];
rz(2.1169128) q[1];
sx q[1];
rz(-1.6352446) q[1];
sx q[1];
rz(1.0120846) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1597848) q[0];
sx q[0];
rz(-0.73043981) q[0];
sx q[0];
rz(-0.83321379) q[0];
x q[1];
rz(-0.32832844) q[2];
sx q[2];
rz(-0.40847455) q[2];
sx q[2];
rz(2.7834053) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.583657) q[1];
sx q[1];
rz(-2.5968938) q[1];
sx q[1];
rz(-0.063636585) q[1];
x q[2];
rz(3.0104962) q[3];
sx q[3];
rz(-0.56250611) q[3];
sx q[3];
rz(-1.2870115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.53081375) q[2];
sx q[2];
rz(-1.6616219) q[2];
sx q[2];
rz(2.7116595) q[2];
rz(2.1271465) q[3];
sx q[3];
rz(-2.7323664) q[3];
sx q[3];
rz(-0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5291418) q[0];
sx q[0];
rz(-2.2021459) q[0];
sx q[0];
rz(0.21417831) q[0];
rz(-2.0902436) q[1];
sx q[1];
rz(-0.21251692) q[1];
sx q[1];
rz(-2.8578551) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8003214) q[0];
sx q[0];
rz(-2.8386136) q[0];
sx q[0];
rz(0.11462258) q[0];
rz(-pi) q[1];
rz(-2.9494638) q[2];
sx q[2];
rz(-1.2461975) q[2];
sx q[2];
rz(-1.0629551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2092065) q[1];
sx q[1];
rz(-1.7665518) q[1];
sx q[1];
rz(2.6372361) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3582718) q[3];
sx q[3];
rz(-2.3187227) q[3];
sx q[3];
rz(2.2058723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(-1.696375) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.701136) q[3];
sx q[3];
rz(0.33932313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.63672367) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(0.069256393) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.2830877) q[1];
sx q[1];
rz(1.5690631) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208961) q[0];
sx q[0];
rz(-2.3224152) q[0];
sx q[0];
rz(-0.92185123) q[0];
x q[1];
rz(-0.32585085) q[2];
sx q[2];
rz(-1.648765) q[2];
sx q[2];
rz(0.98380145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6988277) q[1];
sx q[1];
rz(-2.8932533) q[1];
sx q[1];
rz(-0.34766867) q[1];
rz(-pi) q[2];
rz(0.92215718) q[3];
sx q[3];
rz(-1.4033917) q[3];
sx q[3];
rz(2.5006014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9562324) q[2];
sx q[2];
rz(-0.23024836) q[2];
sx q[2];
rz(-3.0017079) q[2];
rz(0.36758962) q[3];
sx q[3];
rz(-1.1871754) q[3];
sx q[3];
rz(-2.1504413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96520987) q[0];
sx q[0];
rz(-2.7503224) q[0];
sx q[0];
rz(-0.64176732) q[0];
rz(-1.2311252) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7627015) q[0];
sx q[0];
rz(-2.2305616) q[0];
sx q[0];
rz(0.99878175) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3875302) q[2];
sx q[2];
rz(-1.4472618) q[2];
sx q[2];
rz(-2.0545517) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.99991998) q[1];
sx q[1];
rz(-1.8598286) q[1];
sx q[1];
rz(0.33860597) q[1];
rz(-0.28585163) q[3];
sx q[3];
rz(-1.0470069) q[3];
sx q[3];
rz(1.7459735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81007593) q[2];
sx q[2];
rz(-1.5248359) q[2];
sx q[2];
rz(1.6646741) q[2];
rz(-2.8752575) q[3];
sx q[3];
rz(-2.895152) q[3];
sx q[3];
rz(2.5951071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.1521679) q[2];
sx q[2];
rz(-0.94716723) q[2];
sx q[2];
rz(-0.92826044) q[2];
rz(-0.6775425) q[3];
sx q[3];
rz(-1.7384221) q[3];
sx q[3];
rz(1.3330028) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];