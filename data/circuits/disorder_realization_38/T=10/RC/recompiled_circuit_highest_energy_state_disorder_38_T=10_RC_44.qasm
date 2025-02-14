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
rz(2.8528557) q[0];
sx q[0];
rz(5.5878162) q[0];
sx q[0];
rz(6.5509808) q[0];
rz(0.42203045) q[1];
sx q[1];
rz(-2.2208417) q[1];
sx q[1];
rz(1.2738127) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52142329) q[0];
sx q[0];
rz(-1.6156989) q[0];
sx q[0];
rz(1.8123167) q[0];
rz(-pi) q[1];
rz(2.2620413) q[2];
sx q[2];
rz(-0.82720199) q[2];
sx q[2];
rz(0.66397053) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.85292683) q[1];
sx q[1];
rz(-1.5024606) q[1];
sx q[1];
rz(-0.6612551) q[1];
rz(1.3269365) q[3];
sx q[3];
rz(-2.993846) q[3];
sx q[3];
rz(-2.6105931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5219118) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(-0.042595159) q[2];
rz(0.28111449) q[3];
sx q[3];
rz(-1.5646076) q[3];
sx q[3];
rz(0.59578305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5113145) q[0];
sx q[0];
rz(-1.9673286) q[0];
sx q[0];
rz(1.0761155) q[0];
rz(-1.0307182) q[1];
sx q[1];
rz(-2.0298256) q[1];
sx q[1];
rz(-1.3105185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.977518) q[0];
sx q[0];
rz(-0.65854544) q[0];
sx q[0];
rz(1.6048628) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83514799) q[2];
sx q[2];
rz(-1.877376) q[2];
sx q[2];
rz(0.94294244) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0243822) q[1];
sx q[1];
rz(-1.6892994) q[1];
sx q[1];
rz(0.16853965) q[1];
x q[2];
rz(1.1635246) q[3];
sx q[3];
rz(-1.0735452) q[3];
sx q[3];
rz(2.6866792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92404667) q[2];
sx q[2];
rz(-2.3895538) q[2];
sx q[2];
rz(-0.95620608) q[2];
rz(1.8337967) q[3];
sx q[3];
rz(-2.3266413) q[3];
sx q[3];
rz(3.0202878) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2995375) q[0];
sx q[0];
rz(-0.50899035) q[0];
sx q[0];
rz(-0.036238413) q[0];
rz(-2.3402975) q[1];
sx q[1];
rz(-1.5472629) q[1];
sx q[1];
rz(-1.8410199) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7177082) q[0];
sx q[0];
rz(-2.9489655) q[0];
sx q[0];
rz(0.9813375) q[0];
rz(-pi) q[1];
rz(-1.6341798) q[2];
sx q[2];
rz(-2.5434539) q[2];
sx q[2];
rz(1.6063521) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3276931) q[1];
sx q[1];
rz(-2.3997953) q[1];
sx q[1];
rz(-0.90403907) q[1];
rz(-pi) q[2];
rz(1.9375422) q[3];
sx q[3];
rz(-2.2270508) q[3];
sx q[3];
rz(1.7451857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3761882) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(-2.9236531) q[2];
rz(0.91207063) q[3];
sx q[3];
rz(-1.4927161) q[3];
sx q[3];
rz(0.85165858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7243778) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(2.3260314) q[0];
rz(2.0888445) q[1];
sx q[1];
rz(-2.2595854) q[1];
sx q[1];
rz(-1.7900593) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4690875) q[0];
sx q[0];
rz(-0.47576213) q[0];
sx q[0];
rz(-0.16938727) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4250373) q[2];
sx q[2];
rz(-1.5045696) q[2];
sx q[2];
rz(-0.19342455) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0114944) q[1];
sx q[1];
rz(-1.5421499) q[1];
sx q[1];
rz(-1.0953893) q[1];
rz(-pi) q[2];
rz(-1.0686605) q[3];
sx q[3];
rz(-2.2229513) q[3];
sx q[3];
rz(2.7606635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1286596) q[2];
sx q[2];
rz(-1.9482875) q[2];
sx q[2];
rz(-0.43221691) q[2];
rz(-2.2991119) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(0.99561083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1603482) q[0];
sx q[0];
rz(-2.6762185) q[0];
sx q[0];
rz(-0.55111849) q[0];
rz(-2.5235858) q[1];
sx q[1];
rz(-1.8564686) q[1];
sx q[1];
rz(-1.0689703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9286802) q[0];
sx q[0];
rz(-1.3569662) q[0];
sx q[0];
rz(-0.68709157) q[0];
x q[1];
rz(2.6591797) q[2];
sx q[2];
rz(-1.97792) q[2];
sx q[2];
rz(2.3633752) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7643508) q[1];
sx q[1];
rz(-0.83989401) q[1];
sx q[1];
rz(2.9725171) q[1];
rz(-pi) q[2];
rz(-0.25311796) q[3];
sx q[3];
rz(-2.6147644) q[3];
sx q[3];
rz(1.7071068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3507639) q[2];
sx q[2];
rz(-0.5841693) q[2];
sx q[2];
rz(0.9551777) q[2];
rz(2.5480934) q[3];
sx q[3];
rz(-0.94624001) q[3];
sx q[3];
rz(-1.745863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3747568) q[0];
sx q[0];
rz(-0.042348472) q[0];
sx q[0];
rz(1.7497077) q[0];
rz(-1.998924) q[1];
sx q[1];
rz(-1.7997768) q[1];
sx q[1];
rz(-1.3124189) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4785193) q[0];
sx q[0];
rz(-2.0475302) q[0];
sx q[0];
rz(2.3011424) q[0];
x q[1];
rz(0.80760132) q[2];
sx q[2];
rz(-0.86059216) q[2];
sx q[2];
rz(-1.0654895) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23199546) q[1];
sx q[1];
rz(-1.0568406) q[1];
sx q[1];
rz(-1.6327052) q[1];
rz(-pi) q[2];
rz(-3.083701) q[3];
sx q[3];
rz(-1.2558008) q[3];
sx q[3];
rz(-1.0703027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.352508) q[2];
sx q[2];
rz(-2.3859873) q[2];
sx q[2];
rz(1.7279203) q[2];
rz(0.46755725) q[3];
sx q[3];
rz(-2.4831725) q[3];
sx q[3];
rz(-2.5889034) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5107875) q[0];
sx q[0];
rz(-0.063022114) q[0];
sx q[0];
rz(0.68156534) q[0];
rz(-1.1850146) q[1];
sx q[1];
rz(-2.3763035) q[1];
sx q[1];
rz(2.3550745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1543321) q[0];
sx q[0];
rz(-1.2367147) q[0];
sx q[0];
rz(-2.7272237) q[0];
rz(-2.6441387) q[2];
sx q[2];
rz(-1.0062394) q[2];
sx q[2];
rz(-1.4184784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9481755) q[1];
sx q[1];
rz(-1.6778291) q[1];
sx q[1];
rz(2.2501037) q[1];
rz(-0.041140838) q[3];
sx q[3];
rz(-0.84008559) q[3];
sx q[3];
rz(1.5276437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6098392) q[2];
sx q[2];
rz(-2.9015151) q[2];
sx q[2];
rz(1.5698203) q[2];
rz(1.2872559) q[3];
sx q[3];
rz(-1.6183034) q[3];
sx q[3];
rz(0.94304812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.542273) q[0];
sx q[0];
rz(-2.4241408) q[0];
sx q[0];
rz(-1.9532816) q[0];
rz(0.020847281) q[1];
sx q[1];
rz(-2.3884845) q[1];
sx q[1];
rz(0.59897024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99020236) q[0];
sx q[0];
rz(-1.6693475) q[0];
sx q[0];
rz(2.0771785) q[0];
x q[1];
rz(-2.2856667) q[2];
sx q[2];
rz(-1.8740665) q[2];
sx q[2];
rz(1.0452458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5237732) q[1];
sx q[1];
rz(-2.0460772) q[1];
sx q[1];
rz(2.3889524) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3009572) q[3];
sx q[3];
rz(-1.2328237) q[3];
sx q[3];
rz(2.4455051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.282436) q[2];
sx q[2];
rz(-1.891529) q[2];
sx q[2];
rz(0.51521987) q[2];
rz(-0.53269261) q[3];
sx q[3];
rz(-1.0849413) q[3];
sx q[3];
rz(2.2044619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135947) q[0];
sx q[0];
rz(-2.4585215) q[0];
sx q[0];
rz(-0.37044507) q[0];
rz(-2.0619552) q[1];
sx q[1];
rz(-0.41079435) q[1];
sx q[1];
rz(0.058578514) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97998842) q[0];
sx q[0];
rz(-1.4128886) q[0];
sx q[0];
rz(-0.61451332) q[0];
rz(-pi) q[1];
rz(-3.09715) q[2];
sx q[2];
rz(-1.7979597) q[2];
sx q[2];
rz(0.73611605) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9875659) q[1];
sx q[1];
rz(-2.3515764) q[1];
sx q[1];
rz(1.7055737) q[1];
rz(-pi) q[2];
rz(0.83567668) q[3];
sx q[3];
rz(-2.0234851) q[3];
sx q[3];
rz(1.9122461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8687245) q[2];
sx q[2];
rz(-2.3385907) q[2];
sx q[2];
rz(2.8105984) q[2];
rz(3.1281779) q[3];
sx q[3];
rz(-2.1881073) q[3];
sx q[3];
rz(2.0851871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5633504) q[0];
sx q[0];
rz(-0.42911068) q[0];
sx q[0];
rz(0.44813928) q[0];
rz(-3.0478802) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(1.2154481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5872204) q[0];
sx q[0];
rz(-2.0036812) q[0];
sx q[0];
rz(0.032399633) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3289069) q[2];
sx q[2];
rz(-1.2764837) q[2];
sx q[2];
rz(0.049969604) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3443904) q[1];
sx q[1];
rz(-2.0995525) q[1];
sx q[1];
rz(-0.31086287) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9691422) q[3];
sx q[3];
rz(-1.6961548) q[3];
sx q[3];
rz(1.0925869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1756246) q[2];
sx q[2];
rz(-1.5529239) q[2];
sx q[2];
rz(2.442181) q[2];
rz(1.8661963) q[3];
sx q[3];
rz(-1.9857152) q[3];
sx q[3];
rz(1.8705961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4324343) q[0];
sx q[0];
rz(-0.86031886) q[0];
sx q[0];
rz(1.1501089) q[0];
rz(1.3954096) q[1];
sx q[1];
rz(-1.9231053) q[1];
sx q[1];
rz(-1.3833192) q[1];
rz(3.0396661) q[2];
sx q[2];
rz(-0.66833767) q[2];
sx q[2];
rz(2.3705033) q[2];
rz(-0.46635177) q[3];
sx q[3];
rz(-1.6771355) q[3];
sx q[3];
rz(1.9775978) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
