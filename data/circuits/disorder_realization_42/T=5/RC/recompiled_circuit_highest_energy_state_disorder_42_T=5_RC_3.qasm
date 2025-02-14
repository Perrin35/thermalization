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
rz(0.49405721) q[0];
sx q[0];
rz(-0.20209514) q[0];
sx q[0];
rz(0.2986412) q[0];
rz(2.5477297) q[1];
sx q[1];
rz(-1.9057823) q[1];
sx q[1];
rz(-1.2716582) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0814046) q[0];
sx q[0];
rz(-1.5861298) q[0];
sx q[0];
rz(0.091214692) q[0];
x q[1];
rz(0.80940445) q[2];
sx q[2];
rz(-0.17091852) q[2];
sx q[2];
rz(-1.6910291) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9720481) q[1];
sx q[1];
rz(-1.8410334) q[1];
sx q[1];
rz(1.0098021) q[1];
x q[2];
rz(-0.27404578) q[3];
sx q[3];
rz(-1.9290302) q[3];
sx q[3];
rz(-0.90858118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9876081) q[2];
sx q[2];
rz(-0.5927023) q[2];
sx q[2];
rz(0.1565557) q[2];
rz(0.17078677) q[3];
sx q[3];
rz(-0.6905061) q[3];
sx q[3];
rz(1.4120215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60468948) q[0];
sx q[0];
rz(-1.5770183) q[0];
sx q[0];
rz(0.94737303) q[0];
rz(1.0700215) q[1];
sx q[1];
rz(-1.0173631) q[1];
sx q[1];
rz(2.1224461) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3257449) q[0];
sx q[0];
rz(-1.7293892) q[0];
sx q[0];
rz(2.1116381) q[0];
rz(0.48430125) q[2];
sx q[2];
rz(-1.9308085) q[2];
sx q[2];
rz(1.1126435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6883057) q[1];
sx q[1];
rz(-0.46436342) q[1];
sx q[1];
rz(-2.7771441) q[1];
rz(-0.89313497) q[3];
sx q[3];
rz(-2.1669029) q[3];
sx q[3];
rz(0.79465468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.77000874) q[2];
sx q[2];
rz(-1.2086478) q[2];
sx q[2];
rz(0.51119512) q[2];
rz(2.8080071) q[3];
sx q[3];
rz(-0.89190069) q[3];
sx q[3];
rz(0.79264486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12673512) q[0];
sx q[0];
rz(-0.52913409) q[0];
sx q[0];
rz(1.0134617) q[0];
rz(-0.86499372) q[1];
sx q[1];
rz(-1.8987055) q[1];
sx q[1];
rz(1.5369044) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0292708) q[0];
sx q[0];
rz(-1.6942612) q[0];
sx q[0];
rz(1.6804041) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2584707) q[2];
sx q[2];
rz(-0.97849023) q[2];
sx q[2];
rz(-1.2186714) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2257682) q[1];
sx q[1];
rz(-1.1686004) q[1];
sx q[1];
rz(-2.5798227) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99076267) q[3];
sx q[3];
rz(-1.2427075) q[3];
sx q[3];
rz(2.6716731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65261877) q[2];
sx q[2];
rz(-1.7964615) q[2];
sx q[2];
rz(2.6441003) q[2];
rz(0.8688212) q[3];
sx q[3];
rz(-0.38058623) q[3];
sx q[3];
rz(3.0697921) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3196165) q[0];
sx q[0];
rz(-0.0030586547) q[0];
sx q[0];
rz(0.36913607) q[0];
rz(2.4311851) q[1];
sx q[1];
rz(-1.0514759) q[1];
sx q[1];
rz(-1.2999473) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1064736) q[0];
sx q[0];
rz(-1.4958315) q[0];
sx q[0];
rz(2.3142561) q[0];
rz(-pi) q[1];
rz(-2.8244888) q[2];
sx q[2];
rz(-0.46296994) q[2];
sx q[2];
rz(-0.77808873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0365606) q[1];
sx q[1];
rz(-1.2249882) q[1];
sx q[1];
rz(-1.9329835) q[1];
rz(-pi) q[2];
rz(1.7930166) q[3];
sx q[3];
rz(-2.684703) q[3];
sx q[3];
rz(-2.8710988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15061441) q[2];
sx q[2];
rz(-1.4563478) q[2];
sx q[2];
rz(0.72188226) q[2];
rz(2.7041096) q[3];
sx q[3];
rz(-0.36967725) q[3];
sx q[3];
rz(0.30445254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337912) q[0];
sx q[0];
rz(-2.3423539) q[0];
sx q[0];
rz(0.21402982) q[0];
rz(-0.80263823) q[1];
sx q[1];
rz(-1.1624348) q[1];
sx q[1];
rz(-0.42073694) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1871757) q[0];
sx q[0];
rz(-0.19190776) q[0];
sx q[0];
rz(1.9918066) q[0];
rz(-pi) q[1];
rz(-2.7710469) q[2];
sx q[2];
rz(-2.5442313) q[2];
sx q[2];
rz(-1.0295402) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32728095) q[1];
sx q[1];
rz(-0.13209681) q[1];
sx q[1];
rz(-2.0943453) q[1];
rz(0.84263148) q[3];
sx q[3];
rz(-1.6602605) q[3];
sx q[3];
rz(-2.7842229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.78571) q[2];
sx q[2];
rz(-2.7648401) q[2];
sx q[2];
rz(2.9316087) q[2];
rz(-2.0436132) q[3];
sx q[3];
rz(-0.72303253) q[3];
sx q[3];
rz(0.22600225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2881154) q[0];
sx q[0];
rz(-2.2324039) q[0];
sx q[0];
rz(-0.01532456) q[0];
rz(2.3769936) q[1];
sx q[1];
rz(-1.8761643) q[1];
sx q[1];
rz(-0.24078029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6103196) q[0];
sx q[0];
rz(-1.1582147) q[0];
sx q[0];
rz(2.812435) q[0];
x q[1];
rz(1.316761) q[2];
sx q[2];
rz(-2.0587631) q[2];
sx q[2];
rz(0.72869638) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.0097373936) q[1];
sx q[1];
rz(-2.9564361) q[1];
sx q[1];
rz(1.2082165) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0475637) q[3];
sx q[3];
rz(-1.8935455) q[3];
sx q[3];
rz(0.18206319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8814055) q[2];
sx q[2];
rz(-2.0892961) q[2];
sx q[2];
rz(-0.50797272) q[2];
rz(-0.14902614) q[3];
sx q[3];
rz(-0.38145724) q[3];
sx q[3];
rz(1.4582567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7070865) q[0];
sx q[0];
rz(-0.91489804) q[0];
sx q[0];
rz(1.9515422) q[0];
rz(0.03935628) q[1];
sx q[1];
rz(-2.1085565) q[1];
sx q[1];
rz(0.2619704) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9761115) q[0];
sx q[0];
rz(-1.2369177) q[0];
sx q[0];
rz(-1.6209442) q[0];
x q[1];
rz(2.8112721) q[2];
sx q[2];
rz(-1.2093423) q[2];
sx q[2];
rz(-2.9320331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2319164) q[1];
sx q[1];
rz(-1.8662098) q[1];
sx q[1];
rz(-0.40628237) q[1];
x q[2];
rz(2.5983635) q[3];
sx q[3];
rz(-1.3564988) q[3];
sx q[3];
rz(-0.33583904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.94595853) q[2];
sx q[2];
rz(-0.28886637) q[2];
sx q[2];
rz(-1.3836916) q[2];
rz(-1.9687999) q[3];
sx q[3];
rz(-1.5044418) q[3];
sx q[3];
rz(-0.25243944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4381572) q[0];
sx q[0];
rz(-0.96980888) q[0];
sx q[0];
rz(0.54617149) q[0];
rz(2.3225885) q[1];
sx q[1];
rz(-1.1823697) q[1];
sx q[1];
rz(0.46772734) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2692102) q[0];
sx q[0];
rz(-0.2650241) q[0];
sx q[0];
rz(-1.7635959) q[0];
rz(-2.7502656) q[2];
sx q[2];
rz(-1.1839917) q[2];
sx q[2];
rz(-0.70287537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.88087) q[1];
sx q[1];
rz(-2.4308009) q[1];
sx q[1];
rz(-0.81501643) q[1];
rz(2.9727861) q[3];
sx q[3];
rz(-1.8106109) q[3];
sx q[3];
rz(0.60926357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8972682) q[2];
sx q[2];
rz(-1.5398328) q[2];
sx q[2];
rz(2.3648105) q[2];
rz(-1.4463834) q[3];
sx q[3];
rz(-1.8518238) q[3];
sx q[3];
rz(0.18073925) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.986213) q[0];
sx q[0];
rz(-1.4028832) q[0];
sx q[0];
rz(-3.0631995) q[0];
rz(0.29131237) q[1];
sx q[1];
rz(-0.71543175) q[1];
sx q[1];
rz(1.8142456) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9175668) q[0];
sx q[0];
rz(-0.97962681) q[0];
sx q[0];
rz(-0.81637189) q[0];
x q[1];
rz(0.38703309) q[2];
sx q[2];
rz(-3.0245028) q[2];
sx q[2];
rz(-0.69293222) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3811031) q[1];
sx q[1];
rz(-2.2020643) q[1];
sx q[1];
rz(-1.9105511) q[1];
x q[2];
rz(-1.4287656) q[3];
sx q[3];
rz(-1.8544594) q[3];
sx q[3];
rz(-1.7856959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3751117) q[2];
sx q[2];
rz(-1.8152619) q[2];
sx q[2];
rz(-0.43409902) q[2];
rz(0.33657524) q[3];
sx q[3];
rz(-1.429129) q[3];
sx q[3];
rz(-2.7460639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8074938) q[0];
sx q[0];
rz(-2.8039126) q[0];
sx q[0];
rz(0.39921528) q[0];
rz(-2.4576064) q[1];
sx q[1];
rz(-1.3147048) q[1];
sx q[1];
rz(1.7302053) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0610962) q[0];
sx q[0];
rz(-0.37758025) q[0];
sx q[0];
rz(2.6444673) q[0];
rz(-0.67839854) q[2];
sx q[2];
rz(-1.4015371) q[2];
sx q[2];
rz(0.20941833) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4343202) q[1];
sx q[1];
rz(-2.1967683) q[1];
sx q[1];
rz(2.0213338) q[1];
rz(-0.65220368) q[3];
sx q[3];
rz(-0.89487984) q[3];
sx q[3];
rz(1.9151198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4513678) q[2];
sx q[2];
rz(-1.5294315) q[2];
sx q[2];
rz(0.31768793) q[2];
rz(-2.132527) q[3];
sx q[3];
rz(-1.7694446) q[3];
sx q[3];
rz(0.38869977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7379363) q[0];
sx q[0];
rz(-1.4690348) q[0];
sx q[0];
rz(0.22672155) q[0];
rz(-0.33568385) q[1];
sx q[1];
rz(-2.2102368) q[1];
sx q[1];
rz(-0.23063901) q[1];
rz(0.027770859) q[2];
sx q[2];
rz(-0.66412974) q[2];
sx q[2];
rz(0.017505125) q[2];
rz(2.5981591) q[3];
sx q[3];
rz(-2.4590878) q[3];
sx q[3];
rz(1.2868846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
