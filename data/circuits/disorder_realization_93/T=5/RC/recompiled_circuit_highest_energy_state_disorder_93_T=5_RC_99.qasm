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
rz(-0.45577249) q[0];
sx q[0];
rz(-1.1205641) q[0];
sx q[0];
rz(1.4135345) q[0];
rz(0.36355525) q[1];
sx q[1];
rz(4.9699291) q[1];
sx q[1];
rz(15.416754) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7257627) q[0];
sx q[0];
rz(-0.32074577) q[0];
sx q[0];
rz(1.6345535) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9102664) q[2];
sx q[2];
rz(-1.6874773) q[2];
sx q[2];
rz(1.9872023) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26303459) q[1];
sx q[1];
rz(-0.59459762) q[1];
sx q[1];
rz(-0.27185284) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91787429) q[3];
sx q[3];
rz(-2.2736227) q[3];
sx q[3];
rz(0.70706409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1353961) q[2];
sx q[2];
rz(-1.5781032) q[2];
sx q[2];
rz(-0.023999365) q[2];
rz(-2.7855347) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(-0.02478987) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559219) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(0.63496494) q[0];
rz(0.7629281) q[1];
sx q[1];
rz(-1.678391) q[1];
sx q[1];
rz(0.91711226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9207431) q[0];
sx q[0];
rz(-1.5891599) q[0];
sx q[0];
rz(-1.5361274) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6659662) q[2];
sx q[2];
rz(-0.94671072) q[2];
sx q[2];
rz(1.9204572) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.13027993) q[1];
sx q[1];
rz(-0.44751274) q[1];
sx q[1];
rz(1.7496197) q[1];
rz(-pi) q[2];
rz(-2.1907552) q[3];
sx q[3];
rz(-1.7449494) q[3];
sx q[3];
rz(0.0096498409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1941173) q[2];
sx q[2];
rz(-2.6185991) q[2];
sx q[2];
rz(2.4083162) q[2];
rz(1.0670079) q[3];
sx q[3];
rz(-1.4209483) q[3];
sx q[3];
rz(-2.9616621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(2.5718403) q[0];
sx q[0];
rz(-1.1886007) q[0];
sx q[0];
rz(0.52823129) q[0];
rz(0.86604467) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(-0.51308647) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24241867) q[0];
sx q[0];
rz(-2.0089825) q[0];
sx q[0];
rz(0.04695453) q[0];
rz(-pi) q[1];
rz(-2.2202291) q[2];
sx q[2];
rz(-1.8214392) q[2];
sx q[2];
rz(-1.373137) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1248438) q[1];
sx q[1];
rz(-1.4433772) q[1];
sx q[1];
rz(1.9215073) q[1];
rz(-2.641898) q[3];
sx q[3];
rz(-1.5958188) q[3];
sx q[3];
rz(-0.5145038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.036006007) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(3.1028683) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-2.3364412) q[3];
sx q[3];
rz(1.800644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7078581) q[0];
sx q[0];
rz(-0.41446328) q[0];
sx q[0];
rz(1.1210972) q[0];
rz(2.3838249) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(0.68414348) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8389177) q[0];
sx q[0];
rz(-2.3771299) q[0];
sx q[0];
rz(1.27728) q[0];
rz(2.5252417) q[2];
sx q[2];
rz(-3.1292097) q[2];
sx q[2];
rz(-1.3146666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0348548) q[1];
sx q[1];
rz(-1.8632006) q[1];
sx q[1];
rz(-2.8234944) q[1];
x q[2];
rz(-2.6287967) q[3];
sx q[3];
rz(-1.0362451) q[3];
sx q[3];
rz(0.057138047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6571558) q[2];
sx q[2];
rz(-1.6146722) q[2];
sx q[2];
rz(-2.8810775) q[2];
rz(-2.303463) q[3];
sx q[3];
rz(-1.1476436) q[3];
sx q[3];
rz(0.014613541) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218006) q[0];
sx q[0];
rz(-1.7638693) q[0];
sx q[0];
rz(-2.6309784) q[0];
rz(-1.249373) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(1.6872663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4076865) q[0];
sx q[0];
rz(-2.9260194) q[0];
sx q[0];
rz(3.0093497) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32938214) q[2];
sx q[2];
rz(-0.14757809) q[2];
sx q[2];
rz(-1.3775228) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92125398) q[1];
sx q[1];
rz(-0.61437139) q[1];
sx q[1];
rz(-2.6831616) q[1];
rz(-pi) q[2];
rz(0.84564836) q[3];
sx q[3];
rz(-0.21580869) q[3];
sx q[3];
rz(-2.8874021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5877567) q[2];
sx q[2];
rz(-2.7601056) q[2];
sx q[2];
rz(0.39056632) q[2];
rz(-0.25911123) q[3];
sx q[3];
rz(-1.5026389) q[3];
sx q[3];
rz(-2.794877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.7873586) q[0];
sx q[0];
rz(-2.7095095) q[0];
sx q[0];
rz(0.64467347) q[0];
rz(1.8395754) q[1];
sx q[1];
rz(-1.527486) q[1];
sx q[1];
rz(-0.34861809) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020418305) q[0];
sx q[0];
rz(-2.0266337) q[0];
sx q[0];
rz(2.0860079) q[0];
rz(-0.63747823) q[2];
sx q[2];
rz(-2.5036252) q[2];
sx q[2];
rz(1.5359985) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7182838) q[1];
sx q[1];
rz(-1.2759142) q[1];
sx q[1];
rz(2.0096094) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8556267) q[3];
sx q[3];
rz(-3.0197201) q[3];
sx q[3];
rz(1.7274477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1031441) q[2];
sx q[2];
rz(-2.1861173) q[2];
sx q[2];
rz(3.1032584) q[2];
rz(0.96758715) q[3];
sx q[3];
rz(-1.7318232) q[3];
sx q[3];
rz(-0.35965317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86647812) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(3.0679829) q[0];
rz(-1.4069125) q[1];
sx q[1];
rz(-0.55714566) q[1];
sx q[1];
rz(-2.4085192) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382815) q[0];
sx q[0];
rz(-1.065863) q[0];
sx q[0];
rz(2.5236593) q[0];
rz(0.91821155) q[2];
sx q[2];
rz(-0.84647734) q[2];
sx q[2];
rz(-2.0744963) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16790529) q[1];
sx q[1];
rz(-1.9009155) q[1];
sx q[1];
rz(0.37574212) q[1];
x q[2];
rz(-1.8215976) q[3];
sx q[3];
rz(-1.3426675) q[3];
sx q[3];
rz(3.1056014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0343895) q[2];
sx q[2];
rz(-1.4185536) q[2];
sx q[2];
rz(3.1305195) q[2];
rz(-0.30558807) q[3];
sx q[3];
rz(-1.1253076) q[3];
sx q[3];
rz(-1.2529469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3333862) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(-2.407684) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-2.1323233) q[1];
sx q[1];
rz(0.098800585) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0694612) q[0];
sx q[0];
rz(-2.0939576) q[0];
sx q[0];
rz(0.69055542) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96510624) q[2];
sx q[2];
rz(-2.1557376) q[2];
sx q[2];
rz(-3.0966058) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1360201) q[1];
sx q[1];
rz(-0.51678777) q[1];
sx q[1];
rz(-0.52424519) q[1];
x q[2];
rz(-0.55136724) q[3];
sx q[3];
rz(-0.15227642) q[3];
sx q[3];
rz(-0.65413977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29844478) q[2];
sx q[2];
rz(-1.8567905) q[2];
sx q[2];
rz(-1.0996381) q[2];
rz(-0.034865033) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(2.5449424) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027503969) q[0];
sx q[0];
rz(-1.9845668) q[0];
sx q[0];
rz(-1.9644568) q[0];
rz(1.6674532) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(1.6966381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5367573) q[0];
sx q[0];
rz(-2.6492181) q[0];
sx q[0];
rz(2.5144666) q[0];
rz(-pi) q[1];
rz(2.1423526) q[2];
sx q[2];
rz(-2.2543467) q[2];
sx q[2];
rz(0.37887606) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9517773) q[1];
sx q[1];
rz(-0.51945247) q[1];
sx q[1];
rz(1.431926) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2463465) q[3];
sx q[3];
rz(-1.7887247) q[3];
sx q[3];
rz(-2.982614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16704796) q[2];
sx q[2];
rz(-1.5571152) q[2];
sx q[2];
rz(0.17189279) q[2];
rz(0.78957549) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(-1.8062228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852785) q[0];
sx q[0];
rz(-0.21509898) q[0];
sx q[0];
rz(0.54157448) q[0];
rz(1.1594634) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(0.83311876) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3726927) q[0];
sx q[0];
rz(-1.295255) q[0];
sx q[0];
rz(-2.5252456) q[0];
rz(1.9160742) q[2];
sx q[2];
rz(-1.0547027) q[2];
sx q[2];
rz(-2.6654625) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1320859) q[1];
sx q[1];
rz(-2.2530547) q[1];
sx q[1];
rz(-0.29341523) q[1];
x q[2];
rz(2.2984357) q[3];
sx q[3];
rz(-1.5693446) q[3];
sx q[3];
rz(0.42556083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.021576015) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(0.94338256) q[2];
rz(1.0776445) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(2.8130892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42644603) q[0];
sx q[0];
rz(-2.1052512) q[0];
sx q[0];
rz(-0.25682009) q[0];
rz(0.23964755) q[1];
sx q[1];
rz(-2.2365204) q[1];
sx q[1];
rz(2.0047275) q[1];
rz(-0.69375615) q[2];
sx q[2];
rz(-1.0508595) q[2];
sx q[2];
rz(-1.6169242) q[2];
rz(-2.8311493) q[3];
sx q[3];
rz(-1.8127999) q[3];
sx q[3];
rz(-1.8565069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
