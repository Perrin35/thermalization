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
rz(5.1626212) q[0];
sx q[0];
rz(10.838312) q[0];
rz(0.36355525) q[1];
sx q[1];
rz(-1.3132562) q[1];
sx q[1];
rz(-0.29120905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6585892) q[0];
sx q[0];
rz(-1.2507255) q[0];
sx q[0];
rz(-0.021163737) q[0];
rz(1.2313263) q[2];
sx q[2];
rz(-1.4541153) q[2];
sx q[2];
rz(-1.1543903) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.061569842) q[1];
sx q[1];
rz(-2.1407619) q[1];
sx q[1];
rz(-1.391173) q[1];
rz(-2.324027) q[3];
sx q[3];
rz(-1.0888087) q[3];
sx q[3];
rz(2.7369902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0061965813) q[2];
sx q[2];
rz(-1.5781032) q[2];
sx q[2];
rz(0.023999365) q[2];
rz(0.35605797) q[3];
sx q[3];
rz(-2.4672716) q[3];
sx q[3];
rz(-0.02478987) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0559219) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(0.63496494) q[0];
rz(2.3786646) q[1];
sx q[1];
rz(-1.4632016) q[1];
sx q[1];
rz(0.91711226) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.791009) q[0];
sx q[0];
rz(-1.6054594) q[0];
sx q[0];
rz(3.123218) q[0];
rz(-pi) q[1];
rz(0.88998743) q[2];
sx q[2];
rz(-1.9515079) q[2];
sx q[2];
rz(-0.057304545) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2789816) q[1];
sx q[1];
rz(-1.4937506) q[1];
sx q[1];
rz(-1.129523) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1907552) q[3];
sx q[3];
rz(-1.7449494) q[3];
sx q[3];
rz(-0.0096498409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94747535) q[2];
sx q[2];
rz(-0.52299356) q[2];
sx q[2];
rz(-0.73327649) q[2];
rz(-1.0670079) q[3];
sx q[3];
rz(-1.7206444) q[3];
sx q[3];
rz(0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.5718403) q[0];
sx q[0];
rz(-1.1886007) q[0];
sx q[0];
rz(2.6133614) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.3287611) q[1];
sx q[1];
rz(-0.51308647) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7932803) q[0];
sx q[0];
rz(-1.6133119) q[0];
sx q[0];
rz(-1.1321863) q[0];
x q[1];
rz(-0.31103525) q[2];
sx q[2];
rz(-0.94488178) q[2];
sx q[2];
rz(-3.1300822) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5410903) q[1];
sx q[1];
rz(-1.2230495) q[1];
sx q[1];
rz(-3.0060124) q[1];
rz(1.5422899) q[3];
sx q[3];
rz(-1.0712726) q[3];
sx q[3];
rz(-1.0699501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1055866) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(-3.1028683) q[2];
rz(2.683908) q[3];
sx q[3];
rz(-2.3364412) q[3];
sx q[3];
rz(-1.800644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4337346) q[0];
sx q[0];
rz(-0.41446328) q[0];
sx q[0];
rz(-1.1210972) q[0];
rz(0.7577678) q[1];
sx q[1];
rz(-1.6242177) q[1];
sx q[1];
rz(2.4574492) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0476889) q[0];
sx q[0];
rz(-2.2950116) q[0];
sx q[0];
rz(0.27064497) q[0];
rz(-pi) q[1];
rz(0.61635095) q[2];
sx q[2];
rz(-0.01238298) q[2];
sx q[2];
rz(-1.3146666) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0348548) q[1];
sx q[1];
rz(-1.2783921) q[1];
sx q[1];
rz(-2.8234944) q[1];
rz(-0.97400333) q[3];
sx q[3];
rz(-1.1349548) q[3];
sx q[3];
rz(-1.3485933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4844369) q[2];
sx q[2];
rz(-1.6146722) q[2];
sx q[2];
rz(-2.8810775) q[2];
rz(-0.83812964) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(-3.1269791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.8922197) q[1];
sx q[1];
rz(-1.9793972) q[1];
sx q[1];
rz(-1.6872663) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.175486) q[0];
sx q[0];
rz(-1.5425872) q[0];
sx q[0];
rz(2.9278446) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6188443) q[2];
sx q[2];
rz(-1.4312051) q[2];
sx q[2];
rz(-2.0968116) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.37794681) q[1];
sx q[1];
rz(-1.0275405) q[1];
sx q[1];
rz(-1.8733979) q[1];
rz(-pi) q[2];
rz(-2.9972059) q[3];
sx q[3];
rz(-1.4098415) q[3];
sx q[3];
rz(0.99100366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.55383596) q[2];
sx q[2];
rz(-2.7601056) q[2];
sx q[2];
rz(-0.39056632) q[2];
rz(-0.25911123) q[3];
sx q[3];
rz(-1.5026389) q[3];
sx q[3];
rz(-2.794877) q[3];
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
rz(-pi/2) q[0];
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
rz(-2.7873586) q[0];
sx q[0];
rz(-0.43208313) q[0];
sx q[0];
rz(2.4969192) q[0];
rz(1.8395754) q[1];
sx q[1];
rz(-1.527486) q[1];
sx q[1];
rz(-0.34861809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520599) q[0];
sx q[0];
rz(-2.4675998) q[0];
sx q[0];
rz(0.78788449) q[0];
rz(-0.63747823) q[2];
sx q[2];
rz(-0.63796746) q[2];
sx q[2];
rz(-1.5359985) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4069566) q[1];
sx q[1];
rz(-0.52328456) q[1];
sx q[1];
rz(0.95013817) q[1];
x q[2];
rz(0.28596596) q[3];
sx q[3];
rz(-0.12187258) q[3];
sx q[3];
rz(1.7274477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0384486) q[2];
sx q[2];
rz(-2.1861173) q[2];
sx q[2];
rz(0.038334282) q[2];
rz(0.96758715) q[3];
sx q[3];
rz(-1.4097694) q[3];
sx q[3];
rz(0.35965317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86647812) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(0.073609784) q[0];
rz(1.4069125) q[1];
sx q[1];
rz(-0.55714566) q[1];
sx q[1];
rz(-0.73307347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033111) q[0];
sx q[0];
rz(-2.0757296) q[0];
sx q[0];
rz(0.61793332) q[0];
rz(-pi) q[1];
rz(2.5400852) q[2];
sx q[2];
rz(-2.2081293) q[2];
sx q[2];
rz(-1.21797) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0907369) q[1];
sx q[1];
rz(-2.6466718) q[1];
sx q[1];
rz(-0.75116091) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2352318) q[3];
sx q[3];
rz(-1.8149655) q[3];
sx q[3];
rz(-1.5489123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0343895) q[2];
sx q[2];
rz(-1.4185536) q[2];
sx q[2];
rz(3.1305195) q[2];
rz(2.8360046) q[3];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3333862) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(2.407684) q[0];
rz(0.58087307) q[1];
sx q[1];
rz(-2.1323233) q[1];
sx q[1];
rz(-0.098800585) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0990813) q[0];
sx q[0];
rz(-0.83957273) q[0];
sx q[0];
rz(-0.73584105) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4634741) q[2];
sx q[2];
rz(-2.0653915) q[2];
sx q[2];
rz(1.8910318) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0055726) q[1];
sx q[1];
rz(-2.6248049) q[1];
sx q[1];
rz(0.52424519) q[1];
rz(-1.6510165) q[3];
sx q[3];
rz(-1.7003683) q[3];
sx q[3];
rz(-0.097565325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29844478) q[2];
sx q[2];
rz(-1.8567905) q[2];
sx q[2];
rz(2.0419545) q[2];
rz(-0.034865033) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(-0.59665027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1140887) q[0];
sx q[0];
rz(-1.1570258) q[0];
sx q[0];
rz(-1.1771359) q[0];
rz(1.6674532) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(1.6966381) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5367573) q[0];
sx q[0];
rz(-2.6492181) q[0];
sx q[0];
rz(0.62712609) q[0];
rz(-pi) q[1];
rz(0.58622245) q[2];
sx q[2];
rz(-0.86044035) q[2];
sx q[2];
rz(-1.9682056) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.260238) q[1];
sx q[1];
rz(-1.6395651) q[1];
sx q[1];
rz(2.0860904) q[1];
x q[2];
rz(-0.27650303) q[3];
sx q[3];
rz(-0.91405896) q[3];
sx q[3];
rz(-1.9013251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16704796) q[2];
sx q[2];
rz(-1.5571152) q[2];
sx q[2];
rz(2.9696999) q[2];
rz(2.3520172) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(-1.3353698) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852785) q[0];
sx q[0];
rz(-0.21509898) q[0];
sx q[0];
rz(0.54157448) q[0];
rz(-1.1594634) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(2.3084739) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.706382) q[0];
sx q[0];
rz(-0.66775371) q[0];
sx q[0];
rz(2.6866962) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9160742) q[2];
sx q[2];
rz(-2.08689) q[2];
sx q[2];
rz(0.47613019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.009506735) q[1];
sx q[1];
rz(-2.2530547) q[1];
sx q[1];
rz(-0.29341523) q[1];
rz(-1.5686137) q[3];
sx q[3];
rz(-2.4139521) q[3];
sx q[3];
rz(-1.997987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1200166) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(-2.1982101) q[2];
rz(-1.0776445) q[3];
sx q[3];
rz(-1.5543944) q[3];
sx q[3];
rz(0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7151466) q[0];
sx q[0];
rz(-1.0363415) q[0];
sx q[0];
rz(2.8847726) q[0];
rz(-2.9019451) q[1];
sx q[1];
rz(-2.2365204) q[1];
sx q[1];
rz(2.0047275) q[1];
rz(2.4478365) q[2];
sx q[2];
rz(-1.0508595) q[2];
sx q[2];
rz(-1.6169242) q[2];
rz(0.31044337) q[3];
sx q[3];
rz(-1.8127999) q[3];
sx q[3];
rz(-1.8565069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
