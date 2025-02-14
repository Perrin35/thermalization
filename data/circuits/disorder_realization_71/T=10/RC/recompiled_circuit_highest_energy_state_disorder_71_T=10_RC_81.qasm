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
rz(-0.072862536) q[0];
sx q[0];
rz(-2.1003567) q[0];
sx q[0];
rz(0.6676724) q[0];
rz(-0.11198894) q[1];
sx q[1];
rz(-1.656124) q[1];
sx q[1];
rz(-2.8356584) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65016215) q[0];
sx q[0];
rz(-1.1592277) q[0];
sx q[0];
rz(1.6603907) q[0];
rz(0.602416) q[2];
sx q[2];
rz(-0.71068804) q[2];
sx q[2];
rz(1.3606233) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8851807) q[1];
sx q[1];
rz(-0.37799126) q[1];
sx q[1];
rz(2.9609326) q[1];
rz(-pi) q[2];
rz(1.63941) q[3];
sx q[3];
rz(-0.8655394) q[3];
sx q[3];
rz(0.88348639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2970807) q[2];
sx q[2];
rz(-1.4163372) q[2];
sx q[2];
rz(-2.1616006) q[2];
rz(-0.32198191) q[3];
sx q[3];
rz(-1.201509) q[3];
sx q[3];
rz(-2.9856248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62671536) q[0];
sx q[0];
rz(-1.6139655) q[0];
sx q[0];
rz(2.2514586) q[0];
rz(1.1164249) q[1];
sx q[1];
rz(-0.85464388) q[1];
sx q[1];
rz(1.9568303) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0380533) q[0];
sx q[0];
rz(-2.6067043) q[0];
sx q[0];
rz(0.37910342) q[0];
rz(-2.1380139) q[2];
sx q[2];
rz(-2.5126656) q[2];
sx q[2];
rz(2.8975482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5819232) q[1];
sx q[1];
rz(-1.6861177) q[1];
sx q[1];
rz(-1.4588463) q[1];
x q[2];
rz(0.036582935) q[3];
sx q[3];
rz(-2.7162093) q[3];
sx q[3];
rz(0.98664397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4778135) q[2];
sx q[2];
rz(-2.9728643) q[2];
sx q[2];
rz(-1.7185877) q[2];
rz(-2.6723828) q[3];
sx q[3];
rz(-1.4956632) q[3];
sx q[3];
rz(2.7928228) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.19942) q[0];
sx q[0];
rz(-1.9439789) q[0];
sx q[0];
rz(2.0592101) q[0];
rz(2.9669145) q[1];
sx q[1];
rz(-1.382261) q[1];
sx q[1];
rz(-2.6686525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2286005) q[0];
sx q[0];
rz(-0.45683041) q[0];
sx q[0];
rz(-1.940371) q[0];
rz(-pi) q[1];
rz(2.2720415) q[2];
sx q[2];
rz(-1.2954748) q[2];
sx q[2];
rz(2.4710831) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14607695) q[1];
sx q[1];
rz(-2.5716166) q[1];
sx q[1];
rz(1.6126627) q[1];
rz(-pi) q[2];
rz(0.49854203) q[3];
sx q[3];
rz(-0.99851552) q[3];
sx q[3];
rz(1.7576493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6579154) q[2];
sx q[2];
rz(-1.4853073) q[2];
sx q[2];
rz(-1.8243054) q[2];
rz(2.6910736) q[3];
sx q[3];
rz(-2.2299288) q[3];
sx q[3];
rz(-1.7821504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6555742) q[0];
sx q[0];
rz(-1.065932) q[0];
sx q[0];
rz(0.67147142) q[0];
rz(2.2010522) q[1];
sx q[1];
rz(-0.42010072) q[1];
sx q[1];
rz(-2.1817575) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80653866) q[0];
sx q[0];
rz(-1.1919406) q[0];
sx q[0];
rz(-0.23083487) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7966389) q[2];
sx q[2];
rz(-2.2246771) q[2];
sx q[2];
rz(-2.8381062) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6558608) q[1];
sx q[1];
rz(-1.277024) q[1];
sx q[1];
rz(1.345977) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6790326) q[3];
sx q[3];
rz(-2.0104791) q[3];
sx q[3];
rz(-2.3236583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.81022108) q[2];
sx q[2];
rz(-1.081531) q[2];
sx q[2];
rz(0.01037154) q[2];
rz(1.8593908) q[3];
sx q[3];
rz(-0.60716647) q[3];
sx q[3];
rz(-2.8289294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4717167) q[0];
sx q[0];
rz(-0.80729055) q[0];
sx q[0];
rz(0.5759936) q[0];
rz(-2.5895789) q[1];
sx q[1];
rz(-1.8355398) q[1];
sx q[1];
rz(-2.2214417) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33951633) q[0];
sx q[0];
rz(-0.9528725) q[0];
sx q[0];
rz(-2.3137388) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29187496) q[2];
sx q[2];
rz(-1.5988962) q[2];
sx q[2];
rz(2.4296938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.53729) q[1];
sx q[1];
rz(-1.3370945) q[1];
sx q[1];
rz(-1.9275194) q[1];
rz(3.017265) q[3];
sx q[3];
rz(-1.6724068) q[3];
sx q[3];
rz(-1.4426444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0803904) q[2];
sx q[2];
rz(-1.7549425) q[2];
sx q[2];
rz(2.6913255) q[2];
rz(0.14084147) q[3];
sx q[3];
rz(-1.2248657) q[3];
sx q[3];
rz(-2.5950477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30215728) q[0];
sx q[0];
rz(-1.7385229) q[0];
sx q[0];
rz(2.9490525) q[0];
rz(-0.64388609) q[1];
sx q[1];
rz(-1.5637249) q[1];
sx q[1];
rz(0.083316915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24593943) q[0];
sx q[0];
rz(-0.97165474) q[0];
sx q[0];
rz(2.9532632) q[0];
x q[1];
rz(-2.1716921) q[2];
sx q[2];
rz(-2.9978581) q[2];
sx q[2];
rz(-0.78605762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51074886) q[1];
sx q[1];
rz(-0.35707316) q[1];
sx q[1];
rz(-1.0624753) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7658773) q[3];
sx q[3];
rz(-1.8109115) q[3];
sx q[3];
rz(2.1841336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4087499) q[2];
sx q[2];
rz(-2.7538444) q[2];
sx q[2];
rz(-2.5852618) q[2];
rz(-2.2522816) q[3];
sx q[3];
rz(-1.6272864) q[3];
sx q[3];
rz(0.48243943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0675875) q[0];
sx q[0];
rz(-0.50453019) q[0];
sx q[0];
rz(2.0489847) q[0];
rz(-0.69674528) q[1];
sx q[1];
rz(-2.4431591) q[1];
sx q[1];
rz(2.4901857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073805907) q[0];
sx q[0];
rz(-1.7346364) q[0];
sx q[0];
rz(2.4685173) q[0];
rz(-pi) q[1];
rz(2.7138804) q[2];
sx q[2];
rz(-2.1943478) q[2];
sx q[2];
rz(-1.9643188) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5696545) q[1];
sx q[1];
rz(-1.2982315) q[1];
sx q[1];
rz(1.0054719) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1262195) q[3];
sx q[3];
rz(-1.6241933) q[3];
sx q[3];
rz(-1.3863409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50218987) q[2];
sx q[2];
rz(-0.71523372) q[2];
sx q[2];
rz(-2.5779842) q[2];
rz(0.11858502) q[3];
sx q[3];
rz(-1.010681) q[3];
sx q[3];
rz(-2.3732869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80657715) q[0];
sx q[0];
rz(-1.7105569) q[0];
sx q[0];
rz(3.0203876) q[0];
rz(-2.99446) q[1];
sx q[1];
rz(-1.9857261) q[1];
sx q[1];
rz(-2.3187231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74350243) q[0];
sx q[0];
rz(-1.8202298) q[0];
sx q[0];
rz(1.6642847) q[0];
rz(1.5414833) q[2];
sx q[2];
rz(-2.1357059) q[2];
sx q[2];
rz(-0.82290111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40892021) q[1];
sx q[1];
rz(-2.3358445) q[1];
sx q[1];
rz(-2.9288956) q[1];
rz(3.0599252) q[3];
sx q[3];
rz(-1.3586391) q[3];
sx q[3];
rz(1.798925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0269028) q[2];
sx q[2];
rz(-0.59640408) q[2];
sx q[2];
rz(-0.58843311) q[2];
rz(2.9450997) q[3];
sx q[3];
rz(-1.4339707) q[3];
sx q[3];
rz(0.65832037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7115985) q[0];
sx q[0];
rz(-0.64084941) q[0];
sx q[0];
rz(-2.1852921) q[0];
rz(-0.20005964) q[1];
sx q[1];
rz(-1.450489) q[1];
sx q[1];
rz(-0.53093451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5945608) q[0];
sx q[0];
rz(-0.49358019) q[0];
sx q[0];
rz(-1.7112687) q[0];
rz(2.7199581) q[2];
sx q[2];
rz(-2.4608825) q[2];
sx q[2];
rz(-1.9599078) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0772452) q[1];
sx q[1];
rz(-0.56659154) q[1];
sx q[1];
rz(-2.8342325) q[1];
x q[2];
rz(2.8538029) q[3];
sx q[3];
rz(-2.6734201) q[3];
sx q[3];
rz(-0.60540199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42813534) q[2];
sx q[2];
rz(-2.0811446) q[2];
sx q[2];
rz(-1.1327845) q[2];
rz(-0.97635859) q[3];
sx q[3];
rz(-2.4553757) q[3];
sx q[3];
rz(1.2767701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79162947) q[0];
sx q[0];
rz(-3.0296453) q[0];
sx q[0];
rz(-1.5891225) q[0];
rz(2.7516229) q[1];
sx q[1];
rz(-1.5733122) q[1];
sx q[1];
rz(2.3812297) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0438662) q[0];
sx q[0];
rz(-1.1047433) q[0];
sx q[0];
rz(1.4938484) q[0];
rz(-pi) q[1];
rz(-0.46192034) q[2];
sx q[2];
rz(-1.5452146) q[2];
sx q[2];
rz(1.8476576) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8514357) q[1];
sx q[1];
rz(-0.37237114) q[1];
sx q[1];
rz(-1.464616) q[1];
x q[2];
rz(-2.1093587) q[3];
sx q[3];
rz(-1.8377152) q[3];
sx q[3];
rz(2.0570786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4752263) q[2];
sx q[2];
rz(-2.1658289) q[2];
sx q[2];
rz(-0.59398061) q[2];
rz(-2.5178759) q[3];
sx q[3];
rz(-1.2902322) q[3];
sx q[3];
rz(-0.82998961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8754616) q[0];
sx q[0];
rz(-3.0103191) q[0];
sx q[0];
rz(-1.3651146) q[0];
rz(0.020137067) q[1];
sx q[1];
rz(-0.29300856) q[1];
sx q[1];
rz(-1.551052) q[1];
rz(3.1333078) q[2];
sx q[2];
rz(-1.9556254) q[2];
sx q[2];
rz(0.34224184) q[2];
rz(2.9659553) q[3];
sx q[3];
rz(-1.5541811) q[3];
sx q[3];
rz(0.27754996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
