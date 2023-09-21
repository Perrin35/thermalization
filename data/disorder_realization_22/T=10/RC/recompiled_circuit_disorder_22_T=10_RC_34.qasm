OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(-2.9971478) q[0];
rz(-2.5748409) q[1];
sx q[1];
rz(-2.6161939) q[1];
sx q[1];
rz(2.1638343) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9508096) q[0];
sx q[0];
rz(-1.8291744) q[0];
sx q[0];
rz(1.6590614) q[0];
rz(-2.0852282) q[2];
sx q[2];
rz(-1.3008413) q[2];
sx q[2];
rz(-1.1411238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3111167) q[1];
sx q[1];
rz(-1.7101354) q[1];
sx q[1];
rz(1.7140825) q[1];
rz(1.3917189) q[3];
sx q[3];
rz(-0.79154166) q[3];
sx q[3];
rz(0.66034987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-0.93227512) q[2];
rz(-0.19876924) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(-0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7070049) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(2.7804651) q[0];
rz(1.7065642) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(0.8180058) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60048238) q[0];
sx q[0];
rz(-0.4821018) q[0];
sx q[0];
rz(2.5151398) q[0];
rz(-pi) q[1];
rz(1.1992707) q[2];
sx q[2];
rz(-1.6852334) q[2];
sx q[2];
rz(-0.83100806) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86707838) q[1];
sx q[1];
rz(-2.3550905) q[1];
sx q[1];
rz(-2.2134476) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4217671) q[3];
sx q[3];
rz(-0.36747284) q[3];
sx q[3];
rz(3.0847103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8890185) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(1.4206295) q[2];
rz(-1.8255) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(-0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(1.7398359) q[0];
sx q[0];
rz(-2.6494884) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(0.31618205) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(-0.2972163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0270099) q[0];
sx q[0];
rz(-1.5073338) q[0];
sx q[0];
rz(2.7364536) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3702413) q[2];
sx q[2];
rz(-1.0372247) q[2];
sx q[2];
rz(-0.16650621) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9416562) q[1];
sx q[1];
rz(-1.7574649) q[1];
sx q[1];
rz(-3.093064) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3171595) q[3];
sx q[3];
rz(-1.0472877) q[3];
sx q[3];
rz(-1.5200966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3729942) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-0.42759839) q[2];
rz(1.9528495) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(-0.52743131) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4363842) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(0.71290839) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(2.4598222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9457152) q[0];
sx q[0];
rz(-1.3099075) q[0];
sx q[0];
rz(2.9257141) q[0];
x q[1];
rz(-0.21708023) q[2];
sx q[2];
rz(-2.1639369) q[2];
sx q[2];
rz(2.8341688) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53182787) q[1];
sx q[1];
rz(-2.1064639) q[1];
sx q[1];
rz(-0.51171724) q[1];
rz(-1.2069615) q[3];
sx q[3];
rz(-1.893265) q[3];
sx q[3];
rz(2.913167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(2.183389) q[2];
rz(2.0751674) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(-2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(0.56030309) q[0];
rz(-0.99984461) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(-1.6220185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0409659) q[0];
sx q[0];
rz(-1.1275516) q[0];
sx q[0];
rz(-2.6020781) q[0];
x q[1];
rz(-2.1278473) q[2];
sx q[2];
rz(-1.0336913) q[2];
sx q[2];
rz(0.86442664) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49415627) q[1];
sx q[1];
rz(-1.372822) q[1];
sx q[1];
rz(-2.3784749) q[1];
rz(-pi) q[2];
rz(2.1836957) q[3];
sx q[3];
rz(-1.4688204) q[3];
sx q[3];
rz(-2.0436055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4524298) q[2];
sx q[2];
rz(-2.5897557) q[2];
sx q[2];
rz(2.9186644) q[2];
rz(-0.034742268) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8054304) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(-1.0700595) q[0];
rz(-1.3609715) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(1.8575352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491935) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(-2.4248289) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9610923) q[2];
sx q[2];
rz(-1.830606) q[2];
sx q[2];
rz(2.5026472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2388873) q[1];
sx q[1];
rz(-2.0327912) q[1];
sx q[1];
rz(-1.8932896) q[1];
rz(2.0606023) q[3];
sx q[3];
rz(-0.59270699) q[3];
sx q[3];
rz(0.61471516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47508919) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(2.3049138) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(-1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5752983) q[0];
sx q[0];
rz(-0.42995444) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(-0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(-0.93820757) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.339401) q[0];
sx q[0];
rz(-1.2656478) q[0];
sx q[0];
rz(1.9992273) q[0];
rz(1.0023408) q[2];
sx q[2];
rz(-0.92407862) q[2];
sx q[2];
rz(-1.2127753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92749121) q[1];
sx q[1];
rz(-2.2117105) q[1];
sx q[1];
rz(0.44968857) q[1];
x q[2];
rz(0.70721831) q[3];
sx q[3];
rz(-1.3361317) q[3];
sx q[3];
rz(-2.9747687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2618835) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(1.3860469) q[2];
rz(1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(1.8677615) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(-0.83126718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76038218) q[0];
sx q[0];
rz(-1.4236649) q[0];
sx q[0];
rz(0.40572625) q[0];
x q[1];
rz(1.3966884) q[2];
sx q[2];
rz(-0.6422407) q[2];
sx q[2];
rz(-2.1070534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.23849328) q[1];
sx q[1];
rz(-1.0817173) q[1];
sx q[1];
rz(0.78519435) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7940815) q[3];
sx q[3];
rz(-1.3459473) q[3];
sx q[3];
rz(2.4527578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5635809) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(-1.9780805) q[3];
sx q[3];
rz(-0.61146277) q[3];
sx q[3];
rz(0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.7678541) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(0.97737616) q[0];
rz(-1.7550229) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(1.1057373) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2290105) q[0];
sx q[0];
rz(-1.6874896) q[0];
sx q[0];
rz(0.44793655) q[0];
rz(2.2374723) q[2];
sx q[2];
rz(-1.4881926) q[2];
sx q[2];
rz(-2.3711575) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92377201) q[1];
sx q[1];
rz(-0.2982699) q[1];
sx q[1];
rz(0.89426269) q[1];
rz(-2.7902778) q[3];
sx q[3];
rz(-1.8422541) q[3];
sx q[3];
rz(2.9815477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.608312) q[2];
sx q[2];
rz(-0.38485843) q[2];
sx q[2];
rz(-2.9917955) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(-3.0143484) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-0.36247411) q[1];
sx q[1];
rz(2.9737934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0906592) q[0];
sx q[0];
rz(-1.5516084) q[0];
sx q[0];
rz(0.030254342) q[0];
rz(-pi) q[1];
rz(2.9211505) q[2];
sx q[2];
rz(-0.61969212) q[2];
sx q[2];
rz(2.9881791) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4976616) q[1];
sx q[1];
rz(-1.5639468) q[1];
sx q[1];
rz(-2.5531406) q[1];
x q[2];
rz(0.12331788) q[3];
sx q[3];
rz(-1.4535558) q[3];
sx q[3];
rz(0.58837147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(-3.051565) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54820838) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(-0.38800115) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(-2.491647) q[2];
sx q[2];
rz(-0.92234281) q[2];
sx q[2];
rz(-2.6364706) q[2];
rz(-1.0317867) q[3];
sx q[3];
rz(-1.7179334) q[3];
sx q[3];
rz(-2.4739305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
