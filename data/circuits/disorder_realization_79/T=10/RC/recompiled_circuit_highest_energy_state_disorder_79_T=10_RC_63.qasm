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
rz(0.41843721) q[0];
sx q[0];
rz(-0.96324459) q[0];
sx q[0];
rz(0.20382398) q[0];
rz(-2.0889497) q[1];
sx q[1];
rz(-0.79336762) q[1];
sx q[1];
rz(1.6066983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8730072) q[0];
sx q[0];
rz(-2.632336) q[0];
sx q[0];
rz(2.068822) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34687545) q[2];
sx q[2];
rz(-2.2636608) q[2];
sx q[2];
rz(1.3077867) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71434778) q[1];
sx q[1];
rz(-0.74810076) q[1];
sx q[1];
rz(2.1022878) q[1];
x q[2];
rz(0.20333692) q[3];
sx q[3];
rz(-1.7032188) q[3];
sx q[3];
rz(1.9258969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7242929) q[2];
sx q[2];
rz(-2.5827926) q[2];
sx q[2];
rz(0.99754769) q[2];
rz(-0.7558465) q[3];
sx q[3];
rz(-1.7852802) q[3];
sx q[3];
rz(1.0409748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34870979) q[0];
sx q[0];
rz(-0.097276874) q[0];
sx q[0];
rz(1.2874228) q[0];
rz(2.6584794) q[1];
sx q[1];
rz(-1.0419507) q[1];
sx q[1];
rz(-1.9226673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3588803) q[0];
sx q[0];
rz(-0.25280372) q[0];
sx q[0];
rz(0.88835277) q[0];
rz(-0.27147175) q[2];
sx q[2];
rz(-2.4515332) q[2];
sx q[2];
rz(0.59218237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87172958) q[1];
sx q[1];
rz(-0.31530373) q[1];
sx q[1];
rz(2.053715) q[1];
x q[2];
rz(-2.2628701) q[3];
sx q[3];
rz(-1.847953) q[3];
sx q[3];
rz(-0.73776484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37811849) q[2];
sx q[2];
rz(-3.0749622) q[2];
sx q[2];
rz(2.5197869) q[2];
rz(-2.5032737) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(2.2973255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062375) q[0];
sx q[0];
rz(-0.64202809) q[0];
sx q[0];
rz(0.21632347) q[0];
rz(2.5430039) q[1];
sx q[1];
rz(-1.8363771) q[1];
sx q[1];
rz(-0.52282202) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.110958) q[0];
sx q[0];
rz(-1.2485663) q[0];
sx q[0];
rz(0.3229753) q[0];
x q[1];
rz(-1.4141358) q[2];
sx q[2];
rz(-1.1839424) q[2];
sx q[2];
rz(2.000282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93980689) q[1];
sx q[1];
rz(-2.5720189) q[1];
sx q[1];
rz(-0.42167432) q[1];
x q[2];
rz(2.3159695) q[3];
sx q[3];
rz(-0.43681991) q[3];
sx q[3];
rz(-2.2568011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9598976) q[2];
sx q[2];
rz(-1.8671702) q[2];
sx q[2];
rz(-1.6240906) q[2];
rz(2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(-1.6200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14438039) q[0];
sx q[0];
rz(-0.388044) q[0];
sx q[0];
rz(1.5390747) q[0];
rz(-1.0333565) q[1];
sx q[1];
rz(-0.76834232) q[1];
sx q[1];
rz(-1.296952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.813445) q[0];
sx q[0];
rz(-1.7792276) q[0];
sx q[0];
rz(0.84020241) q[0];
rz(-pi) q[1];
rz(-0.06242604) q[2];
sx q[2];
rz(-0.87597825) q[2];
sx q[2];
rz(2.7165987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.32758157) q[1];
sx q[1];
rz(-1.9821229) q[1];
sx q[1];
rz(-2.7279305) q[1];
rz(-pi) q[2];
rz(1.2554712) q[3];
sx q[3];
rz(-1.0389757) q[3];
sx q[3];
rz(-1.9813615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(2.0959334) q[2];
rz(0.21444923) q[3];
sx q[3];
rz(-0.72434536) q[3];
sx q[3];
rz(1.1469871) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3359208) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(2.6222498) q[0];
rz(2.1995811) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(-1.6090144) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2359152) q[0];
sx q[0];
rz(-2.2122243) q[0];
sx q[0];
rz(2.0302613) q[0];
rz(-pi) q[1];
rz(-0.66560676) q[2];
sx q[2];
rz(-0.22826787) q[2];
sx q[2];
rz(0.98787243) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0253801) q[1];
sx q[1];
rz(-2.0713701) q[1];
sx q[1];
rz(0.0059555014) q[1];
rz(-pi) q[2];
rz(-0.71096731) q[3];
sx q[3];
rz(-2.9971854) q[3];
sx q[3];
rz(-2.292423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.64451009) q[2];
sx q[2];
rz(-0.77748674) q[2];
sx q[2];
rz(2.8572148) q[2];
rz(2.583368) q[3];
sx q[3];
rz(-1.7773209) q[3];
sx q[3];
rz(0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068950653) q[0];
sx q[0];
rz(-0.41256368) q[0];
sx q[0];
rz(-0.64055842) q[0];
rz(1.6259646) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(-2.9277149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17992526) q[0];
sx q[0];
rz(-2.7583987) q[0];
sx q[0];
rz(-0.92143329) q[0];
x q[1];
rz(0.40753813) q[2];
sx q[2];
rz(-2.354051) q[2];
sx q[2];
rz(1.0223197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2474413) q[1];
sx q[1];
rz(-0.94235984) q[1];
sx q[1];
rz(3.0559866) q[1];
rz(-pi) q[2];
rz(-2.0908666) q[3];
sx q[3];
rz(-1.5087391) q[3];
sx q[3];
rz(2.7620535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-0.48206097) q[2];
sx q[2];
rz(0.17769979) q[2];
rz(1.7443582) q[3];
sx q[3];
rz(-1.1763562) q[3];
sx q[3];
rz(-0.41745225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419256) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(-2.0360816) q[0];
rz(0.28327495) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(1.4615321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7532336) q[0];
sx q[0];
rz(-0.93030158) q[0];
sx q[0];
rz(0.70029135) q[0];
x q[1];
rz(-2.2143557) q[2];
sx q[2];
rz(-2.3935648) q[2];
sx q[2];
rz(3.1035556) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0431598) q[1];
sx q[1];
rz(-0.57505703) q[1];
sx q[1];
rz(2.9270323) q[1];
x q[2];
rz(0.47847943) q[3];
sx q[3];
rz(-2.5095486) q[3];
sx q[3];
rz(1.5173591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0561515) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(-0.20480569) q[2];
rz(-1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7751223) q[0];
sx q[0];
rz(-2.6752052) q[0];
sx q[0];
rz(3.108016) q[0];
rz(-1.6665943) q[1];
sx q[1];
rz(-2.3971403) q[1];
sx q[1];
rz(1.144145) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5789216) q[0];
sx q[0];
rz(-1.5961338) q[0];
sx q[0];
rz(2.6549005) q[0];
rz(1.8587684) q[2];
sx q[2];
rz(-2.8045553) q[2];
sx q[2];
rz(-2.2971643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.856196) q[1];
sx q[1];
rz(-2.2335529) q[1];
sx q[1];
rz(-1.0876571) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89544483) q[3];
sx q[3];
rz(-0.43206462) q[3];
sx q[3];
rz(0.91401446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7555776) q[2];
sx q[2];
rz(-0.80499804) q[2];
sx q[2];
rz(-1.2516652) q[2];
rz(-1.3517316) q[3];
sx q[3];
rz(-2.2224865) q[3];
sx q[3];
rz(-2.1613817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1397454) q[0];
sx q[0];
rz(-0.63848764) q[0];
sx q[0];
rz(-1.8852604) q[0];
rz(3.046335) q[1];
sx q[1];
rz(-0.88905159) q[1];
sx q[1];
rz(-2.6108066) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19671104) q[0];
sx q[0];
rz(-1.3232348) q[0];
sx q[0];
rz(2.976368) q[0];
rz(-pi) q[1];
rz(-1.6094535) q[2];
sx q[2];
rz(-1.2529904) q[2];
sx q[2];
rz(-3.0927049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6775178) q[1];
sx q[1];
rz(-1.6987659) q[1];
sx q[1];
rz(-0.24125464) q[1];
rz(-0.23654273) q[3];
sx q[3];
rz(-0.45748392) q[3];
sx q[3];
rz(1.903423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(-2.6007593) q[2];
rz(2.3234308) q[3];
sx q[3];
rz(-1.4427789) q[3];
sx q[3];
rz(-2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.624991) q[0];
sx q[0];
rz(-1.3636959) q[0];
sx q[0];
rz(0.3987819) q[0];
rz(-1.3840236) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(-3.0922281) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.787034) q[0];
sx q[0];
rz(-1.6794717) q[0];
sx q[0];
rz(-3.0804481) q[0];
rz(-pi) q[1];
rz(1.9513543) q[2];
sx q[2];
rz(-0.28186381) q[2];
sx q[2];
rz(-2.6388002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4409742) q[1];
sx q[1];
rz(-1.6153533) q[1];
sx q[1];
rz(-0.99127165) q[1];
x q[2];
rz(-2.5096941) q[3];
sx q[3];
rz(-1.4057341) q[3];
sx q[3];
rz(-2.2483605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12386879) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(1.5105985) q[2];
rz(0.92783582) q[3];
sx q[3];
rz(-1.3414914) q[3];
sx q[3];
rz(-1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31502003) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(1.0646959) q[1];
sx q[1];
rz(-2.6152492) q[1];
sx q[1];
rz(-1.1660887) q[1];
rz(-0.76416107) q[2];
sx q[2];
rz(-2.4780826) q[2];
sx q[2];
rz(1.6586951) q[2];
rz(2.1178653) q[3];
sx q[3];
rz(-1.8454875) q[3];
sx q[3];
rz(-1.9133205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
