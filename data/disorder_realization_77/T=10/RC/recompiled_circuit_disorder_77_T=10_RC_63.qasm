OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4958772) q[0];
sx q[0];
rz(-1.0951395) q[0];
sx q[0];
rz(2.9021184) q[0];
rz(-pi) q[1];
rz(0.92476966) q[2];
sx q[2];
rz(-1.3999108) q[2];
sx q[2];
rz(0.31121635) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9819298) q[1];
sx q[1];
rz(-2.5291981) q[1];
sx q[1];
rz(1.0931404) q[1];
x q[2];
rz(1.6151186) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(-1.1322024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(1.8288076) q[0];
rz(2.9361172) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(-1.9899433) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32423702) q[0];
sx q[0];
rz(-1.865987) q[0];
sx q[0];
rz(-2.2833707) q[0];
rz(-pi) q[1];
rz(-2.5580514) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(-1.5577424) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0484867) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(2.6622245) q[1];
rz(-0.50206708) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(-1.7934007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0097222086) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(-2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7650334) q[0];
sx q[0];
rz(-0.74900904) q[0];
sx q[0];
rz(2.2545933) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25643202) q[2];
sx q[2];
rz(-1.7000546) q[2];
sx q[2];
rz(1.749922) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1287071) q[1];
sx q[1];
rz(-1.7696847) q[1];
sx q[1];
rz(-0.30602869) q[1];
rz(-pi) q[2];
rz(0.64485456) q[3];
sx q[3];
rz(-1.8861024) q[3];
sx q[3];
rz(0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(0.55666322) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-1.042897) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2264003) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-2.6779968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.050042) q[0];
sx q[0];
rz(-2.3766962) q[0];
sx q[0];
rz(2.0043623) q[0];
rz(1.5092588) q[2];
sx q[2];
rz(-1.2386285) q[2];
sx q[2];
rz(-0.57519826) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3568748) q[1];
sx q[1];
rz(-0.80662913) q[1];
sx q[1];
rz(2.9033317) q[1];
x q[2];
rz(-0.87426825) q[3];
sx q[3];
rz(-0.6797176) q[3];
sx q[3];
rz(-0.10248871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(-0.4822576) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-0.94435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807242) q[0];
sx q[0];
rz(-1.7797911) q[0];
sx q[0];
rz(3.0955549) q[0];
rz(2.7371251) q[2];
sx q[2];
rz(-2.3239115) q[2];
sx q[2];
rz(0.58194619) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22630616) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(1.5096942) q[1];
x q[2];
rz(2.7086908) q[3];
sx q[3];
rz(-1.8042943) q[3];
sx q[3];
rz(1.3732861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.258761) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-0.8862409) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-3.086673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9960105) q[0];
sx q[0];
rz(-1.3025563) q[0];
sx q[0];
rz(-0.085573816) q[0];
rz(-pi) q[1];
rz(-3.1228742) q[2];
sx q[2];
rz(-1.8939549) q[2];
sx q[2];
rz(-1.9402372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1850486) q[1];
sx q[1];
rz(-1.1951606) q[1];
sx q[1];
rz(0.59021414) q[1];
rz(-pi) q[2];
rz(2.6782126) q[3];
sx q[3];
rz(-2.1043092) q[3];
sx q[3];
rz(2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-2.1248655) q[2];
rz(-2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.8090766) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(2.231266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.401424) q[0];
sx q[0];
rz(-1.6356042) q[0];
sx q[0];
rz(-0.054697371) q[0];
x q[1];
rz(-2.2482713) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(-2.233778) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42156223) q[1];
sx q[1];
rz(-1.3898464) q[1];
sx q[1];
rz(0.49444316) q[1];
rz(-pi) q[2];
rz(0.72083731) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(0.37973675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(-2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(-1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-0.67665726) q[0];
sx q[0];
rz(-3.0122053) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-2.8410889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26586543) q[0];
sx q[0];
rz(-2.2458796) q[0];
sx q[0];
rz(0.48203326) q[0];
x q[1];
rz(-0.76086107) q[2];
sx q[2];
rz(-1.0957452) q[2];
sx q[2];
rz(-2.8617815) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.064425163) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(0.47972958) q[1];
rz(-pi) q[2];
x q[2];
rz(2.825533) q[3];
sx q[3];
rz(-2.3028767) q[3];
sx q[3];
rz(-0.72533208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(2.3596181) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(0.12776275) q[0];
rz(0.54221517) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-0.75884563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7153873) q[0];
sx q[0];
rz(-1.9949159) q[0];
sx q[0];
rz(-3.12294) q[0];
rz(-1.7636289) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(-0.45229518) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5938877) q[1];
sx q[1];
rz(-0.87915671) q[1];
sx q[1];
rz(-1.7734852) q[1];
rz(-pi) q[2];
rz(-3.127029) q[3];
sx q[3];
rz(-2.2363538) q[3];
sx q[3];
rz(-1.2138106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(1.7193433) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35995099) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(3.066257) q[0];
rz(0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(2.5316701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1760575) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(0.39359351) q[0];
rz(2.7650325) q[2];
sx q[2];
rz(-1.3137523) q[2];
sx q[2];
rz(-3.0388447) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1249485) q[1];
sx q[1];
rz(-2.2955107) q[1];
sx q[1];
rz(-1.2490586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53819733) q[3];
sx q[3];
rz(-2.3231069) q[3];
sx q[3];
rz(1.8173816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(1.7238293) q[2];
sx q[2];
rz(-0.19822181) q[2];
sx q[2];
rz(-0.76186686) q[2];
rz(-2.0416904) q[3];
sx q[3];
rz(-1.3445911) q[3];
sx q[3];
rz(0.54683987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
