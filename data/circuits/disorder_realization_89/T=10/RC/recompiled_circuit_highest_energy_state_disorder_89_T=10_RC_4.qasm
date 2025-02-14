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
rz(-3.0638679) q[0];
sx q[0];
rz(-1.024615) q[0];
sx q[0];
rz(1.2999363) q[0];
rz(2.6990702) q[1];
sx q[1];
rz(4.2808851) q[1];
sx q[1];
rz(12.130375) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090442) q[0];
sx q[0];
rz(-1.6761177) q[0];
sx q[0];
rz(1.9645343) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7828658) q[2];
sx q[2];
rz(-0.84237885) q[2];
sx q[2];
rz(-1.6747507) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30937815) q[1];
sx q[1];
rz(-1.9736027) q[1];
sx q[1];
rz(1.554053) q[1];
rz(-pi) q[2];
rz(0.12650872) q[3];
sx q[3];
rz(-0.96130575) q[3];
sx q[3];
rz(-1.3550657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.062833) q[2];
sx q[2];
rz(-1.6597513) q[2];
sx q[2];
rz(0.01595846) q[2];
rz(1.6503085) q[3];
sx q[3];
rz(-1.242638) q[3];
sx q[3];
rz(-2.7119467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780739) q[0];
sx q[0];
rz(-2.0500545) q[0];
sx q[0];
rz(1.3462521) q[0];
rz(-1.9740055) q[1];
sx q[1];
rz(-1.0941894) q[1];
sx q[1];
rz(1.8928554) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561476) q[0];
sx q[0];
rz(-2.1753575) q[0];
sx q[0];
rz(-2.6660835) q[0];
x q[1];
rz(-0.65950583) q[2];
sx q[2];
rz(-2.9444866) q[2];
sx q[2];
rz(1.3172305) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4521422) q[1];
sx q[1];
rz(-0.4045338) q[1];
sx q[1];
rz(-1.4305315) q[1];
x q[2];
rz(-0.98442673) q[3];
sx q[3];
rz(-1.8674486) q[3];
sx q[3];
rz(-0.082060952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1355797) q[2];
sx q[2];
rz(-2.4485782) q[2];
sx q[2];
rz(-2.6683624) q[2];
rz(-3.0114975) q[3];
sx q[3];
rz(-1.3869163) q[3];
sx q[3];
rz(0.010802833) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8148282) q[0];
sx q[0];
rz(-0.15693754) q[0];
sx q[0];
rz(2.9275295) q[0];
rz(-1.3940943) q[1];
sx q[1];
rz(-0.97624818) q[1];
sx q[1];
rz(-3.0551547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6753767) q[0];
sx q[0];
rz(-1.4450184) q[0];
sx q[0];
rz(1.2745538) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8930703) q[2];
sx q[2];
rz(-2.0166335) q[2];
sx q[2];
rz(-2.0634212) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6896539) q[1];
sx q[1];
rz(-1.9898333) q[1];
sx q[1];
rz(-1.7614014) q[1];
x q[2];
rz(-2.2312715) q[3];
sx q[3];
rz(-2.3757138) q[3];
sx q[3];
rz(-0.41492763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.67132407) q[2];
sx q[2];
rz(-0.9321804) q[2];
sx q[2];
rz(-2.0051125) q[2];
rz(1.350435) q[3];
sx q[3];
rz(-1.330749) q[3];
sx q[3];
rz(-2.7854846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39373028) q[0];
sx q[0];
rz(-0.50734729) q[0];
sx q[0];
rz(0.90079975) q[0];
rz(-1.2334709) q[1];
sx q[1];
rz(-0.8546468) q[1];
sx q[1];
rz(2.5305117) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7584832) q[0];
sx q[0];
rz(-2.0462667) q[0];
sx q[0];
rz(0.0080541797) q[0];
rz(-0.22902352) q[2];
sx q[2];
rz(-1.5364858) q[2];
sx q[2];
rz(-1.9051617) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.722695) q[1];
sx q[1];
rz(-1.1686106) q[1];
sx q[1];
rz(0.45940347) q[1];
rz(-pi) q[2];
rz(2.7902725) q[3];
sx q[3];
rz(-1.6804983) q[3];
sx q[3];
rz(0.20449311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5078807) q[2];
sx q[2];
rz(-0.66358006) q[2];
sx q[2];
rz(-0.12082417) q[2];
rz(3.1145596) q[3];
sx q[3];
rz(-0.20172541) q[3];
sx q[3];
rz(0.87593186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9487069) q[0];
sx q[0];
rz(-0.8256194) q[0];
sx q[0];
rz(1.3407619) q[0];
rz(1.5277398) q[1];
sx q[1];
rz(-1.3132881) q[1];
sx q[1];
rz(-0.54940474) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2159113) q[0];
sx q[0];
rz(-1.4936555) q[0];
sx q[0];
rz(-0.74621426) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1099186) q[2];
sx q[2];
rz(-2.1018545) q[2];
sx q[2];
rz(2.2344294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5633302) q[1];
sx q[1];
rz(-1.6206998) q[1];
sx q[1];
rz(1.084855) q[1];
rz(-2.0321376) q[3];
sx q[3];
rz(-1.6008988) q[3];
sx q[3];
rz(2.1738659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.01650979) q[2];
sx q[2];
rz(-1.5089401) q[2];
sx q[2];
rz(2.5588918) q[2];
rz(1.6216283) q[3];
sx q[3];
rz(-2.426332) q[3];
sx q[3];
rz(1.3165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.3358066) q[0];
sx q[0];
rz(-1.0842706) q[0];
sx q[0];
rz(-1.8871319) q[0];
rz(1.1955903) q[1];
sx q[1];
rz(-1.4072199) q[1];
sx q[1];
rz(1.5171299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82255581) q[0];
sx q[0];
rz(-1.5040888) q[0];
sx q[0];
rz(3.1298679) q[0];
rz(-0.80251383) q[2];
sx q[2];
rz(-1.1468059) q[2];
sx q[2];
rz(2.9579332) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.62134777) q[1];
sx q[1];
rz(-0.98437446) q[1];
sx q[1];
rz(-1.919073) q[1];
rz(0.020217309) q[3];
sx q[3];
rz(-0.51512326) q[3];
sx q[3];
rz(-2.8216854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.967041) q[2];
sx q[2];
rz(-1.6074564) q[2];
sx q[2];
rz(-3.0957481) q[2];
rz(0.52715078) q[3];
sx q[3];
rz(-1.0322626) q[3];
sx q[3];
rz(0.84120685) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884278) q[0];
sx q[0];
rz(-0.66075745) q[0];
sx q[0];
rz(-0.38594693) q[0];
rz(1.376232) q[1];
sx q[1];
rz(-1.0877345) q[1];
sx q[1];
rz(-1.2765346) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46158591) q[0];
sx q[0];
rz(-1.9421845) q[0];
sx q[0];
rz(0.71457926) q[0];
rz(-pi) q[1];
rz(3.0936095) q[2];
sx q[2];
rz(-1.1986102) q[2];
sx q[2];
rz(1.3822777) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1305728) q[1];
sx q[1];
rz(-1.6064732) q[1];
sx q[1];
rz(-0.0079572268) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5064729) q[3];
sx q[3];
rz(-1.606588) q[3];
sx q[3];
rz(-2.2944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4253) q[2];
sx q[2];
rz(-1.7159117) q[2];
sx q[2];
rz(1.7745793) q[2];
rz(3.0573209) q[3];
sx q[3];
rz(-1.1140946) q[3];
sx q[3];
rz(-0.082898609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0609584) q[0];
sx q[0];
rz(-0.11140379) q[0];
sx q[0];
rz(-1.8633307) q[0];
rz(3.0807965) q[1];
sx q[1];
rz(-0.95525974) q[1];
sx q[1];
rz(-1.0999058) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9644224) q[0];
sx q[0];
rz(-1.1905429) q[0];
sx q[0];
rz(0.1582665) q[0];
x q[1];
rz(0.75662778) q[2];
sx q[2];
rz(-2.1887458) q[2];
sx q[2];
rz(2.6709887) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5085735) q[1];
sx q[1];
rz(-0.46357511) q[1];
sx q[1];
rz(0.68008382) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61366365) q[3];
sx q[3];
rz(-1.9983158) q[3];
sx q[3];
rz(0.18117426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3595769) q[2];
sx q[2];
rz(-0.64543739) q[2];
sx q[2];
rz(-2.6673356) q[2];
rz(2.5931902) q[3];
sx q[3];
rz(-0.67265284) q[3];
sx q[3];
rz(1.7614346) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7062374) q[0];
sx q[0];
rz(-2.7147003) q[0];
sx q[0];
rz(-1.9306345) q[0];
rz(-1.2779166) q[1];
sx q[1];
rz(-1.5958818) q[1];
sx q[1];
rz(-3.1303774) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0408913) q[0];
sx q[0];
rz(-0.56849231) q[0];
sx q[0];
rz(0.16162737) q[0];
x q[1];
rz(1.8008046) q[2];
sx q[2];
rz(-1.2148148) q[2];
sx q[2];
rz(-1.9614416) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1234731) q[1];
sx q[1];
rz(-1.6247066) q[1];
sx q[1];
rz(2.3192295) q[1];
rz(-0.09402676) q[3];
sx q[3];
rz(-1.8593899) q[3];
sx q[3];
rz(1.2666324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8212905) q[2];
sx q[2];
rz(-1.9052637) q[2];
sx q[2];
rz(0.75766364) q[2];
rz(0.6997987) q[3];
sx q[3];
rz(-2.564513) q[3];
sx q[3];
rz(-0.9052161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5504172) q[0];
sx q[0];
rz(-1.9015522) q[0];
sx q[0];
rz(-0.33429876) q[0];
rz(2.7475157) q[1];
sx q[1];
rz(-0.95029345) q[1];
sx q[1];
rz(-1.2324415) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52215965) q[0];
sx q[0];
rz(-1.8364826) q[0];
sx q[0];
rz(-0.22636087) q[0];
rz(-pi) q[1];
rz(1.3703466) q[2];
sx q[2];
rz(-2.2589189) q[2];
sx q[2];
rz(-1.8375979) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0403172) q[1];
sx q[1];
rz(-0.35666944) q[1];
sx q[1];
rz(-2.2471395) q[1];
rz(0.27421342) q[3];
sx q[3];
rz(-1.0875487) q[3];
sx q[3];
rz(-2.8770214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1831827) q[2];
sx q[2];
rz(-0.56768688) q[2];
sx q[2];
rz(3.0779823) q[2];
rz(2.375864) q[3];
sx q[3];
rz(-2.016341) q[3];
sx q[3];
rz(0.061802797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80339377) q[0];
sx q[0];
rz(-0.82427187) q[0];
sx q[0];
rz(-1.6765539) q[0];
rz(-1.6784531) q[1];
sx q[1];
rz(-1.2668162) q[1];
sx q[1];
rz(-0.75513671) q[1];
rz(-2.9207567) q[2];
sx q[2];
rz(-2.4855843) q[2];
sx q[2];
rz(1.9334855) q[2];
rz(2.9764497) q[3];
sx q[3];
rz(-2.2950635) q[3];
sx q[3];
rz(0.78499817) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
