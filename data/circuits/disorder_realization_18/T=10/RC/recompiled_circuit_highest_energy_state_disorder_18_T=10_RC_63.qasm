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
rz(1.7914766) q[0];
sx q[0];
rz(-0.67576367) q[0];
sx q[0];
rz(-0.0078553353) q[0];
rz(0.76771843) q[1];
sx q[1];
rz(-1.5239198) q[1];
sx q[1];
rz(-0.24267264) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9822916) q[0];
sx q[0];
rz(-0.86573273) q[0];
sx q[0];
rz(2.9755735) q[0];
x q[1];
rz(2.484876) q[2];
sx q[2];
rz(-0.55865951) q[2];
sx q[2];
rz(0.61121537) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3585904) q[1];
sx q[1];
rz(-1.5944408) q[1];
sx q[1];
rz(-0.35332291) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8318066) q[3];
sx q[3];
rz(-2.4521146) q[3];
sx q[3];
rz(2.8387217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0365888) q[2];
sx q[2];
rz(-1.236311) q[2];
sx q[2];
rz(2.429602) q[2];
rz(0.30501929) q[3];
sx q[3];
rz(-0.22235338) q[3];
sx q[3];
rz(0.045507889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104093) q[0];
sx q[0];
rz(-1.2774066) q[0];
sx q[0];
rz(-1.7741868) q[0];
rz(-3.101688) q[1];
sx q[1];
rz(-2.3343562) q[1];
sx q[1];
rz(0.59919277) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78135787) q[0];
sx q[0];
rz(-2.0443981) q[0];
sx q[0];
rz(-1.6998769) q[0];
rz(-pi) q[1];
rz(-0.56053253) q[2];
sx q[2];
rz(-0.90122242) q[2];
sx q[2];
rz(0.59434429) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6015905) q[1];
sx q[1];
rz(-2.4036602) q[1];
sx q[1];
rz(-2.6327391) q[1];
rz(2.8017524) q[3];
sx q[3];
rz(-2.426894) q[3];
sx q[3];
rz(3.0106737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0590608) q[2];
sx q[2];
rz(-0.56190562) q[2];
sx q[2];
rz(-1.3085636) q[2];
rz(0.023905309) q[3];
sx q[3];
rz(-1.5289565) q[3];
sx q[3];
rz(1.5628975) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4557274) q[0];
sx q[0];
rz(-0.66615921) q[0];
sx q[0];
rz(0.84079963) q[0];
rz(1.0288382) q[1];
sx q[1];
rz(-2.7327171) q[1];
sx q[1];
rz(-1.4604481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90636364) q[0];
sx q[0];
rz(-1.1113941) q[0];
sx q[0];
rz(1.4109341) q[0];
x q[1];
rz(2.4912848) q[2];
sx q[2];
rz(-0.57194369) q[2];
sx q[2];
rz(0.24947333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3884133) q[1];
sx q[1];
rz(-2.3824661) q[1];
sx q[1];
rz(2.6469346) q[1];
x q[2];
rz(-1.5140686) q[3];
sx q[3];
rz(-2.1773844) q[3];
sx q[3];
rz(-1.1067672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2568405) q[2];
sx q[2];
rz(-1.931087) q[2];
sx q[2];
rz(0.15667008) q[2];
rz(1.6414075) q[3];
sx q[3];
rz(-1.2431966) q[3];
sx q[3];
rz(-0.18178864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665862) q[0];
sx q[0];
rz(-1.2970507) q[0];
sx q[0];
rz(0.080667607) q[0];
rz(2.0196041) q[1];
sx q[1];
rz(-0.64067084) q[1];
sx q[1];
rz(1.5871619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51475731) q[0];
sx q[0];
rz(-0.74585241) q[0];
sx q[0];
rz(-1.4738073) q[0];
rz(-2.8055259) q[2];
sx q[2];
rz(-0.12305752) q[2];
sx q[2];
rz(0.50257909) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39496468) q[1];
sx q[1];
rz(-2.0680799) q[1];
sx q[1];
rz(-2.4056466) q[1];
x q[2];
rz(1.2359592) q[3];
sx q[3];
rz(-0.99042884) q[3];
sx q[3];
rz(2.133647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2598205) q[2];
sx q[2];
rz(-1.0910923) q[2];
sx q[2];
rz(1.0238949) q[2];
rz(0.77670589) q[3];
sx q[3];
rz(-1.9909765) q[3];
sx q[3];
rz(2.1385433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40989947) q[0];
sx q[0];
rz(-0.85407805) q[0];
sx q[0];
rz(0.62752974) q[0];
rz(-2.4420786) q[1];
sx q[1];
rz(-1.0001837) q[1];
sx q[1];
rz(2.2342009) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8193247) q[0];
sx q[0];
rz(-1.5307679) q[0];
sx q[0];
rz(2.3884474) q[0];
rz(-pi) q[1];
rz(0.19030119) q[2];
sx q[2];
rz(-1.4138599) q[2];
sx q[2];
rz(2.4630594) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6945362) q[1];
sx q[1];
rz(-1.2061822) q[1];
sx q[1];
rz(-2.1237808) q[1];
x q[2];
rz(1.0346342) q[3];
sx q[3];
rz(-1.6214633) q[3];
sx q[3];
rz(1.9127653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13713914) q[2];
sx q[2];
rz(-1.8567825) q[2];
sx q[2];
rz(0.84490204) q[2];
rz(-0.6984624) q[3];
sx q[3];
rz(-2.4013077) q[3];
sx q[3];
rz(-0.79160488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7591105) q[0];
sx q[0];
rz(-1.6805205) q[0];
sx q[0];
rz(-1.8735029) q[0];
rz(-1.6146487) q[1];
sx q[1];
rz(-2.0497649) q[1];
sx q[1];
rz(2.7472034) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6170664) q[0];
sx q[0];
rz(-1.6163278) q[0];
sx q[0];
rz(-2.2152053) q[0];
rz(-pi) q[1];
rz(-3.0494681) q[2];
sx q[2];
rz(-0.73706223) q[2];
sx q[2];
rz(-1.2110405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7463267) q[1];
sx q[1];
rz(-1.526064) q[1];
sx q[1];
rz(2.4886063) q[1];
x q[2];
rz(1.4930356) q[3];
sx q[3];
rz(-1.2111411) q[3];
sx q[3];
rz(-0.80870562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7093198) q[2];
sx q[2];
rz(-3.0797112) q[2];
sx q[2];
rz(-2.5800932) q[2];
rz(-2.0105441) q[3];
sx q[3];
rz(-2.3384422) q[3];
sx q[3];
rz(0.48857442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442218) q[0];
sx q[0];
rz(-0.18246305) q[0];
sx q[0];
rz(-0.23319787) q[0];
rz(3.0274262) q[1];
sx q[1];
rz(-2.3515067) q[1];
sx q[1];
rz(-0.17359576) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6703807) q[0];
sx q[0];
rz(-0.91915536) q[0];
sx q[0];
rz(2.6782413) q[0];
rz(-pi) q[1];
rz(-2.0162705) q[2];
sx q[2];
rz(-1.2230754) q[2];
sx q[2];
rz(1.7537774) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50362316) q[1];
sx q[1];
rz(-1.4084647) q[1];
sx q[1];
rz(1.7945402) q[1];
rz(1.5830481) q[3];
sx q[3];
rz(-1.7594595) q[3];
sx q[3];
rz(-1.0808627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0357828) q[2];
sx q[2];
rz(-2.1828987) q[2];
sx q[2];
rz(0.85476533) q[2];
rz(-0.0828951) q[3];
sx q[3];
rz(-1.9708743) q[3];
sx q[3];
rz(-0.054072592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4891124) q[0];
sx q[0];
rz(-1.880045) q[0];
sx q[0];
rz(-1.9626807) q[0];
rz(2.1401801) q[1];
sx q[1];
rz(-1.2636355) q[1];
sx q[1];
rz(-2.9875535) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4074041) q[0];
sx q[0];
rz(-1.2449236) q[0];
sx q[0];
rz(2.4209321) q[0];
x q[1];
rz(-2.270584) q[2];
sx q[2];
rz(-0.99620512) q[2];
sx q[2];
rz(-0.99829295) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6045842) q[1];
sx q[1];
rz(-1.2717383) q[1];
sx q[1];
rz(1.0826151) q[1];
x q[2];
rz(-2.4064111) q[3];
sx q[3];
rz(-0.38215853) q[3];
sx q[3];
rz(2.2459386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.68086326) q[2];
sx q[2];
rz(-2.0622084) q[2];
sx q[2];
rz(0.66780773) q[2];
rz(1.7364511) q[3];
sx q[3];
rz(-1.7738155) q[3];
sx q[3];
rz(0.63647979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77968303) q[0];
sx q[0];
rz(-1.6661665) q[0];
sx q[0];
rz(0.59378004) q[0];
rz(1.2807912) q[1];
sx q[1];
rz(-2.180876) q[1];
sx q[1];
rz(1.8437754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7379388) q[0];
sx q[0];
rz(-0.58890051) q[0];
sx q[0];
rz(-0.1299917) q[0];
rz(2.7251352) q[2];
sx q[2];
rz(-0.67085941) q[2];
sx q[2];
rz(2.351298) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1480963) q[1];
sx q[1];
rz(-1.2132753) q[1];
sx q[1];
rz(1.7054249) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7421977) q[3];
sx q[3];
rz(-0.77091445) q[3];
sx q[3];
rz(-2.7293929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7456776) q[2];
sx q[2];
rz(-0.021947689) q[2];
sx q[2];
rz(-1.5094666) q[2];
rz(0.11219003) q[3];
sx q[3];
rz(-1.0271007) q[3];
sx q[3];
rz(1.3794911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.6701732) q[0];
sx q[0];
rz(-1.0059953) q[0];
sx q[0];
rz(0.26563409) q[0];
rz(0.98948014) q[1];
sx q[1];
rz(-1.8786636) q[1];
sx q[1];
rz(-2.8094453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6270094) q[0];
sx q[0];
rz(-0.3334612) q[0];
sx q[0];
rz(-1.5824806) q[0];
rz(-pi) q[1];
rz(-2.1493456) q[2];
sx q[2];
rz(-1.3253115) q[2];
sx q[2];
rz(3.0089889) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0964342) q[1];
sx q[1];
rz(-1.2306884) q[1];
sx q[1];
rz(0.051326871) q[1];
rz(-3.0933558) q[3];
sx q[3];
rz(-1.3833481) q[3];
sx q[3];
rz(1.5099883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0418479) q[2];
sx q[2];
rz(-1.0724649) q[2];
sx q[2];
rz(2.9912046) q[2];
rz(0.8485052) q[3];
sx q[3];
rz(-0.34981194) q[3];
sx q[3];
rz(1.295804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083241845) q[0];
sx q[0];
rz(-1.8235089) q[0];
sx q[0];
rz(1.9705809) q[0];
rz(1.8104443) q[1];
sx q[1];
rz(-2.34453) q[1];
sx q[1];
rz(-1.1042368) q[1];
rz(2.3103279) q[2];
sx q[2];
rz(-1.6192042) q[2];
sx q[2];
rz(-1.529196) q[2];
rz(1.5965309) q[3];
sx q[3];
rz(-2.2774057) q[3];
sx q[3];
rz(-2.0792815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
