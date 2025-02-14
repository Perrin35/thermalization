OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8389799) q[0];
sx q[0];
rz(-1.7054727) q[0];
sx q[0];
rz(2.9247395) q[0];
rz(2.6537553) q[1];
sx q[1];
rz(-2.2626329) q[1];
sx q[1];
rz(0.15329696) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72399272) q[0];
sx q[0];
rz(-1.2823021) q[0];
sx q[0];
rz(-0.27758502) q[0];
x q[1];
rz(2.6716483) q[2];
sx q[2];
rz(-1.1026088) q[2];
sx q[2];
rz(-0.17344698) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44854201) q[1];
sx q[1];
rz(-1.3227751) q[1];
sx q[1];
rz(-0.51183707) q[1];
x q[2];
rz(-2.9847095) q[3];
sx q[3];
rz(-2.3290754) q[3];
sx q[3];
rz(0.37655991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6903901) q[2];
sx q[2];
rz(-0.58053747) q[2];
sx q[2];
rz(0.85980493) q[2];
rz(-2.9016923) q[3];
sx q[3];
rz(-2.4247215) q[3];
sx q[3];
rz(2.8587604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4848223) q[0];
sx q[0];
rz(-1.7253933) q[0];
sx q[0];
rz(0.91944486) q[0];
rz(-2.2942309) q[1];
sx q[1];
rz(-2.5571926) q[1];
sx q[1];
rz(-1.970361) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6028588) q[0];
sx q[0];
rz(-0.26358381) q[0];
sx q[0];
rz(-1.3766422) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2861346) q[2];
sx q[2];
rz(-1.8252862) q[2];
sx q[2];
rz(-1.364691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1619809) q[1];
sx q[1];
rz(-1.2250326) q[1];
sx q[1];
rz(2.8678721) q[1];
rz(-pi) q[2];
rz(-1.5492577) q[3];
sx q[3];
rz(-2.1088558) q[3];
sx q[3];
rz(-0.027854161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.072711572) q[2];
sx q[2];
rz(-0.88901192) q[2];
sx q[2];
rz(-1.9286801) q[2];
rz(-2.9291901) q[3];
sx q[3];
rz(-2.0272777) q[3];
sx q[3];
rz(-0.67306486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71639672) q[0];
sx q[0];
rz(-1.9028417) q[0];
sx q[0];
rz(0.58309251) q[0];
rz(1.8816226) q[1];
sx q[1];
rz(-0.59695736) q[1];
sx q[1];
rz(-1.8276385) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1309926) q[0];
sx q[0];
rz(-2.6525462) q[0];
sx q[0];
rz(2.2940192) q[0];
rz(-pi) q[1];
rz(1.5384244) q[2];
sx q[2];
rz(-2.3043423) q[2];
sx q[2];
rz(2.5493279) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78216923) q[1];
sx q[1];
rz(-0.12453989) q[1];
sx q[1];
rz(-1.729639) q[1];
rz(-1.9989955) q[3];
sx q[3];
rz(-1.7230265) q[3];
sx q[3];
rz(0.58109944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6446357) q[2];
sx q[2];
rz(-0.76377112) q[2];
sx q[2];
rz(2.606707) q[2];
rz(-2.6584117) q[3];
sx q[3];
rz(-1.6202241) q[3];
sx q[3];
rz(0.11972891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508761) q[0];
sx q[0];
rz(-3.0829939) q[0];
sx q[0];
rz(-1.4709877) q[0];
rz(-1.927467) q[1];
sx q[1];
rz(-1.7045226) q[1];
sx q[1];
rz(-0.4462744) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8295815) q[0];
sx q[0];
rz(-2.952842) q[0];
sx q[0];
rz(1.944456) q[0];
rz(-pi) q[1];
rz(0.032164737) q[2];
sx q[2];
rz(-0.52903658) q[2];
sx q[2];
rz(0.47495237) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41149263) q[1];
sx q[1];
rz(-2.2124083) q[1];
sx q[1];
rz(1.2506676) q[1];
rz(-3.0837184) q[3];
sx q[3];
rz(-0.24316517) q[3];
sx q[3];
rz(-1.1806219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2176167) q[2];
sx q[2];
rz(-2.6805704) q[2];
sx q[2];
rz(-0.32611845) q[2];
rz(2.7255132) q[3];
sx q[3];
rz(-1.6385498) q[3];
sx q[3];
rz(-0.7777586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.07051) q[0];
sx q[0];
rz(-1.5639045) q[0];
sx q[0];
rz(-2.7974906) q[0];
rz(2.7857419) q[1];
sx q[1];
rz(-2.3882723) q[1];
sx q[1];
rz(0.25486249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0094285) q[0];
sx q[0];
rz(-2.4309506) q[0];
sx q[0];
rz(-0.406213) q[0];
rz(1.1084313) q[2];
sx q[2];
rz(-1.4810665) q[2];
sx q[2];
rz(-1.7347908) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8336941) q[1];
sx q[1];
rz(-0.98684084) q[1];
sx q[1];
rz(2.6112982) q[1];
rz(-pi) q[2];
rz(1.9790824) q[3];
sx q[3];
rz(-0.84813839) q[3];
sx q[3];
rz(1.6734139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43763375) q[2];
sx q[2];
rz(-2.2767229) q[2];
sx q[2];
rz(3.051009) q[2];
rz(-0.80638742) q[3];
sx q[3];
rz(-1.3017637) q[3];
sx q[3];
rz(1.8346132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0818455) q[0];
sx q[0];
rz(-1.9192001) q[0];
sx q[0];
rz(-2.5229689) q[0];
rz(0.6126569) q[1];
sx q[1];
rz(-2.5109992) q[1];
sx q[1];
rz(-0.15288606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83251774) q[0];
sx q[0];
rz(-2.8075657) q[0];
sx q[0];
rz(-0.79132737) q[0];
x q[1];
rz(2.9866508) q[2];
sx q[2];
rz(-0.58523387) q[2];
sx q[2];
rz(0.80807782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5093435) q[1];
sx q[1];
rz(-1.1763089) q[1];
sx q[1];
rz(-3.0837653) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99075192) q[3];
sx q[3];
rz(-1.5914122) q[3];
sx q[3];
rz(2.0772658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9360518) q[2];
sx q[2];
rz(-2.3369393) q[2];
sx q[2];
rz(-2.0118227) q[2];
rz(1.0063082) q[3];
sx q[3];
rz(-1.8215424) q[3];
sx q[3];
rz(1.6442851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048700843) q[0];
sx q[0];
rz(-1.0660271) q[0];
sx q[0];
rz(-2.997828) q[0];
rz(-2.9394506) q[1];
sx q[1];
rz(-0.527924) q[1];
sx q[1];
rz(2.0595097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3976404) q[0];
sx q[0];
rz(-2.2778768) q[0];
sx q[0];
rz(-2.6567064) q[0];
x q[1];
rz(-2.4074102) q[2];
sx q[2];
rz(-2.1819644) q[2];
sx q[2];
rz(-2.2847459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52817594) q[1];
sx q[1];
rz(-0.74494637) q[1];
sx q[1];
rz(-1.8862392) q[1];
x q[2];
rz(-2.1752259) q[3];
sx q[3];
rz(-2.2141738) q[3];
sx q[3];
rz(-2.9288187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14273345) q[2];
sx q[2];
rz(-1.961668) q[2];
sx q[2];
rz(0.72706968) q[2];
rz(-1.3149698) q[3];
sx q[3];
rz(-1.8948137) q[3];
sx q[3];
rz(1.4962014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975248) q[0];
sx q[0];
rz(-1.0777363) q[0];
sx q[0];
rz(-1.9986073) q[0];
rz(2.9179528) q[1];
sx q[1];
rz(-2.5032005) q[1];
sx q[1];
rz(-1.9167831) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95668225) q[0];
sx q[0];
rz(-1.3089797) q[0];
sx q[0];
rz(-1.8641406) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46240004) q[2];
sx q[2];
rz(-2.3650596) q[2];
sx q[2];
rz(0.76867104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30149192) q[1];
sx q[1];
rz(-0.21243851) q[1];
sx q[1];
rz(-1.2919687) q[1];
x q[2];
rz(0.56153058) q[3];
sx q[3];
rz(-1.6609471) q[3];
sx q[3];
rz(0.60240373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58330047) q[2];
sx q[2];
rz(-1.4536828) q[2];
sx q[2];
rz(-2.3395786) q[2];
rz(-1.924104) q[3];
sx q[3];
rz(-3.0018482) q[3];
sx q[3];
rz(1.2453992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8121346) q[0];
sx q[0];
rz(-1.7725002) q[0];
sx q[0];
rz(-1.580397) q[0];
rz(2.2691057) q[1];
sx q[1];
rz(-2.8552738) q[1];
sx q[1];
rz(-0.40503851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0319388) q[0];
sx q[0];
rz(-0.73470107) q[0];
sx q[0];
rz(-0.74560179) q[0];
rz(-pi) q[1];
rz(1.4680176) q[2];
sx q[2];
rz(-0.98066247) q[2];
sx q[2];
rz(-1.128732) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2916221) q[1];
sx q[1];
rz(-2.1287781) q[1];
sx q[1];
rz(0.72429742) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9946278) q[3];
sx q[3];
rz(-1.5515447) q[3];
sx q[3];
rz(1.7363154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27434719) q[2];
sx q[2];
rz(-0.87928191) q[2];
sx q[2];
rz(-2.530063) q[2];
rz(1.8379755) q[3];
sx q[3];
rz(-0.36208624) q[3];
sx q[3];
rz(-2.005579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63854727) q[0];
sx q[0];
rz(-1.2940116) q[0];
sx q[0];
rz(-3.0882623) q[0];
rz(2.5697925) q[1];
sx q[1];
rz(-1.6245533) q[1];
sx q[1];
rz(-1.8624051) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8227575) q[0];
sx q[0];
rz(-1.822572) q[0];
sx q[0];
rz(2.84637) q[0];
x q[1];
rz(2.0512852) q[2];
sx q[2];
rz(-1.4391403) q[2];
sx q[2];
rz(2.5327794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6257929) q[1];
sx q[1];
rz(-2.1560921) q[1];
sx q[1];
rz(-3.0477691) q[1];
rz(-pi) q[2];
rz(1.199715) q[3];
sx q[3];
rz(-1.7744007) q[3];
sx q[3];
rz(-0.18818391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6751487) q[2];
sx q[2];
rz(-1.5326591) q[2];
sx q[2];
rz(2.4565728) q[2];
rz(0.78392616) q[3];
sx q[3];
rz(-2.7214366) q[3];
sx q[3];
rz(2.4043731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6017799) q[0];
sx q[0];
rz(-1.2417326) q[0];
sx q[0];
rz(-1.7057521) q[0];
rz(-0.80541366) q[1];
sx q[1];
rz(-0.46331159) q[1];
sx q[1];
rz(0.70809271) q[1];
rz(-2.851839) q[2];
sx q[2];
rz(-1.7782246) q[2];
sx q[2];
rz(-2.092337) q[2];
rz(-0.2652771) q[3];
sx q[3];
rz(-2.7133085) q[3];
sx q[3];
rz(2.6648236) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
