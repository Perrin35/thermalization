OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61385566) q[0];
sx q[0];
rz(4.6392415) q[0];
sx q[0];
rz(8.5949329) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012955879) q[0];
sx q[0];
rz(-0.46015938) q[0];
sx q[0];
rz(-2.602306) q[0];
rz(-pi) q[1];
rz(0.9853739) q[2];
sx q[2];
rz(-1.7314163) q[2];
sx q[2];
rz(-1.0759682) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4727002) q[1];
sx q[1];
rz(-1.9007705) q[1];
sx q[1];
rz(3.1256691) q[1];
rz(-2.797804) q[3];
sx q[3];
rz(-2.7123835) q[3];
sx q[3];
rz(-0.49376282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9709388) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(2.412964) q[2];
rz(2.6206) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(0.20761028) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347297) q[0];
sx q[0];
rz(-1.9711718) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(0.25575486) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(0.87444011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53423184) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(2.1952654) q[0];
rz(-pi) q[1];
rz(-0.15408709) q[2];
sx q[2];
rz(-2.140897) q[2];
sx q[2];
rz(0.64683435) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16972152) q[1];
sx q[1];
rz(-2.019878) q[1];
sx q[1];
rz(-2.4596618) q[1];
rz(1.3014684) q[3];
sx q[3];
rz(-1.519031) q[3];
sx q[3];
rz(2.2477637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4008537) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(0.43593105) q[2];
rz(-2.46051) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044559) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-0.83956051) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(2.7456465) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0229189) q[0];
sx q[0];
rz(-0.6624822) q[0];
sx q[0];
rz(2.2454717) q[0];
x q[1];
rz(1.5244353) q[2];
sx q[2];
rz(-0.29435396) q[2];
sx q[2];
rz(-2.4183395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1173276) q[1];
sx q[1];
rz(-1.5612484) q[1];
sx q[1];
rz(-0.68318232) q[1];
rz(-pi) q[2];
rz(3.0523236) q[3];
sx q[3];
rz(-0.25651989) q[3];
sx q[3];
rz(1.741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(-0.23920693) q[2];
rz(0.075332969) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72162119) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(-0.81992942) q[0];
rz(2.6539102) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(2.908169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4963835) q[0];
sx q[0];
rz(-1.5102981) q[0];
sx q[0];
rz(1.4738826) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31077023) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(3.0096465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9277716) q[1];
sx q[1];
rz(-1.908761) q[1];
sx q[1];
rz(0.73422276) q[1];
x q[2];
rz(2.3573973) q[3];
sx q[3];
rz(-1.2894221) q[3];
sx q[3];
rz(-0.97660645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(-0.74742571) q[2];
rz(0.22339544) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(-0.064237021) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(-0.76104004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8262186) q[0];
sx q[0];
rz(-1.6008908) q[0];
sx q[0];
rz(1.3099567) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3741578) q[2];
sx q[2];
rz(-0.27786294) q[2];
sx q[2];
rz(-2.1674736) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.099247301) q[1];
sx q[1];
rz(-1.5464916) q[1];
sx q[1];
rz(0.32056067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77869271) q[3];
sx q[3];
rz(-0.38749309) q[3];
sx q[3];
rz(0.029822895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(0.09207329) q[2];
rz(-0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(-1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(-2.81566) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8317141) q[0];
sx q[0];
rz(-1.9448115) q[0];
sx q[0];
rz(2.5783587) q[0];
x q[1];
rz(2.6961156) q[2];
sx q[2];
rz(-1.2947086) q[2];
sx q[2];
rz(-0.38976994) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0210452) q[1];
sx q[1];
rz(-1.435934) q[1];
sx q[1];
rz(-0.30121505) q[1];
rz(3.0405007) q[3];
sx q[3];
rz(-1.5540082) q[3];
sx q[3];
rz(-1.8198131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(-0.65417543) q[2];
rz(1.4298965) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(2.4196999) q[0];
rz(-1.7294653) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(0.11925764) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0810453) q[0];
sx q[0];
rz(-1.3067055) q[0];
sx q[0];
rz(-2.413785) q[0];
rz(3.0316333) q[2];
sx q[2];
rz(-1.409515) q[2];
sx q[2];
rz(-1.1683299) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1489361) q[1];
sx q[1];
rz(-1.438237) q[1];
sx q[1];
rz(-2.0723144) q[1];
x q[2];
rz(-1.2001541) q[3];
sx q[3];
rz(-0.89768411) q[3];
sx q[3];
rz(0.17236575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(-2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(-0.27451441) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(-2.2559821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53997707) q[0];
sx q[0];
rz(-0.5479387) q[0];
sx q[0];
rz(-2.6849296) q[0];
x q[1];
rz(-1.5468855) q[2];
sx q[2];
rz(-2.5619321) q[2];
sx q[2];
rz(1.5469345) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0091128) q[1];
sx q[1];
rz(-1.1785893) q[1];
sx q[1];
rz(-2.1177887) q[1];
rz(0.76836821) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(-0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(-2.0098861) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30329147) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-1.0820748) q[0];
rz(1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-1.1358322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866859) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(2.9493939) q[0];
rz(-pi) q[1];
rz(1.5485974) q[2];
sx q[2];
rz(-1.1361085) q[2];
sx q[2];
rz(-1.7322025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.218704) q[1];
sx q[1];
rz(-1.6195546) q[1];
sx q[1];
rz(-2.7152275) q[1];
rz(0.42580749) q[3];
sx q[3];
rz(-0.96447456) q[3];
sx q[3];
rz(3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0231126) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(0.87289587) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(-2.9516454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89235598) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(1.8632442) q[0];
rz(2.1168013) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(1.1970253) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4599633) q[0];
sx q[0];
rz(-0.75461331) q[0];
sx q[0];
rz(-2.2370536) q[0];
rz(1.3329266) q[2];
sx q[2];
rz(-1.4537721) q[2];
sx q[2];
rz(2.7931917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69233209) q[1];
sx q[1];
rz(-0.83220607) q[1];
sx q[1];
rz(-2.2663121) q[1];
x q[2];
rz(-0.61693807) q[3];
sx q[3];
rz(-0.60041282) q[3];
sx q[3];
rz(2.2850349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8355576) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(2.3416134) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4326614) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(-0.52195436) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(-1.21571) q[2];
sx q[2];
rz(-1.2181031) q[2];
sx q[2];
rz(1.9545771) q[2];
rz(0.63511499) q[3];
sx q[3];
rz(-2.1182346) q[3];
sx q[3];
rz(2.3253141) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
