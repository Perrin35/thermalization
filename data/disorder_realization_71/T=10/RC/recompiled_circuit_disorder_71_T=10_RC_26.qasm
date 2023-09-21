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
rz(2.2652664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012955879) q[0];
sx q[0];
rz(-0.46015938) q[0];
sx q[0];
rz(0.53928661) q[0];
rz(1.285577) q[2];
sx q[2];
rz(-0.60456317) q[2];
sx q[2];
rz(2.4100458) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7180011) q[1];
sx q[1];
rz(-0.33034409) q[1];
sx q[1];
rz(-1.524339) q[1];
rz(-pi) q[2];
rz(1.417744) q[3];
sx q[3];
rz(-1.1682086) q[3];
sx q[3];
rz(-3.0229085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17065389) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(0.5209926) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-2.0200502) q[0];
rz(-2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(0.87444011) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.550299) q[0];
sx q[0];
rz(-1.8739788) q[0];
sx q[0];
rz(2.0198054) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1463257) q[2];
sx q[2];
rz(-1.4412291) q[2];
sx q[2];
rz(-2.1339983) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9718711) q[1];
sx q[1];
rz(-2.019878) q[1];
sx q[1];
rz(0.68193087) q[1];
rz(0.053697649) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(-2.4503436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.740739) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(2.7056616) q[2];
rz(-2.46051) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23713672) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-2.3020321) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-2.7456465) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2303282) q[0];
sx q[0];
rz(-1.0697782) q[0];
sx q[0];
rz(-2.688174) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2767407) q[2];
sx q[2];
rz(-1.5842423) q[2];
sx q[2];
rz(2.3384192) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0242651) q[1];
sx q[1];
rz(-1.5612484) q[1];
sx q[1];
rz(-2.4584103) q[1];
rz(-pi) q[2];
rz(-0.08926908) q[3];
sx q[3];
rz(-0.25651989) q[3];
sx q[3];
rz(1.741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(-0.23920693) q[2];
rz(-3.0662597) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(-1.3031561) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72162119) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(2.3216632) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-2.908169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48196402) q[0];
sx q[0];
rz(-3.0273962) q[0];
sx q[0];
rz(-2.1301079) q[0];
rz(-pi) q[1];
rz(0.31077023) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(-0.13194612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64775002) q[1];
sx q[1];
rz(-2.2550681) q[1];
sx q[1];
rz(-2.0129596) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1831746) q[3];
sx q[3];
rz(-0.82510199) q[3];
sx q[3];
rz(2.2771698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(0.74742571) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(-0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(-3.0773556) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(2.3805526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8262186) q[0];
sx q[0];
rz(-1.5407019) q[0];
sx q[0];
rz(1.3099567) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76743482) q[2];
sx q[2];
rz(-2.8637297) q[2];
sx q[2];
rz(-0.9741191) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3984822) q[1];
sx q[1];
rz(-2.8201436) q[1];
sx q[1];
rz(-3.0645963) q[1];
rz(0.77869271) q[3];
sx q[3];
rz(-0.38749309) q[3];
sx q[3];
rz(-0.029822895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80660194) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-3.0495194) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(-1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-0.026542149) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(-0.32593265) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8317141) q[0];
sx q[0];
rz(-1.9448115) q[0];
sx q[0];
rz(2.5783587) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8750149) q[2];
sx q[2];
rz(-1.9982669) q[2];
sx q[2];
rz(-1.3104591) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1833916) q[1];
sx q[1];
rz(-2.8124053) q[1];
sx q[1];
rz(-2.7126461) q[1];
rz(-pi) q[2];
rz(0.10109191) q[3];
sx q[3];
rz(-1.5875845) q[3];
sx q[3];
rz(-1.8198131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5439593) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.7116961) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(0.11925764) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79524604) q[0];
sx q[0];
rz(-2.3756785) q[0];
sx q[0];
rz(-0.38608293) q[0];
rz(-pi) q[1];
rz(-1.4085521) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(0.42019368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5057999) q[1];
sx q[1];
rz(-1.074082) q[1];
sx q[1];
rz(0.15091166) q[1];
rz(-pi) q[2];
rz(1.9414385) q[3];
sx q[3];
rz(-2.2439085) q[3];
sx q[3];
rz(2.9692269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(-0.75891495) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(1.0634364) q[0];
rz(-2.8670782) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(0.88561052) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53997707) q[0];
sx q[0];
rz(-2.5936539) q[0];
sx q[0];
rz(-0.45666306) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9912668) q[2];
sx q[2];
rz(-1.5577003) q[2];
sx q[2];
rz(-3.1377369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0091128) q[1];
sx q[1];
rz(-1.9630034) q[1];
sx q[1];
rz(-2.1177887) q[1];
rz(-2.3732244) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(-0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(-2.0098861) q[2];
rz(1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.926698) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30329147) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(1.0820748) q[0];
rz(1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(1.1358322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95490676) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(-0.19219877) q[0];
rz(1.5485974) q[2];
sx q[2];
rz(-2.0054842) q[2];
sx q[2];
rz(1.7322025) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62577265) q[1];
sx q[1];
rz(-1.1449709) q[1];
sx q[1];
rz(1.6243402) q[1];
rz(-0.42580749) q[3];
sx q[3];
rz(-0.96447456) q[3];
sx q[3];
rz(-3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(-0.84351271) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-1.1268076) q[1];
sx q[1];
rz(1.9445673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7322757) q[0];
sx q[0];
rz(-2.0079552) q[0];
sx q[0];
rz(-2.2073295) q[0];
x q[1];
rz(0.12038259) q[2];
sx q[2];
rz(-1.8070081) q[2];
sx q[2];
rz(1.8908959) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69233209) q[1];
sx q[1];
rz(-0.83220607) q[1];
sx q[1];
rz(2.2663121) q[1];
x q[2];
rz(1.9479806) q[3];
sx q[3];
rz(-2.0497272) q[3];
sx q[3];
rz(-2.9951028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(0.79997921) q[2];
rz(-1.1768613) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(-0.95361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(0.52195436) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-1.9258826) q[2];
sx q[2];
rz(-1.9234895) q[2];
sx q[2];
rz(-1.1870155) q[2];
rz(-0.79896169) q[3];
sx q[3];
rz(-0.81294717) q[3];
sx q[3];
rz(-3.0019928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
