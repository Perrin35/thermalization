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
rz(-1.6439438) q[0];
sx q[0];
rz(-0.82984501) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0656466) q[0];
sx q[0];
rz(-1.3407205) q[0];
sx q[0];
rz(0.40212698) q[0];
rz(-1.8560156) q[2];
sx q[2];
rz(-0.60456317) q[2];
sx q[2];
rz(-0.7315469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7180011) q[1];
sx q[1];
rz(-0.33034409) q[1];
sx q[1];
rz(-1.524339) q[1];
rz(0.40684367) q[3];
sx q[3];
rz(-1.7115271) q[3];
sx q[3];
rz(1.391747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9709388) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(-0.7286287) q[2];
rz(2.6206) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(2.0200502) q[0];
rz(0.25575486) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(-0.87444011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5912936) q[0];
sx q[0];
rz(-1.8739788) q[0];
sx q[0];
rz(2.0198054) q[0];
rz(-0.15408709) q[2];
sx q[2];
rz(-2.140897) q[2];
sx q[2];
rz(-2.4947583) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9718711) q[1];
sx q[1];
rz(-1.1217146) q[1];
sx q[1];
rz(0.68193087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3784834) q[3];
sx q[3];
rz(-0.27413878) q[3];
sx q[3];
rz(2.6499234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.740739) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(-2.7056616) q[2];
rz(2.46051) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(-0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.23713672) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(1.0748192) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-0.39594617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9112644) q[0];
sx q[0];
rz(-2.0718144) q[0];
sx q[0];
rz(-0.45341861) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6171574) q[2];
sx q[2];
rz(-2.8472387) q[2];
sx q[2];
rz(2.4183395) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1173276) q[1];
sx q[1];
rz(-1.5803442) q[1];
sx q[1];
rz(2.4584103) q[1];
x q[2];
rz(1.5941761) q[3];
sx q[3];
rz(-1.3153207) q[3];
sx q[3];
rz(-1.3083096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-2.9023857) q[2];
rz(0.075332969) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(-1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(2.3216632) q[0];
rz(-0.48768249) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(2.908169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6596286) q[0];
sx q[0];
rz(-0.11419645) q[0];
sx q[0];
rz(-2.1301079) q[0];
x q[1];
rz(1.3350305) q[2];
sx q[2];
rz(-2.1996017) q[2];
sx q[2];
rz(-2.6218888) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4938426) q[1];
sx q[1];
rz(-2.2550681) q[1];
sx q[1];
rz(2.0129596) q[1];
rz(-pi) q[2];
rz(-0.38846429) q[3];
sx q[3];
rz(-0.82287517) q[3];
sx q[3];
rz(0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(0.74742571) q[2];
rz(0.22339544) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(-0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(-3.0773556) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(0.76104004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1431664) q[0];
sx q[0];
rz(-2.879062) q[0];
sx q[0];
rz(-1.6869998) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3752851) q[2];
sx q[2];
rz(-1.3720781) q[2];
sx q[2];
rz(-0.18713258) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0423454) q[1];
sx q[1];
rz(-1.595101) q[1];
sx q[1];
rz(-0.32056067) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28273545) q[3];
sx q[3];
rz(-1.3021819) q[3];
sx q[3];
rz(-0.80073592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3349907) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(3.0495194) q[2];
rz(2.4798685) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(-1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46835607) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(-2.81566) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78544261) q[0];
sx q[0];
rz(-0.66473648) q[0];
sx q[0];
rz(0.63389969) q[0];
rz(-pi) q[1];
rz(-1.2665777) q[2];
sx q[2];
rz(-1.9982669) q[2];
sx q[2];
rz(1.8311335) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1205475) q[1];
sx q[1];
rz(-1.7056587) q[1];
sx q[1];
rz(0.30121505) q[1];
rz(1.553922) q[3];
sx q[3];
rz(-1.4697187) q[3];
sx q[3];
rz(0.24731393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(1.7294653) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(3.022335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0810453) q[0];
sx q[0];
rz(-1.8348872) q[0];
sx q[0];
rz(-0.72780769) q[0];
x q[1];
rz(-1.4085521) q[2];
sx q[2];
rz(-1.679323) q[2];
sx q[2];
rz(2.721399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5057999) q[1];
sx q[1];
rz(-1.074082) q[1];
sx q[1];
rz(-0.15091166) q[1];
rz(-pi) q[2];
rz(-0.42640949) q[3];
sx q[3];
rz(-0.75424131) q[3];
sx q[3];
rz(-0.38503669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44935903) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(0.75891495) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83157241) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(1.0634364) q[0];
rz(0.27451441) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(-0.88561052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793593) q[0];
sx q[0];
rz(-1.0842807) q[0];
sx q[0];
rz(-1.3079206) q[0];
rz(-pi) q[1];
rz(-1.5947072) q[2];
sx q[2];
rz(-2.5619321) q[2];
sx q[2];
rz(1.5946582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.209621) q[1];
sx q[1];
rz(-2.0721657) q[1];
sx q[1];
rz(-2.6905836) q[1];
rz(-pi) q[2];
rz(-2.3732244) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(-0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4132335) q[2];
sx q[2];
rz(-2.3791172) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30329147) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-2.0595179) q[0];
rz(-1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-2.0057604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1866859) q[0];
sx q[0];
rz(-1.0404772) q[0];
sx q[0];
rz(-0.19219877) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.047770569) q[2];
sx q[2];
rz(-2.7063745) q[2];
sx q[2];
rz(1.4620632) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75479924) q[1];
sx q[1];
rz(-0.42897412) q[1];
sx q[1];
rz(-0.11744833) q[1];
rz(0.42580749) q[3];
sx q[3];
rz(-2.1771181) q[3];
sx q[3];
rz(-3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(-2.2686968) q[2];
rz(-0.84351271) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(0.18994722) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2492367) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(1.2783485) q[0];
rz(-2.1168013) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(-1.9445673) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9997864) q[0];
sx q[0];
rz(-2.1394661) q[0];
sx q[0];
rz(-2.6151711) q[0];
rz(-pi) q[1];
rz(-0.12038259) q[2];
sx q[2];
rz(-1.8070081) q[2];
sx q[2];
rz(1.2506968) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3665109) q[1];
sx q[1];
rz(-2.0644036) q[1];
sx q[1];
rz(-0.8702741) q[1];
rz(-pi) q[2];
x q[2];
rz(1.193612) q[3];
sx q[3];
rz(-2.0497272) q[3];
sx q[3];
rz(2.9951028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3060351) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(2.3416134) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(-2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70893127) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(-2.3847053) q[2];
sx q[2];
rz(-2.6464528) q[2];
sx q[2];
rz(1.1337627) q[2];
rz(-2.342631) q[3];
sx q[3];
rz(-2.3286455) q[3];
sx q[3];
rz(0.13959985) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];