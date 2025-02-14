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
rz(-2.4565826) q[0];
sx q[0];
rz(-0.69184875) q[0];
sx q[0];
rz(2.703171) q[0];
rz(3.0216079) q[1];
sx q[1];
rz(-2.9000403) q[1];
sx q[1];
rz(-0.99689364) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0783421) q[0];
sx q[0];
rz(-1.41875) q[0];
sx q[0];
rz(-0.087910533) q[0];
rz(-0.27049944) q[2];
sx q[2];
rz(-1.5322313) q[2];
sx q[2];
rz(2.6730516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62694395) q[1];
sx q[1];
rz(-2.6730826) q[1];
sx q[1];
rz(2.1283988) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9044457) q[3];
sx q[3];
rz(-1.3714694) q[3];
sx q[3];
rz(1.3913499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3991656) q[2];
sx q[2];
rz(-0.70968598) q[2];
sx q[2];
rz(2.4281375) q[2];
rz(-0.82792884) q[3];
sx q[3];
rz(-1.6498339) q[3];
sx q[3];
rz(-2.2152065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68726081) q[0];
sx q[0];
rz(-2.5255272) q[0];
sx q[0];
rz(-1.2607505) q[0];
rz(-2.8997391) q[1];
sx q[1];
rz(-1.2956023) q[1];
sx q[1];
rz(0.58580011) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6144852) q[0];
sx q[0];
rz(-1.4117084) q[0];
sx q[0];
rz(1.2522526) q[0];
rz(-2.7136373) q[2];
sx q[2];
rz(-1.578786) q[2];
sx q[2];
rz(0.34763476) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9487544) q[1];
sx q[1];
rz(-0.84872972) q[1];
sx q[1];
rz(-0.39871902) q[1];
rz(-pi) q[2];
rz(1.5942138) q[3];
sx q[3];
rz(-0.011772884) q[3];
sx q[3];
rz(0.6977607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5213617) q[2];
sx q[2];
rz(-0.62180454) q[2];
sx q[2];
rz(2.8561031) q[2];
rz(-2.1776958) q[3];
sx q[3];
rz(-0.6670835) q[3];
sx q[3];
rz(-0.17531659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46850878) q[0];
sx q[0];
rz(-2.5991169) q[0];
sx q[0];
rz(-0.28955224) q[0];
rz(2.8383004) q[1];
sx q[1];
rz(-1.2723609) q[1];
sx q[1];
rz(-1.9151275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68115409) q[0];
sx q[0];
rz(-1.1484206) q[0];
sx q[0];
rz(2.284221) q[0];
rz(-pi) q[1];
rz(-0.3704088) q[2];
sx q[2];
rz(-0.61282149) q[2];
sx q[2];
rz(1.0072264) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37003041) q[1];
sx q[1];
rz(-2.4874333) q[1];
sx q[1];
rz(-1.0942208) q[1];
rz(-1.6885593) q[3];
sx q[3];
rz(-0.78826815) q[3];
sx q[3];
rz(0.39358172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44846416) q[2];
sx q[2];
rz(-2.144564) q[2];
sx q[2];
rz(-0.86359751) q[2];
rz(-3.038285) q[3];
sx q[3];
rz(-1.7938675) q[3];
sx q[3];
rz(-3.0062655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1358262) q[0];
sx q[0];
rz(-2.5209881) q[0];
sx q[0];
rz(3.0964858) q[0];
rz(0.82334423) q[1];
sx q[1];
rz(-0.38213676) q[1];
sx q[1];
rz(-2.8270922) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5161834) q[0];
sx q[0];
rz(-2.4780031) q[0];
sx q[0];
rz(0.13701464) q[0];
x q[1];
rz(-1.034819) q[2];
sx q[2];
rz(-1.935144) q[2];
sx q[2];
rz(-1.2743769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.040995251) q[1];
sx q[1];
rz(-1.081859) q[1];
sx q[1];
rz(-1.7799499) q[1];
rz(-pi) q[2];
rz(-0.32415402) q[3];
sx q[3];
rz(-1.7739033) q[3];
sx q[3];
rz(2.5517011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0663674) q[2];
sx q[2];
rz(-1.6547357) q[2];
sx q[2];
rz(0.60869795) q[2];
rz(-1.6913951) q[3];
sx q[3];
rz(-0.83289731) q[3];
sx q[3];
rz(2.6628185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4920014) q[0];
sx q[0];
rz(-0.86968017) q[0];
sx q[0];
rz(-0.087652303) q[0];
rz(1.2610669) q[1];
sx q[1];
rz(-1.4050211) q[1];
sx q[1];
rz(-2.2671949) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028624) q[0];
sx q[0];
rz(-2.4695354) q[0];
sx q[0];
rz(-1.623339) q[0];
x q[1];
rz(1.6220785) q[2];
sx q[2];
rz(-2.5502594) q[2];
sx q[2];
rz(-1.4695449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4147784) q[1];
sx q[1];
rz(-0.86514478) q[1];
sx q[1];
rz(2.3470641) q[1];
rz(-pi) q[2];
rz(1.9256853) q[3];
sx q[3];
rz(-0.57587934) q[3];
sx q[3];
rz(1.8047386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7586691) q[2];
sx q[2];
rz(-2.1047968) q[2];
sx q[2];
rz(2.5385638) q[2];
rz(2.1488819) q[3];
sx q[3];
rz(-0.38645667) q[3];
sx q[3];
rz(-2.4934798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25931609) q[0];
sx q[0];
rz(-2.3969605) q[0];
sx q[0];
rz(-2.4477006) q[0];
rz(-0.71677417) q[1];
sx q[1];
rz(-2.0718772) q[1];
sx q[1];
rz(-2.1778291) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5187274) q[0];
sx q[0];
rz(-1.4528926) q[0];
sx q[0];
rz(0.43772719) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3974254) q[2];
sx q[2];
rz(-1.7503305) q[2];
sx q[2];
rz(-2.4712359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48719826) q[1];
sx q[1];
rz(-2.4787612) q[1];
sx q[1];
rz(0.93516949) q[1];
x q[2];
rz(-1.5429872) q[3];
sx q[3];
rz(-2.1429792) q[3];
sx q[3];
rz(-1.5861301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11513772) q[2];
sx q[2];
rz(-0.28885767) q[2];
sx q[2];
rz(0.10819437) q[2];
rz(-2.1859956) q[3];
sx q[3];
rz(-3.1028265) q[3];
sx q[3];
rz(-0.55967104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3947802) q[0];
sx q[0];
rz(-0.17212269) q[0];
sx q[0];
rz(-3.0346003) q[0];
rz(2.6394898) q[1];
sx q[1];
rz(-1.5109477) q[1];
sx q[1];
rz(-2.9923901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82115) q[0];
sx q[0];
rz(-1.7217198) q[0];
sx q[0];
rz(1.3942777) q[0];
rz(1.0525429) q[2];
sx q[2];
rz(-2.0796806) q[2];
sx q[2];
rz(-2.756812) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2857692) q[1];
sx q[1];
rz(-2.2659784) q[1];
sx q[1];
rz(-2.7421363) q[1];
rz(-pi) q[2];
rz(-0.42909166) q[3];
sx q[3];
rz(-1.4403995) q[3];
sx q[3];
rz(0.71500378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2627829) q[2];
sx q[2];
rz(-0.97027367) q[2];
sx q[2];
rz(2.523017) q[2];
rz(-0.68459073) q[3];
sx q[3];
rz(-0.22817831) q[3];
sx q[3];
rz(0.98606199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2672284) q[0];
sx q[0];
rz(-1.3716797) q[0];
sx q[0];
rz(2.5688963) q[0];
rz(2.122208) q[1];
sx q[1];
rz(-2.9044594) q[1];
sx q[1];
rz(3.0651029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6394326) q[0];
sx q[0];
rz(-1.0442088) q[0];
sx q[0];
rz(-2.9798085) q[0];
rz(2.8440153) q[2];
sx q[2];
rz(-1.7817162) q[2];
sx q[2];
rz(-0.95266137) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.98497154) q[1];
sx q[1];
rz(-0.67423361) q[1];
sx q[1];
rz(-1.3944425) q[1];
rz(-pi) q[2];
rz(-2.7143257) q[3];
sx q[3];
rz(-2.2133641) q[3];
sx q[3];
rz(1.6258607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0387592) q[2];
sx q[2];
rz(-0.88492727) q[2];
sx q[2];
rz(2.6494675) q[2];
rz(0.37880185) q[3];
sx q[3];
rz(-2.7087961) q[3];
sx q[3];
rz(-2.2545599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73961306) q[0];
sx q[0];
rz(-0.48960632) q[0];
sx q[0];
rz(2.711645) q[0];
rz(-0.49067378) q[1];
sx q[1];
rz(-2.6538167) q[1];
sx q[1];
rz(1.0027764) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7512892) q[0];
sx q[0];
rz(-1.0868158) q[0];
sx q[0];
rz(-2.8439265) q[0];
x q[1];
rz(-0.077024027) q[2];
sx q[2];
rz(-0.90156889) q[2];
sx q[2];
rz(0.86996824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2247315) q[1];
sx q[1];
rz(-1.5229578) q[1];
sx q[1];
rz(-2.6754154) q[1];
x q[2];
rz(-2.6096576) q[3];
sx q[3];
rz(-0.77456512) q[3];
sx q[3];
rz(2.3331785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7542725) q[2];
sx q[2];
rz(-0.60101271) q[2];
sx q[2];
rz(0.69166541) q[2];
rz(0.12868853) q[3];
sx q[3];
rz(-1.5920937) q[3];
sx q[3];
rz(-0.58840978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9733031) q[0];
sx q[0];
rz(-0.099270865) q[0];
sx q[0];
rz(0.6231935) q[0];
rz(-0.19459952) q[1];
sx q[1];
rz(-2.01229) q[1];
sx q[1];
rz(0.87619877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.122418) q[0];
sx q[0];
rz(-3.0184973) q[0];
sx q[0];
rz(1.9341296) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16531971) q[2];
sx q[2];
rz(-0.7484127) q[2];
sx q[2];
rz(-0.95647631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33223486) q[1];
sx q[1];
rz(-1.7064662) q[1];
sx q[1];
rz(2.7103781) q[1];
rz(-pi) q[2];
rz(0.51214062) q[3];
sx q[3];
rz(-2.509553) q[3];
sx q[3];
rz(-1.8919945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1759922) q[2];
sx q[2];
rz(-2.8922562) q[2];
sx q[2];
rz(2.5052137) q[2];
rz(-1.5198358) q[3];
sx q[3];
rz(-0.88080019) q[3];
sx q[3];
rz(-0.22708587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31430055) q[0];
sx q[0];
rz(-1.5999595) q[0];
sx q[0];
rz(-0.40524361) q[0];
rz(-1.7397407) q[1];
sx q[1];
rz(-1.4973462) q[1];
sx q[1];
rz(-3.0037465) q[1];
rz(0.62217106) q[2];
sx q[2];
rz(-1.6710186) q[2];
sx q[2];
rz(0.98115151) q[2];
rz(1.6062395) q[3];
sx q[3];
rz(-1.0349904) q[3];
sx q[3];
rz(-2.4258202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
