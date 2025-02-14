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
rz(0.74547493) q[0];
sx q[0];
rz(-0.59362721) q[0];
sx q[0];
rz(-2.797085) q[0];
rz(1.1747107) q[1];
sx q[1];
rz(-2.7732958) q[1];
sx q[1];
rz(-2.5370497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060319107) q[0];
sx q[0];
rz(-0.7501157) q[0];
sx q[0];
rz(-1.4846205) q[0];
rz(-pi) q[1];
rz(0.069434631) q[2];
sx q[2];
rz(-1.3926818) q[2];
sx q[2];
rz(-1.3645862) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74605265) q[1];
sx q[1];
rz(-2.722858) q[1];
sx q[1];
rz(0.0050425649) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56960241) q[3];
sx q[3];
rz(-2.0853373) q[3];
sx q[3];
rz(-1.2698184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1623666) q[2];
sx q[2];
rz(-1.6567076) q[2];
sx q[2];
rz(1.7171198) q[2];
rz(3.034397) q[3];
sx q[3];
rz(-2.9469979) q[3];
sx q[3];
rz(-2.181982) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87218881) q[0];
sx q[0];
rz(-0.93282455) q[0];
sx q[0];
rz(1.3099571) q[0];
rz(0.45920363) q[1];
sx q[1];
rz(-1.2745067) q[1];
sx q[1];
rz(-2.8244663) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32401356) q[0];
sx q[0];
rz(-1.5234103) q[0];
sx q[0];
rz(-0.13806483) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0917815) q[2];
sx q[2];
rz(-1.5684109) q[2];
sx q[2];
rz(2.002169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.210798) q[1];
sx q[1];
rz(-2.5902777) q[1];
sx q[1];
rz(-3.0153794) q[1];
x q[2];
rz(2.9212115) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(3.0241682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7121048) q[2];
sx q[2];
rz(-1.8450582) q[2];
sx q[2];
rz(-0.8075766) q[2];
rz(2.5782222) q[3];
sx q[3];
rz(-0.49401504) q[3];
sx q[3];
rz(2.0886683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.5092369) q[0];
sx q[0];
rz(-2.95166) q[0];
sx q[0];
rz(0.83455363) q[0];
rz(2.6368311) q[1];
sx q[1];
rz(-2.7752462) q[1];
sx q[1];
rz(1.4588446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014011009) q[0];
sx q[0];
rz(-2.547285) q[0];
sx q[0];
rz(2.7796714) q[0];
x q[1];
rz(0.81291109) q[2];
sx q[2];
rz(-2.2208919) q[2];
sx q[2];
rz(1.4073828) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9533938) q[1];
sx q[1];
rz(-2.2911304) q[1];
sx q[1];
rz(1.6745643) q[1];
rz(2.1670413) q[3];
sx q[3];
rz(-0.85874346) q[3];
sx q[3];
rz(0.94030583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0316281) q[2];
sx q[2];
rz(-1.6692903) q[2];
sx q[2];
rz(0.36299452) q[2];
rz(1.1686769) q[3];
sx q[3];
rz(-0.76255885) q[3];
sx q[3];
rz(-2.3762083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4204243) q[0];
sx q[0];
rz(-2.6154501) q[0];
sx q[0];
rz(-2.5508733) q[0];
rz(-2.4144454) q[1];
sx q[1];
rz(-0.37494451) q[1];
sx q[1];
rz(-2.4730543) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0787857) q[0];
sx q[0];
rz(-1.5292087) q[0];
sx q[0];
rz(1.6522626) q[0];
rz(0.48769571) q[2];
sx q[2];
rz(-2.3137011) q[2];
sx q[2];
rz(-2.6726892) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6741028) q[1];
sx q[1];
rz(-0.88574848) q[1];
sx q[1];
rz(2.5193307) q[1];
rz(-pi) q[2];
rz(1.0428069) q[3];
sx q[3];
rz(-0.90195459) q[3];
sx q[3];
rz(-2.5537415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0380402) q[2];
sx q[2];
rz(-1.870564) q[2];
sx q[2];
rz(-1.6345778) q[2];
rz(-0.42190894) q[3];
sx q[3];
rz(-0.76015893) q[3];
sx q[3];
rz(2.4222477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6307395) q[0];
sx q[0];
rz(-0.35950867) q[0];
sx q[0];
rz(2.0810293) q[0];
rz(3.0781436) q[1];
sx q[1];
rz(-0.63065204) q[1];
sx q[1];
rz(-2.5877156) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4872525) q[0];
sx q[0];
rz(-2.3759288) q[0];
sx q[0];
rz(2.2310217) q[0];
x q[1];
rz(-0.62184288) q[2];
sx q[2];
rz(-2.1178195) q[2];
sx q[2];
rz(1.1012929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6468473) q[1];
sx q[1];
rz(-1.8738197) q[1];
sx q[1];
rz(-1.5290497) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2451671) q[3];
sx q[3];
rz(-1.4417329) q[3];
sx q[3];
rz(1.3323402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8744897) q[2];
sx q[2];
rz(-0.59659448) q[2];
sx q[2];
rz(2.5388517) q[2];
rz(3.0722669) q[3];
sx q[3];
rz(-1.5452789) q[3];
sx q[3];
rz(0.16665211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9423187) q[0];
sx q[0];
rz(-2.5397904) q[0];
sx q[0];
rz(0.45403516) q[0];
rz(-1.2557238) q[1];
sx q[1];
rz(-0.95971003) q[1];
sx q[1];
rz(1.7032636) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3011264) q[0];
sx q[0];
rz(-0.92341237) q[0];
sx q[0];
rz(2.6324138) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8642162) q[2];
sx q[2];
rz(-1.2834594) q[2];
sx q[2];
rz(-0.6544906) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.033129902) q[1];
sx q[1];
rz(-2.3537043) q[1];
sx q[1];
rz(2.2910396) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71261974) q[3];
sx q[3];
rz(-1.9639059) q[3];
sx q[3];
rz(0.46526965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2890275) q[2];
sx q[2];
rz(-1.270741) q[2];
sx q[2];
rz(1.3555869) q[2];
rz(-1.2191314) q[3];
sx q[3];
rz(-1.6617323) q[3];
sx q[3];
rz(1.5752972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3179625) q[0];
sx q[0];
rz(-0.83331236) q[0];
sx q[0];
rz(-1.3939567) q[0];
rz(0.45669237) q[1];
sx q[1];
rz(-2.2629181) q[1];
sx q[1];
rz(1.7535694) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0972892) q[0];
sx q[0];
rz(-2.7345022) q[0];
sx q[0];
rz(1.9096229) q[0];
rz(2.9406383) q[2];
sx q[2];
rz(-1.4575053) q[2];
sx q[2];
rz(-2.410881) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2051831) q[1];
sx q[1];
rz(-0.35297063) q[1];
sx q[1];
rz(0.38683968) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5230266) q[3];
sx q[3];
rz(-1.9962709) q[3];
sx q[3];
rz(-0.24263231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2726511) q[2];
sx q[2];
rz(-2.36918) q[2];
sx q[2];
rz(-2.9504377) q[2];
rz(-1.8027421) q[3];
sx q[3];
rz(-0.50722417) q[3];
sx q[3];
rz(-2.6195781) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1538447) q[0];
sx q[0];
rz(-2.4697883) q[0];
sx q[0];
rz(-1.2220569) q[0];
rz(0.98474312) q[1];
sx q[1];
rz(-0.9877111) q[1];
sx q[1];
rz(2.9827859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5903871) q[0];
sx q[0];
rz(-0.90999167) q[0];
sx q[0];
rz(-1.4676276) q[0];
rz(-0.42694636) q[2];
sx q[2];
rz(-1.5265577) q[2];
sx q[2];
rz(-0.51170631) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1013704) q[1];
sx q[1];
rz(-1.985038) q[1];
sx q[1];
rz(-0.48419063) q[1];
rz(0.21130224) q[3];
sx q[3];
rz(-1.7083454) q[3];
sx q[3];
rz(2.7338303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.020393546) q[2];
sx q[2];
rz(-0.32200107) q[2];
sx q[2];
rz(0.89312345) q[2];
rz(-2.761306) q[3];
sx q[3];
rz(-1.2023353) q[3];
sx q[3];
rz(-1.7072385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0685843) q[0];
sx q[0];
rz(-0.46982729) q[0];
sx q[0];
rz(1.6597066) q[0];
rz(0.98663729) q[1];
sx q[1];
rz(-1.6529129) q[1];
sx q[1];
rz(2.5843487) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38735861) q[0];
sx q[0];
rz(-3.0172303) q[0];
sx q[0];
rz(2.6919305) q[0];
rz(-2.8403867) q[2];
sx q[2];
rz(-1.2264892) q[2];
sx q[2];
rz(2.3489281) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7047953) q[1];
sx q[1];
rz(-0.76334243) q[1];
sx q[1];
rz(-2.3176058) q[1];
rz(-pi) q[2];
rz(3.1246947) q[3];
sx q[3];
rz(-1.5819614) q[3];
sx q[3];
rz(-1.4433235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3096932) q[2];
sx q[2];
rz(-2.5688186) q[2];
sx q[2];
rz(1.0581623) q[2];
rz(1.38331) q[3];
sx q[3];
rz(-1.0388831) q[3];
sx q[3];
rz(0.96341187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6365373) q[0];
sx q[0];
rz(-1.9872682) q[0];
sx q[0];
rz(0.43168798) q[0];
rz(0.70026669) q[1];
sx q[1];
rz(-1.9077178) q[1];
sx q[1];
rz(-2.4468927) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2852531) q[0];
sx q[0];
rz(-2.4972557) q[0];
sx q[0];
rz(2.2088693) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3586559) q[2];
sx q[2];
rz(-2.6797469) q[2];
sx q[2];
rz(-1.0657276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6208261) q[1];
sx q[1];
rz(-1.0559901) q[1];
sx q[1];
rz(-0.52052814) q[1];
rz(-pi) q[2];
rz(2.905013) q[3];
sx q[3];
rz(-0.18356174) q[3];
sx q[3];
rz(-2.4381886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.485864) q[2];
sx q[2];
rz(-0.30980095) q[2];
sx q[2];
rz(2.7244205) q[2];
rz(0.74472767) q[3];
sx q[3];
rz(-1.797902) q[3];
sx q[3];
rz(0.97851306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3485296) q[0];
sx q[0];
rz(-1.2611669) q[0];
sx q[0];
rz(1.4766759) q[0];
rz(-1.2186125) q[1];
sx q[1];
rz(-2.0017793) q[1];
sx q[1];
rz(0.58560169) q[1];
rz(0.89305604) q[2];
sx q[2];
rz(-1.1184659) q[2];
sx q[2];
rz(2.804166) q[2];
rz(-2.1914235) q[3];
sx q[3];
rz(-1.2990549) q[3];
sx q[3];
rz(-2.3334734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
