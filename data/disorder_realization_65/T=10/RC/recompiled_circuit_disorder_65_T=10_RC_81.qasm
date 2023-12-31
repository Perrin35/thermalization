OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(-1.4927827) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(-1.0746497) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.933691) q[0];
sx q[0];
rz(-1.5760734) q[0];
sx q[0];
rz(-0.018215608) q[0];
rz(-pi) q[1];
rz(2.3518402) q[2];
sx q[2];
rz(-1.4573922) q[2];
sx q[2];
rz(-1.5838768) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5731116) q[1];
sx q[1];
rz(-2.9955578) q[1];
sx q[1];
rz(2.5558429) q[1];
rz(-pi) q[2];
rz(-0.69016506) q[3];
sx q[3];
rz(-2.4881722) q[3];
sx q[3];
rz(1.1284459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82912123) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(-1.1734022) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(-2.9287958) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56869498) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(-2.9808295) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(1.5011903) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879724) q[0];
sx q[0];
rz(-3.0164218) q[0];
sx q[0];
rz(-0.20010389) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31200774) q[2];
sx q[2];
rz(-1.1216251) q[2];
sx q[2];
rz(-2.1829407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7306819) q[1];
sx q[1];
rz(-0.69697471) q[1];
sx q[1];
rz(2.8721786) q[1];
rz(-pi) q[2];
rz(-3.1355751) q[3];
sx q[3];
rz(-1.179751) q[3];
sx q[3];
rz(-2.1276377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.5141053) q[2];
rz(-1.9747915) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353772) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.5037781) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5119748) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(2.2114121) q[0];
x q[1];
rz(0.47567993) q[2];
sx q[2];
rz(-2.1026346) q[2];
sx q[2];
rz(1.4059517) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7392002) q[1];
sx q[1];
rz(-1.8715197) q[1];
sx q[1];
rz(-1.5261569) q[1];
rz(-3.010473) q[3];
sx q[3];
rz(-0.91851202) q[3];
sx q[3];
rz(-0.91344792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(2.3977996) q[2];
rz(2.8841694) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96823111) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(1.051735) q[0];
rz(-0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.9925041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.944753) q[0];
sx q[0];
rz(-1.4595932) q[0];
sx q[0];
rz(2.8993594) q[0];
rz(-pi) q[1];
rz(2.1707702) q[2];
sx q[2];
rz(-0.68471013) q[2];
sx q[2];
rz(3.0175356) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.541084) q[1];
sx q[1];
rz(-2.2693172) q[1];
sx q[1];
rz(0.38703106) q[1];
x q[2];
rz(0.74242384) q[3];
sx q[3];
rz(-1.8291049) q[3];
sx q[3];
rz(1.0148259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-0.66876283) q[2];
rz(-1.8956005) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(3.1252089) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(1.0523798) q[0];
rz(-1.9309689) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8540871) q[0];
sx q[0];
rz(-0.26212087) q[0];
sx q[0];
rz(0.8192807) q[0];
rz(2.7344633) q[2];
sx q[2];
rz(-1.9756769) q[2];
sx q[2];
rz(1.4332353) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4864038) q[1];
sx q[1];
rz(-2.9826394) q[1];
sx q[1];
rz(-2.7093191) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0377117) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(-1.0853634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(-1.0926251) q[2];
rz(-0.24400273) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(1.8238235) q[0];
rz(2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(0.96907369) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9861525) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(-0.040779671) q[0];
x q[1];
rz(-1.6386119) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(2.3358047) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.59086696) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(2.4025737) q[1];
rz(2.6763775) q[3];
sx q[3];
rz(-1.0930659) q[3];
sx q[3];
rz(0.623869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(-2.7039841) q[2];
rz(-3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(-1.3597885) q[0];
rz(-1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.4356027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5098269) q[0];
sx q[0];
rz(-2.0938211) q[0];
sx q[0];
rz(-0.57720362) q[0];
x q[1];
rz(3.1232386) q[2];
sx q[2];
rz(-1.5510501) q[2];
sx q[2];
rz(-0.28723082) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85305271) q[1];
sx q[1];
rz(-1.2097881) q[1];
sx q[1];
rz(2.8193874) q[1];
rz(-pi) q[2];
rz(0.42922677) q[3];
sx q[3];
rz(-1.1042522) q[3];
sx q[3];
rz(0.45688094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4009565) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(-2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-0.1329578) q[0];
rz(-2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.7763604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62194659) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(1.7072862) q[0];
rz(0.021990701) q[2];
sx q[2];
rz(-1.6985661) q[2];
sx q[2];
rz(2.2760504) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8147162) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(0.53175064) q[1];
x q[2];
rz(-0.39799277) q[3];
sx q[3];
rz(-2.6905305) q[3];
sx q[3];
rz(3.1068902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0662971) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(-0.42262849) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(-2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-0.92393595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5101178) q[0];
sx q[0];
rz(-0.81561136) q[0];
sx q[0];
rz(0.60879137) q[0];
rz(-pi) q[1];
rz(0.48859476) q[2];
sx q[2];
rz(-1.5602419) q[2];
sx q[2];
rz(-0.34305629) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3194771) q[1];
sx q[1];
rz(-1.176782) q[1];
sx q[1];
rz(-2.134321) q[1];
rz(-2.4531035) q[3];
sx q[3];
rz(-2.2985035) q[3];
sx q[3];
rz(1.9635995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(-1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37344638) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(-1.3884397) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-3.1034234) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725806) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(1.5241745) q[0];
rz(-pi) q[1];
rz(-1.795904) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(0.12347808) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.888615) q[1];
sx q[1];
rz(-1.9441609) q[1];
sx q[1];
rz(3.0967767) q[1];
x q[2];
rz(-2.2323654) q[3];
sx q[3];
rz(-2.0709166) q[3];
sx q[3];
rz(1.7928894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(1.4981131) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-3.0315728) q[2];
sx q[2];
rz(-1.6286055) q[2];
sx q[2];
rz(1.5946228) q[2];
rz(0.86919541) q[3];
sx q[3];
rz(-1.8164608) q[3];
sx q[3];
rz(1.7432004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
