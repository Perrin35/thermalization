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
rz(1.5384742) q[0];
sx q[0];
rz(-1.3863907) q[0];
sx q[0];
rz(-1.2586799) q[0];
rz(1.4409244) q[1];
sx q[1];
rz(-2.6135593) q[1];
sx q[1];
rz(0.33906403) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.514076) q[0];
sx q[0];
rz(-1.7154804) q[0];
sx q[0];
rz(1.9067212) q[0];
rz(-pi) q[1];
rz(-2.4615844) q[2];
sx q[2];
rz(-1.8915904) q[2];
sx q[2];
rz(-2.8217905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2175581) q[1];
sx q[1];
rz(-1.6372895) q[1];
sx q[1];
rz(-2.1832863) q[1];
x q[2];
rz(2.4851299) q[3];
sx q[3];
rz(-1.5243666) q[3];
sx q[3];
rz(-1.186478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6429834) q[2];
sx q[2];
rz(-2.6080841) q[2];
sx q[2];
rz(1.741893) q[2];
rz(1.1534322) q[3];
sx q[3];
rz(-1.6250936) q[3];
sx q[3];
rz(-1.4199408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145849) q[0];
sx q[0];
rz(-1.4684533) q[0];
sx q[0];
rz(-0.46838316) q[0];
rz(-2.1539457) q[1];
sx q[1];
rz(-1.7665607) q[1];
sx q[1];
rz(2.1013026) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0752234) q[0];
sx q[0];
rz(-2.1589737) q[0];
sx q[0];
rz(0.69216125) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78590337) q[2];
sx q[2];
rz(-2.6214947) q[2];
sx q[2];
rz(1.6883858) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2548488) q[1];
sx q[1];
rz(-0.4720531) q[1];
sx q[1];
rz(0.2173425) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3181245) q[3];
sx q[3];
rz(-0.57909617) q[3];
sx q[3];
rz(-1.1163348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4686244) q[2];
sx q[2];
rz(-2.0277703) q[2];
sx q[2];
rz(1.4581397) q[2];
rz(-0.80312076) q[3];
sx q[3];
rz(-1.2642658) q[3];
sx q[3];
rz(-1.3706114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2405546) q[0];
sx q[0];
rz(-2.5032208) q[0];
sx q[0];
rz(1.4580131) q[0];
rz(1.1395678) q[1];
sx q[1];
rz(-1.640806) q[1];
sx q[1];
rz(-2.4511888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6656356) q[0];
sx q[0];
rz(-2.1607669) q[0];
sx q[0];
rz(2.0796805) q[0];
x q[1];
rz(-1.5883303) q[2];
sx q[2];
rz(-0.25871535) q[2];
sx q[2];
rz(2.0257547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8048794) q[1];
sx q[1];
rz(-1.1142528) q[1];
sx q[1];
rz(-1.9976569) q[1];
rz(-pi) q[2];
rz(1.6355945) q[3];
sx q[3];
rz(-2.4931723) q[3];
sx q[3];
rz(-2.5823926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0940501) q[2];
sx q[2];
rz(-1.8886781) q[2];
sx q[2];
rz(2.3599153) q[2];
rz(2.9983669) q[3];
sx q[3];
rz(-2.4476624) q[3];
sx q[3];
rz(-3.1128939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52611065) q[0];
sx q[0];
rz(-0.019793864) q[0];
sx q[0];
rz(-0.97446781) q[0];
rz(-3.114585) q[1];
sx q[1];
rz(-1.4911431) q[1];
sx q[1];
rz(0.19198051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5088599) q[0];
sx q[0];
rz(-1.8562537) q[0];
sx q[0];
rz(-2.767059) q[0];
rz(-pi) q[1];
rz(-1.7186195) q[2];
sx q[2];
rz(-1.333892) q[2];
sx q[2];
rz(-0.29659268) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7168107) q[1];
sx q[1];
rz(-1.5871986) q[1];
sx q[1];
rz(1.8540127) q[1];
x q[2];
rz(-2.8497898) q[3];
sx q[3];
rz(-0.58819669) q[3];
sx q[3];
rz(0.10313973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8359022) q[2];
sx q[2];
rz(-1.4985936) q[2];
sx q[2];
rz(-0.72176019) q[2];
rz(0.52779683) q[3];
sx q[3];
rz(-0.58776394) q[3];
sx q[3];
rz(1.9378978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.0498306) q[0];
sx q[0];
rz(-1.993655) q[0];
sx q[0];
rz(0.21035305) q[0];
rz(2.6166088) q[1];
sx q[1];
rz(-2.4139082) q[1];
sx q[1];
rz(-2.4997589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9729298) q[0];
sx q[0];
rz(-2.7073858) q[0];
sx q[0];
rz(-0.28247873) q[0];
rz(-1.3566229) q[2];
sx q[2];
rz(-1.3371981) q[2];
sx q[2];
rz(0.24326277) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0090897) q[1];
sx q[1];
rz(-2.9727474) q[1];
sx q[1];
rz(-1.7573207) q[1];
x q[2];
rz(0.028270311) q[3];
sx q[3];
rz(-1.8688374) q[3];
sx q[3];
rz(0.6606797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.019023808) q[2];
sx q[2];
rz(-1.4539315) q[2];
sx q[2];
rz(2.8008578) q[2];
rz(0.46180284) q[3];
sx q[3];
rz(-1.267649) q[3];
sx q[3];
rz(-0.99892282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13071624) q[0];
sx q[0];
rz(-2.0858481) q[0];
sx q[0];
rz(-0.51489949) q[0];
rz(-1.2841094) q[1];
sx q[1];
rz(-2.1045556) q[1];
sx q[1];
rz(-2.5426755) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0055858) q[0];
sx q[0];
rz(-1.4201627) q[0];
sx q[0];
rz(2.1120911) q[0];
rz(-pi) q[1];
rz(3.0294543) q[2];
sx q[2];
rz(-1.7742998) q[2];
sx q[2];
rz(0.017353915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22344124) q[1];
sx q[1];
rz(-2.591003) q[1];
sx q[1];
rz(-1.7718901) q[1];
rz(-1.9768561) q[3];
sx q[3];
rz(-0.58771509) q[3];
sx q[3];
rz(-0.57460659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14737031) q[2];
sx q[2];
rz(-1.3753336) q[2];
sx q[2];
rz(-2.4293409) q[2];
rz(-2.0272523) q[3];
sx q[3];
rz(-2.2388191) q[3];
sx q[3];
rz(0.72079349) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0459868) q[0];
sx q[0];
rz(-1.1271891) q[0];
sx q[0];
rz(1.4229232) q[0];
rz(-1.0320041) q[1];
sx q[1];
rz(-1.874066) q[1];
sx q[1];
rz(1.2123908) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.732837) q[0];
sx q[0];
rz(-1.8538332) q[0];
sx q[0];
rz(-2.416859) q[0];
rz(-pi) q[1];
rz(2.5676377) q[2];
sx q[2];
rz(-2.7864153) q[2];
sx q[2];
rz(1.4188302) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8960986) q[1];
sx q[1];
rz(-2.3521455) q[1];
sx q[1];
rz(1.1723799) q[1];
x q[2];
rz(-0.79030307) q[3];
sx q[3];
rz(-0.41158446) q[3];
sx q[3];
rz(-1.5800588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57264868) q[2];
sx q[2];
rz(-1.135301) q[2];
sx q[2];
rz(2.9917319) q[2];
rz(-0.0830689) q[3];
sx q[3];
rz(-0.84851256) q[3];
sx q[3];
rz(-1.2868631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5158373) q[0];
sx q[0];
rz(-2.5560162) q[0];
sx q[0];
rz(2.9199912) q[0];
rz(-1.340647) q[1];
sx q[1];
rz(-2.4602175) q[1];
sx q[1];
rz(0.63366205) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0283141) q[0];
sx q[0];
rz(-1.6080556) q[0];
sx q[0];
rz(1.5118044) q[0];
rz(-pi) q[1];
rz(-0.35291283) q[2];
sx q[2];
rz(-0.9448337) q[2];
sx q[2];
rz(-2.7865648) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.934856) q[1];
sx q[1];
rz(-2.2168787) q[1];
sx q[1];
rz(1.6621432) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.996534) q[3];
sx q[3];
rz(-0.71509711) q[3];
sx q[3];
rz(1.7760269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0076912045) q[2];
sx q[2];
rz(-2.2654221) q[2];
sx q[2];
rz(0.42789856) q[2];
rz(-2.2513572) q[3];
sx q[3];
rz(-0.72042239) q[3];
sx q[3];
rz(-3.1261352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4175005) q[0];
sx q[0];
rz(-0.98560792) q[0];
sx q[0];
rz(-0.86522657) q[0];
rz(-0.73156196) q[1];
sx q[1];
rz(-2.1831436) q[1];
sx q[1];
rz(-0.98347121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9051233) q[0];
sx q[0];
rz(-0.60689245) q[0];
sx q[0];
rz(-2.7913362) q[0];
x q[1];
rz(1.3147361) q[2];
sx q[2];
rz(-1.4863925) q[2];
sx q[2];
rz(-2.7483482) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.379885) q[1];
sx q[1];
rz(-1.7049676) q[1];
sx q[1];
rz(-0.11685808) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2526475) q[3];
sx q[3];
rz(-1.3414012) q[3];
sx q[3];
rz(1.655024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5359155) q[2];
sx q[2];
rz(-1.8226049) q[2];
sx q[2];
rz(0.21161045) q[2];
rz(-1.051988) q[3];
sx q[3];
rz(-2.7973723) q[3];
sx q[3];
rz(0.91867623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4695404) q[0];
sx q[0];
rz(-2.7993027) q[0];
sx q[0];
rz(2.4218609) q[0];
rz(1.4620048) q[1];
sx q[1];
rz(-0.77762496) q[1];
sx q[1];
rz(-0.07196149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77047399) q[0];
sx q[0];
rz(-1.8544056) q[0];
sx q[0];
rz(-0.78259612) q[0];
rz(1.8481723) q[2];
sx q[2];
rz(-0.69239834) q[2];
sx q[2];
rz(-1.5497249) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.601974) q[1];
sx q[1];
rz(-1.7672136) q[1];
sx q[1];
rz(0.41069527) q[1];
rz(-pi) q[2];
rz(0.34787223) q[3];
sx q[3];
rz(-1.8765896) q[3];
sx q[3];
rz(-2.3527462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0993242) q[2];
sx q[2];
rz(-0.86809731) q[2];
sx q[2];
rz(-0.68359366) q[2];
rz(-1.7259224) q[3];
sx q[3];
rz(-0.96368027) q[3];
sx q[3];
rz(1.2832618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0409446) q[0];
sx q[0];
rz(-1.9293979) q[0];
sx q[0];
rz(1.9579493) q[0];
rz(2.4965141) q[1];
sx q[1];
rz(-1.3007785) q[1];
sx q[1];
rz(0.29161463) q[1];
rz(-1.8963457) q[2];
sx q[2];
rz(-1.799188) q[2];
sx q[2];
rz(-1.3273018) q[2];
rz(2.162289) q[3];
sx q[3];
rz(-0.42604705) q[3];
sx q[3];
rz(2.7548509) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
