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
rz(2.1287542) q[0];
sx q[0];
rz(-1.1829809) q[0];
sx q[0];
rz(-2.8943789) q[0];
rz(-1.5850868) q[1];
sx q[1];
rz(4.1559846) q[1];
sx q[1];
rz(13.520887) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92468372) q[0];
sx q[0];
rz(-1.1966853) q[0];
sx q[0];
rz(3.1103136) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11516042) q[2];
sx q[2];
rz(-0.45591337) q[2];
sx q[2];
rz(-1.0418089) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87434972) q[1];
sx q[1];
rz(-1.2552208) q[1];
sx q[1];
rz(-2.3467031) q[1];
rz(1.263721) q[3];
sx q[3];
rz(-2.319984) q[3];
sx q[3];
rz(1.6916569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9556094) q[2];
sx q[2];
rz(-0.60638967) q[2];
sx q[2];
rz(0.28140226) q[2];
rz(-1.227281) q[3];
sx q[3];
rz(-1.809092) q[3];
sx q[3];
rz(1.2138155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94673741) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(-0.58709225) q[0];
rz(-0.52014822) q[1];
sx q[1];
rz(-1.2112434) q[1];
sx q[1];
rz(0.21336666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30946975) q[0];
sx q[0];
rz(-1.42064) q[0];
sx q[0];
rz(3.0792159) q[0];
x q[1];
rz(0.92170282) q[2];
sx q[2];
rz(-1.9919167) q[2];
sx q[2];
rz(-2.365961) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0736463) q[1];
sx q[1];
rz(-1.526462) q[1];
sx q[1];
rz(1.4774051) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1150241) q[3];
sx q[3];
rz(-2.5682862) q[3];
sx q[3];
rz(-0.70070964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1063891) q[2];
sx q[2];
rz(-2.8530402) q[2];
sx q[2];
rz(-2.8698548) q[2];
rz(0.64618293) q[3];
sx q[3];
rz(-1.313442) q[3];
sx q[3];
rz(-1.9338231) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3216517) q[0];
sx q[0];
rz(-2.2714748) q[0];
sx q[0];
rz(0.2680378) q[0];
rz(-2.9053814) q[1];
sx q[1];
rz(-1.3583438) q[1];
sx q[1];
rz(1.7265629) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9504323) q[0];
sx q[0];
rz(-1.6163905) q[0];
sx q[0];
rz(-0.10204236) q[0];
rz(-pi) q[1];
rz(-0.36249749) q[2];
sx q[2];
rz(-1.4115184) q[2];
sx q[2];
rz(-0.25937072) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8047976) q[1];
sx q[1];
rz(-2.6141254) q[1];
sx q[1];
rz(2.5673366) q[1];
rz(-pi) q[2];
rz(-2.799166) q[3];
sx q[3];
rz(-2.4391616) q[3];
sx q[3];
rz(-2.680955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.95431027) q[2];
sx q[2];
rz(-2.8128746) q[2];
sx q[2];
rz(1.7598565) q[2];
rz(-0.36661822) q[3];
sx q[3];
rz(-1.4724052) q[3];
sx q[3];
rz(-0.097675145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919375) q[0];
sx q[0];
rz(-0.20579919) q[0];
sx q[0];
rz(3.1254712) q[0];
rz(1.6751809) q[1];
sx q[1];
rz(-1.5212395) q[1];
sx q[1];
rz(-0.88502562) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2766874) q[0];
sx q[0];
rz(-1.7430796) q[0];
sx q[0];
rz(2.7322463) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22815223) q[2];
sx q[2];
rz(-1.9288262) q[2];
sx q[2];
rz(-1.4756853) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4964543) q[1];
sx q[1];
rz(-1.1119324) q[1];
sx q[1];
rz(1.9978844) q[1];
rz(0.4039558) q[3];
sx q[3];
rz(-1.0896297) q[3];
sx q[3];
rz(1.6649401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5625988) q[2];
sx q[2];
rz(-2.6851974) q[2];
sx q[2];
rz(1.558051) q[2];
rz(-1.3715749) q[3];
sx q[3];
rz(-1.2746425) q[3];
sx q[3];
rz(-0.19545999) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564932) q[0];
sx q[0];
rz(-1.3753966) q[0];
sx q[0];
rz(2.9920355) q[0];
rz(2.0984446) q[1];
sx q[1];
rz(-0.97277343) q[1];
sx q[1];
rz(-2.3057888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5476033) q[0];
sx q[0];
rz(-2.7271252) q[0];
sx q[0];
rz(0.91048059) q[0];
rz(0.88045516) q[2];
sx q[2];
rz(-1.5780851) q[2];
sx q[2];
rz(-1.4237504) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0936613) q[1];
sx q[1];
rz(-1.3627628) q[1];
sx q[1];
rz(3.1407977) q[1];
rz(0.18330611) q[3];
sx q[3];
rz(-1.3245031) q[3];
sx q[3];
rz(3.1251647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0189556) q[2];
sx q[2];
rz(-1.6654207) q[2];
sx q[2];
rz(0.41610757) q[2];
rz(0.082402669) q[3];
sx q[3];
rz(-0.96937886) q[3];
sx q[3];
rz(0.20127067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06269726) q[0];
sx q[0];
rz(-1.0286999) q[0];
sx q[0];
rz(-1.0908352) q[0];
rz(-1.1385607) q[1];
sx q[1];
rz(-1.9416315) q[1];
sx q[1];
rz(-2.535215) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17601062) q[0];
sx q[0];
rz(-0.93744288) q[0];
sx q[0];
rz(1.8853582) q[0];
x q[1];
rz(-0.44754812) q[2];
sx q[2];
rz(-2.2879763) q[2];
sx q[2];
rz(0.56499519) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5412288) q[1];
sx q[1];
rz(-2.0906587) q[1];
sx q[1];
rz(-2.4281349) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6929469) q[3];
sx q[3];
rz(-2.2854317) q[3];
sx q[3];
rz(-0.84343101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.31327569) q[2];
sx q[2];
rz(-0.54855359) q[2];
sx q[2];
rz(-0.11087785) q[2];
rz(-1.5634792) q[3];
sx q[3];
rz(-2.9302247) q[3];
sx q[3];
rz(1.8160688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368847) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(-0.65716609) q[0];
rz(-1.3232629) q[1];
sx q[1];
rz(-2.3902049) q[1];
sx q[1];
rz(1.125186) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1745963) q[0];
sx q[0];
rz(-1.5856992) q[0];
sx q[0];
rz(-3.1257939) q[0];
rz(-pi) q[1];
rz(1.6680587) q[2];
sx q[2];
rz(-1.2583557) q[2];
sx q[2];
rz(0.51431235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4954056) q[1];
sx q[1];
rz(-2.0584724) q[1];
sx q[1];
rz(-1.133092) q[1];
rz(-1.8784889) q[3];
sx q[3];
rz(-2.2224226) q[3];
sx q[3];
rz(-3.0857112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7323759) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(-0.96119514) q[2];
rz(-0.0088648908) q[3];
sx q[3];
rz(-2.2564087) q[3];
sx q[3];
rz(-1.9302906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0601785) q[0];
sx q[0];
rz(-0.1135122) q[0];
sx q[0];
rz(2.7139582) q[0];
rz(1.112452) q[1];
sx q[1];
rz(-1.3149657) q[1];
sx q[1];
rz(2.1449259) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450981) q[0];
sx q[0];
rz(-1.0834435) q[0];
sx q[0];
rz(2.4250406) q[0];
x q[1];
rz(-1.8420503) q[2];
sx q[2];
rz(-0.82640663) q[2];
sx q[2];
rz(2.8959993) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8045051) q[1];
sx q[1];
rz(-1.8356016) q[1];
sx q[1];
rz(-0.9881641) q[1];
x q[2];
rz(-0.31672041) q[3];
sx q[3];
rz(-0.59419981) q[3];
sx q[3];
rz(2.7946928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4619649) q[2];
sx q[2];
rz(-0.5373911) q[2];
sx q[2];
rz(-2.4930387) q[2];
rz(0.31351659) q[3];
sx q[3];
rz(-2.2532012) q[3];
sx q[3];
rz(-3.0891109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.40470966) q[0];
sx q[0];
rz(-0.22060224) q[0];
sx q[0];
rz(2.1746461) q[0];
rz(2.1838358) q[1];
sx q[1];
rz(-1.0916595) q[1];
sx q[1];
rz(-2.0584094) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5388018) q[0];
sx q[0];
rz(-1.1447009) q[0];
sx q[0];
rz(1.8417131) q[0];
rz(0.90980977) q[2];
sx q[2];
rz(-0.95558117) q[2];
sx q[2];
rz(-1.3843975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5450648) q[1];
sx q[1];
rz(-1.9633496) q[1];
sx q[1];
rz(2.3190581) q[1];
rz(-pi) q[2];
rz(3.1414491) q[3];
sx q[3];
rz(-0.4793491) q[3];
sx q[3];
rz(1.0841313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0208685) q[2];
sx q[2];
rz(-0.82896295) q[2];
sx q[2];
rz(-1.6797569) q[2];
rz(-2.7802137) q[3];
sx q[3];
rz(-1.8166108) q[3];
sx q[3];
rz(-0.98226515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35230377) q[0];
sx q[0];
rz(-2.4109349) q[0];
sx q[0];
rz(-0.58255449) q[0];
rz(1.3652868) q[1];
sx q[1];
rz(-1.0368232) q[1];
sx q[1];
rz(-2.9127311) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.992975) q[0];
sx q[0];
rz(-1.9662736) q[0];
sx q[0];
rz(1.3735885) q[0];
x q[1];
rz(2.9624356) q[2];
sx q[2];
rz(-1.3115885) q[2];
sx q[2];
rz(1.955206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6534253) q[1];
sx q[1];
rz(-1.5060802) q[1];
sx q[1];
rz(1.4117227) q[1];
rz(-pi) q[2];
rz(1.1353605) q[3];
sx q[3];
rz(-1.9276016) q[3];
sx q[3];
rz(-2.8716068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87054306) q[2];
sx q[2];
rz(-2.5898263) q[2];
sx q[2];
rz(-1.40847) q[2];
rz(-2.4208141) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(-1.6225947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085717399) q[0];
sx q[0];
rz(-1.7039104) q[0];
sx q[0];
rz(-1.3088551) q[0];
rz(-1.5695288) q[1];
sx q[1];
rz(-0.61620284) q[1];
sx q[1];
rz(-2.9175704) q[1];
rz(-1.9523377) q[2];
sx q[2];
rz(-2.5425662) q[2];
sx q[2];
rz(-2.7903916) q[2];
rz(0.15275501) q[3];
sx q[3];
rz(-1.412231) q[3];
sx q[3];
rz(-2.9291736) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
