OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.32970348) q[0];
sx q[0];
rz(-2.831037) q[0];
sx q[0];
rz(-1.1986873) q[0];
rz(-2.0074453) q[1];
sx q[1];
rz(-0.79678798) q[1];
sx q[1];
rz(-2.8400583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.806592) q[0];
sx q[0];
rz(-1.3863239) q[0];
sx q[0];
rz(-2.9770699) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8704988) q[2];
sx q[2];
rz(-1.3915919) q[2];
sx q[2];
rz(-1.186893) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5079699) q[1];
sx q[1];
rz(-2.3862641) q[1];
sx q[1];
rz(1.4719187) q[1];
rz(-1.1490378) q[3];
sx q[3];
rz(-0.43614498) q[3];
sx q[3];
rz(2.8835322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6887168) q[2];
sx q[2];
rz(-1.0998925) q[2];
sx q[2];
rz(-2.0067046) q[2];
rz(-3.0196043) q[3];
sx q[3];
rz(-2.205409) q[3];
sx q[3];
rz(-3.0854935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6720471) q[0];
sx q[0];
rz(-2.8944954) q[0];
sx q[0];
rz(0.0023512996) q[0];
rz(2.7408842) q[1];
sx q[1];
rz(-1.7727163) q[1];
sx q[1];
rz(2.2854663) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5120289) q[0];
sx q[0];
rz(-2.1855547) q[0];
sx q[0];
rz(2.828086) q[0];
x q[1];
rz(2.8555774) q[2];
sx q[2];
rz(-1.0723503) q[2];
sx q[2];
rz(-1.1666544) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5634671) q[1];
sx q[1];
rz(-1.5208066) q[1];
sx q[1];
rz(2.0983138) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5991531) q[3];
sx q[3];
rz(-0.63839165) q[3];
sx q[3];
rz(-2.8505489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1043642) q[2];
sx q[2];
rz(-1.9375485) q[2];
sx q[2];
rz(-0.47002235) q[2];
rz(0.59703279) q[3];
sx q[3];
rz(-3.0266914) q[3];
sx q[3];
rz(-1.443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700579) q[0];
sx q[0];
rz(-0.4929339) q[0];
sx q[0];
rz(2.9135627) q[0];
rz(2.6920964) q[1];
sx q[1];
rz(-2.1411965) q[1];
sx q[1];
rz(-2.1953348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11433521) q[0];
sx q[0];
rz(-1.9693621) q[0];
sx q[0];
rz(-0.27862676) q[0];
x q[1];
rz(0.7545289) q[2];
sx q[2];
rz(-0.81475329) q[2];
sx q[2];
rz(0.94902869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11383948) q[1];
sx q[1];
rz(-1.7530754) q[1];
sx q[1];
rz(-2.7992642) q[1];
rz(-1.2872996) q[3];
sx q[3];
rz(-1.722966) q[3];
sx q[3];
rz(2.7622472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9049282) q[2];
sx q[2];
rz(-1.4651352) q[2];
sx q[2];
rz(1.6620592) q[2];
rz(2.9605401) q[3];
sx q[3];
rz(-2.3513887) q[3];
sx q[3];
rz(-0.886206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7774696) q[0];
sx q[0];
rz(-1.5793707) q[0];
sx q[0];
rz(-0.89853483) q[0];
rz(-0.50615519) q[1];
sx q[1];
rz(-1.8471142) q[1];
sx q[1];
rz(1.6815965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130044) q[0];
sx q[0];
rz(-0.53697907) q[0];
sx q[0];
rz(-3.0130638) q[0];
rz(-pi) q[1];
rz(-2.0115667) q[2];
sx q[2];
rz(-0.87723644) q[2];
sx q[2];
rz(-2.9206985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4913702) q[1];
sx q[1];
rz(-1.4563602) q[1];
sx q[1];
rz(3.070862) q[1];
rz(3.0811062) q[3];
sx q[3];
rz(-0.91676676) q[3];
sx q[3];
rz(0.56246434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26607457) q[2];
sx q[2];
rz(-0.47553277) q[2];
sx q[2];
rz(1.5024705) q[2];
rz(1.8858887) q[3];
sx q[3];
rz(-1.4016822) q[3];
sx q[3];
rz(-0.046253117) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083549) q[0];
sx q[0];
rz(-0.73035705) q[0];
sx q[0];
rz(-2.9610942) q[0];
rz(-1.3795229) q[1];
sx q[1];
rz(-0.5677529) q[1];
sx q[1];
rz(1.0951805) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8081854) q[0];
sx q[0];
rz(-1.6327724) q[0];
sx q[0];
rz(-0.13570762) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0600431) q[2];
sx q[2];
rz(-2.3429541) q[2];
sx q[2];
rz(2.7872686) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3928581) q[1];
sx q[1];
rz(-1.431141) q[1];
sx q[1];
rz(-0.77844324) q[1];
rz(-pi) q[2];
rz(-0.13325444) q[3];
sx q[3];
rz(-2.2550003) q[3];
sx q[3];
rz(1.9397936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77202648) q[2];
sx q[2];
rz(-2.2942746) q[2];
sx q[2];
rz(-1.0687211) q[2];
rz(0.42701834) q[3];
sx q[3];
rz(-0.87978274) q[3];
sx q[3];
rz(-1.8900185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1123947) q[0];
sx q[0];
rz(-0.91570941) q[0];
sx q[0];
rz(-0.29119626) q[0];
rz(1.4578106) q[1];
sx q[1];
rz(-2.1144512) q[1];
sx q[1];
rz(2.2529032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55808961) q[0];
sx q[0];
rz(-2.0348685) q[0];
sx q[0];
rz(0.91832535) q[0];
rz(3.1142919) q[2];
sx q[2];
rz(-2.3841287) q[2];
sx q[2];
rz(-2.2006019) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86119628) q[1];
sx q[1];
rz(-0.79842192) q[1];
sx q[1];
rz(0.75099985) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8055775) q[3];
sx q[3];
rz(-0.42308261) q[3];
sx q[3];
rz(-2.7132636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39151829) q[2];
sx q[2];
rz(-1.3041648) q[2];
sx q[2];
rz(-0.94775003) q[2];
rz(-2.9446972) q[3];
sx q[3];
rz(-0.24053776) q[3];
sx q[3];
rz(2.8314364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.8295558) q[0];
sx q[0];
rz(-1.1126248) q[0];
sx q[0];
rz(2.6716676) q[0];
rz(1.6795109) q[1];
sx q[1];
rz(-1.8406248) q[1];
sx q[1];
rz(1.7376815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4257098) q[0];
sx q[0];
rz(-2.4615767) q[0];
sx q[0];
rz(0.99773565) q[0];
rz(-1.092488) q[2];
sx q[2];
rz(-1.385139) q[2];
sx q[2];
rz(2.3379679) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0216071) q[1];
sx q[1];
rz(-2.48263) q[1];
sx q[1];
rz(0.95665415) q[1];
rz(-1.2265731) q[3];
sx q[3];
rz(-2.1267517) q[3];
sx q[3];
rz(-1.8575625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.64138428) q[2];
sx q[2];
rz(-2.310014) q[2];
sx q[2];
rz(-3.0339962) q[2];
rz(-0.61947668) q[3];
sx q[3];
rz(-1.8457125) q[3];
sx q[3];
rz(1.3639601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3951931) q[0];
sx q[0];
rz(-1.0522319) q[0];
sx q[0];
rz(-2.8885544) q[0];
rz(3.053275) q[1];
sx q[1];
rz(-2.2679236) q[1];
sx q[1];
rz(-0.94892445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.388577) q[0];
sx q[0];
rz(-2.8070794) q[0];
sx q[0];
rz(-1.0981513) q[0];
x q[1];
rz(0.19472943) q[2];
sx q[2];
rz(-1.1695332) q[2];
sx q[2];
rz(-1.8068318) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.526405) q[1];
sx q[1];
rz(-1.624409) q[1];
sx q[1];
rz(1.0907021) q[1];
rz(-pi) q[2];
rz(-2.3036495) q[3];
sx q[3];
rz(-2.4334014) q[3];
sx q[3];
rz(1.4670682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39763149) q[2];
sx q[2];
rz(-2.2546015) q[2];
sx q[2];
rz(0.31069791) q[2];
rz(0.24142309) q[3];
sx q[3];
rz(-1.9390691) q[3];
sx q[3];
rz(-2.0588622) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13984118) q[0];
sx q[0];
rz(-2.7353291) q[0];
sx q[0];
rz(-1.9061506) q[0];
rz(0.48108092) q[1];
sx q[1];
rz(-2.1886539) q[1];
sx q[1];
rz(-1.4280041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3612551) q[0];
sx q[0];
rz(-1.5596103) q[0];
sx q[0];
rz(-0.7594603) q[0];
x q[1];
rz(-2.9129006) q[2];
sx q[2];
rz(-1.2989825) q[2];
sx q[2];
rz(-0.52955176) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.42983152) q[1];
sx q[1];
rz(-1.3138384) q[1];
sx q[1];
rz(0.25196948) q[1];
x q[2];
rz(1.8974182) q[3];
sx q[3];
rz(-0.87509586) q[3];
sx q[3];
rz(0.43750924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1741751) q[2];
sx q[2];
rz(-0.82505161) q[2];
sx q[2];
rz(0.027912557) q[2];
rz(0.86999718) q[3];
sx q[3];
rz(-2.3614707) q[3];
sx q[3];
rz(1.1869259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-3.1160527) q[0];
sx q[0];
rz(-1.5929796) q[0];
sx q[0];
rz(1.0593587) q[0];
rz(-1.7268044) q[1];
sx q[1];
rz(-1.6635514) q[1];
sx q[1];
rz(-0.47763225) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1519417) q[0];
sx q[0];
rz(-0.7583337) q[0];
sx q[0];
rz(0.21679057) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5651836) q[2];
sx q[2];
rz(-0.33585784) q[2];
sx q[2];
rz(-0.72959057) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2207909) q[1];
sx q[1];
rz(-3.0169636) q[1];
sx q[1];
rz(1.1225379) q[1];
rz(-pi) q[2];
rz(-0.82179339) q[3];
sx q[3];
rz(-2.0374277) q[3];
sx q[3];
rz(0.94576437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.16964218) q[2];
sx q[2];
rz(-2.1596238) q[2];
sx q[2];
rz(-2.7726445) q[2];
rz(-1.87489) q[3];
sx q[3];
rz(-2.4167175) q[3];
sx q[3];
rz(-0.96316159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13851588) q[0];
sx q[0];
rz(-1.930548) q[0];
sx q[0];
rz(1.7532274) q[0];
rz(2.6970462) q[1];
sx q[1];
rz(-0.81041705) q[1];
sx q[1];
rz(3.0141426) q[1];
rz(2.6981392) q[2];
sx q[2];
rz(-1.2724066) q[2];
sx q[2];
rz(2.3497697) q[2];
rz(-3.1034341) q[3];
sx q[3];
rz(-0.46783075) q[3];
sx q[3];
rz(-0.086188407) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
