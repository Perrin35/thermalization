OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6634231) q[0];
sx q[0];
rz(-2.2278251) q[0];
sx q[0];
rz(-1.0650286) q[0];
rz(-1.5364667) q[1];
sx q[1];
rz(4.2470266) q[1];
sx q[1];
rz(12.469835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4783515) q[0];
sx q[0];
rz(-2.6697201) q[0];
sx q[0];
rz(1.1058411) q[0];
rz(-pi) q[1];
rz(2.7235254) q[2];
sx q[2];
rz(-0.56116784) q[2];
sx q[2];
rz(-0.18218606) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21743821) q[1];
sx q[1];
rz(-2.4352142) q[1];
sx q[1];
rz(-2.6732754) q[1];
rz(3.0482951) q[3];
sx q[3];
rz(-0.99280706) q[3];
sx q[3];
rz(-1.597412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1221293) q[2];
sx q[2];
rz(-0.18522842) q[2];
sx q[2];
rz(-1.8786001) q[2];
rz(2.7428135) q[3];
sx q[3];
rz(-1.140927) q[3];
sx q[3];
rz(-2.1122475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(2.1102092) q[0];
sx q[0];
rz(-0.6321913) q[0];
sx q[0];
rz(-1.4612009) q[0];
rz(-0.070143135) q[1];
sx q[1];
rz(-1.5488011) q[1];
sx q[1];
rz(-0.4313012) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2099144) q[0];
sx q[0];
rz(-0.062206833) q[0];
sx q[0];
rz(-1.4181523) q[0];
x q[1];
rz(0.29745086) q[2];
sx q[2];
rz(-1.97142) q[2];
sx q[2];
rz(1.9971776) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0924393) q[1];
sx q[1];
rz(-0.91040694) q[1];
sx q[1];
rz(0.67164849) q[1];
rz(-pi) q[2];
x q[2];
rz(0.044593354) q[3];
sx q[3];
rz(-1.591288) q[3];
sx q[3];
rz(0.12497917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.18616072) q[2];
sx q[2];
rz(-1.7308851) q[2];
sx q[2];
rz(-0.53039941) q[2];
rz(-1.0839328) q[3];
sx q[3];
rz(-1.089774) q[3];
sx q[3];
rz(1.8918022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0386061) q[0];
sx q[0];
rz(-2.5654721) q[0];
sx q[0];
rz(-1.2023793) q[0];
rz(2.1309958) q[1];
sx q[1];
rz(-1.3253515) q[1];
sx q[1];
rz(3.1114263) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8081431) q[0];
sx q[0];
rz(-1.5238058) q[0];
sx q[0];
rz(0.065568173) q[0];
rz(1.2717275) q[2];
sx q[2];
rz(-0.054882955) q[2];
sx q[2];
rz(1.043821) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.846993) q[1];
sx q[1];
rz(-1.2275808) q[1];
sx q[1];
rz(-1.4354475) q[1];
x q[2];
rz(-0.10134931) q[3];
sx q[3];
rz(-1.0004442) q[3];
sx q[3];
rz(2.1150908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61115894) q[2];
sx q[2];
rz(-2.9035089) q[2];
sx q[2];
rz(-0.37453026) q[2];
rz(-2.020906) q[3];
sx q[3];
rz(-1.6569542) q[3];
sx q[3];
rz(-2.4499031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.197914) q[0];
sx q[0];
rz(-1.2926084) q[0];
sx q[0];
rz(1.9669272) q[0];
rz(-1.8265751) q[1];
sx q[1];
rz(-1.1066133) q[1];
sx q[1];
rz(-0.16474251) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024285) q[0];
sx q[0];
rz(-1.0902349) q[0];
sx q[0];
rz(-0.64569999) q[0];
rz(-pi) q[1];
rz(-1.7354749) q[2];
sx q[2];
rz(-1.0633755) q[2];
sx q[2];
rz(-2.8230482) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2811198) q[1];
sx q[1];
rz(-1.9083284) q[1];
sx q[1];
rz(-2.0915477) q[1];
rz(-pi) q[2];
rz(1.4993787) q[3];
sx q[3];
rz(-1.5319699) q[3];
sx q[3];
rz(-0.13664548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4508007) q[2];
sx q[2];
rz(-1.6054634) q[2];
sx q[2];
rz(-0.50626051) q[2];
rz(2.6311503) q[3];
sx q[3];
rz(-0.78845316) q[3];
sx q[3];
rz(-2.3300664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4617758) q[0];
sx q[0];
rz(-2.8051069) q[0];
sx q[0];
rz(-0.80999723) q[0];
rz(-3.0432155) q[1];
sx q[1];
rz(-2.7521303) q[1];
sx q[1];
rz(0.32442763) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9085037) q[0];
sx q[0];
rz(-1.5688741) q[0];
sx q[0];
rz(-1.4483099) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7067028) q[2];
sx q[2];
rz(-1.3506442) q[2];
sx q[2];
rz(2.1957601) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.95448976) q[1];
sx q[1];
rz(-1.3924161) q[1];
sx q[1];
rz(-2.5599077) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47629997) q[3];
sx q[3];
rz(-1.3838963) q[3];
sx q[3];
rz(2.4812061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35679945) q[2];
sx q[2];
rz(-1.5609799) q[2];
sx q[2];
rz(0.55850935) q[2];
rz(-1.2950581) q[3];
sx q[3];
rz(-1.7241155) q[3];
sx q[3];
rz(-2.7891125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9117987) q[0];
sx q[0];
rz(-2.7138382) q[0];
sx q[0];
rz(1.3707004) q[0];
rz(1.3937344) q[1];
sx q[1];
rz(-1.0153898) q[1];
sx q[1];
rz(2.9975991) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32329473) q[0];
sx q[0];
rz(-3.0157964) q[0];
sx q[0];
rz(-2.1677141) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0387467) q[2];
sx q[2];
rz(-2.6399595) q[2];
sx q[2];
rz(-2.9930121) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3578971) q[1];
sx q[1];
rz(-1.6829957) q[1];
sx q[1];
rz(2.1782801) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9861476) q[3];
sx q[3];
rz(-2.1710733) q[3];
sx q[3];
rz(2.5714175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7250942) q[2];
sx q[2];
rz(-2.7855253) q[2];
sx q[2];
rz(0.98102513) q[2];
rz(2.6128926) q[3];
sx q[3];
rz(-1.3542391) q[3];
sx q[3];
rz(-0.56021148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8037146) q[0];
sx q[0];
rz(-0.80334544) q[0];
sx q[0];
rz(-1.9377608) q[0];
rz(-0.0040668049) q[1];
sx q[1];
rz(-1.4264359) q[1];
sx q[1];
rz(-2.8003661) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2409397) q[0];
sx q[0];
rz(-1.1497069) q[0];
sx q[0];
rz(0.76068855) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1401743) q[2];
sx q[2];
rz(-2.2621415) q[2];
sx q[2];
rz(-1.6306016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3918889) q[1];
sx q[1];
rz(-0.14924696) q[1];
sx q[1];
rz(3.0796389) q[1];
rz(-0.34269603) q[3];
sx q[3];
rz(-1.7083941) q[3];
sx q[3];
rz(0.15145853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61871201) q[2];
sx q[2];
rz(-1.2423542) q[2];
sx q[2];
rz(-2.8704571) q[2];
rz(-1.6297657) q[3];
sx q[3];
rz(-2.1745493) q[3];
sx q[3];
rz(2.7395774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66255823) q[0];
sx q[0];
rz(-2.5738578) q[0];
sx q[0];
rz(0.16831368) q[0];
rz(-2.1315101) q[1];
sx q[1];
rz(-1.964566) q[1];
sx q[1];
rz(-3.041306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4254978) q[0];
sx q[0];
rz(-1.8779813) q[0];
sx q[0];
rz(1.2745122) q[0];
rz(-pi) q[1];
rz(1.3919207) q[2];
sx q[2];
rz(-0.84576195) q[2];
sx q[2];
rz(-2.7357312) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0834962) q[1];
sx q[1];
rz(-1.5462977) q[1];
sx q[1];
rz(-0.20170881) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0151305) q[3];
sx q[3];
rz(-2.0499886) q[3];
sx q[3];
rz(-1.5492619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4526796) q[2];
sx q[2];
rz(-0.49892384) q[2];
sx q[2];
rz(1.4208687) q[2];
rz(-1.0639327) q[3];
sx q[3];
rz(-1.7199793) q[3];
sx q[3];
rz(-3.0498144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38561472) q[0];
sx q[0];
rz(-1.4571964) q[0];
sx q[0];
rz(0.036855999) q[0];
rz(1.8846177) q[1];
sx q[1];
rz(-1.5024065) q[1];
sx q[1];
rz(-1.3704376) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12808558) q[0];
sx q[0];
rz(-2.8261746) q[0];
sx q[0];
rz(2.0924139) q[0];
rz(-pi) q[1];
rz(1.7212935) q[2];
sx q[2];
rz(-1.2062819) q[2];
sx q[2];
rz(-2.8340659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3664497) q[1];
sx q[1];
rz(-2.4058172) q[1];
sx q[1];
rz(-1.1476357) q[1];
x q[2];
rz(1.1537136) q[3];
sx q[3];
rz(-1.6110824) q[3];
sx q[3];
rz(-1.2800486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86159414) q[2];
sx q[2];
rz(-1.7656606) q[2];
sx q[2];
rz(-0.17246788) q[2];
rz(-0.55781588) q[3];
sx q[3];
rz(-2.6004801) q[3];
sx q[3];
rz(1.4403758) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7114792) q[0];
sx q[0];
rz(-2.5028296) q[0];
sx q[0];
rz(2.8016256) q[0];
rz(2.4760447) q[1];
sx q[1];
rz(-1.7219209) q[1];
sx q[1];
rz(-1.3131712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7198044) q[0];
sx q[0];
rz(-1.7504842) q[0];
sx q[0];
rz(-1.522904) q[0];
rz(-pi) q[1];
rz(1.4943141) q[2];
sx q[2];
rz(-1.8308365) q[2];
sx q[2];
rz(1.2225012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.03913) q[1];
sx q[1];
rz(-0.38929554) q[1];
sx q[1];
rz(-2.2401458) q[1];
rz(-pi) q[2];
rz(0.40295593) q[3];
sx q[3];
rz(-1.4266485) q[3];
sx q[3];
rz(-0.47553167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7790935) q[2];
sx q[2];
rz(-2.7069147) q[2];
sx q[2];
rz(-2.4707826) q[2];
rz(0.32495156) q[3];
sx q[3];
rz(-0.81348014) q[3];
sx q[3];
rz(-1.0130829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.308607) q[0];
sx q[0];
rz(-1.2815463) q[0];
sx q[0];
rz(1.4776342) q[0];
rz(-2.1113405) q[1];
sx q[1];
rz(-1.6696842) q[1];
sx q[1];
rz(-1.3546863) q[1];
rz(2.6162888) q[2];
sx q[2];
rz(-1.1266194) q[2];
sx q[2];
rz(2.5818679) q[2];
rz(0.58086953) q[3];
sx q[3];
rz(-1.8884251) q[3];
sx q[3];
rz(-2.6407218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
