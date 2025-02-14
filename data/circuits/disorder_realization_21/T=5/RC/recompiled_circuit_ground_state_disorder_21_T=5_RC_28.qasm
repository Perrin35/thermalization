OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2668827) q[0];
sx q[0];
rz(-1.0615791) q[0];
sx q[0];
rz(2.6258186) q[0];
rz(0.90142673) q[1];
sx q[1];
rz(-0.20144784) q[1];
sx q[1];
rz(2.1405061) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62681055) q[0];
sx q[0];
rz(-1.0812461) q[0];
sx q[0];
rz(0.63386699) q[0];
rz(-2.7269709) q[2];
sx q[2];
rz(-1.5765363) q[2];
sx q[2];
rz(-1.8513377) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6382003) q[1];
sx q[1];
rz(-0.93825785) q[1];
sx q[1];
rz(2.1293541) q[1];
rz(0.68925011) q[3];
sx q[3];
rz(-2.9514276) q[3];
sx q[3];
rz(2.3306757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72656816) q[2];
sx q[2];
rz(-0.51076204) q[2];
sx q[2];
rz(-1.6815574) q[2];
rz(0.36871746) q[3];
sx q[3];
rz(-1.5548778) q[3];
sx q[3];
rz(-1.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52861315) q[0];
sx q[0];
rz(-3.0280085) q[0];
sx q[0];
rz(-3.0533277) q[0];
rz(0.9322235) q[1];
sx q[1];
rz(-2.014522) q[1];
sx q[1];
rz(0.50311911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7403443) q[0];
sx q[0];
rz(-0.49763864) q[0];
sx q[0];
rz(-0.21792441) q[0];
x q[1];
rz(-1.8064151) q[2];
sx q[2];
rz(-2.4580147) q[2];
sx q[2];
rz(2.8238817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14993653) q[1];
sx q[1];
rz(-2.7767065) q[1];
sx q[1];
rz(-2.8349692) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6121268) q[3];
sx q[3];
rz(-2.240431) q[3];
sx q[3];
rz(2.1807699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3162389) q[2];
sx q[2];
rz(-2.1672858) q[2];
sx q[2];
rz(-0.45315722) q[2];
rz(0.00099269021) q[3];
sx q[3];
rz(-1.5627292) q[3];
sx q[3];
rz(2.4003975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2082625) q[0];
sx q[0];
rz(-0.18284155) q[0];
sx q[0];
rz(0.46874794) q[0];
rz(-2.5543429) q[1];
sx q[1];
rz(-1.8533555) q[1];
sx q[1];
rz(-2.8900878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6856988) q[0];
sx q[0];
rz(-0.60916466) q[0];
sx q[0];
rz(-2.3745911) q[0];
x q[1];
rz(1.3590888) q[2];
sx q[2];
rz(-1.9518753) q[2];
sx q[2];
rz(2.7104988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1826102) q[1];
sx q[1];
rz(-1.6291367) q[1];
sx q[1];
rz(-1.6946409) q[1];
rz(0.65966123) q[3];
sx q[3];
rz(-0.53852496) q[3];
sx q[3];
rz(2.7646567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.59715366) q[2];
sx q[2];
rz(-0.73828283) q[2];
sx q[2];
rz(0.047133751) q[2];
rz(2.2733287) q[3];
sx q[3];
rz(-1.2406113) q[3];
sx q[3];
rz(0.49938437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78859162) q[0];
sx q[0];
rz(-2.7356739) q[0];
sx q[0];
rz(1.3141919) q[0];
rz(0.81002533) q[1];
sx q[1];
rz(-1.516073) q[1];
sx q[1];
rz(2.8257418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0994713) q[0];
sx q[0];
rz(-2.1232623) q[0];
sx q[0];
rz(2.5312159) q[0];
rz(-pi) q[1];
rz(-1.3988109) q[2];
sx q[2];
rz(-1.2575527) q[2];
sx q[2];
rz(-1.9364374) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1882945) q[1];
sx q[1];
rz(-0.85894924) q[1];
sx q[1];
rz(-1.9602592) q[1];
x q[2];
rz(-2.9258435) q[3];
sx q[3];
rz(-2.5255754) q[3];
sx q[3];
rz(-2.8293138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6541859) q[2];
sx q[2];
rz(-2.5155289) q[2];
sx q[2];
rz(-2.249304) q[2];
rz(-1.4687294) q[3];
sx q[3];
rz(-2.0817616) q[3];
sx q[3];
rz(-2.0784126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47657087) q[0];
sx q[0];
rz(-0.63711089) q[0];
sx q[0];
rz(-2.2549905) q[0];
rz(0.84238148) q[1];
sx q[1];
rz(-1.3096755) q[1];
sx q[1];
rz(-1.1096035) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32217596) q[0];
sx q[0];
rz(-2.3482394) q[0];
sx q[0];
rz(0.94950907) q[0];
x q[1];
rz(-0.72317041) q[2];
sx q[2];
rz(-1.7077272) q[2];
sx q[2];
rz(-1.8214769) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.230727) q[1];
sx q[1];
rz(-2.334667) q[1];
sx q[1];
rz(0.74300503) q[1];
rz(-2.4585633) q[3];
sx q[3];
rz(-1.08324) q[3];
sx q[3];
rz(1.725793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2716918) q[2];
sx q[2];
rz(-0.39125189) q[2];
sx q[2];
rz(-0.58763495) q[2];
rz(-2.9247126) q[3];
sx q[3];
rz(-1.2809332) q[3];
sx q[3];
rz(1.3542401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.13123913) q[0];
sx q[0];
rz(-2.6709747) q[0];
sx q[0];
rz(1.3838029) q[0];
rz(0.17419392) q[1];
sx q[1];
rz(-1.1015588) q[1];
sx q[1];
rz(-1.9536288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096322358) q[0];
sx q[0];
rz(-2.850527) q[0];
sx q[0];
rz(-1.8610659) q[0];
x q[1];
rz(2.8415478) q[2];
sx q[2];
rz(-2.6159673) q[2];
sx q[2];
rz(0.92074652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7158051) q[1];
sx q[1];
rz(-1.6015698) q[1];
sx q[1];
rz(-1.4620859) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8781885) q[3];
sx q[3];
rz(-2.3114021) q[3];
sx q[3];
rz(2.8678558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69372988) q[2];
sx q[2];
rz(-1.8665946) q[2];
sx q[2];
rz(-1.0542144) q[2];
rz(2.5888455) q[3];
sx q[3];
rz(-0.25944969) q[3];
sx q[3];
rz(-2.0235846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066910557) q[0];
sx q[0];
rz(-2.5801165) q[0];
sx q[0];
rz(-0.018360227) q[0];
rz(1.8439937) q[1];
sx q[1];
rz(-1.3341787) q[1];
sx q[1];
rz(-1.1150572) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0783065) q[0];
sx q[0];
rz(-3.0111599) q[0];
sx q[0];
rz(0.13551573) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31324306) q[2];
sx q[2];
rz(-2.0038824) q[2];
sx q[2];
rz(-1.05561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12488481) q[1];
sx q[1];
rz(-0.77335268) q[1];
sx q[1];
rz(2.4340637) q[1];
rz(-pi) q[2];
rz(-2.2998718) q[3];
sx q[3];
rz(-2.5091565) q[3];
sx q[3];
rz(-0.034733437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4516419) q[2];
sx q[2];
rz(-0.74863282) q[2];
sx q[2];
rz(-0.79317036) q[2];
rz(0.66148174) q[3];
sx q[3];
rz(-1.7659148) q[3];
sx q[3];
rz(0.64600265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9328705) q[0];
sx q[0];
rz(-1.1247617) q[0];
sx q[0];
rz(0.17054184) q[0];
rz(1.0199245) q[1];
sx q[1];
rz(-0.41671697) q[1];
sx q[1];
rz(-2.4822617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.764544) q[0];
sx q[0];
rz(-1.555849) q[0];
sx q[0];
rz(0.055146761) q[0];
rz(-pi) q[1];
rz(-0.77240981) q[2];
sx q[2];
rz(-2.4669224) q[2];
sx q[2];
rz(-1.2090558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48261343) q[1];
sx q[1];
rz(-0.70090997) q[1];
sx q[1];
rz(-2.6644876) q[1];
x q[2];
rz(-2.5468656) q[3];
sx q[3];
rz(-0.85053125) q[3];
sx q[3];
rz(-1.8320465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3063804) q[2];
sx q[2];
rz(-2.7472718) q[2];
sx q[2];
rz(2.1739056) q[2];
rz(-2.5505998) q[3];
sx q[3];
rz(-1.3700181) q[3];
sx q[3];
rz(-1.437457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10385926) q[0];
sx q[0];
rz(-2.2028706) q[0];
sx q[0];
rz(-0.18873225) q[0];
rz(2.0813148) q[1];
sx q[1];
rz(-0.66245586) q[1];
sx q[1];
rz(2.2417384) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3093787) q[0];
sx q[0];
rz(-0.77118528) q[0];
sx q[0];
rz(0.82892771) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84328358) q[2];
sx q[2];
rz(-1.5688063) q[2];
sx q[2];
rz(0.28283027) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91119345) q[1];
sx q[1];
rz(-2.491103) q[1];
sx q[1];
rz(0.002928811) q[1];
rz(-0.84964801) q[3];
sx q[3];
rz(-1.7150522) q[3];
sx q[3];
rz(0.80286959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73682827) q[2];
sx q[2];
rz(-1.9177723) q[2];
sx q[2];
rz(-1.9723816) q[2];
rz(-0.87559187) q[3];
sx q[3];
rz(-2.4647522) q[3];
sx q[3];
rz(-3.0524047) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.635023) q[0];
sx q[0];
rz(-1.4590141) q[0];
sx q[0];
rz(2.7154679) q[0];
rz(-1.2395073) q[1];
sx q[1];
rz(-1.0303717) q[1];
sx q[1];
rz(-0.32136163) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5581767) q[0];
sx q[0];
rz(-0.9984278) q[0];
sx q[0];
rz(3.0757927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.05409492) q[2];
sx q[2];
rz(-2.0129497) q[2];
sx q[2];
rz(1.7345605) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4453955) q[1];
sx q[1];
rz(-1.9085669) q[1];
sx q[1];
rz(-2.9902451) q[1];
x q[2];
rz(1.0635183) q[3];
sx q[3];
rz(-2.2174156) q[3];
sx q[3];
rz(0.77669981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9876124) q[2];
sx q[2];
rz(-0.93903792) q[2];
sx q[2];
rz(-3.0737446) q[2];
rz(0.96505729) q[3];
sx q[3];
rz(-0.32876757) q[3];
sx q[3];
rz(-1.4062101) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52473749) q[0];
sx q[0];
rz(-2.6829834) q[0];
sx q[0];
rz(0.99880698) q[0];
rz(0.85159341) q[1];
sx q[1];
rz(-1.0706182) q[1];
sx q[1];
rz(0.65382438) q[1];
rz(-1.5077672) q[2];
sx q[2];
rz(-1.0152643) q[2];
sx q[2];
rz(2.0282657) q[2];
rz(1.5953993) q[3];
sx q[3];
rz(-2.2886124) q[3];
sx q[3];
rz(-0.99289865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
