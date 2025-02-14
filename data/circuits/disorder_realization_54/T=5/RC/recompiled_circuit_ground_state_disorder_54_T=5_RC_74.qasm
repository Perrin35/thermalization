OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7689826) q[0];
sx q[0];
rz(3.1867653) q[0];
sx q[0];
rz(10.09633) q[0];
rz(2.1454732) q[1];
sx q[1];
rz(5.4944333) q[1];
sx q[1];
rz(6.5715437) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8771633) q[0];
sx q[0];
rz(-2.9403939) q[0];
sx q[0];
rz(1.7539133) q[0];
rz(-2.4774083) q[2];
sx q[2];
rz(-1.786288) q[2];
sx q[2];
rz(-1.8813949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25569281) q[1];
sx q[1];
rz(-1.5330557) q[1];
sx q[1];
rz(-2.6934212) q[1];
rz(-pi) q[2];
rz(2.8513808) q[3];
sx q[3];
rz(-2.6701616) q[3];
sx q[3];
rz(-0.71347845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4396189) q[2];
sx q[2];
rz(-2.7957323) q[2];
sx q[2];
rz(-0.71887476) q[2];
rz(-1.6950722) q[3];
sx q[3];
rz(-1.4868163) q[3];
sx q[3];
rz(2.1046624) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.059939) q[0];
sx q[0];
rz(-2.7103598) q[0];
sx q[0];
rz(0.28847873) q[0];
rz(2.5892995) q[1];
sx q[1];
rz(-1.0918795) q[1];
sx q[1];
rz(-1.9658032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78647731) q[0];
sx q[0];
rz(-1.8821431) q[0];
sx q[0];
rz(-0.15736736) q[0];
rz(-1.7442877) q[2];
sx q[2];
rz(-1.6855006) q[2];
sx q[2];
rz(-0.8300654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.50880269) q[1];
sx q[1];
rz(-2.1450217) q[1];
sx q[1];
rz(-2.9935163) q[1];
x q[2];
rz(-0.001143841) q[3];
sx q[3];
rz(-1.5387156) q[3];
sx q[3];
rz(-0.45878057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4484078) q[2];
sx q[2];
rz(-1.8790481) q[2];
sx q[2];
rz(0.022424879) q[2];
rz(1.4261931) q[3];
sx q[3];
rz(-1.8636999) q[3];
sx q[3];
rz(-0.9001596) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65396032) q[0];
sx q[0];
rz(-1.2314236) q[0];
sx q[0];
rz(2.3768429) q[0];
rz(-1.7817616) q[1];
sx q[1];
rz(-1.9197074) q[1];
sx q[1];
rz(-1.2545895) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0645197) q[0];
sx q[0];
rz(-1.4409587) q[0];
sx q[0];
rz(-0.10420756) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8171982) q[2];
sx q[2];
rz(-1.430871) q[2];
sx q[2];
rz(-0.10332271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1327127) q[1];
sx q[1];
rz(-0.55435743) q[1];
sx q[1];
rz(2.5704513) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34864254) q[3];
sx q[3];
rz(-2.7789634) q[3];
sx q[3];
rz(1.2606106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4059056) q[2];
sx q[2];
rz(-0.76852208) q[2];
sx q[2];
rz(0.020966919) q[2];
rz(1.5603125) q[3];
sx q[3];
rz(-2.0798648) q[3];
sx q[3];
rz(-0.90644065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2727994) q[0];
sx q[0];
rz(-0.11196207) q[0];
sx q[0];
rz(1.772076) q[0];
rz(-0.70961332) q[1];
sx q[1];
rz(-1.2779002) q[1];
sx q[1];
rz(0.69724625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0727834) q[0];
sx q[0];
rz(-1.8801196) q[0];
sx q[0];
rz(2.7129422) q[0];
rz(2.4575649) q[2];
sx q[2];
rz(-0.942217) q[2];
sx q[2];
rz(-2.7492439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5768421) q[1];
sx q[1];
rz(-2.4240342) q[1];
sx q[1];
rz(-2.0400042) q[1];
rz(-pi) q[2];
rz(1.3342821) q[3];
sx q[3];
rz(-1.7380889) q[3];
sx q[3];
rz(0.084137045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9559481) q[2];
sx q[2];
rz(-0.94719013) q[2];
sx q[2];
rz(2.4626125) q[2];
rz(-1.125157) q[3];
sx q[3];
rz(-1.9966639) q[3];
sx q[3];
rz(-0.2230491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090002447) q[0];
sx q[0];
rz(-1.2035878) q[0];
sx q[0];
rz(-3.0644655) q[0];
rz(1.1513101) q[1];
sx q[1];
rz(-1.5733893) q[1];
sx q[1];
rz(-0.77879771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1098233) q[0];
sx q[0];
rz(-1.5228049) q[0];
sx q[0];
rz(0.035295156) q[0];
x q[1];
rz(-0.16806099) q[2];
sx q[2];
rz(-1.8521143) q[2];
sx q[2];
rz(-2.9945053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38074782) q[1];
sx q[1];
rz(-1.8151974) q[1];
sx q[1];
rz(1.9158855) q[1];
rz(-pi) q[2];
rz(0.058249931) q[3];
sx q[3];
rz(-0.16213972) q[3];
sx q[3];
rz(-1.3221962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.516958) q[2];
sx q[2];
rz(-1.4411074) q[2];
sx q[2];
rz(1.9514294) q[2];
rz(-0.55365753) q[3];
sx q[3];
rz(-0.62479574) q[3];
sx q[3];
rz(0.73738086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6051642) q[0];
sx q[0];
rz(-1.1506511) q[0];
sx q[0];
rz(2.8705257) q[0];
rz(2.1412663) q[1];
sx q[1];
rz(-1.2157636) q[1];
sx q[1];
rz(2.2727374) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2216744) q[0];
sx q[0];
rz(-2.1561047) q[0];
sx q[0];
rz(-2.1409537) q[0];
rz(-pi) q[1];
rz(-1.4914091) q[2];
sx q[2];
rz(-2.304791) q[2];
sx q[2];
rz(-0.47576093) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9166475) q[1];
sx q[1];
rz(-0.38991726) q[1];
sx q[1];
rz(2.030158) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6713761) q[3];
sx q[3];
rz(-0.84132552) q[3];
sx q[3];
rz(2.7470392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.038593682) q[2];
sx q[2];
rz(-0.30240348) q[2];
sx q[2];
rz(-1.137286) q[2];
rz(2.8310827) q[3];
sx q[3];
rz(-0.91028428) q[3];
sx q[3];
rz(-0.92946068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0248658) q[0];
sx q[0];
rz(-1.7127345) q[0];
sx q[0];
rz(-1.0799991) q[0];
rz(-2.179821) q[1];
sx q[1];
rz(-0.56517833) q[1];
sx q[1];
rz(0.22294179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9586102) q[0];
sx q[0];
rz(-1.5633564) q[0];
sx q[0];
rz(-1.6056662) q[0];
rz(-pi) q[1];
rz(-1.9688897) q[2];
sx q[2];
rz(-0.68813339) q[2];
sx q[2];
rz(-1.11042) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0698439) q[1];
sx q[1];
rz(-0.78673601) q[1];
sx q[1];
rz(1.1757502) q[1];
rz(-pi) q[2];
rz(-0.97496521) q[3];
sx q[3];
rz(-1.5585853) q[3];
sx q[3];
rz(1.7705767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.96233931) q[2];
sx q[2];
rz(-0.27955678) q[2];
sx q[2];
rz(-1.1780098) q[2];
rz(-2.8152605) q[3];
sx q[3];
rz(-0.81063619) q[3];
sx q[3];
rz(-1.1403722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2278263) q[0];
sx q[0];
rz(-0.95781177) q[0];
sx q[0];
rz(-0.80818278) q[0];
rz(-0.35762865) q[1];
sx q[1];
rz(-1.6817776) q[1];
sx q[1];
rz(2.0590032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9900695) q[0];
sx q[0];
rz(-2.4381579) q[0];
sx q[0];
rz(1.3057054) q[0];
rz(-1.0183187) q[2];
sx q[2];
rz(-1.3541156) q[2];
sx q[2];
rz(1.073217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4112473) q[1];
sx q[1];
rz(-2.3180006) q[1];
sx q[1];
rz(-0.27539416) q[1];
rz(-1.0937505) q[3];
sx q[3];
rz(-1.9821321) q[3];
sx q[3];
rz(2.3871832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.064934405) q[2];
sx q[2];
rz(-2.6532463) q[2];
sx q[2];
rz(1.2154382) q[2];
rz(2.2293034) q[3];
sx q[3];
rz(-2.2513697) q[3];
sx q[3];
rz(2.598855) q[3];
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
rz(1.2927581) q[0];
sx q[0];
rz(-2.500535) q[0];
sx q[0];
rz(-0.42375281) q[0];
rz(1.1774225) q[1];
sx q[1];
rz(-2.2260428) q[1];
sx q[1];
rz(0.37568572) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59058023) q[0];
sx q[0];
rz(-0.19119054) q[0];
sx q[0];
rz(-1.952233) q[0];
x q[1];
rz(1.5440953) q[2];
sx q[2];
rz(-1.3320247) q[2];
sx q[2];
rz(1.8516061) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5392307) q[1];
sx q[1];
rz(-1.5054387) q[1];
sx q[1];
rz(1.0687625) q[1];
rz(-pi) q[2];
rz(-1.2311835) q[3];
sx q[3];
rz(-2.2582128) q[3];
sx q[3];
rz(-0.31654762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11289135) q[2];
sx q[2];
rz(-1.6501004) q[2];
sx q[2];
rz(-0.69226199) q[2];
rz(2.4568457) q[3];
sx q[3];
rz(-2.1919577) q[3];
sx q[3];
rz(3.0013066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5570062) q[0];
sx q[0];
rz(-0.96730119) q[0];
sx q[0];
rz(1.517357) q[0];
rz(-1.5380305) q[1];
sx q[1];
rz(-1.3471194) q[1];
sx q[1];
rz(-1.2204407) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081471215) q[0];
sx q[0];
rz(-0.94917008) q[0];
sx q[0];
rz(0.073304852) q[0];
x q[1];
rz(2.6399355) q[2];
sx q[2];
rz(-1.3634472) q[2];
sx q[2];
rz(-0.32688552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5904078) q[1];
sx q[1];
rz(-2.3687994) q[1];
sx q[1];
rz(0.40829746) q[1];
rz(-pi) q[2];
rz(-2.4046005) q[3];
sx q[3];
rz(-1.7919645) q[3];
sx q[3];
rz(-1.8870518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5662235) q[2];
sx q[2];
rz(-1.0074002) q[2];
sx q[2];
rz(-2.2929906) q[2];
rz(2.1620915) q[3];
sx q[3];
rz(-1.5749911) q[3];
sx q[3];
rz(-3.1041253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10850567) q[0];
sx q[0];
rz(-1.4773693) q[0];
sx q[0];
rz(0.69837511) q[0];
rz(-0.55950821) q[1];
sx q[1];
rz(-0.63897501) q[1];
sx q[1];
rz(-0.41313304) q[1];
rz(2.8200061) q[2];
sx q[2];
rz(-1.814331) q[2];
sx q[2];
rz(-0.41650256) q[2];
rz(-3.0499023) q[3];
sx q[3];
rz(-0.58135645) q[3];
sx q[3];
rz(-2.9945798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
