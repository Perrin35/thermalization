OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.96707764) q[0];
sx q[0];
rz(-1.4580026) q[0];
sx q[0];
rz(-2.8414677) q[0];
rz(1.1974273) q[1];
sx q[1];
rz(4.6680968) q[1];
sx q[1];
rz(11.065281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5434721) q[0];
sx q[0];
rz(-1.5859005) q[0];
sx q[0];
rz(1.5669109) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3387942) q[2];
sx q[2];
rz(-1.8539696) q[2];
sx q[2];
rz(-0.51529037) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7850656) q[1];
sx q[1];
rz(-1.7887113) q[1];
sx q[1];
rz(2.480605) q[1];
rz(0.0048801076) q[3];
sx q[3];
rz(-0.80880755) q[3];
sx q[3];
rz(-0.67833662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40452051) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(-0.70297757) q[2];
rz(0.56973714) q[3];
sx q[3];
rz(-2.2629786) q[3];
sx q[3];
rz(1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0074145929) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(1.9084357) q[0];
rz(0.59457072) q[1];
sx q[1];
rz(-1.0392799) q[1];
sx q[1];
rz(-2.0436683) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5054975) q[0];
sx q[0];
rz(-1.5042802) q[0];
sx q[0];
rz(1.3696704) q[0];
rz(-pi) q[1];
rz(-2.6498943) q[2];
sx q[2];
rz(-1.9145116) q[2];
sx q[2];
rz(0.51174405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8658376) q[1];
sx q[1];
rz(-3.0601353) q[1];
sx q[1];
rz(-1.2553196) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0497562) q[3];
sx q[3];
rz(-0.79376924) q[3];
sx q[3];
rz(-0.61749279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7132831) q[2];
sx q[2];
rz(-1.2998591) q[2];
sx q[2];
rz(1.606288) q[2];
rz(-2.4226268) q[3];
sx q[3];
rz(-0.9459559) q[3];
sx q[3];
rz(0.21397056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035901) q[0];
sx q[0];
rz(-0.8135697) q[0];
sx q[0];
rz(2.549262) q[0];
rz(-1.2447641) q[1];
sx q[1];
rz(-1.1400305) q[1];
sx q[1];
rz(-0.14561428) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62632769) q[0];
sx q[0];
rz(-1.087731) q[0];
sx q[0];
rz(-1.417571) q[0];
rz(-2.2592741) q[2];
sx q[2];
rz(-1.6711418) q[2];
sx q[2];
rz(2.0813297) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5779931) q[1];
sx q[1];
rz(-1.0219928) q[1];
sx q[1];
rz(-0.33434681) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26992195) q[3];
sx q[3];
rz(-1.0974435) q[3];
sx q[3];
rz(-0.68655754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40272063) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(1.2325475) q[2];
rz(2.9018719) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(-2.6565552) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0204912) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(-3.1109911) q[0];
rz(-0.55514151) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(-2.0487002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3045438) q[0];
sx q[0];
rz(-0.92307675) q[0];
sx q[0];
rz(2.442028) q[0];
rz(-3.080064) q[2];
sx q[2];
rz(-0.76632753) q[2];
sx q[2];
rz(-1.7249223) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3446185) q[1];
sx q[1];
rz(-0.94616468) q[1];
sx q[1];
rz(2.888604) q[1];
rz(0.67375848) q[3];
sx q[3];
rz(-2.1762037) q[3];
sx q[3];
rz(-0.76402367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0577724) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(0.27285451) q[2];
rz(-1.7504292) q[3];
sx q[3];
rz(-2.302156) q[3];
sx q[3];
rz(-0.951989) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6735753) q[0];
sx q[0];
rz(-1.3382358) q[0];
sx q[0];
rz(-1.2982298) q[0];
rz(-0.6908373) q[1];
sx q[1];
rz(-1.9419443) q[1];
sx q[1];
rz(-1.615049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1197136) q[0];
sx q[0];
rz(-2.102431) q[0];
sx q[0];
rz(2.3945532) q[0];
rz(-0.66730209) q[2];
sx q[2];
rz(-0.4385674) q[2];
sx q[2];
rz(-1.9090609) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6479128) q[1];
sx q[1];
rz(-0.70131749) q[1];
sx q[1];
rz(1.5307242) q[1];
rz(2.4590309) q[3];
sx q[3];
rz(-0.45387156) q[3];
sx q[3];
rz(-1.0938494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4096628) q[2];
sx q[2];
rz(-2.3356428) q[2];
sx q[2];
rz(0.16732495) q[2];
rz(1.3903728) q[3];
sx q[3];
rz(-0.53283397) q[3];
sx q[3];
rz(0.65756857) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.14199) q[0];
sx q[0];
rz(-2.7775192) q[0];
sx q[0];
rz(-1.863119) q[0];
rz(1.4356042) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(-2.7659168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784793) q[0];
sx q[0];
rz(-2.7857384) q[0];
sx q[0];
rz(-2.2694017) q[0];
rz(2.2907545) q[2];
sx q[2];
rz(-0.22432835) q[2];
sx q[2];
rz(1.4590603) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3842114) q[1];
sx q[1];
rz(-0.39296341) q[1];
sx q[1];
rz(-2.1860366) q[1];
rz(-pi) q[2];
rz(-2.385684) q[3];
sx q[3];
rz(-1.7497853) q[3];
sx q[3];
rz(-2.8132015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0756691) q[2];
sx q[2];
rz(-1.064294) q[2];
sx q[2];
rz(-2.2806878) q[2];
rz(1.9979477) q[3];
sx q[3];
rz(-2.836561) q[3];
sx q[3];
rz(2.8633269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1160195) q[0];
sx q[0];
rz(-1.3207859) q[0];
sx q[0];
rz(-2.386911) q[0];
rz(-0.93983752) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(-1.8584724) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5425576) q[0];
sx q[0];
rz(-1.5640266) q[0];
sx q[0];
rz(-1.8775131) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2290984) q[2];
sx q[2];
rz(-2.5844838) q[2];
sx q[2];
rz(1.0583641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2312113) q[1];
sx q[1];
rz(-1.2802135) q[1];
sx q[1];
rz(2.742275) q[1];
rz(-pi) q[2];
rz(0.96638443) q[3];
sx q[3];
rz(-2.3570286) q[3];
sx q[3];
rz(1.094629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54520404) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(1.0835353) q[2];
rz(-2.3986473) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(-1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525986) q[0];
sx q[0];
rz(-2.3567663) q[0];
sx q[0];
rz(-2.3866744) q[0];
rz(-0.49916357) q[1];
sx q[1];
rz(-0.9639591) q[1];
sx q[1];
rz(-2.4430433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23125739) q[0];
sx q[0];
rz(-0.6930348) q[0];
sx q[0];
rz(0.49654754) q[0];
rz(-pi) q[1];
rz(2.8590747) q[2];
sx q[2];
rz(-0.48136863) q[2];
sx q[2];
rz(3.0309739) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4515775) q[1];
sx q[1];
rz(-1.5507586) q[1];
sx q[1];
rz(2.103014) q[1];
x q[2];
rz(1.5683062) q[3];
sx q[3];
rz(-2.7225113) q[3];
sx q[3];
rz(-0.97657874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90427202) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(-2.3392056) q[2];
rz(-0.55195105) q[3];
sx q[3];
rz(-0.80058432) q[3];
sx q[3];
rz(2.5210181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17860831) q[0];
sx q[0];
rz(-1.4405788) q[0];
sx q[0];
rz(-2.0340023) q[0];
rz(-1.7604609) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(1.6345056) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12023036) q[0];
sx q[0];
rz(-1.430738) q[0];
sx q[0];
rz(-2.7441478) q[0];
rz(-pi) q[1];
rz(1.1610031) q[2];
sx q[2];
rz(-1.2772182) q[2];
sx q[2];
rz(-2.8382206) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23907614) q[1];
sx q[1];
rz(-2.7488144) q[1];
sx q[1];
rz(-2.7536489) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3146995) q[3];
sx q[3];
rz(-1.6312459) q[3];
sx q[3];
rz(-1.353144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0570809) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(-1.8692807) q[2];
rz(-1.48014) q[3];
sx q[3];
rz(-1.1353761) q[3];
sx q[3];
rz(-2.541466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965974) q[0];
sx q[0];
rz(-2.5801881) q[0];
sx q[0];
rz(-1.2835314) q[0];
rz(0.01783477) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(-3.0160115) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3572306) q[0];
sx q[0];
rz(-1.0166369) q[0];
sx q[0];
rz(3.1249376) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0415885) q[2];
sx q[2];
rz(-0.94816899) q[2];
sx q[2];
rz(2.560844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8931742) q[1];
sx q[1];
rz(-1.9918973) q[1];
sx q[1];
rz(2.0401272) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9405648) q[3];
sx q[3];
rz(-1.8327692) q[3];
sx q[3];
rz(0.51823814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77832001) q[2];
sx q[2];
rz(-1.8724226) q[2];
sx q[2];
rz(-0.8030836) q[2];
rz(-2.965029) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(0.77809063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9618027) q[0];
sx q[0];
rz(-0.50146865) q[0];
sx q[0];
rz(1.5487221) q[0];
rz(-1.8190307) q[1];
sx q[1];
rz(-1.3643199) q[1];
sx q[1];
rz(1.551569) q[1];
rz(-1.5079458) q[2];
sx q[2];
rz(-1.3286776) q[2];
sx q[2];
rz(1.5059289) q[2];
rz(-1.5104891) q[3];
sx q[3];
rz(-1.9077488) q[3];
sx q[3];
rz(0.13691402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
