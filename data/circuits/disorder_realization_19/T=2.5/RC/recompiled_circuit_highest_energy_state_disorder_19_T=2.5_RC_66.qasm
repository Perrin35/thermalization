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
rz(2.1228696) q[0];
sx q[0];
rz(-2.2824204) q[0];
sx q[0];
rz(0.81508842) q[0];
rz(-1.1852784) q[1];
sx q[1];
rz(-1.4108682) q[1];
sx q[1];
rz(1.0676395) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51781228) q[0];
sx q[0];
rz(-1.5470439) q[0];
sx q[0];
rz(-1.1393113) q[0];
x q[1];
rz(2.5065866) q[2];
sx q[2];
rz(-2.5189812) q[2];
sx q[2];
rz(-0.16281637) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63032615) q[1];
sx q[1];
rz(-2.8095803) q[1];
sx q[1];
rz(1.270833) q[1];
rz(-0.26431636) q[3];
sx q[3];
rz(-2.2393763) q[3];
sx q[3];
rz(0.17449915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4472569) q[2];
sx q[2];
rz(-2.3981636) q[2];
sx q[2];
rz(1.2747964) q[2];
rz(-2.6796807) q[3];
sx q[3];
rz(-0.67449823) q[3];
sx q[3];
rz(-2.0089202) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79134113) q[0];
sx q[0];
rz(-0.30838648) q[0];
sx q[0];
rz(-1.2530918) q[0];
rz(2.99627) q[1];
sx q[1];
rz(-1.7456313) q[1];
sx q[1];
rz(1.0911509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017830124) q[0];
sx q[0];
rz(-1.995867) q[0];
sx q[0];
rz(-0.2353038) q[0];
rz(-pi) q[1];
rz(0.65204377) q[2];
sx q[2];
rz(-2.5863159) q[2];
sx q[2];
rz(2.9342143) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5877643) q[1];
sx q[1];
rz(-1.3979646) q[1];
sx q[1];
rz(0.21765222) q[1];
rz(-pi) q[2];
rz(-1.9698148) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(-2.6103013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48065177) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(-2.6118028) q[2];
rz(-2.3482813) q[3];
sx q[3];
rz(-1.5373693) q[3];
sx q[3];
rz(2.5180499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48524258) q[0];
sx q[0];
rz(-0.95004496) q[0];
sx q[0];
rz(2.6128838) q[0];
rz(-0.54620019) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(2.8012457) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0673235) q[0];
sx q[0];
rz(-2.7087697) q[0];
sx q[0];
rz(-0.69822635) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5916948) q[2];
sx q[2];
rz(-0.27355121) q[2];
sx q[2];
rz(-3.105148) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3888549) q[1];
sx q[1];
rz(-1.5675263) q[1];
sx q[1];
rz(-3.0470303) q[1];
rz(-1.9779376) q[3];
sx q[3];
rz(-0.7837067) q[3];
sx q[3];
rz(-0.78898417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9693552) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(-2.7117512) q[2];
rz(-2.3853081) q[3];
sx q[3];
rz(-2.9679306) q[3];
sx q[3];
rz(-0.82591301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38465685) q[0];
sx q[0];
rz(-1.6357559) q[0];
sx q[0];
rz(1.6480308) q[0];
rz(-0.76796302) q[1];
sx q[1];
rz(-0.47052828) q[1];
sx q[1];
rz(0.54642645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38619216) q[0];
sx q[0];
rz(-1.6419171) q[0];
sx q[0];
rz(3.0767308) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52164061) q[2];
sx q[2];
rz(-2.0211136) q[2];
sx q[2];
rz(-0.7618103) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5810952) q[1];
sx q[1];
rz(-0.17737922) q[1];
sx q[1];
rz(-2.5590798) q[1];
x q[2];
rz(-2.4599471) q[3];
sx q[3];
rz(-2.1053616) q[3];
sx q[3];
rz(-0.93972423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4297318) q[2];
sx q[2];
rz(-1.4873361) q[2];
sx q[2];
rz(-2.8374953) q[2];
rz(1.3652623) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(2.064866) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51233184) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(-3.1410134) q[0];
rz(-2.0594788) q[1];
sx q[1];
rz(-0.52434701) q[1];
sx q[1];
rz(0.93926114) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1056553) q[0];
sx q[0];
rz(-2.2024184) q[0];
sx q[0];
rz(-1.0211591) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1823857) q[2];
sx q[2];
rz(-1.8610753) q[2];
sx q[2];
rz(-0.91735754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56783695) q[1];
sx q[1];
rz(-1.7845961) q[1];
sx q[1];
rz(-0.9937728) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22721283) q[3];
sx q[3];
rz(-1.876653) q[3];
sx q[3];
rz(1.7744886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.038736343) q[2];
sx q[2];
rz(-1.880371) q[2];
sx q[2];
rz(2.8803414) q[2];
rz(-2.2767565) q[3];
sx q[3];
rz(-1.4120925) q[3];
sx q[3];
rz(-0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84457266) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(-2.936506) q[0];
rz(-1.0467485) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(2.7632025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9721824) q[0];
sx q[0];
rz(-0.43146389) q[0];
sx q[0];
rz(-1.1110825) q[0];
rz(2.66983) q[2];
sx q[2];
rz(-1.7182351) q[2];
sx q[2];
rz(-0.23213895) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7415315) q[1];
sx q[1];
rz(-2.0963256) q[1];
sx q[1];
rz(-2.7353213) q[1];
rz(1.1388402) q[3];
sx q[3];
rz(-2.216385) q[3];
sx q[3];
rz(2.1172252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4947074) q[2];
sx q[2];
rz(-1.1581706) q[2];
sx q[2];
rz(2.0873439) q[2];
rz(2.4833637) q[3];
sx q[3];
rz(-0.64526486) q[3];
sx q[3];
rz(2.6575507) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14195104) q[0];
sx q[0];
rz(-1.9331837) q[0];
sx q[0];
rz(0.64055881) q[0];
rz(2.7236252) q[1];
sx q[1];
rz(-1.9211946) q[1];
sx q[1];
rz(-2.3366065) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5403596) q[0];
sx q[0];
rz(-0.86358294) q[0];
sx q[0];
rz(-1.5401031) q[0];
rz(-pi) q[1];
rz(-1.1518794) q[2];
sx q[2];
rz(-2.3275314) q[2];
sx q[2];
rz(-0.69413041) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0143483) q[1];
sx q[1];
rz(-1.5509203) q[1];
sx q[1];
rz(1.5991421) q[1];
rz(0.42607362) q[3];
sx q[3];
rz(-0.63562993) q[3];
sx q[3];
rz(-2.9866708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54887041) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(2.3804046) q[2];
rz(-2.393764) q[3];
sx q[3];
rz(-1.3176094) q[3];
sx q[3];
rz(-0.67659155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51650301) q[0];
sx q[0];
rz(-2.857132) q[0];
sx q[0];
rz(-1.6647343) q[0];
rz(-0.45627108) q[1];
sx q[1];
rz(-1.7218593) q[1];
sx q[1];
rz(-2.2705073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7719771) q[0];
sx q[0];
rz(-0.61423683) q[0];
sx q[0];
rz(-0.02158879) q[0];
rz(-pi) q[1];
rz(-0.41042491) q[2];
sx q[2];
rz(-2.7685809) q[2];
sx q[2];
rz(-0.57508627) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5567315) q[1];
sx q[1];
rz(-1.8980527) q[1];
sx q[1];
rz(1.8021939) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7026448) q[3];
sx q[3];
rz(-2.7119497) q[3];
sx q[3];
rz(2.1983918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2386834) q[2];
sx q[2];
rz(-2.6148655) q[2];
sx q[2];
rz(-2.5229559) q[2];
rz(3.1213308) q[3];
sx q[3];
rz(-2.198115) q[3];
sx q[3];
rz(-0.77997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063754931) q[0];
sx q[0];
rz(-1.5564593) q[0];
sx q[0];
rz(0.42386398) q[0];
rz(0.18691143) q[1];
sx q[1];
rz(-2.321545) q[1];
sx q[1];
rz(-1.4580457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20987147) q[0];
sx q[0];
rz(-0.11882028) q[0];
sx q[0];
rz(-0.11943905) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0571613) q[2];
sx q[2];
rz(-1.3328526) q[2];
sx q[2];
rz(-2.6702211) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8193389) q[1];
sx q[1];
rz(-1.6221433) q[1];
sx q[1];
rz(-0.13887413) q[1];
rz(-pi) q[2];
rz(0.23641674) q[3];
sx q[3];
rz(-0.36521736) q[3];
sx q[3];
rz(-1.055298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3905048) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(-1.0941774) q[2];
rz(-1.68082) q[3];
sx q[3];
rz(-1.3760309) q[3];
sx q[3];
rz(-1.338753) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3222892) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(-1.639701) q[0];
rz(-0.27885258) q[1];
sx q[1];
rz(-0.96062213) q[1];
sx q[1];
rz(1.0241114) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8968643) q[0];
sx q[0];
rz(-2.6435268) q[0];
sx q[0];
rz(-2.9480272) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0908998) q[2];
sx q[2];
rz(-1.2963352) q[2];
sx q[2];
rz(-1.7982184) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0054202) q[1];
sx q[1];
rz(-1.7917624) q[1];
sx q[1];
rz(1.7762089) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9387705) q[3];
sx q[3];
rz(-2.362613) q[3];
sx q[3];
rz(2.413903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9378822) q[2];
sx q[2];
rz(-1.2869765) q[2];
sx q[2];
rz(2.6046275) q[2];
rz(1.2282061) q[3];
sx q[3];
rz(-2.1824586) q[3];
sx q[3];
rz(-0.29135191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0252329) q[0];
sx q[0];
rz(-1.3675084) q[0];
sx q[0];
rz(1.0687923) q[0];
rz(0.7069201) q[1];
sx q[1];
rz(-2.0126577) q[1];
sx q[1];
rz(-0.001002034) q[1];
rz(0.42264414) q[2];
sx q[2];
rz(-2.7012237) q[2];
sx q[2];
rz(1.2958432) q[2];
rz(-2.1643987) q[3];
sx q[3];
rz(-2.8799812) q[3];
sx q[3];
rz(2.9577586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
