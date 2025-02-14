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
rz(0.81646252) q[0];
sx q[0];
rz(3.2433885) q[0];
sx q[0];
rz(9.9643702) q[0];
rz(-2.645283) q[1];
sx q[1];
rz(-2.8318383) q[1];
sx q[1];
rz(2.6024979) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91796275) q[0];
sx q[0];
rz(-1.5453891) q[0];
sx q[0];
rz(-0.44141234) q[0];
rz(-pi) q[1];
rz(-0.38549785) q[2];
sx q[2];
rz(-0.38542569) q[2];
sx q[2];
rz(-0.62776923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4161168) q[1];
sx q[1];
rz(-1.892964) q[1];
sx q[1];
rz(-2.9941818) q[1];
x q[2];
rz(1.4629389) q[3];
sx q[3];
rz(-2.4898006) q[3];
sx q[3];
rz(-2.6848313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76217905) q[2];
sx q[2];
rz(-0.87612408) q[2];
sx q[2];
rz(-1.1890821) q[2];
rz(1.9573697) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(-1.5949465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.14520833) q[0];
sx q[0];
rz(-1.5897911) q[0];
sx q[0];
rz(1.8183964) q[0];
rz(0.48201758) q[1];
sx q[1];
rz(-0.91164416) q[1];
sx q[1];
rz(-0.97420305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27069399) q[0];
sx q[0];
rz(-0.93241954) q[0];
sx q[0];
rz(2.9186072) q[0];
rz(0.09928273) q[2];
sx q[2];
rz(-1.6834604) q[2];
sx q[2];
rz(-2.7275865) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1434553) q[1];
sx q[1];
rz(-1.0364729) q[1];
sx q[1];
rz(-1.8443395) q[1];
rz(1.810435) q[3];
sx q[3];
rz(-1.0175034) q[3];
sx q[3];
rz(3.1146793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1366068) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(-2.3919487) q[2];
rz(0.43198112) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(1.8384793) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5169446) q[0];
sx q[0];
rz(-2.8693146) q[0];
sx q[0];
rz(2.5352449) q[0];
rz(1.3336522) q[1];
sx q[1];
rz(-1.2140112) q[1];
sx q[1];
rz(2.8820754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57588803) q[0];
sx q[0];
rz(-1.7262926) q[0];
sx q[0];
rz(1.8207654) q[0];
rz(-0.86028966) q[2];
sx q[2];
rz(-2.8767071) q[2];
sx q[2];
rz(0.67240326) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1705679) q[1];
sx q[1];
rz(-1.1491421) q[1];
sx q[1];
rz(0.41280156) q[1];
rz(-pi) q[2];
rz(0.43478888) q[3];
sx q[3];
rz(-2.090095) q[3];
sx q[3];
rz(-3.0360707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8848662) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(1.5941031) q[2];
rz(0.78440845) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(-1.5756395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.49482685) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(2.8821017) q[0];
rz(-1.9208113) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(1.4422013) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3078559) q[0];
sx q[0];
rz(-2.4015732) q[0];
sx q[0];
rz(-2.3038175) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7624315) q[2];
sx q[2];
rz(-0.60376061) q[2];
sx q[2];
rz(1.9866643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4220548) q[1];
sx q[1];
rz(-1.2201078) q[1];
sx q[1];
rz(-1.3077875) q[1];
x q[2];
rz(-0.070017858) q[3];
sx q[3];
rz(-2.2517831) q[3];
sx q[3];
rz(-2.3511166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1059025) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(-0.40994677) q[2];
rz(1.2116872) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(1.80779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164417) q[0];
sx q[0];
rz(-2.4254159) q[0];
sx q[0];
rz(-0.27405611) q[0];
rz(-2.7033499) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(2.4058707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0817524) q[0];
sx q[0];
rz(-2.3222089) q[0];
sx q[0];
rz(-1.451534) q[0];
rz(-0.76916285) q[2];
sx q[2];
rz(-0.20649466) q[2];
sx q[2];
rz(1.002305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2949038) q[1];
sx q[1];
rz(-1.8195933) q[1];
sx q[1];
rz(1.873068) q[1];
rz(0.93345668) q[3];
sx q[3];
rz(-0.75691716) q[3];
sx q[3];
rz(0.33801038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.938544) q[2];
sx q[2];
rz(-1.6837348) q[2];
sx q[2];
rz(-0.61384821) q[2];
rz(0.40890536) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(-0.74234211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78855377) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(1.4452112) q[1];
sx q[1];
rz(-0.45672363) q[1];
sx q[1];
rz(-2.9663185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69578457) q[0];
sx q[0];
rz(-2.9668167) q[0];
sx q[0];
rz(-1.4265027) q[0];
rz(-2.1138328) q[2];
sx q[2];
rz(-0.98688417) q[2];
sx q[2];
rz(-1.2500545) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1239667) q[1];
sx q[1];
rz(-1.062289) q[1];
sx q[1];
rz(1.0715436) q[1];
x q[2];
rz(-2.184559) q[3];
sx q[3];
rz(-1.0816469) q[3];
sx q[3];
rz(-0.027994284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.96482977) q[2];
sx q[2];
rz(-1.333933) q[2];
sx q[2];
rz(0.97243398) q[2];
rz(-1.4771627) q[3];
sx q[3];
rz(-1.1494613) q[3];
sx q[3];
rz(-2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2870188) q[0];
sx q[0];
rz(-3.119097) q[0];
sx q[0];
rz(0.35312411) q[0];
rz(-1.0278206) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(-1.0110528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7142657) q[0];
sx q[0];
rz(-1.9280701) q[0];
sx q[0];
rz(1.5861804) q[0];
rz(-pi) q[1];
rz(-0.6017466) q[2];
sx q[2];
rz(-0.95843747) q[2];
sx q[2];
rz(-1.6342721) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2680929) q[1];
sx q[1];
rz(-2.756739) q[1];
sx q[1];
rz(1.1420239) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2561463) q[3];
sx q[3];
rz(-1.8487747) q[3];
sx q[3];
rz(-1.4643471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8280243) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(-1.3671406) q[2];
rz(-1.8732871) q[3];
sx q[3];
rz(-1.4788078) q[3];
sx q[3];
rz(0.84793276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0276412) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(2.0116346) q[0];
rz(-2.9226411) q[1];
sx q[1];
rz(-1.4066701) q[1];
sx q[1];
rz(2.2311282) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51800358) q[0];
sx q[0];
rz(-1.2388889) q[0];
sx q[0];
rz(2.5123572) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47825079) q[2];
sx q[2];
rz(-1.4946117) q[2];
sx q[2];
rz(-2.8752799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1840802) q[1];
sx q[1];
rz(-1.0614479) q[1];
sx q[1];
rz(-2.8345076) q[1];
rz(2.4164356) q[3];
sx q[3];
rz(-1.8817543) q[3];
sx q[3];
rz(-0.5288611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3160481) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(-2.3160589) q[2];
rz(2.6801706) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(-0.44574827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(0.6398741) q[0];
sx q[0];
rz(-1.8032782) q[0];
sx q[0];
rz(-0.83183944) q[0];
rz(-2.4141451) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(0.096253455) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98672685) q[0];
sx q[0];
rz(-0.71457246) q[0];
sx q[0];
rz(-2.3113219) q[0];
rz(-pi) q[1];
rz(0.23230884) q[2];
sx q[2];
rz(-0.26089982) q[2];
sx q[2];
rz(1.0123073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5288594) q[1];
sx q[1];
rz(-1.7816356) q[1];
sx q[1];
rz(0.092540578) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89266606) q[3];
sx q[3];
rz(-2.8916807) q[3];
sx q[3];
rz(-1.753861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40621296) q[2];
sx q[2];
rz(-1.4344183) q[2];
sx q[2];
rz(-1.5926788) q[2];
rz(-0.63273543) q[3];
sx q[3];
rz(-1.2318719) q[3];
sx q[3];
rz(0.24937853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525472) q[0];
sx q[0];
rz(-2.1010375) q[0];
sx q[0];
rz(-0.35300514) q[0];
rz(2.8189335) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(2.7696612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73814017) q[0];
sx q[0];
rz(-0.36309013) q[0];
sx q[0];
rz(0.9842666) q[0];
rz(-0.48666059) q[2];
sx q[2];
rz(-2.7722205) q[2];
sx q[2];
rz(-0.83115679) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5185753) q[1];
sx q[1];
rz(-1.9037582) q[1];
sx q[1];
rz(1.4013616) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2953565) q[3];
sx q[3];
rz(-1.5414943) q[3];
sx q[3];
rz(-1.3361564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9670664) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(-1.2971499) q[2];
rz(0.8477115) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(-2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2321155) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(-1.3879981) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(0.83609381) q[2];
sx q[2];
rz(-2.2420364) q[2];
sx q[2];
rz(-1.4902761) q[2];
rz(1.150849) q[3];
sx q[3];
rz(-2.3227878) q[3];
sx q[3];
rz(-0.25521758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
