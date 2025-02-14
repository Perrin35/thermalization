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
rz(1.6871356) q[0];
sx q[0];
rz(1.9782826) q[0];
sx q[0];
rz(10.692698) q[0];
rz(1.6549702) q[1];
sx q[1];
rz(4.4284664) q[1];
sx q[1];
rz(9.5944302) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21714144) q[0];
sx q[0];
rz(-2.0471153) q[0];
sx q[0];
rz(0.40556559) q[0];
rz(-2.0889766) q[2];
sx q[2];
rz(-1.7917601) q[2];
sx q[2];
rz(0.25928869) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4028266) q[1];
sx q[1];
rz(-1.9662274) q[1];
sx q[1];
rz(0.02301245) q[1];
rz(-2.9291487) q[3];
sx q[3];
rz(-1.2646958) q[3];
sx q[3];
rz(2.3116632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22592813) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(0.020326745) q[2];
rz(-2.9018371) q[3];
sx q[3];
rz(-2.8862947) q[3];
sx q[3];
rz(2.3026626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13424419) q[0];
sx q[0];
rz(-2.5237995) q[0];
sx q[0];
rz(2.0452621) q[0];
rz(-1.4758551) q[1];
sx q[1];
rz(-1.9225537) q[1];
sx q[1];
rz(-0.79442564) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83704575) q[0];
sx q[0];
rz(-1.908806) q[0];
sx q[0];
rz(1.2370626) q[0];
rz(-0.99296221) q[2];
sx q[2];
rz(-0.90415274) q[2];
sx q[2];
rz(-1.7844328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.55358821) q[1];
sx q[1];
rz(-0.4376227) q[1];
sx q[1];
rz(2.8904084) q[1];
rz(2.9555213) q[3];
sx q[3];
rz(-2.3181462) q[3];
sx q[3];
rz(-0.88445437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.085792556) q[2];
sx q[2];
rz(-1.5689359) q[2];
sx q[2];
rz(2.9158578) q[2];
rz(2.8977532) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(-2.9624511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2589515) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(2.7146345) q[0];
rz(-0.34307617) q[1];
sx q[1];
rz(-1.1767358) q[1];
sx q[1];
rz(-0.25161904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7386603) q[0];
sx q[0];
rz(-2.1485747) q[0];
sx q[0];
rz(-1.8389804) q[0];
rz(-1.8015091) q[2];
sx q[2];
rz(-0.78015155) q[2];
sx q[2];
rz(-3.0456269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3494663) q[1];
sx q[1];
rz(-1.2361965) q[1];
sx q[1];
rz(-1.7484596) q[1];
rz(-pi) q[2];
rz(0.050046845) q[3];
sx q[3];
rz(-2.280683) q[3];
sx q[3];
rz(2.0662226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.130927) q[2];
sx q[2];
rz(-1.6341354) q[2];
sx q[2];
rz(-0.44535401) q[2];
rz(1.2238067) q[3];
sx q[3];
rz(-1.2659975) q[3];
sx q[3];
rz(-0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7447516) q[0];
sx q[0];
rz(-2.3362609) q[0];
sx q[0];
rz(-1.9824363) q[0];
rz(-1.9388439) q[1];
sx q[1];
rz(-1.1801327) q[1];
sx q[1];
rz(2.4045827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7173968) q[0];
sx q[0];
rz(-1.5079466) q[0];
sx q[0];
rz(-2.8934188) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1221755) q[2];
sx q[2];
rz(-0.72690847) q[2];
sx q[2];
rz(-0.87269937) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9056978) q[1];
sx q[1];
rz(-1.3970084) q[1];
sx q[1];
rz(-2.5008965) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6287732) q[3];
sx q[3];
rz(-1.7372469) q[3];
sx q[3];
rz(-0.49515192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0143339) q[2];
sx q[2];
rz(-2.2517683) q[2];
sx q[2];
rz(-2.3298402) q[2];
rz(2.3938866) q[3];
sx q[3];
rz(-0.85993189) q[3];
sx q[3];
rz(-0.060613304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6478445) q[0];
sx q[0];
rz(-2.0525377) q[0];
sx q[0];
rz(-1.7895948) q[0];
rz(-0.79212517) q[1];
sx q[1];
rz(-0.89901662) q[1];
sx q[1];
rz(-2.1314714) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5631752) q[0];
sx q[0];
rz(-2.0978598) q[0];
sx q[0];
rz(-0.38946797) q[0];
x q[1];
rz(0.13212684) q[2];
sx q[2];
rz(-2.7458522) q[2];
sx q[2];
rz(2.6504315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6163075) q[1];
sx q[1];
rz(-2.1263564) q[1];
sx q[1];
rz(1.6643334) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95577464) q[3];
sx q[3];
rz(-2.9412651) q[3];
sx q[3];
rz(2.1708787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.036756214) q[2];
sx q[2];
rz(-0.75883055) q[2];
sx q[2];
rz(1.3834312) q[2];
rz(3.1116327) q[3];
sx q[3];
rz(-2.1432917) q[3];
sx q[3];
rz(1.8647319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6303915) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(-3.0101486) q[0];
rz(2.1427872) q[1];
sx q[1];
rz(-1.1591594) q[1];
sx q[1];
rz(-1.3712032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5509439) q[0];
sx q[0];
rz(-1.9134132) q[0];
sx q[0];
rz(-0.37775535) q[0];
rz(-pi) q[1];
rz(1.7741469) q[2];
sx q[2];
rz(-2.2468745) q[2];
sx q[2];
rz(1.0695367) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7853493) q[1];
sx q[1];
rz(-1.8438854) q[1];
sx q[1];
rz(2.9778984) q[1];
rz(-3.1375299) q[3];
sx q[3];
rz(-2.4306979) q[3];
sx q[3];
rz(-1.9181741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3580631) q[2];
sx q[2];
rz(-2.4130267) q[2];
sx q[2];
rz(-2.0520468) q[2];
rz(-0.36695668) q[3];
sx q[3];
rz(-2.7373382) q[3];
sx q[3];
rz(-2.3914316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51339665) q[0];
sx q[0];
rz(-2.3371526) q[0];
sx q[0];
rz(-0.87271571) q[0];
rz(1.577042) q[1];
sx q[1];
rz(-1.6460452) q[1];
sx q[1];
rz(0.73211342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67938618) q[0];
sx q[0];
rz(-2.4024978) q[0];
sx q[0];
rz(0.30090354) q[0];
rz(-1.6953515) q[2];
sx q[2];
rz(-0.8890748) q[2];
sx q[2];
rz(-2.9602697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70591761) q[1];
sx q[1];
rz(-0.49723782) q[1];
sx q[1];
rz(1.9051208) q[1];
x q[2];
rz(-1.9663345) q[3];
sx q[3];
rz(-2.8962781) q[3];
sx q[3];
rz(-0.79110629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5319891) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(0.58094376) q[2];
rz(-1.1254492) q[3];
sx q[3];
rz(-1.9476798) q[3];
sx q[3];
rz(1.9862991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693212) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(-0.17053764) q[0];
rz(-1.2455617) q[1];
sx q[1];
rz(-1.5252557) q[1];
sx q[1];
rz(-2.4571498) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.754622) q[0];
sx q[0];
rz(-2.5477161) q[0];
sx q[0];
rz(-2.519849) q[0];
rz(-pi) q[1];
rz(-2.6587517) q[2];
sx q[2];
rz(-1.4582658) q[2];
sx q[2];
rz(-1.4275328) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0741006) q[1];
sx q[1];
rz(-1.3252186) q[1];
sx q[1];
rz(-0.26439338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4926458) q[3];
sx q[3];
rz(-1.5228378) q[3];
sx q[3];
rz(-3.0779882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.838852) q[2];
sx q[2];
rz(-1.3417696) q[2];
sx q[2];
rz(-1.1203241) q[2];
rz(2.4773347) q[3];
sx q[3];
rz(-0.40226007) q[3];
sx q[3];
rz(-0.50814381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6054194) q[0];
sx q[0];
rz(-0.37261951) q[0];
sx q[0];
rz(2.5033409) q[0];
rz(0.0087180184) q[1];
sx q[1];
rz(-1.1161085) q[1];
sx q[1];
rz(-2.7353824) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531949) q[0];
sx q[0];
rz(-2.8450845) q[0];
sx q[0];
rz(-1.1738846) q[0];
rz(-pi) q[1];
rz(3.0893095) q[2];
sx q[2];
rz(-2.0039275) q[2];
sx q[2];
rz(-1.0811102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9487755) q[1];
sx q[1];
rz(-1.684976) q[1];
sx q[1];
rz(0.78928708) q[1];
rz(-0.073447179) q[3];
sx q[3];
rz(-1.5961507) q[3];
sx q[3];
rz(3.1318093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0476039) q[2];
sx q[2];
rz(-2.0145907) q[2];
sx q[2];
rz(2.7809533) q[2];
rz(-2.4141198) q[3];
sx q[3];
rz(-0.68220264) q[3];
sx q[3];
rz(0.31807652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.3280846) q[0];
sx q[0];
rz(-2.1961975) q[0];
sx q[0];
rz(-0.16950053) q[0];
rz(2.4318579) q[1];
sx q[1];
rz(-1.7252445) q[1];
sx q[1];
rz(-2.0555326) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5616578) q[0];
sx q[0];
rz(-0.57889639) q[0];
sx q[0];
rz(-2.7282342) q[0];
x q[1];
rz(-2.0073838) q[2];
sx q[2];
rz(-1.2870169) q[2];
sx q[2];
rz(-1.6095609) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5063613) q[1];
sx q[1];
rz(-1.6690134) q[1];
sx q[1];
rz(2.1327095) q[1];
x q[2];
rz(2.707241) q[3];
sx q[3];
rz(-0.41623033) q[3];
sx q[3];
rz(-0.21823635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3544932) q[2];
sx q[2];
rz(-0.084986173) q[2];
sx q[2];
rz(-0.3698012) q[2];
rz(-1.9434816) q[3];
sx q[3];
rz(-1.5503649) q[3];
sx q[3];
rz(-0.91355356) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13854606) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(1.0176324) q[1];
sx q[1];
rz(-1.3980649) q[1];
sx q[1];
rz(-1.3921888) q[1];
rz(-0.77059435) q[2];
sx q[2];
rz(-3.040585) q[2];
sx q[2];
rz(0.2069006) q[2];
rz(-0.16416488) q[3];
sx q[3];
rz(-1.6995242) q[3];
sx q[3];
rz(-1.8229998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
