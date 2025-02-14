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
rz(-2.2588377) q[0];
sx q[0];
rz(-0.14482276) q[0];
sx q[0];
rz(2.3902399) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(5.3137988) q[1];
sx q[1];
rz(9.245524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0796666) q[0];
sx q[0];
rz(-2.5653337) q[0];
sx q[0];
rz(2.1012596) q[0];
rz(-pi) q[1];
rz(0.037161552) q[2];
sx q[2];
rz(-2.4183309) q[2];
sx q[2];
rz(0.0497555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1826659) q[1];
sx q[1];
rz(-1.8799056) q[1];
sx q[1];
rz(-1.4953509) q[1];
rz(-pi) q[2];
rz(3.0823067) q[3];
sx q[3];
rz(-2.2036457) q[3];
sx q[3];
rz(2.2135753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2396607) q[2];
sx q[2];
rz(-1.2942945) q[2];
sx q[2];
rz(-0.61061668) q[2];
rz(-0.88879746) q[3];
sx q[3];
rz(-0.68032467) q[3];
sx q[3];
rz(2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4438542) q[0];
sx q[0];
rz(-0.30174169) q[0];
sx q[0];
rz(-1.8119716) q[0];
rz(-1.4854206) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(-0.86404538) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5764424) q[0];
sx q[0];
rz(-2.8684542) q[0];
sx q[0];
rz(-1.9594203) q[0];
rz(-pi) q[1];
rz(-2.1435166) q[2];
sx q[2];
rz(-0.53600271) q[2];
sx q[2];
rz(-2.0886476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5031918) q[1];
sx q[1];
rz(-1.4309037) q[1];
sx q[1];
rz(-1.7452471) q[1];
rz(-pi) q[2];
rz(2.3584189) q[3];
sx q[3];
rz(-2.0370954) q[3];
sx q[3];
rz(-0.1649905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0835421) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(2.9502499) q[2];
rz(2.8355016) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(0.93650854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8323583) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(2.7171296) q[0];
rz(-1.3407432) q[1];
sx q[1];
rz(-2.2258591) q[1];
sx q[1];
rz(-0.65139604) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49418816) q[0];
sx q[0];
rz(-1.5779102) q[0];
sx q[0];
rz(1.4071541) q[0];
rz(-pi) q[1];
rz(1.3878787) q[2];
sx q[2];
rz(-1.1929034) q[2];
sx q[2];
rz(1.479508) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0733903) q[1];
sx q[1];
rz(-1.7531839) q[1];
sx q[1];
rz(0.50478151) q[1];
rz(2.1742854) q[3];
sx q[3];
rz(-0.27316948) q[3];
sx q[3];
rz(-0.73707132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.035630781) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(0.6967217) q[2];
rz(-0.68909711) q[3];
sx q[3];
rz(-1.9870116) q[3];
sx q[3];
rz(2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7032787) q[0];
sx q[0];
rz(-2.8499481) q[0];
sx q[0];
rz(-2.1412204) q[0];
rz(1.6814303) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(-1.9283074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.723169) q[0];
sx q[0];
rz(-0.30659404) q[0];
sx q[0];
rz(-1.9830389) q[0];
rz(-pi) q[1];
rz(2.6268509) q[2];
sx q[2];
rz(-0.79823433) q[2];
sx q[2];
rz(-1.9412184) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69566947) q[1];
sx q[1];
rz(-1.6918105) q[1];
sx q[1];
rz(-2.8112429) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7217851) q[3];
sx q[3];
rz(-2.5979492) q[3];
sx q[3];
rz(2.061894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0393684) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(-3.0774806) q[2];
rz(0.72470775) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(-2.7418315) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940014) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(-0.79175788) q[0];
rz(2.0296312) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(-1.3349104) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70671088) q[0];
sx q[0];
rz(-2.0405053) q[0];
sx q[0];
rz(2.241894) q[0];
rz(-pi) q[1];
rz(2.9826775) q[2];
sx q[2];
rz(-0.58664413) q[2];
sx q[2];
rz(1.9288837) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68685442) q[1];
sx q[1];
rz(-1.6461186) q[1];
sx q[1];
rz(-1.7416864) q[1];
rz(-pi) q[2];
rz(2.9369257) q[3];
sx q[3];
rz(-2.137326) q[3];
sx q[3];
rz(0.03290225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.84805924) q[2];
sx q[2];
rz(-2.2424825) q[2];
sx q[2];
rz(1.1406356) q[2];
rz(3.0349777) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(-0.30985668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.8150197) q[0];
sx q[0];
rz(-0.44875479) q[0];
sx q[0];
rz(-0.003224592) q[0];
rz(3.0942753) q[1];
sx q[1];
rz(-1.4553757) q[1];
sx q[1];
rz(-1.3501732) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88104507) q[0];
sx q[0];
rz(-0.25940093) q[0];
sx q[0];
rz(1.3400643) q[0];
rz(-pi) q[1];
rz(2.240594) q[2];
sx q[2];
rz(-1.0987079) q[2];
sx q[2];
rz(2.6028518) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9612757) q[1];
sx q[1];
rz(-2.2928228) q[1];
sx q[1];
rz(2.3555593) q[1];
rz(-pi) q[2];
x q[2];
rz(0.030478625) q[3];
sx q[3];
rz(-1.7913831) q[3];
sx q[3];
rz(-2.7762846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1515767) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(-3.0604176) q[2];
rz(-2.2422527) q[3];
sx q[3];
rz(-2.3266561) q[3];
sx q[3];
rz(1.8544633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7154295) q[0];
sx q[0];
rz(-1.8303215) q[0];
sx q[0];
rz(-2.6370866) q[0];
rz(2.1465178) q[1];
sx q[1];
rz(-2.1213396) q[1];
sx q[1];
rz(3.1212433) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58344719) q[0];
sx q[0];
rz(-2.7648395) q[0];
sx q[0];
rz(1.8280297) q[0];
rz(-pi) q[1];
rz(-1.3624914) q[2];
sx q[2];
rz(-1.6951188) q[2];
sx q[2];
rz(0.97036874) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6183747) q[1];
sx q[1];
rz(-1.8357539) q[1];
sx q[1];
rz(2.9441903) q[1];
rz(-pi) q[2];
x q[2];
rz(0.048681569) q[3];
sx q[3];
rz(-1.4845856) q[3];
sx q[3];
rz(1.8415368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4703579) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(1.3746369) q[2];
rz(-1.6124407) q[3];
sx q[3];
rz(-2.0917454) q[3];
sx q[3];
rz(3.0488739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2328211) q[0];
sx q[0];
rz(-3.0290373) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(0.96083653) q[1];
sx q[1];
rz(-1.6583534) q[1];
sx q[1];
rz(0.94295162) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43555194) q[0];
sx q[0];
rz(-1.0953951) q[0];
sx q[0];
rz(0.33669223) q[0];
rz(-0.93859886) q[2];
sx q[2];
rz(-1.7603121) q[2];
sx q[2];
rz(-2.0507484) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3642204) q[1];
sx q[1];
rz(-2.1465214) q[1];
sx q[1];
rz(0.74933021) q[1];
rz(-pi) q[2];
rz(1.2632964) q[3];
sx q[3];
rz(-1.6700744) q[3];
sx q[3];
rz(-0.57534522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97041398) q[2];
sx q[2];
rz(-1.3033988) q[2];
sx q[2];
rz(-1.1478434) q[2];
rz(2.7796699) q[3];
sx q[3];
rz(-1.3779093) q[3];
sx q[3];
rz(1.743478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4107133) q[0];
sx q[0];
rz(-0.35922265) q[0];
sx q[0];
rz(-3.0391589) q[0];
rz(-0.61406413) q[1];
sx q[1];
rz(-1.026261) q[1];
sx q[1];
rz(-2.3238497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2376643) q[0];
sx q[0];
rz(-0.92887628) q[0];
sx q[0];
rz(-0.68591811) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7162343) q[2];
sx q[2];
rz(-2.465473) q[2];
sx q[2];
rz(0.506625) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35303822) q[1];
sx q[1];
rz(-0.8416881) q[1];
sx q[1];
rz(-0.025351449) q[1];
rz(0.46157591) q[3];
sx q[3];
rz(-2.1510923) q[3];
sx q[3];
rz(3.1212774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3488591) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(-0.31361541) q[2];
rz(0.46401986) q[3];
sx q[3];
rz(-0.84657621) q[3];
sx q[3];
rz(2.2217506) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87026507) q[0];
sx q[0];
rz(-0.79065228) q[0];
sx q[0];
rz(2.9050997) q[0];
rz(0.21367167) q[1];
sx q[1];
rz(-0.90286076) q[1];
sx q[1];
rz(2.9170091) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7296627) q[0];
sx q[0];
rz(-0.87422919) q[0];
sx q[0];
rz(2.4921472) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5245536) q[2];
sx q[2];
rz(-0.65215014) q[2];
sx q[2];
rz(2.1422276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1285291) q[1];
sx q[1];
rz(-2.6010102) q[1];
sx q[1];
rz(2.9180727) q[1];
rz(-0.8524695) q[3];
sx q[3];
rz(-1.895088) q[3];
sx q[3];
rz(1.7773903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1647722) q[2];
sx q[2];
rz(-2.2701264) q[2];
sx q[2];
rz(2.9898047) q[2];
rz(-3.1101036) q[3];
sx q[3];
rz(-1.3969235) q[3];
sx q[3];
rz(3.013124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(-2.0634154) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(1.0833441) q[1];
sx q[1];
rz(-1.5768408) q[1];
sx q[1];
rz(1.5595938) q[1];
rz(-2.940098) q[2];
sx q[2];
rz(-1.6215743) q[2];
sx q[2];
rz(-1.8334186) q[2];
rz(-0.29240378) q[3];
sx q[3];
rz(-0.3429827) q[3];
sx q[3];
rz(-0.61957785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
