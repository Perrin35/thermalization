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
rz(0.39245519) q[0];
sx q[0];
rz(-0.28372228) q[0];
sx q[0];
rz(-1.9537881) q[0];
rz(2.6999733) q[1];
sx q[1];
rz(-1.266357) q[1];
sx q[1];
rz(-1.1511572) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1063163) q[0];
sx q[0];
rz(-1.6589249) q[0];
sx q[0];
rz(0.022038925) q[0];
rz(2.5106836) q[2];
sx q[2];
rz(-2.563347) q[2];
sx q[2];
rz(-1.9672245) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3656039) q[1];
sx q[1];
rz(-1.1465441) q[1];
sx q[1];
rz(-1.7107023) q[1];
rz(-pi) q[2];
rz(1.0813484) q[3];
sx q[3];
rz(-0.57751211) q[3];
sx q[3];
rz(-2.5546527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2426131) q[2];
sx q[2];
rz(-2.2448764) q[2];
sx q[2];
rz(2.8420281) q[2];
rz(2.7859531) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(-0.43660823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3886609) q[0];
sx q[0];
rz(-2.5932073) q[0];
sx q[0];
rz(2.8652628) q[0];
rz(-3.062124) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(-0.4963378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68688289) q[0];
sx q[0];
rz(-1.0609896) q[0];
sx q[0];
rz(-2.8532993) q[0];
rz(0.21644105) q[2];
sx q[2];
rz(-2.7958044) q[2];
sx q[2];
rz(-1.8255289) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73442921) q[1];
sx q[1];
rz(-1.2447272) q[1];
sx q[1];
rz(0.77496518) q[1];
rz(-0.0032418781) q[3];
sx q[3];
rz(-1.4088944) q[3];
sx q[3];
rz(-1.9954132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3695099) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(-0.17862865) q[2];
rz(-1.6102839) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(-2.3791651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4302706) q[0];
sx q[0];
rz(-1.0866168) q[0];
sx q[0];
rz(-1.7311199) q[0];
rz(-1.3871644) q[1];
sx q[1];
rz(-1.6491363) q[1];
sx q[1];
rz(-1.656146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5532865) q[0];
sx q[0];
rz(-1.5380412) q[0];
sx q[0];
rz(2.4577228) q[0];
rz(-pi) q[1];
rz(2.335821) q[2];
sx q[2];
rz(-1.8745595) q[2];
sx q[2];
rz(0.31651503) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8285813) q[1];
sx q[1];
rz(-1.269956) q[1];
sx q[1];
rz(-1.1238696) q[1];
rz(-2.2353506) q[3];
sx q[3];
rz(-1.8633047) q[3];
sx q[3];
rz(-0.14636235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2014655) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(-0.49261576) q[2];
rz(2.7228739) q[3];
sx q[3];
rz(-2.1188354) q[3];
sx q[3];
rz(-0.59188265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017460499) q[0];
sx q[0];
rz(-3.0191665) q[0];
sx q[0];
rz(-0.19805743) q[0];
rz(2.4824704) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(0.19703279) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27749946) q[0];
sx q[0];
rz(-2.6098893) q[0];
sx q[0];
rz(-0.066389991) q[0];
rz(-1.8577099) q[2];
sx q[2];
rz(-1.3409145) q[2];
sx q[2];
rz(0.89864602) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9024762) q[1];
sx q[1];
rz(-0.65985876) q[1];
sx q[1];
rz(-0.59999864) q[1];
x q[2];
rz(1.5102618) q[3];
sx q[3];
rz(-0.54577561) q[3];
sx q[3];
rz(2.7013403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3531602) q[2];
sx q[2];
rz(-1.8400695) q[2];
sx q[2];
rz(1.2737459) q[2];
rz(1.5431917) q[3];
sx q[3];
rz(-1.9485995) q[3];
sx q[3];
rz(-0.25632349) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0850328) q[0];
sx q[0];
rz(-1.4075449) q[0];
sx q[0];
rz(-3.0388167) q[0];
rz(-0.13725266) q[1];
sx q[1];
rz(-2.6922701) q[1];
sx q[1];
rz(1.7030254) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53046561) q[0];
sx q[0];
rz(-1.0004064) q[0];
sx q[0];
rz(2.3378951) q[0];
x q[1];
rz(-1.5248143) q[2];
sx q[2];
rz(-2.0334938) q[2];
sx q[2];
rz(1.3563987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4754191) q[1];
sx q[1];
rz(-1.3120781) q[1];
sx q[1];
rz(2.9889876) q[1];
rz(-pi) q[2];
rz(-0.73682555) q[3];
sx q[3];
rz(-2.9758436) q[3];
sx q[3];
rz(2.152657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30456257) q[2];
sx q[2];
rz(-1.4745159) q[2];
sx q[2];
rz(-2.4414818) q[2];
rz(2.2434798) q[3];
sx q[3];
rz(-0.63106314) q[3];
sx q[3];
rz(1.496544) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9333618) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(1.5775648) q[0];
rz(0.6598407) q[1];
sx q[1];
rz(-2.3520062) q[1];
sx q[1];
rz(0.97445625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0147043) q[0];
sx q[0];
rz(-3.1278353) q[0];
sx q[0];
rz(2.5123054) q[0];
rz(-1.3979369) q[2];
sx q[2];
rz(-1.0893911) q[2];
sx q[2];
rz(2.0916033) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9012857) q[1];
sx q[1];
rz(-2.0723596) q[1];
sx q[1];
rz(1.5830513) q[1];
rz(1.5100689) q[3];
sx q[3];
rz(-2.1083045) q[3];
sx q[3];
rz(2.2030061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4727891) q[2];
sx q[2];
rz(-2.4432224) q[2];
sx q[2];
rz(-0.93639708) q[2];
rz(-0.89013731) q[3];
sx q[3];
rz(-1.7809159) q[3];
sx q[3];
rz(-0.20588017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6137921) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(3.0278681) q[0];
rz(-1.744489) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(-0.024959175) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67201383) q[0];
sx q[0];
rz(-0.22525283) q[0];
sx q[0];
rz(1.3242871) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8309085) q[2];
sx q[2];
rz(-0.94365135) q[2];
sx q[2];
rz(2.8994034) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0164707) q[1];
sx q[1];
rz(-0.85713398) q[1];
sx q[1];
rz(3.0550856) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9437441) q[3];
sx q[3];
rz(-1.6258844) q[3];
sx q[3];
rz(2.4428989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52594841) q[2];
sx q[2];
rz(-1.7082461) q[2];
sx q[2];
rz(1.0053267) q[2];
rz(-2.2510236) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(0.65517202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5263851) q[0];
sx q[0];
rz(-2.7387775) q[0];
sx q[0];
rz(2.030754) q[0];
rz(0.088518294) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(1.7875338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18830472) q[0];
sx q[0];
rz(-0.79995709) q[0];
sx q[0];
rz(-0.15733457) q[0];
rz(2.0769139) q[2];
sx q[2];
rz(-0.30914834) q[2];
sx q[2];
rz(0.34999166) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7343261) q[1];
sx q[1];
rz(-1.3248492) q[1];
sx q[1];
rz(1.6489563) q[1];
x q[2];
rz(0.99127524) q[3];
sx q[3];
rz(-2.2439828) q[3];
sx q[3];
rz(1.1330825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2985349) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(2.7561772) q[2];
rz(-2.391583) q[3];
sx q[3];
rz(-2.1842726) q[3];
sx q[3];
rz(-3.0756522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176004) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(-0.83220926) q[0];
rz(2.0596793) q[1];
sx q[1];
rz(-0.81580201) q[1];
sx q[1];
rz(-1.5218511) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31829054) q[0];
sx q[0];
rz(-2.5099237) q[0];
sx q[0];
rz(-2.2715827) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84110259) q[2];
sx q[2];
rz(-2.7698467) q[2];
sx q[2];
rz(-1.8018307) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9389006) q[1];
sx q[1];
rz(-2.1317516) q[1];
sx q[1];
rz(2.2748018) q[1];
rz(-pi) q[2];
rz(2.9004292) q[3];
sx q[3];
rz(-2.2134288) q[3];
sx q[3];
rz(-1.3295028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99721432) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.2388371) q[2];
rz(-1.8501806) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(-0.13135697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111897) q[0];
sx q[0];
rz(-0.953453) q[0];
sx q[0];
rz(-0.74380547) q[0];
rz(-0.047608308) q[1];
sx q[1];
rz(-1.1046537) q[1];
sx q[1];
rz(-1.1415175) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7615712) q[0];
sx q[0];
rz(-2.0788045) q[0];
sx q[0];
rz(0.53701274) q[0];
x q[1];
rz(1.9594934) q[2];
sx q[2];
rz(-2.1517589) q[2];
sx q[2];
rz(-2.3240391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7773173) q[1];
sx q[1];
rz(-2.2688443) q[1];
sx q[1];
rz(2.1080522) q[1];
rz(-pi) q[2];
rz(-1.3895681) q[3];
sx q[3];
rz(-2.2076105) q[3];
sx q[3];
rz(-0.94094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6514674) q[2];
sx q[2];
rz(-2.9059124) q[2];
sx q[2];
rz(2.8362595) q[2];
rz(1.4281645) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(1.292424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3931428) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(-2.582386) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(-2.7087173) q[2];
sx q[2];
rz(-0.71472157) q[2];
sx q[2];
rz(1.0403487) q[2];
rz(0.84973002) q[3];
sx q[3];
rz(-0.17912514) q[3];
sx q[3];
rz(2.9368243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
