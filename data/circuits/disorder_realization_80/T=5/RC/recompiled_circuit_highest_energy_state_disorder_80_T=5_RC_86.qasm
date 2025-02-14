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
rz(1.1878045) q[0];
rz(-0.44161931) q[1];
sx q[1];
rz(-1.8752357) q[1];
sx q[1];
rz(1.1511572) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53357994) q[0];
sx q[0];
rz(-1.5927497) q[0];
sx q[0];
rz(-1.4826464) q[0];
x q[1];
rz(-1.9382989) q[2];
sx q[2];
rz(-1.1137059) q[2];
sx q[2];
rz(-2.6845396) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85271231) q[1];
sx q[1];
rz(-1.4433644) q[1];
sx q[1];
rz(0.42794771) q[1];
x q[2];
rz(-2.0602442) q[3];
sx q[3];
rz(-0.57751211) q[3];
sx q[3];
rz(-2.5546527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2426131) q[2];
sx q[2];
rz(-0.8967163) q[2];
sx q[2];
rz(-0.29956451) q[2];
rz(2.7859531) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(-0.43660823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75293175) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(-2.8652628) q[0];
rz(0.07946864) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(2.6452549) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68688289) q[0];
sx q[0];
rz(-2.080603) q[0];
sx q[0];
rz(-2.8532993) q[0];
x q[1];
rz(2.8032643) q[2];
sx q[2];
rz(-1.4979432) q[2];
sx q[2];
rz(-0.05073994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4071634) q[1];
sx q[1];
rz(-1.8968655) q[1];
sx q[1];
rz(0.77496518) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7326991) q[3];
sx q[3];
rz(-1.5739958) q[3];
sx q[3];
rz(0.42409431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77208272) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(-0.17862865) q[2];
rz(1.6102839) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(2.3791651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302706) q[0];
sx q[0];
rz(-2.0549759) q[0];
sx q[0];
rz(-1.4104728) q[0];
rz(1.3871644) q[1];
sx q[1];
rz(-1.4924563) q[1];
sx q[1];
rz(-1.656146) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1992599) q[0];
sx q[0];
rz(-2.457065) q[0];
sx q[0];
rz(-0.051816414) q[0];
rz(-pi) q[1];
rz(-1.1457655) q[2];
sx q[2];
rz(-0.81159893) q[2];
sx q[2];
rz(-1.5852864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.3988603) q[1];
sx q[1];
rz(-1.1452864) q[1];
sx q[1];
rz(0.33136435) q[1];
rz(-pi) q[2];
rz(1.1165421) q[3];
sx q[3];
rz(-2.4245533) q[3];
sx q[3];
rz(1.7771135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2014655) q[2];
sx q[2];
rz(-1.0041693) q[2];
sx q[2];
rz(2.6489769) q[2];
rz(-0.41871873) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(-2.54971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017460499) q[0];
sx q[0];
rz(-5/(13*pi)) q[0];
sx q[0];
rz(2.9435352) q[0];
rz(0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-0.19703279) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9055419) q[0];
sx q[0];
rz(-1.5371548) q[0];
sx q[0];
rz(2.6108526) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2617499) q[2];
sx q[2];
rz(-2.7759399) q[2];
sx q[2];
rz(-1.8118389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9024762) q[1];
sx q[1];
rz(-2.4817339) q[1];
sx q[1];
rz(2.541594) q[1];
rz(-pi) q[2];
rz(1.0258338) q[3];
sx q[3];
rz(-1.6022046) q[3];
sx q[3];
rz(1.9592913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3531602) q[2];
sx q[2];
rz(-1.3015231) q[2];
sx q[2];
rz(1.2737459) q[2];
rz(-1.5431917) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(-0.25632349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0850328) q[0];
sx q[0];
rz(-1.7340478) q[0];
sx q[0];
rz(0.10277596) q[0];
rz(3.00434) q[1];
sx q[1];
rz(-0.44932258) q[1];
sx q[1];
rz(1.4385673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.611127) q[0];
sx q[0];
rz(-1.0004064) q[0];
sx q[0];
rz(-2.3378951) q[0];
rz(-pi) q[1];
rz(-1.6167783) q[2];
sx q[2];
rz(-2.0334938) q[2];
sx q[2];
rz(1.785194) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2762982) q[1];
sx q[1];
rz(-1.423308) q[1];
sx q[1];
rz(-1.8324204) q[1];
rz(-pi) q[2];
rz(1.4588615) q[3];
sx q[3];
rz(-1.6932958) q[3];
sx q[3];
rz(-1.4089597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.30456257) q[2];
sx q[2];
rz(-1.4745159) q[2];
sx q[2];
rz(0.70011082) q[2];
rz(-2.2434798) q[3];
sx q[3];
rz(-2.5105295) q[3];
sx q[3];
rz(1.496544) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20823088) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(1.5640278) q[0];
rz(2.481752) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(0.97445625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.385372) q[0];
sx q[0];
rz(-1.5819183) q[0];
sx q[0];
rz(-1.5788938) q[0];
x q[1];
rz(1.3979369) q[2];
sx q[2];
rz(-2.0522016) q[2];
sx q[2];
rz(-1.0499894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32459697) q[1];
sx q[1];
rz(-1.5815418) q[1];
sx q[1];
rz(0.50159494) q[1];
x q[2];
rz(-3.0401214) q[3];
sx q[3];
rz(-0.54059282) q[3];
sx q[3];
rz(-2.0848047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4727891) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(0.93639708) q[2];
rz(-2.2514553) q[3];
sx q[3];
rz(-1.7809159) q[3];
sx q[3];
rz(0.20588017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(1.6137921) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(-3.0278681) q[0];
rz(1.744489) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(0.024959175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41938101) q[0];
sx q[0];
rz(-1.352465) q[0];
sx q[0];
rz(0.055856987) q[0];
rz(2.8005902) q[2];
sx q[2];
rz(-2.4694168) q[2];
sx q[2];
rz(-2.9580046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.125122) q[1];
sx q[1];
rz(-0.85713398) q[1];
sx q[1];
rz(-0.086507052) q[1];
rz(-1.5146144) q[3];
sx q[3];
rz(-1.3732519) q[3];
sx q[3];
rz(0.88314013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6156442) q[2];
sx q[2];
rz(-1.7082461) q[2];
sx q[2];
rz(2.136266) q[2];
rz(-2.2510236) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(-2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152076) q[0];
sx q[0];
rz(-2.7387775) q[0];
sx q[0];
rz(2.030754) q[0];
rz(-0.088518294) q[1];
sx q[1];
rz(-0.42461494) q[1];
sx q[1];
rz(1.7875338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.272404) q[0];
sx q[0];
rz(-1.6834295) q[0];
sx q[0];
rz(-2.3478481) q[0];
rz(-pi) q[1];
rz(2.9879775) q[2];
sx q[2];
rz(-1.3014463) q[2];
sx q[2];
rz(-0.17698032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14446324) q[1];
sx q[1];
rz(-1.494993) q[1];
sx q[1];
rz(-0.24667011) q[1];
rz(2.380221) q[3];
sx q[3];
rz(-2.0132445) q[3];
sx q[3];
rz(0.82514742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2985349) q[2];
sx q[2];
rz(-2.3044846) q[2];
sx q[2];
rz(-2.7561772) q[2];
rz(2.391583) q[3];
sx q[3];
rz(-0.95732006) q[3];
sx q[3];
rz(0.065940417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1239922) q[0];
sx q[0];
rz(-1.0753068) q[0];
sx q[0];
rz(-0.83220926) q[0];
rz(2.0596793) q[1];
sx q[1];
rz(-0.81580201) q[1];
sx q[1];
rz(-1.5218511) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31829054) q[0];
sx q[0];
rz(-0.63166895) q[0];
sx q[0];
rz(-0.87000991) q[0];
x q[1];
rz(1.85361) q[2];
sx q[2];
rz(-1.8153802) q[2];
sx q[2];
rz(-2.2158538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9389006) q[1];
sx q[1];
rz(-1.0098411) q[1];
sx q[1];
rz(-0.86679082) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2275739) q[3];
sx q[3];
rz(-1.3784215) q[3];
sx q[3];
rz(3.0466444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99721432) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.9027556) q[2];
rz(1.8501806) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(-3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111897) q[0];
sx q[0];
rz(-0.953453) q[0];
sx q[0];
rz(-2.3977872) q[0];
rz(-0.047608308) q[1];
sx q[1];
rz(-1.1046537) q[1];
sx q[1];
rz(2.0000752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38002148) q[0];
sx q[0];
rz(-2.0788045) q[0];
sx q[0];
rz(-0.53701274) q[0];
x q[1];
rz(0.52351482) q[2];
sx q[2];
rz(-0.68632978) q[2];
sx q[2];
rz(-1.6829784) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1591349) q[1];
sx q[1];
rz(-1.1679113) q[1];
sx q[1];
rz(-2.3680658) q[1];
x q[2];
rz(-2.4968653) q[3];
sx q[3];
rz(-1.7162185) q[3];
sx q[3];
rz(2.4032088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6514674) q[2];
sx q[2];
rz(-2.9059124) q[2];
sx q[2];
rz(2.8362595) q[2];
rz(1.7134282) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(-1.292424) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74844985) q[0];
sx q[0];
rz(-0.23101692) q[0];
sx q[0];
rz(1.1363181) q[0];
rz(-0.5592067) q[1];
sx q[1];
rz(-1.3312085) q[1];
sx q[1];
rz(-2.8937173) q[1];
rz(0.43287533) q[2];
sx q[2];
rz(-0.71472157) q[2];
sx q[2];
rz(1.0403487) q[2];
rz(-0.1189726) q[3];
sx q[3];
rz(-1.705022) q[3];
sx q[3];
rz(-2.6172887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
