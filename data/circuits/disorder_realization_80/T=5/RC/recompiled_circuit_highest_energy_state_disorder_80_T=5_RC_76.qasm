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
rz(-2.7491375) q[0];
sx q[0];
rz(-2.8578704) q[0];
sx q[0];
rz(-1.1878045) q[0];
rz(2.6999733) q[1];
sx q[1];
rz(-1.266357) q[1];
sx q[1];
rz(-1.1511572) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6080127) q[0];
sx q[0];
rz(-1.5488429) q[0];
sx q[0];
rz(1.4826464) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9382989) q[2];
sx q[2];
rz(-2.0278868) q[2];
sx q[2];
rz(2.6845396) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85271231) q[1];
sx q[1];
rz(-1.6982283) q[1];
sx q[1];
rz(2.7136449) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8443195) q[3];
sx q[3];
rz(-1.0680388) q[3];
sx q[3];
rz(-1.1532602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2426131) q[2];
sx q[2];
rz(-0.8967163) q[2];
sx q[2];
rz(2.8420281) q[2];
rz(0.35563955) q[3];
sx q[3];
rz(-0.99323291) q[3];
sx q[3];
rz(2.7049844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.75293175) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(-0.27632982) q[0];
rz(0.07946864) q[1];
sx q[1];
rz(-0.53229585) q[1];
sx q[1];
rz(-0.4963378) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68688289) q[0];
sx q[0];
rz(-1.0609896) q[0];
sx q[0];
rz(0.28829337) q[0];
rz(1.6480108) q[2];
sx q[2];
rz(-1.2334012) q[2];
sx q[2];
rz(-1.5959306) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6092094) q[1];
sx q[1];
rz(-0.84621284) q[1];
sx q[1];
rz(-2.0128472) q[1];
x q[2];
rz(-1.5906421) q[3];
sx q[3];
rz(-2.9796585) q[3];
sx q[3];
rz(1.1662883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77208272) q[2];
sx q[2];
rz(-1.3630867) q[2];
sx q[2];
rz(-2.962964) q[2];
rz(1.5313088) q[3];
sx q[3];
rz(-0.9261927) q[3];
sx q[3];
rz(0.76242751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71132201) q[0];
sx q[0];
rz(-2.0549759) q[0];
sx q[0];
rz(1.7311199) q[0];
rz(1.7544282) q[1];
sx q[1];
rz(-1.4924563) q[1];
sx q[1];
rz(-14/(3*pi)) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58830611) q[0];
sx q[0];
rz(-1.6035514) q[0];
sx q[0];
rz(0.68386987) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1457655) q[2];
sx q[2];
rz(-0.81159893) q[2];
sx q[2];
rz(-1.5563063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8456257) q[1];
sx q[1];
rz(-0.53301552) q[1];
sx q[1];
rz(2.19341) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90624202) q[3];
sx q[3];
rz(-1.2782879) q[3];
sx q[3];
rz(0.14636235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9401271) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(-0.49261576) q[2];
rz(0.41871873) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(2.54971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640932) q[0];
sx q[0];
rz(-2.6098893) q[0];
sx q[0];
rz(-3.0752027) q[0];
rz(2.9022759) q[2];
sx q[2];
rz(-1.2916358) q[2];
sx q[2];
rz(0.60501947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1888194) q[1];
sx q[1];
rz(-1.0403301) q[1];
sx q[1];
rz(1.9837045) q[1];
x q[2];
rz(1.0258338) q[3];
sx q[3];
rz(-1.539388) q[3];
sx q[3];
rz(1.1823013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3531602) q[2];
sx q[2];
rz(-1.8400695) q[2];
sx q[2];
rz(1.8678467) q[2];
rz(1.5431917) q[3];
sx q[3];
rz(-1.9485995) q[3];
sx q[3];
rz(2.8852692) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0850328) q[0];
sx q[0];
rz(-1.7340478) q[0];
sx q[0];
rz(0.10277596) q[0];
rz(-3.00434) q[1];
sx q[1];
rz(-0.44932258) q[1];
sx q[1];
rz(1.7030254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5202839) q[0];
sx q[0];
rz(-0.94693333) q[0];
sx q[0];
rz(2.413723) q[0];
rz(0.091890865) q[2];
sx q[2];
rz(-2.6767807) q[2];
sx q[2];
rz(1.6824695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2074429) q[1];
sx q[1];
rz(-2.8420919) q[1];
sx q[1];
rz(1.0494166) q[1];
rz(-pi) q[2];
rz(-2.4047671) q[3];
sx q[3];
rz(-0.16574905) q[3];
sx q[3];
rz(-0.98893569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30456257) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(-2.4414818) q[2];
rz(-0.89811283) q[3];
sx q[3];
rz(-2.5105295) q[3];
sx q[3];
rz(-1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9333618) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(1.5775648) q[0];
rz(-2.481752) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(2.1671364) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0147043) q[0];
sx q[0];
rz(-0.013757324) q[0];
sx q[0];
rz(-0.62928726) q[0];
rz(-1.7436557) q[2];
sx q[2];
rz(-2.0522016) q[2];
sx q[2];
rz(2.0916033) q[2];
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
rz(0.53831999) q[3];
sx q[3];
rz(-1.5186406) q[3];
sx q[3];
rz(2.5405034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66880354) q[2];
sx q[2];
rz(-0.69837022) q[2];
sx q[2];
rz(-2.2051956) q[2];
rz(0.89013731) q[3];
sx q[3];
rz(-1.7809159) q[3];
sx q[3];
rz(-2.9357125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6137921) q[0];
sx q[0];
rz(-2.9449936) q[0];
sx q[0];
rz(3.0278681) q[0];
rz(1.3971036) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(3.1166335) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4695788) q[0];
sx q[0];
rz(-2.9163398) q[0];
sx q[0];
rz(1.8173056) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34100248) q[2];
sx q[2];
rz(-2.4694168) q[2];
sx q[2];
rz(2.9580046) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2568367) q[1];
sx q[1];
rz(-0.71797007) q[1];
sx q[1];
rz(-1.6702449) q[1];
rz(1.5146144) q[3];
sx q[3];
rz(-1.3732519) q[3];
sx q[3];
rz(2.2584525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52594841) q[2];
sx q[2];
rz(-1.7082461) q[2];
sx q[2];
rz(-1.0053267) q[2];
rz(2.2510236) q[3];
sx q[3];
rz(-1.4054479) q[3];
sx q[3];
rz(2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152076) q[0];
sx q[0];
rz(-2.7387775) q[0];
sx q[0];
rz(-1.1108387) q[0];
rz(-0.088518294) q[1];
sx q[1];
rz(-0.42461494) q[1];
sx q[1];
rz(-1.3540589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18830472) q[0];
sx q[0];
rz(-2.3416356) q[0];
sx q[0];
rz(-2.9842581) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1536151) q[2];
sx q[2];
rz(-1.8401463) q[2];
sx q[2];
rz(0.17698032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.718487) q[1];
sx q[1];
rz(-0.25782789) q[1];
sx q[1];
rz(2.8400374) q[1];
x q[2];
rz(-2.1503174) q[3];
sx q[3];
rz(-2.2439828) q[3];
sx q[3];
rz(1.1330825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2985349) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(-2.7561772) q[2];
rz(-2.391583) q[3];
sx q[3];
rz(-0.95732006) q[3];
sx q[3];
rz(3.0756522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0176004) q[0];
sx q[0];
rz(-1.0753068) q[0];
sx q[0];
rz(2.3093834) q[0];
rz(-1.0819134) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(-1.6197416) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523217) q[0];
sx q[0];
rz(-1.1025363) q[0];
sx q[0];
rz(0.44083198) q[0];
rz(-0.84110259) q[2];
sx q[2];
rz(-0.37174598) q[2];
sx q[2];
rz(-1.8018307) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9389006) q[1];
sx q[1];
rz(-2.1317516) q[1];
sx q[1];
rz(2.2748018) q[1];
rz(-pi) q[2];
rz(1.8796104) q[3];
sx q[3];
rz(-2.4612456) q[3];
sx q[3];
rz(1.4226567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1443783) q[2];
sx q[2];
rz(-1.8961366) q[2];
sx q[2];
rz(1.2388371) q[2];
rz(-1.8501806) q[3];
sx q[3];
rz(-1.8700799) q[3];
sx q[3];
rz(-3.0102357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111897) q[0];
sx q[0];
rz(-0.953453) q[0];
sx q[0];
rz(-2.3977872) q[0];
rz(3.0939843) q[1];
sx q[1];
rz(-2.036939) q[1];
sx q[1];
rz(1.1415175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50573802) q[0];
sx q[0];
rz(-2.4199652) q[0];
sx q[0];
rz(0.82765915) q[0];
x q[1];
rz(2.6180778) q[2];
sx q[2];
rz(-0.68632978) q[2];
sx q[2];
rz(-1.4586142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9824578) q[1];
sx q[1];
rz(-1.1679113) q[1];
sx q[1];
rz(0.77352685) q[1];
rz(-pi) q[2];
rz(0.23903592) q[3];
sx q[3];
rz(-0.65863684) q[3];
sx q[3];
rz(-0.64202587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6514674) q[2];
sx q[2];
rz(-0.23568025) q[2];
sx q[2];
rz(2.8362595) q[2];
rz(1.4281645) q[3];
sx q[3];
rz(-1.3729264) q[3];
sx q[3];
rz(-1.8491687) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
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
rz(-2.582386) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(-0.43287533) q[2];
sx q[2];
rz(-2.4268711) q[2];
sx q[2];
rz(-2.101244) q[2];
rz(-1.4356267) q[3];
sx q[3];
rz(-1.6886938) q[3];
sx q[3];
rz(2.0791048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
