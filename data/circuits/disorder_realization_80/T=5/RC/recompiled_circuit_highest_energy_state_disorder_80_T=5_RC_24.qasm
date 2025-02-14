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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1063163) q[0];
sx q[0];
rz(-1.6589249) q[0];
sx q[0];
rz(0.022038925) q[0];
rz(-1.2032937) q[2];
sx q[2];
rz(-1.1137059) q[2];
sx q[2];
rz(-0.45705308) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3656039) q[1];
sx q[1];
rz(-1.9950486) q[1];
sx q[1];
rz(1.4308903) q[1];
x q[2];
rz(2.8443195) q[3];
sx q[3];
rz(-1.0680388) q[3];
sx q[3];
rz(1.9883324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89897951) q[2];
sx q[2];
rz(-0.8967163) q[2];
sx q[2];
rz(2.8420281) q[2];
rz(-2.7859531) q[3];
sx q[3];
rz(-2.1483597) q[3];
sx q[3];
rz(-2.7049844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3886609) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(-0.27632982) q[0];
rz(-0.07946864) q[1];
sx q[1];
rz(-2.6092968) q[1];
sx q[1];
rz(-0.4963378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4014028) q[0];
sx q[0];
rz(-1.3200134) q[0];
sx q[0];
rz(1.0428421) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4935819) q[2];
sx q[2];
rz(-1.9081915) q[2];
sx q[2];
rz(1.545662) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53238324) q[1];
sx q[1];
rz(-0.84621284) q[1];
sx q[1];
rz(2.0128472) q[1];
rz(-pi) q[2];
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
rz(-0.77208272) q[2];
sx q[2];
rz(-1.778506) q[2];
sx q[2];
rz(2.962964) q[2];
rz(1.5313088) q[3];
sx q[3];
rz(-2.2154) q[3];
sx q[3];
rz(-0.76242751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4302706) q[0];
sx q[0];
rz(-1.0866168) q[0];
sx q[0];
rz(-1.4104728) q[0];
rz(1.3871644) q[1];
sx q[1];
rz(-1.4924563) q[1];
sx q[1];
rz(14/(3*pi)) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091772) q[0];
sx q[0];
rz(-0.88736358) q[0];
sx q[0];
rz(-1.6130436) q[0];
x q[1];
rz(0.40992592) q[2];
sx q[2];
rz(-2.2927444) q[2];
sx q[2];
rz(0.97463911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8456257) q[1];
sx q[1];
rz(-0.53301552) q[1];
sx q[1];
rz(-0.94818268) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7762122) q[3];
sx q[3];
rz(-2.2024904) q[3];
sx q[3];
rz(-1.2022579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9401271) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(0.49261576) q[2];
rz(0.41871873) q[3];
sx q[3];
rz(-2.1188354) q[3];
sx q[3];
rz(-2.54971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.017460499) q[0];
sx q[0];
rz(-3.0191665) q[0];
sx q[0];
rz(2.9435352) q[0];
rz(-0.65912229) q[1];
sx q[1];
rz(-0.56050572) q[1];
sx q[1];
rz(2.9445599) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448387) q[0];
sx q[0];
rz(-1.0403883) q[0];
sx q[0];
rz(-1.5317937) q[0];
x q[1];
rz(0.23931673) q[2];
sx q[2];
rz(-1.8499568) q[2];
sx q[2];
rz(-2.5365732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1888194) q[1];
sx q[1];
rz(-1.0403301) q[1];
sx q[1];
rz(1.1578881) q[1];
rz(3.1048685) q[3];
sx q[3];
rz(-2.1154599) q[3];
sx q[3];
rz(2.7721318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.78843242) q[2];
sx q[2];
rz(-1.3015231) q[2];
sx q[2];
rz(1.8678467) q[2];
rz(-1.5984009) q[3];
sx q[3];
rz(-1.1929932) q[3];
sx q[3];
rz(-2.8852692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0850328) q[0];
sx q[0];
rz(-1.7340478) q[0];
sx q[0];
rz(-0.10277596) q[0];
rz(0.13725266) q[1];
sx q[1];
rz(-0.44932258) q[1];
sx q[1];
rz(-1.4385673) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.611821) q[0];
sx q[0];
rz(-0.91982929) q[0];
sx q[0];
rz(-2.3168827) q[0];
rz(3.0497018) q[2];
sx q[2];
rz(-0.46481195) q[2];
sx q[2];
rz(1.6824695) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4754191) q[1];
sx q[1];
rz(-1.8295145) q[1];
sx q[1];
rz(0.15260507) q[1];
rz(3.0183295) q[3];
sx q[3];
rz(-1.4597037) q[3];
sx q[3];
rz(0.14810239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8370301) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(2.4414818) q[2];
rz(0.89811283) q[3];
sx q[3];
rz(-0.63106314) q[3];
sx q[3];
rz(1.6450487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20823088) q[0];
sx q[0];
rz(-0.58657402) q[0];
sx q[0];
rz(1.5775648) q[0];
rz(0.6598407) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(-0.97445625) q[1];
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
rz(1.5626989) q[0];
rz(-pi) q[1];
rz(0.48759385) q[2];
sx q[2];
rz(-1.7238443) q[2];
sx q[2];
rz(-2.7014554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8169957) q[1];
sx q[1];
rz(-1.5600509) q[1];
sx q[1];
rz(2.6399977) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0401214) q[3];
sx q[3];
rz(-2.6009998) q[3];
sx q[3];
rz(1.056788) q[3];
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
rz(0.93639708) q[2];
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
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6137921) q[0];
sx q[0];
rz(-0.19659909) q[0];
sx q[0];
rz(0.11372456) q[0];
rz(-1.3971036) q[1];
sx q[1];
rz(-2.0753658) q[1];
sx q[1];
rz(-3.1166335) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.002288) q[0];
sx q[0];
rz(-1.5162667) q[0];
sx q[0];
rz(1.3521347) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34100248) q[2];
sx q[2];
rz(-2.4694168) q[2];
sx q[2];
rz(-2.9580046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88475591) q[1];
sx q[1];
rz(-2.4236226) q[1];
sx q[1];
rz(1.6702449) q[1];
x q[2];
rz(-1.6269782) q[3];
sx q[3];
rz(-1.7683408) q[3];
sx q[3];
rz(-2.2584525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152076) q[0];
sx q[0];
rz(-2.7387775) q[0];
sx q[0];
rz(2.030754) q[0];
rz(3.0530744) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(-1.7875338) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41218573) q[0];
sx q[0];
rz(-0.78347396) q[0];
sx q[0];
rz(-1.7307348) q[0];
x q[1];
rz(-1.8432003) q[2];
sx q[2];
rz(-1.4227616) q[2];
sx q[2];
rz(1.7065976) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7343261) q[1];
sx q[1];
rz(-1.3248492) q[1];
sx q[1];
rz(-1.6489563) q[1];
rz(0.60175394) q[3];
sx q[3];
rz(-2.2838998) q[3];
sx q[3];
rz(-2.8181638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84305772) q[2];
sx q[2];
rz(-2.3044846) q[2];
sx q[2];
rz(0.38541547) q[2];
rz(-2.391583) q[3];
sx q[3];
rz(-2.1842726) q[3];
sx q[3];
rz(-3.0756522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239922) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(-2.3093834) q[0];
rz(-2.0596793) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(1.6197416) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2913301) q[0];
sx q[0];
rz(-1.9614152) q[0];
sx q[0];
rz(-1.0608835) q[0];
rz(-pi) q[1];
rz(-0.84110259) q[2];
sx q[2];
rz(-0.37174598) q[2];
sx q[2];
rz(-1.8018307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79241133) q[1];
sx q[1];
rz(-2.1508625) q[1];
sx q[1];
rz(-2.4522454) q[1];
x q[2];
rz(-2.2275739) q[3];
sx q[3];
rz(-1.7631712) q[3];
sx q[3];
rz(-0.094948204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1443783) q[2];
sx q[2];
rz(-1.245456) q[2];
sx q[2];
rz(1.9027556) q[2];
rz(-1.291412) q[3];
sx q[3];
rz(-1.2715128) q[3];
sx q[3];
rz(-3.0102357) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111897) q[0];
sx q[0];
rz(-2.1881396) q[0];
sx q[0];
rz(0.74380547) q[0];
rz(3.0939843) q[1];
sx q[1];
rz(-1.1046537) q[1];
sx q[1];
rz(-1.1415175) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7615712) q[0];
sx q[0];
rz(-1.0627882) q[0];
sx q[0];
rz(-2.6045799) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9594934) q[2];
sx q[2];
rz(-0.98983374) q[2];
sx q[2];
rz(2.3240391) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9824578) q[1];
sx q[1];
rz(-1.1679113) q[1];
sx q[1];
rz(2.3680658) q[1];
rz(-pi) q[2];
rz(2.4968653) q[3];
sx q[3];
rz(-1.4253741) q[3];
sx q[3];
rz(2.4032088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49012524) q[2];
sx q[2];
rz(-2.9059124) q[2];
sx q[2];
rz(2.8362595) q[2];
rz(-1.4281645) q[3];
sx q[3];
rz(-1.7686663) q[3];
sx q[3];
rz(-1.8491687) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(0.66721532) q[2];
sx q[2];
rz(-1.2922774) q[2];
sx q[2];
rz(-0.86624672) q[2];
rz(1.4356267) q[3];
sx q[3];
rz(-1.4528989) q[3];
sx q[3];
rz(-1.0624878) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
