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
rz(1.2394387) q[0];
sx q[0];
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48403063) q[0];
sx q[0];
rz(-3.1320509) q[0];
sx q[0];
rz(1.5242759) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9340408) q[2];
sx q[2];
rz(-2.2170728) q[2];
sx q[2];
rz(2.5763047) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7211044) q[1];
sx q[1];
rz(-1.8454362) q[1];
sx q[1];
rz(-0.65914776) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3232949) q[3];
sx q[3];
rz(-2.1967271) q[3];
sx q[3];
rz(1.4656386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2892896) q[2];
sx q[2];
rz(-1.0785582) q[2];
sx q[2];
rz(3.0618073) q[2];
rz(-0.96528178) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(-0.7307581) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6642283) q[0];
sx q[0];
rz(-2.1381162) q[0];
sx q[0];
rz(-0.088951237) q[0];
rz(-1.3720007) q[1];
sx q[1];
rz(-1.7141432) q[1];
sx q[1];
rz(-1.1044097) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1386344) q[0];
sx q[0];
rz(-1.8002274) q[0];
sx q[0];
rz(-2.5087439) q[0];
rz(-pi) q[1];
rz(-1.0816108) q[2];
sx q[2];
rz(-0.85027611) q[2];
sx q[2];
rz(0.94948506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30518499) q[1];
sx q[1];
rz(-0.47023222) q[1];
sx q[1];
rz(-1.2364619) q[1];
rz(-pi) q[2];
rz(1.6703963) q[3];
sx q[3];
rz(-1.0336813) q[3];
sx q[3];
rz(-2.296534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.264512) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(-1.4667) q[2];
rz(-0.029021164) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(-0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90917176) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(2.8636041) q[0];
rz(1.4311283) q[1];
sx q[1];
rz(-0.93228308) q[1];
sx q[1];
rz(2.7412282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5059197) q[0];
sx q[0];
rz(-2.189496) q[0];
sx q[0];
rz(-2.1257504) q[0];
x q[1];
rz(1.8258926) q[2];
sx q[2];
rz(-2.2410903) q[2];
sx q[2];
rz(-2.3885661) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0644933) q[1];
sx q[1];
rz(-2.0410186) q[1];
sx q[1];
rz(-1.3920067) q[1];
x q[2];
rz(0.34137643) q[3];
sx q[3];
rz(-1.2814643) q[3];
sx q[3];
rz(-3.0738664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.532734) q[2];
sx q[2];
rz(-2.686783) q[2];
rz(1.8999892) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(-2.4212867) q[3];
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
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95551816) q[0];
sx q[0];
rz(-1.8121413) q[0];
sx q[0];
rz(-3.1410826) q[0];
rz(0.60091248) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(0.14437637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6116981) q[0];
sx q[0];
rz(-1.3641883) q[0];
sx q[0];
rz(2.1312461) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8881489) q[2];
sx q[2];
rz(-1.4254693) q[2];
sx q[2];
rz(-3.1325454) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3444654) q[1];
sx q[1];
rz(-0.68736156) q[1];
sx q[1];
rz(-0.3631773) q[1];
rz(2.9668429) q[3];
sx q[3];
rz(-2.5876382) q[3];
sx q[3];
rz(-2.9308211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.942975) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(0.96251881) q[2];
rz(-1.6019542) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(2.8008154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9050423) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(1.9566253) q[0];
rz(0.22625893) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(-2.102898) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.561956) q[0];
sx q[0];
rz(-2.547894) q[0];
sx q[0];
rz(0.14621347) q[0];
x q[1];
rz(2.4657927) q[2];
sx q[2];
rz(-0.86188176) q[2];
sx q[2];
rz(-2.076864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2782368) q[1];
sx q[1];
rz(-2.0788631) q[1];
sx q[1];
rz(-0.79973508) q[1];
rz(2.0108972) q[3];
sx q[3];
rz(-1.6301042) q[3];
sx q[3];
rz(-1.9594994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7096536) q[2];
sx q[2];
rz(-0.69810549) q[2];
sx q[2];
rz(2.8254438) q[2];
rz(1.6759253) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(-1.3795615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50452152) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(-1.1035408) q[0];
rz(-0.48031131) q[1];
sx q[1];
rz(-2.4791398) q[1];
sx q[1];
rz(0.59741098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6107268) q[0];
sx q[0];
rz(-1.7590176) q[0];
sx q[0];
rz(1.3590616) q[0];
x q[1];
rz(-0.74540794) q[2];
sx q[2];
rz(-1.7495724) q[2];
sx q[2];
rz(2.488236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.061042) q[1];
sx q[1];
rz(-1.4054728) q[1];
sx q[1];
rz(-2.8870261) q[1];
rz(0.091986309) q[3];
sx q[3];
rz(-1.9402656) q[3];
sx q[3];
rz(1.2928158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.20208134) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-2.3484774) q[3];
sx q[3];
rz(0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(-0.038473815) q[0];
rz(3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(-2.9077392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9486987) q[0];
sx q[0];
rz(-1.3854376) q[0];
sx q[0];
rz(-0.22484397) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87217561) q[2];
sx q[2];
rz(-2.6769014) q[2];
sx q[2];
rz(1.2810436) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0044627) q[1];
sx q[1];
rz(-0.65441583) q[1];
sx q[1];
rz(-2.9942102) q[1];
rz(-pi) q[2];
rz(1.5130299) q[3];
sx q[3];
rz(-1.4411297) q[3];
sx q[3];
rz(-0.71021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48876277) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-0.59824198) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(-0.84754506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365731) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(-0.2463499) q[0];
rz(-1.8440638) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(-2.6447703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5733872) q[0];
sx q[0];
rz(-2.3620689) q[0];
sx q[0];
rz(-0.60978344) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4364528) q[2];
sx q[2];
rz(-1.6821096) q[2];
sx q[2];
rz(1.4347347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0223479) q[1];
sx q[1];
rz(-1.3356707) q[1];
sx q[1];
rz(-1.0757331) q[1];
rz(-1.2874576) q[3];
sx q[3];
rz(-0.50931286) q[3];
sx q[3];
rz(-2.4278502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7184427) q[2];
sx q[2];
rz(-0.53360525) q[2];
sx q[2];
rz(1.9971087) q[2];
rz(-1.1540958) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(2.9437039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941187) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(0.06614729) q[0];
rz(1.9937531) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(0.60417169) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19949958) q[0];
sx q[0];
rz(-1.2954933) q[0];
sx q[0];
rz(-1.3574187) q[0];
rz(-pi) q[1];
rz(-1.4140903) q[2];
sx q[2];
rz(-2.092166) q[2];
sx q[2];
rz(-0.40796134) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38018885) q[1];
sx q[1];
rz(-1.0127002) q[1];
sx q[1];
rz(1.3437043) q[1];
x q[2];
rz(-0.87852134) q[3];
sx q[3];
rz(-2.1067348) q[3];
sx q[3];
rz(-1.5978447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26312795) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(2.7086332) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1935254) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(1.2731592) q[0];
rz(0.46514568) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(-0.21496162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33782712) q[0];
sx q[0];
rz(-1.8233607) q[0];
sx q[0];
rz(2.4830677) q[0];
rz(-pi) q[1];
rz(2.6097492) q[2];
sx q[2];
rz(-1.1924679) q[2];
sx q[2];
rz(-1.7500306) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7219639) q[1];
sx q[1];
rz(-2.2441494) q[1];
sx q[1];
rz(2.2339905) q[1];
rz(-pi) q[2];
rz(1.8278024) q[3];
sx q[3];
rz(-1.0238092) q[3];
sx q[3];
rz(-1.6099324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2994069) q[2];
sx q[2];
rz(-0.81373787) q[2];
sx q[2];
rz(2.1827533) q[2];
rz(-1.9793319) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(1.4452665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5398298) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(2.4330347) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(-2.6513349) q[2];
sx q[2];
rz(-2.5210862) q[2];
sx q[2];
rz(1.6747337) q[2];
rz(0.21227588) q[3];
sx q[3];
rz(-0.73321453) q[3];
sx q[3];
rz(-1.1238255) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
