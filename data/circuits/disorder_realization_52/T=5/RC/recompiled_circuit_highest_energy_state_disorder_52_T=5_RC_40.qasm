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
rz(1.7088543) q[0];
sx q[0];
rz(-2.813485) q[0];
sx q[0];
rz(-0.83972591) q[0];
rz(1.4108763) q[1];
sx q[1];
rz(1.6003992) q[1];
sx q[1];
rz(8.6551854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107723) q[0];
sx q[0];
rz(-1.4913173) q[0];
sx q[0];
rz(3.1067418) q[0];
x q[1];
rz(1.5430449) q[2];
sx q[2];
rz(-2.395438) q[2];
sx q[2];
rz(-2.1389824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2857267) q[1];
sx q[1];
rz(-1.9225504) q[1];
sx q[1];
rz(0.00067623059) q[1];
x q[2];
rz(-2.0887718) q[3];
sx q[3];
rz(-2.149579) q[3];
sx q[3];
rz(-0.56182789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2644234) q[2];
sx q[2];
rz(-1.8921655) q[2];
sx q[2];
rz(1.9682311) q[2];
rz(-0.42826432) q[3];
sx q[3];
rz(-0.56848017) q[3];
sx q[3];
rz(-2.4885524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0083171) q[0];
sx q[0];
rz(-0.44047099) q[0];
sx q[0];
rz(1.0622729) q[0];
rz(1.5509037) q[1];
sx q[1];
rz(-0.56113243) q[1];
sx q[1];
rz(-2.0468457) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8991722) q[0];
sx q[0];
rz(-1.8255485) q[0];
sx q[0];
rz(-0.052847742) q[0];
x q[1];
rz(-2.1065478) q[2];
sx q[2];
rz(-0.75281116) q[2];
sx q[2];
rz(2.8967146) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9506663) q[1];
sx q[1];
rz(-2.6537173) q[1];
sx q[1];
rz(0.64779727) q[1];
rz(-pi) q[2];
rz(2.3660819) q[3];
sx q[3];
rz(-1.417096) q[3];
sx q[3];
rz(-2.6284144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.67794472) q[2];
sx q[2];
rz(-0.28767583) q[2];
sx q[2];
rz(-2.2410683) q[2];
rz(-1.0362222) q[3];
sx q[3];
rz(-3.0039054) q[3];
sx q[3];
rz(-0.010995939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66505945) q[0];
sx q[0];
rz(-1.1256951) q[0];
sx q[0];
rz(1.6735459) q[0];
rz(2.2131069) q[1];
sx q[1];
rz(-2.1378345) q[1];
sx q[1];
rz(1.7195864) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30303653) q[0];
sx q[0];
rz(-1.8183129) q[0];
sx q[0];
rz(-1.4776138) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2879351) q[2];
sx q[2];
rz(-0.89353925) q[2];
sx q[2];
rz(3.1129865) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.068068935) q[1];
sx q[1];
rz(-1.0367107) q[1];
sx q[1];
rz(1.1492386) q[1];
rz(-0.57392759) q[3];
sx q[3];
rz(-0.49387723) q[3];
sx q[3];
rz(-0.18791325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9906087) q[2];
sx q[2];
rz(-1.432212) q[2];
sx q[2];
rz(-0.81986156) q[2];
rz(3.0153583) q[3];
sx q[3];
rz(-1.9641967) q[3];
sx q[3];
rz(-1.646515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6821297) q[0];
sx q[0];
rz(-1.7403025) q[0];
sx q[0];
rz(1.5392186) q[0];
rz(1.0914717) q[1];
sx q[1];
rz(-2.3075054) q[1];
sx q[1];
rz(1.1702671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71890408) q[0];
sx q[0];
rz(-0.50245521) q[0];
sx q[0];
rz(-0.86645856) q[0];
rz(-pi) q[1];
x q[1];
rz(2.290904) q[2];
sx q[2];
rz(-1.7682835) q[2];
sx q[2];
rz(-0.46839223) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2728426) q[1];
sx q[1];
rz(-3.0210872) q[1];
sx q[1];
rz(2.4122448) q[1];
rz(-pi) q[2];
rz(1.7863196) q[3];
sx q[3];
rz(-1.2053262) q[3];
sx q[3];
rz(-2.6579223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9200661) q[2];
sx q[2];
rz(-2.2130794) q[2];
sx q[2];
rz(-0.083219223) q[2];
rz(2.4890238) q[3];
sx q[3];
rz(-0.021952732) q[3];
sx q[3];
rz(0.86161247) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6457152) q[0];
sx q[0];
rz(-2.8577514) q[0];
sx q[0];
rz(2.9448217) q[0];
rz(-1.1605877) q[1];
sx q[1];
rz(-0.44191688) q[1];
sx q[1];
rz(-1.388185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8854839) q[0];
sx q[0];
rz(-1.5692595) q[0];
sx q[0];
rz(-2.9716757) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9391727) q[2];
sx q[2];
rz(-1.4223863) q[2];
sx q[2];
rz(-0.14904101) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1886471) q[1];
sx q[1];
rz(-1.398723) q[1];
sx q[1];
rz(-2.2285372) q[1];
x q[2];
rz(-2.7775473) q[3];
sx q[3];
rz(-0.87556404) q[3];
sx q[3];
rz(2.1882868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56472003) q[2];
sx q[2];
rz(-1.9331965) q[2];
sx q[2];
rz(-1.8604856) q[2];
rz(-2.7992904) q[3];
sx q[3];
rz(-0.050693158) q[3];
sx q[3];
rz(2.2127693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36011919) q[0];
sx q[0];
rz(-2.0442648) q[0];
sx q[0];
rz(-2.1088364) q[0];
rz(0.67689854) q[1];
sx q[1];
rz(-2.758226) q[1];
sx q[1];
rz(-2.3597609) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14836687) q[0];
sx q[0];
rz(-1.2162104) q[0];
sx q[0];
rz(0.605159) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1091613) q[2];
sx q[2];
rz(-1.6219181) q[2];
sx q[2];
rz(-2.0699208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1345152) q[1];
sx q[1];
rz(-0.65943968) q[1];
sx q[1];
rz(-2.6263528) q[1];
x q[2];
rz(1.5277465) q[3];
sx q[3];
rz(-1.2331748) q[3];
sx q[3];
rz(1.4211602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3013762) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(2.2678579) q[2];
rz(-0.78911632) q[3];
sx q[3];
rz(-1.5849761) q[3];
sx q[3];
rz(2.6553787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9588722) q[0];
sx q[0];
rz(-0.046165753) q[0];
sx q[0];
rz(-3.0507372) q[0];
rz(0.60984045) q[1];
sx q[1];
rz(-1.5515386) q[1];
sx q[1];
rz(-0.59744936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4700902) q[0];
sx q[0];
rz(-1.4250942) q[0];
sx q[0];
rz(-1.7007366) q[0];
rz(-pi) q[1];
rz(-0.33281271) q[2];
sx q[2];
rz(-1.2166585) q[2];
sx q[2];
rz(0.86610438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56126334) q[1];
sx q[1];
rz(-0.65960192) q[1];
sx q[1];
rz(-2.329201) q[1];
rz(-pi) q[2];
rz(1.6555637) q[3];
sx q[3];
rz(-0.42094195) q[3];
sx q[3];
rz(-1.8404901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7647543) q[2];
sx q[2];
rz(-2.59616) q[2];
sx q[2];
rz(0.1405912) q[2];
rz(-1.8600672) q[3];
sx q[3];
rz(-1.8257273) q[3];
sx q[3];
rz(-0.52711058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.918688) q[0];
sx q[0];
rz(-1.0162901) q[0];
sx q[0];
rz(1.5284398) q[0];
rz(-0.81360045) q[1];
sx q[1];
rz(-0.10491144) q[1];
sx q[1];
rz(0.80748564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9514526) q[0];
sx q[0];
rz(-1.373111) q[0];
sx q[0];
rz(-2.9485767) q[0];
rz(-pi) q[1];
rz(-3.0119704) q[2];
sx q[2];
rz(-1.1309012) q[2];
sx q[2];
rz(0.18378809) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2827816) q[1];
sx q[1];
rz(-2.3859897) q[1];
sx q[1];
rz(1.0048466) q[1];
x q[2];
rz(-2.2178879) q[3];
sx q[3];
rz(-0.64416235) q[3];
sx q[3];
rz(-1.3180863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3537102) q[2];
sx q[2];
rz(-1.7870125) q[2];
sx q[2];
rz(2.9969969) q[2];
rz(2.0374129) q[3];
sx q[3];
rz(-2.927533) q[3];
sx q[3];
rz(2.1000699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5398194) q[0];
sx q[0];
rz(-2.0105392) q[0];
sx q[0];
rz(-2.5850776) q[0];
rz(-2.3717608) q[1];
sx q[1];
rz(-3.0172805) q[1];
sx q[1];
rz(-2.6803023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3614013) q[0];
sx q[0];
rz(-2.5690329) q[0];
sx q[0];
rz(1.6587692) q[0];
rz(-0.84750743) q[2];
sx q[2];
rz(-0.88107785) q[2];
sx q[2];
rz(0.72874069) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53992311) q[1];
sx q[1];
rz(-0.8569255) q[1];
sx q[1];
rz(-1.9636092) q[1];
rz(2.5784745) q[3];
sx q[3];
rz(-2.4909935) q[3];
sx q[3];
rz(1.4103149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7028659) q[2];
sx q[2];
rz(-0.30539572) q[2];
sx q[2];
rz(-2.3766282) q[2];
rz(-1.9276098) q[3];
sx q[3];
rz(-2.5524804) q[3];
sx q[3];
rz(0.37863076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42534378) q[0];
sx q[0];
rz(-0.049976293) q[0];
sx q[0];
rz(-0.36323994) q[0];
rz(-1.4092457) q[1];
sx q[1];
rz(-0.98594085) q[1];
sx q[1];
rz(-0.22663103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0186414) q[0];
sx q[0];
rz(-1.5737857) q[0];
sx q[0];
rz(-0.13521533) q[0];
x q[1];
rz(-2.2025467) q[2];
sx q[2];
rz(-0.9264731) q[2];
sx q[2];
rz(0.34560386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2001463) q[1];
sx q[1];
rz(-2.2392803) q[1];
sx q[1];
rz(0.41819765) q[1];
rz(2.8037386) q[3];
sx q[3];
rz(-1.5808592) q[3];
sx q[3];
rz(-2.4110766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5994485) q[2];
sx q[2];
rz(-0.52079529) q[2];
sx q[2];
rz(0.64860541) q[2];
rz(-0.058622807) q[3];
sx q[3];
rz(-2.5128745) q[3];
sx q[3];
rz(2.4962943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36515737) q[0];
sx q[0];
rz(-1.9226274) q[0];
sx q[0];
rz(2.6489039) q[0];
rz(-1.0538712) q[1];
sx q[1];
rz(-0.55640472) q[1];
sx q[1];
rz(0.41592204) q[1];
rz(1.8070167) q[2];
sx q[2];
rz(-1.5470355) q[2];
sx q[2];
rz(2.7337337) q[2];
rz(-0.33526996) q[3];
sx q[3];
rz(-0.90349586) q[3];
sx q[3];
rz(-0.60140453) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
