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
rz(-0.47973862) q[0];
sx q[0];
rz(-2.0070183) q[0];
sx q[0];
rz(-1.467508) q[0];
rz(-1.462734) q[1];
sx q[1];
rz(-0.94574133) q[1];
sx q[1];
rz(2.6374964) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60668463) q[0];
sx q[0];
rz(-2.1143171) q[0];
sx q[0];
rz(-2.1786819) q[0];
x q[1];
rz(1.7376704) q[2];
sx q[2];
rz(-1.9144399) q[2];
sx q[2];
rz(-2.0640896) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2368184) q[1];
sx q[1];
rz(-2.1874551) q[1];
sx q[1];
rz(-0.84994933) q[1];
rz(-pi) q[2];
rz(1.5220187) q[3];
sx q[3];
rz(-1.2360473) q[3];
sx q[3];
rz(-0.77222811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7291339) q[2];
sx q[2];
rz(-1.8987741) q[2];
sx q[2];
rz(-2.9553555) q[2];
rz(1.7765744) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(1.9369102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.81011009) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(-0.052074281) q[0];
rz(-1.0366108) q[1];
sx q[1];
rz(-2.5317445) q[1];
sx q[1];
rz(2.5263272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66917574) q[0];
sx q[0];
rz(-0.76509464) q[0];
sx q[0];
rz(1.0932176) q[0];
rz(2.6411166) q[2];
sx q[2];
rz(-0.26929528) q[2];
sx q[2];
rz(-3.0840906) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76814684) q[1];
sx q[1];
rz(-2.2229338) q[1];
sx q[1];
rz(-3.0708205) q[1];
rz(-pi) q[2];
rz(1.3619945) q[3];
sx q[3];
rz(-0.80392814) q[3];
sx q[3];
rz(-1.3458136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.44876862) q[2];
sx q[2];
rz(-1.8961467) q[2];
sx q[2];
rz(1.2358865) q[2];
rz(-2.0125194) q[3];
sx q[3];
rz(-0.90714199) q[3];
sx q[3];
rz(-2.3045519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6579984) q[0];
sx q[0];
rz(-0.75355419) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(-2.6013069) q[1];
sx q[1];
rz(-2.1018201) q[1];
sx q[1];
rz(1.4556063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7175563) q[0];
sx q[0];
rz(-0.40923318) q[0];
sx q[0];
rz(-0.17828973) q[0];
rz(-pi) q[1];
rz(2.7305538) q[2];
sx q[2];
rz(-0.22145311) q[2];
sx q[2];
rz(-2.7483181) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1108577) q[1];
sx q[1];
rz(-0.95672119) q[1];
sx q[1];
rz(-2.7053737) q[1];
rz(-pi) q[2];
rz(1.0946526) q[3];
sx q[3];
rz(-0.75679251) q[3];
sx q[3];
rz(-1.688129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49440631) q[2];
sx q[2];
rz(-0.33438412) q[2];
sx q[2];
rz(1.4462659) q[2];
rz(1.9485731) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0541662) q[0];
sx q[0];
rz(-2.8717201) q[0];
sx q[0];
rz(2.7179981) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.7901763) q[1];
sx q[1];
rz(0.097600309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0711084) q[0];
sx q[0];
rz(-1.1997265) q[0];
sx q[0];
rz(2.5778997) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3988858) q[2];
sx q[2];
rz(-0.70016501) q[2];
sx q[2];
rz(1.3106048) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1244036) q[1];
sx q[1];
rz(-1.9162971) q[1];
sx q[1];
rz(0.33331897) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9097435) q[3];
sx q[3];
rz(-2.2710861) q[3];
sx q[3];
rz(-0.46745121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61509722) q[2];
sx q[2];
rz(-2.6690833) q[2];
sx q[2];
rz(0.16461593) q[2];
rz(-0.047867157) q[3];
sx q[3];
rz(-1.7367626) q[3];
sx q[3];
rz(-0.78711787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.98016244) q[0];
sx q[0];
rz(-2.6949368) q[0];
sx q[0];
rz(-1.5572146) q[0];
rz(-0.36852512) q[1];
sx q[1];
rz(-1.3216113) q[1];
sx q[1];
rz(-1.6065067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53705789) q[0];
sx q[0];
rz(-0.7365444) q[0];
sx q[0];
rz(-0.81088541) q[0];
x q[1];
rz(2.755411) q[2];
sx q[2];
rz(-2.2791822) q[2];
sx q[2];
rz(-2.2798373) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6905744) q[1];
sx q[1];
rz(-1.8201882) q[1];
sx q[1];
rz(-2.0687769) q[1];
rz(-pi) q[2];
rz(-3.0626397) q[3];
sx q[3];
rz(-1.8584849) q[3];
sx q[3];
rz(2.1514078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9102455) q[2];
sx q[2];
rz(-2.7635837) q[2];
sx q[2];
rz(-1.0450012) q[2];
rz(-0.16799489) q[3];
sx q[3];
rz(-1.1863656) q[3];
sx q[3];
rz(0.35429889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198233) q[0];
sx q[0];
rz(-1.0291809) q[0];
sx q[0];
rz(-1.9371012) q[0];
rz(-0.80351859) q[1];
sx q[1];
rz(-0.81356994) q[1];
sx q[1];
rz(-0.71570754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2432118) q[0];
sx q[0];
rz(-1.0168494) q[0];
sx q[0];
rz(2.5701163) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44797051) q[2];
sx q[2];
rz(-1.7501737) q[2];
sx q[2];
rz(0.57283869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0625851) q[1];
sx q[1];
rz(-2.2299096) q[1];
sx q[1];
rz(-1.4277894) q[1];
rz(-2.8479667) q[3];
sx q[3];
rz(-1.1291227) q[3];
sx q[3];
rz(-2.3229118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41010007) q[2];
sx q[2];
rz(-0.64971739) q[2];
sx q[2];
rz(-2.0984207) q[2];
rz(-3.0336174) q[3];
sx q[3];
rz(-2.2004746) q[3];
sx q[3];
rz(-2.030355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9589979) q[0];
sx q[0];
rz(-1.3866871) q[0];
sx q[0];
rz(-1.9166272) q[0];
rz(-0.23172465) q[1];
sx q[1];
rz(-0.8871612) q[1];
sx q[1];
rz(-1.8909594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7608322) q[0];
sx q[0];
rz(-2.3835813) q[0];
sx q[0];
rz(-2.9635628) q[0];
x q[1];
rz(-1.3725946) q[2];
sx q[2];
rz(-1.5924675) q[2];
sx q[2];
rz(2.596851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9385443) q[1];
sx q[1];
rz(-2.7640759) q[1];
sx q[1];
rz(-2.50752) q[1];
rz(-2.5851303) q[3];
sx q[3];
rz(-1.5598522) q[3];
sx q[3];
rz(1.8045604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0687678) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(-3.0094299) q[2];
rz(0.69563785) q[3];
sx q[3];
rz(-1.9591103) q[3];
sx q[3];
rz(-2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19145963) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(2.4323442) q[0];
rz(-1.5059772) q[1];
sx q[1];
rz(-1.154107) q[1];
sx q[1];
rz(2.0567315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1902311) q[0];
sx q[0];
rz(-1.6282646) q[0];
sx q[0];
rz(-2.5404853) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3368767) q[2];
sx q[2];
rz(-1.448481) q[2];
sx q[2];
rz(1.027077) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0811396) q[1];
sx q[1];
rz(-0.53489242) q[1];
sx q[1];
rz(-0.51523955) q[1];
rz(1.3832757) q[3];
sx q[3];
rz(-0.061358364) q[3];
sx q[3];
rz(0.96422577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52013493) q[2];
sx q[2];
rz(-2.1389213) q[2];
sx q[2];
rz(1.8247068) q[2];
rz(-0.78553158) q[3];
sx q[3];
rz(-1.1979878) q[3];
sx q[3];
rz(0.42640105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.921628) q[0];
sx q[0];
rz(-1.4861318) q[0];
sx q[0];
rz(-2.7332136) q[0];
rz(2.1829103) q[1];
sx q[1];
rz(-0.32902333) q[1];
sx q[1];
rz(1.5214517) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1216089) q[0];
sx q[0];
rz(-0.43916288) q[0];
sx q[0];
rz(-0.99187054) q[0];
rz(-2.5720308) q[2];
sx q[2];
rz(-1.6120835) q[2];
sx q[2];
rz(2.2531367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7359594) q[1];
sx q[1];
rz(-2.209747) q[1];
sx q[1];
rz(2.6239441) q[1];
x q[2];
rz(1.0318087) q[3];
sx q[3];
rz(-2.781311) q[3];
sx q[3];
rz(2.173645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.3756322) q[2];
sx q[2];
rz(-2.0291294) q[2];
sx q[2];
rz(0.56606236) q[2];
rz(1.3426956) q[3];
sx q[3];
rz(-2.7261966) q[3];
sx q[3];
rz(0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.600243) q[0];
sx q[0];
rz(-0.039027795) q[0];
sx q[0];
rz(-1.4254697) q[0];
rz(1.0936945) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(-0.38965449) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3411639) q[0];
sx q[0];
rz(-0.50398705) q[0];
sx q[0];
rz(-1.0188854) q[0];
rz(-pi) q[1];
rz(-0.25954982) q[2];
sx q[2];
rz(-0.96298157) q[2];
sx q[2];
rz(1.9454625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.006778) q[1];
sx q[1];
rz(-1.14535) q[1];
sx q[1];
rz(-0.025789217) q[1];
rz(-pi) q[2];
rz(0.49764244) q[3];
sx q[3];
rz(-1.3013869) q[3];
sx q[3];
rz(0.26925081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9044372) q[2];
sx q[2];
rz(-2.6025786) q[2];
sx q[2];
rz(-1.944444) q[2];
rz(-2.9947179) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(2.1541514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.420153) q[0];
sx q[0];
rz(-2.1099821) q[0];
sx q[0];
rz(-1.4954062) q[0];
rz(2.738476) q[1];
sx q[1];
rz(-1.3251726) q[1];
sx q[1];
rz(-1.6641738) q[1];
rz(-2.5667122) q[2];
sx q[2];
rz(-1.6745865) q[2];
sx q[2];
rz(-2.9048017) q[2];
rz(1.758214) q[3];
sx q[3];
rz(-1.9500749) q[3];
sx q[3];
rz(-0.73350026) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
