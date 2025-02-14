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
rz(1.6788586) q[1];
sx q[1];
rz(0.94574133) q[1];
sx q[1];
rz(8.9206817) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5228793) q[0];
sx q[0];
rz(-2.0815432) q[0];
sx q[0];
rz(0.63453959) q[0];
rz(-pi) q[1];
rz(1.7376704) q[2];
sx q[2];
rz(-1.2271527) q[2];
sx q[2];
rz(2.0640896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2368184) q[1];
sx q[1];
rz(-2.1874551) q[1];
sx q[1];
rz(-0.84994933) q[1];
rz(-1.5220187) q[3];
sx q[3];
rz(-1.2360473) q[3];
sx q[3];
rz(-2.3693645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41245875) q[2];
sx q[2];
rz(-1.8987741) q[2];
sx q[2];
rz(-2.9553555) q[2];
rz(1.3650182) q[3];
sx q[3];
rz(-0.61187196) q[3];
sx q[3];
rz(-1.2046825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314826) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(-3.0895184) q[0];
rz(-1.0366108) q[1];
sx q[1];
rz(-2.5317445) q[1];
sx q[1];
rz(-0.61526543) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66917574) q[0];
sx q[0];
rz(-2.376498) q[0];
sx q[0];
rz(1.0932176) q[0];
rz(-0.5004761) q[2];
sx q[2];
rz(-2.8722974) q[2];
sx q[2];
rz(-0.057502086) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4897291) q[1];
sx q[1];
rz(-2.4861838) q[1];
sx q[1];
rz(-1.4784527) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9297057) q[3];
sx q[3];
rz(-2.35244) q[3];
sx q[3];
rz(1.0494389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.692824) q[2];
sx q[2];
rz(-1.8961467) q[2];
sx q[2];
rz(1.9057062) q[2];
rz(-1.1290733) q[3];
sx q[3];
rz(-2.2344507) q[3];
sx q[3];
rz(0.83704078) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4835943) q[0];
sx q[0];
rz(-0.75355419) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(-0.54028571) q[1];
sx q[1];
rz(-1.0397725) q[1];
sx q[1];
rz(-1.6859863) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7175563) q[0];
sx q[0];
rz(-0.40923318) q[0];
sx q[0];
rz(2.9633029) q[0];
rz(-pi) q[1];
rz(-1.6605145) q[2];
sx q[2];
rz(-1.3680581) q[2];
sx q[2];
rz(2.3281472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7109622) q[1];
sx q[1];
rz(-2.4050131) q[1];
sx q[1];
rz(1.0308835) q[1];
x q[2];
rz(-2.2690587) q[3];
sx q[3];
rz(-1.8909406) q[3];
sx q[3];
rz(2.9001989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6471863) q[2];
sx q[2];
rz(-2.8072085) q[2];
sx q[2];
rz(1.4462659) q[2];
rz(-1.9485731) q[3];
sx q[3];
rz(-0.97508591) q[3];
sx q[3];
rz(1.4421991) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0541662) q[0];
sx q[0];
rz(-0.26987258) q[0];
sx q[0];
rz(2.7179981) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.3514163) q[1];
sx q[1];
rz(-0.097600309) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711084) q[0];
sx q[0];
rz(-1.1997265) q[0];
sx q[0];
rz(-2.5778997) q[0];
rz(1.7427069) q[2];
sx q[2];
rz(-2.4414276) q[2];
sx q[2];
rz(-1.3106048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3279539) q[1];
sx q[1];
rz(-2.6662152) q[1];
sx q[1];
rz(-0.83303501) q[1];
rz(-pi) q[2];
rz(1.2318491) q[3];
sx q[3];
rz(-0.87050658) q[3];
sx q[3];
rz(0.46745121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5264954) q[2];
sx q[2];
rz(-2.6690833) q[2];
sx q[2];
rz(2.9769767) q[2];
rz(-3.0937255) q[3];
sx q[3];
rz(-1.7367626) q[3];
sx q[3];
rz(-2.3544748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614302) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(-1.5572146) q[0];
rz(2.7730675) q[1];
sx q[1];
rz(-1.3216113) q[1];
sx q[1];
rz(-1.6065067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212352) q[0];
sx q[0];
rz(-2.0518655) q[0];
sx q[0];
rz(-2.1522983) q[0];
rz(-2.755411) q[2];
sx q[2];
rz(-0.86241041) q[2];
sx q[2];
rz(-2.2798373) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8884224) q[1];
sx q[1];
rz(-2.0520323) q[1];
sx q[1];
rz(-0.28216823) q[1];
x q[2];
rz(-0.078952958) q[3];
sx q[3];
rz(-1.8584849) q[3];
sx q[3];
rz(-2.1514078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9102455) q[2];
sx q[2];
rz(-0.37800899) q[2];
sx q[2];
rz(1.0450012) q[2];
rz(-2.9735978) q[3];
sx q[3];
rz(-1.1863656) q[3];
sx q[3];
rz(-0.35429889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.3280227) q[1];
sx q[1];
rz(0.71570754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3580644) q[0];
sx q[0];
rz(-0.77371374) q[0];
sx q[0];
rz(0.85229243) q[0];
rz(0.44797051) q[2];
sx q[2];
rz(-1.7501737) q[2];
sx q[2];
rz(-0.57283869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0790076) q[1];
sx q[1];
rz(-0.91168303) q[1];
sx q[1];
rz(-1.4277894) q[1];
rz(-pi) q[2];
rz(2.8479667) q[3];
sx q[3];
rz(-2.01247) q[3];
sx q[3];
rz(-2.3229118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7314926) q[2];
sx q[2];
rz(-2.4918753) q[2];
sx q[2];
rz(-2.0984207) q[2];
rz(0.10797524) q[3];
sx q[3];
rz(-0.94111809) q[3];
sx q[3];
rz(2.030355) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1825948) q[0];
sx q[0];
rz(-1.3866871) q[0];
sx q[0];
rz(1.2249655) q[0];
rz(0.23172465) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(-1.8909594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38076048) q[0];
sx q[0];
rz(-2.3835813) q[0];
sx q[0];
rz(0.17802989) q[0];
x q[1];
rz(-1.3725946) q[2];
sx q[2];
rz(-1.5491252) q[2];
sx q[2];
rz(0.54474165) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9385443) q[1];
sx q[1];
rz(-0.37751679) q[1];
sx q[1];
rz(0.63407268) q[1];
x q[2];
rz(0.55646236) q[3];
sx q[3];
rz(-1.5817405) q[3];
sx q[3];
rz(1.3370322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0687678) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(2.4459548) q[3];
sx q[3];
rz(-1.9591103) q[3];
sx q[3];
rz(-0.80719358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.950133) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(2.4323442) q[0];
rz(-1.5059772) q[1];
sx q[1];
rz(-1.9874856) q[1];
sx q[1];
rz(-2.0567315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4385972) q[0];
sx q[0];
rz(-0.60351047) q[0];
sx q[0];
rz(-3.0402157) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9726344) q[2];
sx q[2];
rz(-0.8118793) q[2];
sx q[2];
rz(0.66056992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0811396) q[1];
sx q[1];
rz(-0.53489242) q[1];
sx q[1];
rz(-0.51523955) q[1];
rz(1.510511) q[3];
sx q[3];
rz(-1.5593646) q[3];
sx q[3];
rz(-0.41939467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.52013493) q[2];
sx q[2];
rz(-2.1389213) q[2];
sx q[2];
rz(1.3168859) q[2];
rz(-0.78553158) q[3];
sx q[3];
rz(-1.1979878) q[3];
sx q[3];
rz(0.42640105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.921628) q[0];
sx q[0];
rz(-1.4861318) q[0];
sx q[0];
rz(0.40837902) q[0];
rz(0.95868239) q[1];
sx q[1];
rz(-0.32902333) q[1];
sx q[1];
rz(1.6201409) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7470541) q[0];
sx q[0];
rz(-1.9346721) q[0];
sx q[0];
rz(2.8900212) q[0];
rz(1.6198115) q[2];
sx q[2];
rz(-1.00178) q[2];
sx q[2];
rz(-2.4328277) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4056333) q[1];
sx q[1];
rz(-0.93184566) q[1];
sx q[1];
rz(-0.51764857) q[1];
x q[2];
rz(-1.0318087) q[3];
sx q[3];
rz(-2.781311) q[3];
sx q[3];
rz(-2.173645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.3756322) q[2];
sx q[2];
rz(-1.1124632) q[2];
sx q[2];
rz(0.56606236) q[2];
rz(-1.3426956) q[3];
sx q[3];
rz(-2.7261966) q[3];
sx q[3];
rz(2.8528163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.600243) q[0];
sx q[0];
rz(-3.1025649) q[0];
sx q[0];
rz(1.716123) q[0];
rz(-1.0936945) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(-2.7519382) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.728316) q[0];
sx q[0];
rz(-1.1470058) q[0];
sx q[0];
rz(-2.8601147) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94691838) q[2];
sx q[2];
rz(-1.7830666) q[2];
sx q[2];
rz(-0.52516261) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.13481465) q[1];
sx q[1];
rz(-1.14535) q[1];
sx q[1];
rz(-3.1158034) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6171754) q[3];
sx q[3];
rz(-0.56045415) q[3];
sx q[3];
rz(-0.84597142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23715544) q[2];
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
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.420153) q[0];
sx q[0];
rz(-2.1099821) q[0];
sx q[0];
rz(-1.4954062) q[0];
rz(-2.738476) q[1];
sx q[1];
rz(-1.8164201) q[1];
sx q[1];
rz(1.4774189) q[1];
rz(1.6942799) q[2];
sx q[2];
rz(-0.99939838) q[2];
sx q[2];
rz(1.8746092) q[2];
rz(2.7561989) q[3];
sx q[3];
rz(-1.7447532) q[3];
sx q[3];
rz(0.90739653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
