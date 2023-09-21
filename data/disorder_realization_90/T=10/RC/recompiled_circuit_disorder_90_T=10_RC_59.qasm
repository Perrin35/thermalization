OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4164299) q[0];
sx q[0];
rz(-0.13983146) q[0];
sx q[0];
rz(-2.5319985) q[0];
rz(-2.4729589) q[1];
sx q[1];
rz(-0.86548391) q[1];
sx q[1];
rz(-3.0545711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10843033) q[0];
sx q[0];
rz(-2.172029) q[0];
sx q[0];
rz(2.6866459) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1331698) q[2];
sx q[2];
rz(-0.77169092) q[2];
sx q[2];
rz(1.1222249) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.062292369) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(2.9855707) q[1];
rz(-0.10591412) q[3];
sx q[3];
rz(-0.4142136) q[3];
sx q[3];
rz(-0.29990444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13880754) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(-0.3368245) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.2574629) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.194838) q[0];
rz(3.0589814) q[1];
sx q[1];
rz(-1.1673085) q[1];
sx q[1];
rz(-3.1412178) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99129358) q[0];
sx q[0];
rz(-0.55643686) q[0];
sx q[0];
rz(-0.87470212) q[0];
rz(-2.2194901) q[2];
sx q[2];
rz(-0.32425913) q[2];
sx q[2];
rz(2.1622554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8723215) q[1];
sx q[1];
rz(-1.806013) q[1];
sx q[1];
rz(0.61551952) q[1];
rz(-pi) q[2];
rz(-1.5156636) q[3];
sx q[3];
rz(-2.030636) q[3];
sx q[3];
rz(-0.22051792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(-0.17671281) q[2];
rz(0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(-2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394543) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(2.7368271) q[0];
rz(-1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.3084897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015941042) q[0];
sx q[0];
rz(-1.1106655) q[0];
sx q[0];
rz(2.5815651) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3345509) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(2.6785148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1723459) q[1];
sx q[1];
rz(-1.8241276) q[1];
sx q[1];
rz(-2.2971056) q[1];
x q[2];
rz(-1.4806467) q[3];
sx q[3];
rz(-2.0079068) q[3];
sx q[3];
rz(-2.3428832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(0.02040872) q[2];
rz(2.0698047) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(-2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(2.2256057) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-1.0659734) q[1];
sx q[1];
rz(1.0571009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14295828) q[0];
sx q[0];
rz(-0.43804533) q[0];
sx q[0];
rz(-0.060158718) q[0];
rz(-0.57025036) q[2];
sx q[2];
rz(-1.0409365) q[2];
sx q[2];
rz(2.3392764) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6009439) q[1];
sx q[1];
rz(-1.800866) q[1];
sx q[1];
rz(3.1275216) q[1];
x q[2];
rz(-2.9017157) q[3];
sx q[3];
rz(-0.7719709) q[3];
sx q[3];
rz(-3.1023657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(0.34269732) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(-2.4263583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(-2.6729029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86779224) q[0];
sx q[0];
rz(-0.81875728) q[0];
sx q[0];
rz(1.0206945) q[0];
rz(-pi) q[1];
rz(-0.52207559) q[2];
sx q[2];
rz(-2.7895088) q[2];
sx q[2];
rz(-0.075369518) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9710755) q[1];
sx q[1];
rz(-1.7726388) q[1];
sx q[1];
rz(1.7032196) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81434806) q[3];
sx q[3];
rz(-1.6296248) q[3];
sx q[3];
rz(-1.7562263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9995352) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7140759) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(-2.9445904) q[0];
rz(1.7794094) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(2.8053455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67664458) q[0];
sx q[0];
rz(-2.0083545) q[0];
sx q[0];
rz(2.4462571) q[0];
rz(-0.7746081) q[2];
sx q[2];
rz(-0.75658549) q[2];
sx q[2];
rz(2.2677383) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.058192249) q[1];
sx q[1];
rz(-0.5914878) q[1];
sx q[1];
rz(3.0565492) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9365262) q[3];
sx q[3];
rz(-2.4347217) q[3];
sx q[3];
rz(-0.52458483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(-1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-0.0011750778) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7476615) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.7566453) q[1];
sx q[1];
rz(1.0940201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8261995) q[0];
sx q[0];
rz(-2.0059735) q[0];
sx q[0];
rz(0.3776334) q[0];
rz(-pi) q[1];
rz(1.6254243) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(-0.54578997) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5309696) q[1];
sx q[1];
rz(-1.5986643) q[1];
sx q[1];
rz(-1.665297) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.592698) q[3];
sx q[3];
rz(-1.6688445) q[3];
sx q[3];
rz(2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(0.32067498) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86673474) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(-0.59016219) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(-2.337713) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3748976) q[0];
sx q[0];
rz(-2.156213) q[0];
sx q[0];
rz(-1.6181519) q[0];
x q[1];
rz(2.5875823) q[2];
sx q[2];
rz(-1.7286574) q[2];
sx q[2];
rz(-2.6487034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67919532) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(1.4560844) q[1];
rz(0.11532468) q[3];
sx q[3];
rz(-0.66676312) q[3];
sx q[3];
rz(-2.525884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-1.0104898) q[2];
rz(-2.5381952) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(-0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-2.7340775) q[0];
rz(-2.852476) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(-2.3908652) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4780316) q[0];
sx q[0];
rz(-2.7285517) q[0];
sx q[0];
rz(-2.3703299) q[0];
rz(-pi) q[1];
rz(0.27955555) q[2];
sx q[2];
rz(-1.6492372) q[2];
sx q[2];
rz(-2.560175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1474485) q[1];
sx q[1];
rz(-1.3683043) q[1];
sx q[1];
rz(1.2575498) q[1];
rz(-0.30131807) q[3];
sx q[3];
rz(-1.5958061) q[3];
sx q[3];
rz(0.1015639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0687381) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(-1.1661952) q[2];
rz(1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5831379) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(1.7512084) q[0];
rz(0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(1.4601382) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351785) q[0];
sx q[0];
rz(-1.1513452) q[0];
sx q[0];
rz(1.0181396) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0748126) q[2];
sx q[2];
rz(-1.5911615) q[2];
sx q[2];
rz(2.6838944) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9863661) q[1];
sx q[1];
rz(-0.44253293) q[1];
sx q[1];
rz(-2.7000484) q[1];
x q[2];
rz(-0.7695997) q[3];
sx q[3];
rz(-2.1160612) q[3];
sx q[3];
rz(-1.8850957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(-1.6798518) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6256975) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(-1.760578) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(2.7881651) q[2];
sx q[2];
rz(-0.8197439) q[2];
sx q[2];
rz(2.843429) q[2];
rz(-2.4685728) q[3];
sx q[3];
rz(-1.8802079) q[3];
sx q[3];
rz(0.019284266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
