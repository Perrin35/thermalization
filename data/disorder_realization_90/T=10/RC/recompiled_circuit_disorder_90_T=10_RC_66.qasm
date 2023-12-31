OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(-0.60959417) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(-0.087021526) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10843033) q[0];
sx q[0];
rz(-0.96956367) q[0];
sx q[0];
rz(-2.6866459) q[0];
x q[1];
rz(-1.6992703) q[2];
sx q[2];
rz(-2.3339084) q[2];
sx q[2];
rz(1.3070004) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39057186) q[1];
sx q[1];
rz(-2.6387072) q[1];
sx q[1];
rz(1.8608171) q[1];
rz(-pi) q[2];
rz(2.7294455) q[3];
sx q[3];
rz(-1.5282359) q[3];
sx q[3];
rz(1.3679078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0027851) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(-2.7677317) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(1.8841298) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.9467547) q[0];
rz(0.082611235) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(0.00037489051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9276792) q[0];
sx q[0];
rz(-1.1535026) q[0];
sx q[0];
rz(-0.37950619) q[0];
rz(-0.20034321) q[2];
sx q[2];
rz(-1.8274954) q[2];
sx q[2];
rz(-0.30470195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5217168) q[1];
sx q[1];
rz(-2.4881425) q[1];
sx q[1];
rz(0.39342777) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46044465) q[3];
sx q[3];
rz(-1.5213955) q[3];
sx q[3];
rz(-1.3747665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3327545) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(2.3475032) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(-2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394543) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(-1.2813214) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(-1.8331029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015941042) q[0];
sx q[0];
rz(-1.1106655) q[0];
sx q[0];
rz(0.56002754) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53738014) q[2];
sx q[2];
rz(-1.3668622) q[2];
sx q[2];
rz(2.1539719) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8150427) q[1];
sx q[1];
rz(-0.76154852) q[1];
sx q[1];
rz(1.9425068) q[1];
rz(-1.6609459) q[3];
sx q[3];
rz(-1.1336859) q[3];
sx q[3];
rz(0.79870943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9222766) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(2.0698047) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(2.2256057) q[0];
rz(2.6782716) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(1.0571009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6592641) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(-2.7042424) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1787203) q[2];
sx q[2];
rz(-2.0553556) q[2];
sx q[2];
rz(1.0819266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6625729) q[1];
sx q[1];
rz(-2.9111007) q[1];
sx q[1];
rz(-1.6307994) q[1];
rz(1.7980868) q[3];
sx q[3];
rz(-0.82633457) q[3];
sx q[3];
rz(-2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(1.6977067) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(-0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-0.086439565) q[0];
rz(-1.7516288) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(-2.6729029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8353032) q[0];
sx q[0];
rz(-1.9625184) q[0];
sx q[0];
rz(-2.3098371) q[0];
rz(-0.30829633) q[2];
sx q[2];
rz(-1.7436276) q[2];
sx q[2];
rz(-2.1413213) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9710755) q[1];
sx q[1];
rz(-1.7726388) q[1];
sx q[1];
rz(1.7032196) q[1];
rz(-1.4851941) q[3];
sx q[3];
rz(-2.3833131) q[3];
sx q[3];
rz(-3.0183834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1420574) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(-2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(-2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(-2.9445904) q[0];
rz(1.7794094) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(-2.8053455) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4649481) q[0];
sx q[0];
rz(-1.1332382) q[0];
sx q[0];
rz(-0.69533555) q[0];
rz(-pi) q[1];
rz(2.5480812) q[2];
sx q[2];
rz(-1.0700018) q[2];
sx q[2];
rz(1.3154495) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.058192249) q[1];
sx q[1];
rz(-2.5501049) q[1];
sx q[1];
rz(3.0565492) q[1];
x q[2];
rz(2.9365262) q[3];
sx q[3];
rz(-2.4347217) q[3];
sx q[3];
rz(0.52458483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53529915) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(-0.25892648) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-1.0940201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0703697) q[0];
sx q[0];
rz(-2.5734512) q[0];
sx q[0];
rz(-0.90026654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6254243) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(0.54578997) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.5309696) q[1];
sx q[1];
rz(-1.5986643) q[1];
sx q[1];
rz(-1.4762957) q[1];
rz(-pi) q[2];
rz(-2.9225227) q[3];
sx q[3];
rz(-0.1004569) q[3];
sx q[3];
rz(0.62517525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9672433) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(-0.32067498) q[2];
rz(2.705412) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2748579) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(1.0048237) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(-2.337713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.852254) q[0];
sx q[0];
rz(-2.5544871) q[0];
sx q[0];
rz(-3.0703074) q[0];
rz(0.55401037) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(0.49288921) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0058606) q[1];
sx q[1];
rz(-1.5809959) q[1];
sx q[1];
rz(-1.482153) q[1];
x q[2];
rz(-0.11532468) q[3];
sx q[3];
rz(-2.4748295) q[3];
sx q[3];
rz(-2.525884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(1.0104898) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(-2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5378961) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-2.7340775) q[0];
rz(2.852476) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(2.3908652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84839612) q[0];
sx q[0];
rz(-1.2788532) q[0];
sx q[0];
rz(-1.8672724) q[0];
rz(-2.8620371) q[2];
sx q[2];
rz(-1.6492372) q[2];
sx q[2];
rz(-2.560175) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0093065) q[1];
sx q[1];
rz(-2.7704151) q[1];
sx q[1];
rz(-2.1585141) q[1];
rz(0.084089355) q[3];
sx q[3];
rz(-0.30232271) q[3];
sx q[3];
rz(-1.38894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0728545) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-3.1197746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5584548) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(1.7512084) q[0];
rz(0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.6814544) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18171039) q[0];
sx q[0];
rz(-1.0707756) q[0];
sx q[0];
rz(0.48258968) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5286469) q[2];
sx q[2];
rz(-2.6372006) q[2];
sx q[2];
rz(1.1500037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9863661) q[1];
sx q[1];
rz(-2.6990597) q[1];
sx q[1];
rz(-0.44154422) q[1];
rz(-pi) q[2];
rz(-2.4246033) q[3];
sx q[3];
rz(-2.2319712) q[3];
sx q[3];
rz(2.9644074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7918487) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.4617408) q[2];
rz(2.0215624) q[3];
sx q[3];
rz(-2.5199065) q[3];
sx q[3];
rz(0.45564836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.35342758) q[2];
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
