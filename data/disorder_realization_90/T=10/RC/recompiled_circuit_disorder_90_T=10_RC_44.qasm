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
rz(-2.4729589) q[1];
sx q[1];
rz(-0.86548391) q[1];
sx q[1];
rz(-3.0545711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60458175) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(-1.0010615) q[0];
rz(-3.0084228) q[2];
sx q[2];
rz(-2.3699017) q[2];
sx q[2];
rz(1.1222249) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0793003) q[1];
sx q[1];
rz(-1.0907409) q[1];
sx q[1];
rz(0.15602195) q[1];
x q[2];
rz(1.617241) q[3];
sx q[3];
rz(-1.159045) q[3];
sx q[3];
rz(2.9573033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0027851) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(2.9132304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(-1.9467547) q[0];
rz(3.0589814) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(-0.00037489051) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99129358) q[0];
sx q[0];
rz(-0.55643686) q[0];
sx q[0];
rz(-2.2668905) q[0];
x q[1];
rz(1.3090918) q[2];
sx q[2];
rz(-1.7644901) q[2];
sx q[2];
rz(-1.2145834) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0034069) q[1];
sx q[1];
rz(-2.1669743) q[1];
sx q[1];
rz(1.8562993) q[1];
x q[2];
rz(0.46044465) q[3];
sx q[3];
rz(-1.5213955) q[3];
sx q[3];
rz(1.3747665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3327545) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(0.17671281) q[2];
rz(-2.3475032) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(-2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394543) q[0];
sx q[0];
rz(-0.70394009) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.8331029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93861261) q[0];
sx q[0];
rz(-2.4327607) q[0];
sx q[0];
rz(-0.75074408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3345509) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(2.6785148) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3209826) q[1];
sx q[1];
rz(-2.2690986) q[1];
sx q[1];
rz(-0.33336158) q[1];
rz(-1.4806467) q[3];
sx q[3];
rz(-2.0079068) q[3];
sx q[3];
rz(0.79870943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(-1.071788) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(-0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9999009) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14295828) q[0];
sx q[0];
rz(-0.43804533) q[0];
sx q[0];
rz(-0.060158718) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96287231) q[2];
sx q[2];
rz(-1.086237) q[2];
sx q[2];
rz(-1.0819266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4790198) q[1];
sx q[1];
rz(-2.9111007) q[1];
sx q[1];
rz(1.5107933) q[1];
rz(-pi) q[2];
rz(-0.75745849) q[3];
sx q[3];
rz(-1.7372903) q[3];
sx q[3];
rz(-1.7050626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2477734) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(-0.34269732) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(-2.4263583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(-0.086439565) q[0];
rz(-1.7516288) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(0.46868971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3062895) q[0];
sx q[0];
rz(-1.9625184) q[0];
sx q[0];
rz(0.83175559) q[0];
rz(-0.52207559) q[2];
sx q[2];
rz(-0.35208382) q[2];
sx q[2];
rz(-3.0662231) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9710755) q[1];
sx q[1];
rz(-1.3689539) q[1];
sx q[1];
rz(1.7032196) q[1];
rz(-pi) q[2];
rz(3.0607871) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9995352) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-0.86432499) q[2];
rz(2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(0.19700225) q[0];
rz(1.7794094) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(0.33624712) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2339736) q[0];
sx q[0];
rz(-2.189878) q[0];
sx q[0];
rz(-2.1179849) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3669846) q[2];
sx q[2];
rz(-2.3850072) q[2];
sx q[2];
rz(-0.87385439) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5583401) q[1];
sx q[1];
rz(-1.5234158) q[1];
sx q[1];
rz(2.5517795) q[1];
rz(2.9365262) q[3];
sx q[3];
rz(-0.70687095) q[3];
sx q[3];
rz(-0.52458483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.53529915) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(0.84632787) q[2];
rz(1.9019295) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7476615) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(0.25892648) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(2.0475725) q[1];
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
rz(-2.7639593) q[0];
rz(-pi) q[1];
rz(-1.5161683) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(-0.54578997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5309696) q[1];
sx q[1];
rz(-1.5986643) q[1];
sx q[1];
rz(1.665297) q[1];
x q[2];
rz(3.0435211) q[3];
sx q[3];
rz(-1.5925928) q[3];
sx q[3];
rz(1.9779713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9672433) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(2.8209177) q[2];
rz(-0.43618068) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.86673474) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(-1.0048237) q[0];
rz(-0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(-0.80387962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893387) q[0];
sx q[0];
rz(-0.58710557) q[0];
sx q[0];
rz(3.0703074) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5875823) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(-2.6487034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4623973) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(-1.4560844) q[1];
rz(1.4804833) q[3];
sx q[3];
rz(-0.90925018) q[3];
sx q[3];
rz(0.46935287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(1.0104898) q[2];
rz(0.60339749) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(-0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-0.4075152) q[0];
rz(2.852476) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-2.3908652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2931965) q[0];
sx q[0];
rz(-1.8627394) q[0];
sx q[0];
rz(-1.2743203) q[0];
x q[1];
rz(-0.27751343) q[2];
sx q[2];
rz(-0.29007426) q[2];
sx q[2];
rz(1.8857423) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99414413) q[1];
sx q[1];
rz(-1.3683043) q[1];
sx q[1];
rz(-1.2575498) q[1];
x q[2];
rz(-1.5446072) q[3];
sx q[3];
rz(-1.8720172) q[3];
sx q[3];
rz(1.6645886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0687381) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(1.1661952) q[2];
rz(-1.3646305) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.5831379) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(-1.7512084) q[0];
rz(2.8109) q[1];
sx q[1];
rz(-0.76534098) q[1];
sx q[1];
rz(1.4601382) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64797348) q[0];
sx q[0];
rz(-2.4612392) q[0];
sx q[0];
rz(0.86662678) q[0];
rz(-pi) q[1];
rz(3.1183364) q[2];
sx q[2];
rz(-1.0668945) q[2];
sx q[2];
rz(-2.0397253) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8148988) q[1];
sx q[1];
rz(-1.1732475) q[1];
sx q[1];
rz(1.7705998) q[1];
rz(-pi) q[2];
rz(-0.86942418) q[3];
sx q[3];
rz(-0.93360177) q[3];
sx q[3];
rz(-2.3616392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7918487) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.4617408) q[2];
rz(-1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5158952) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(1.760578) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(-1.9258326) q[2];
sx q[2];
rz(-2.3264865) q[2];
sx q[2];
rz(-2.9441499) q[2];
rz(-1.1827042) q[3];
sx q[3];
rz(-0.93508616) q[3];
sx q[3];
rz(-1.3133776) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
