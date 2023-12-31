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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0331623) q[0];
sx q[0];
rz(-2.172029) q[0];
sx q[0];
rz(-0.45494672) q[0];
rz(3.0084228) q[2];
sx q[2];
rz(-2.3699017) q[2];
sx q[2];
rz(2.0193677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4359856) q[1];
sx q[1];
rz(-1.7090624) q[1];
sx q[1];
rz(-1.0857173) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0356785) q[3];
sx q[3];
rz(-2.7273791) q[3];
sx q[3];
rz(2.8416882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8841298) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.9467547) q[0];
rz(-0.082611235) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(3.1412178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1502991) q[0];
sx q[0];
rz(-2.5851558) q[0];
sx q[0];
rz(-2.2668905) q[0];
rz(-pi) q[1];
rz(-1.3090918) q[2];
sx q[2];
rz(-1.3771025) q[2];
sx q[2];
rz(-1.2145834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8723215) q[1];
sx q[1];
rz(-1.3355796) q[1];
sx q[1];
rz(2.5260731) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11081103) q[3];
sx q[3];
rz(-0.46289819) q[3];
sx q[3];
rz(-0.09679951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3327545) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(0.79408944) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394543) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(0.40476558) q[0];
rz(-1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(-1.3084897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8264016) q[0];
sx q[0];
rz(-1.0747402) q[0];
sx q[0];
rz(-2.1000923) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3345509) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(-0.46307785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96924671) q[1];
sx q[1];
rz(-1.8241276) q[1];
sx q[1];
rz(-2.2971056) q[1];
rz(-1.6609459) q[3];
sx q[3];
rz(-1.1336859) q[3];
sx q[3];
rz(0.79870943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(-3.1211839) q[2];
rz(2.0698047) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(-0.91598696) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(2.0844918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0650478) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(1.5426427) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1787203) q[2];
sx q[2];
rz(-2.0553556) q[2];
sx q[2];
rz(-1.0819266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.033356655) q[1];
sx q[1];
rz(-1.557096) q[1];
sx q[1];
rz(1.800888) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3841342) q[3];
sx q[3];
rz(-1.4043024) q[3];
sx q[3];
rz(1.43653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(2.7988953) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-0.87564898) q[3];
sx q[3];
rz(-0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.7856359) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(-0.086439565) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(-0.46868971) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5422573) q[0];
sx q[0];
rz(-0.89878966) q[0];
sx q[0];
rz(0.50962781) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8332963) q[2];
sx q[2];
rz(-1.7436276) q[2];
sx q[2];
rz(1.0002713) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75697452) q[1];
sx q[1];
rz(-0.24090919) q[1];
sx q[1];
rz(0.5730281) q[1];
x q[2];
rz(-0.81434806) q[3];
sx q[3];
rz(-1.6296248) q[3];
sx q[3];
rz(1.7562263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1420574) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(0.21720973) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(2.9445904) q[0];
rz(1.7794094) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(-0.33624712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67664458) q[0];
sx q[0];
rz(-1.1332382) q[0];
sx q[0];
rz(2.4462571) q[0];
rz(-pi) q[1];
rz(0.7746081) q[2];
sx q[2];
rz(-2.3850072) q[2];
sx q[2];
rz(2.2677383) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.058192249) q[1];
sx q[1];
rz(-0.5914878) q[1];
sx q[1];
rz(3.0565492) q[1];
rz(-1.7429966) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53529915) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(-2.2952648) q[2];
rz(1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3939312) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(-1.0940201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31539311) q[0];
sx q[0];
rz(-1.1356192) q[0];
sx q[0];
rz(0.3776334) q[0];
rz(-pi) q[1];
rz(-1.3482434) q[2];
sx q[2];
rz(-1.5828653) q[2];
sx q[2];
rz(1.0782858) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5309696) q[1];
sx q[1];
rz(-1.5986643) q[1];
sx q[1];
rz(-1.4762957) q[1];
rz(-pi) q[2];
rz(3.0435211) q[3];
sx q[3];
rz(-1.5925928) q[3];
sx q[3];
rz(1.9779713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.1743494) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-0.32067498) q[2];
rz(2.705412) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86673474) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(2.136769) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-0.90548038) q[1];
sx q[1];
rz(-0.80387962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.852254) q[0];
sx q[0];
rz(-0.58710557) q[0];
sx q[0];
rz(-0.071285204) q[0];
x q[1];
rz(2.8477746) q[2];
sx q[2];
rz(-0.57379469) q[2];
sx q[2];
rz(0.82905713) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4623973) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(1.4560844) q[1];
rz(0.6635267) q[3];
sx q[3];
rz(-1.6420206) q[3];
sx q[3];
rz(-2.0957259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(-1.0104898) q[2];
rz(-0.60339749) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-0.4075152) q[0];
rz(0.28911668) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(-2.3908652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84839612) q[0];
sx q[0];
rz(-1.2788532) q[0];
sx q[0];
rz(1.2743203) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8640792) q[2];
sx q[2];
rz(-2.8515184) q[2];
sx q[2];
rz(-1.8857423) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.99414413) q[1];
sx q[1];
rz(-1.7732883) q[1];
sx q[1];
rz(1.8840428) q[1];
rz(-2.8402746) q[3];
sx q[3];
rz(-1.5958061) q[3];
sx q[3];
rz(3.0400288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0728545) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(1.9753974) q[2];
rz(1.3646305) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5831379) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(1.7512084) q[0];
rz(-2.8109) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(1.4601382) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351785) q[0];
sx q[0];
rz(-1.9902475) q[0];
sx q[0];
rz(2.123453) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.023256217) q[2];
sx q[2];
rz(-2.0746982) q[2];
sx q[2];
rz(-1.1018673) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3223443) q[1];
sx q[1];
rz(-1.3867612) q[1];
sx q[1];
rz(-0.40477246) q[1];
rz(-pi) q[2];
rz(-0.71698935) q[3];
sx q[3];
rz(-0.90962142) q[3];
sx q[3];
rz(-0.17718525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.6798518) q[2];
rz(-1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6256975) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(1.3810146) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(-1.9258326) q[2];
sx q[2];
rz(-2.3264865) q[2];
sx q[2];
rz(-2.9441499) q[2];
rz(2.667726) q[3];
sx q[3];
rz(-0.73054536) q[3];
sx q[3];
rz(1.2252145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
