OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37501332) q[0];
sx q[0];
rz(4.3809173) q[0];
sx q[0];
rz(8.3349174) q[0];
rz(-2.0774948) q[1];
sx q[1];
rz(-1.2531333) q[1];
sx q[1];
rz(-0.90322948) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21891864) q[0];
sx q[0];
rz(-0.9454596) q[0];
sx q[0];
rz(-0.7801642) q[0];
rz(1.5527234) q[2];
sx q[2];
rz(-2.4583092) q[2];
sx q[2];
rz(2.9073213) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.049202327) q[1];
sx q[1];
rz(-2.4396439) q[1];
sx q[1];
rz(-0.54767056) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0704449) q[3];
sx q[3];
rz(-1.4649434) q[3];
sx q[3];
rz(1.9377886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-3.004965) q[2];
sx q[2];
rz(2.7083) q[2];
rz(-2.8516234) q[3];
sx q[3];
rz(-2.2765997) q[3];
sx q[3];
rz(-2.8210631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1746154) q[0];
sx q[0];
rz(-0.86242914) q[0];
sx q[0];
rz(1.2906661) q[0];
rz(0.67944747) q[1];
sx q[1];
rz(-0.77630711) q[1];
sx q[1];
rz(2.2490833) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7793625) q[0];
sx q[0];
rz(-0.17900544) q[0];
sx q[0];
rz(-1.711861) q[0];
rz(1.9694445) q[2];
sx q[2];
rz(-0.95602334) q[2];
sx q[2];
rz(3.1190256) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1074236) q[1];
sx q[1];
rz(-0.71798199) q[1];
sx q[1];
rz(2.7173244) q[1];
x q[2];
rz(-2.1492747) q[3];
sx q[3];
rz(-0.43799339) q[3];
sx q[3];
rz(2.7051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9925925) q[2];
sx q[2];
rz(-1.0789824) q[2];
sx q[2];
rz(2.0089669) q[2];
rz(2.4752786) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(1.9168436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76729524) q[0];
sx q[0];
rz(-2.7485924) q[0];
sx q[0];
rz(-2.1597916) q[0];
rz(2.1994195) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(0.91032496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86493409) q[0];
sx q[0];
rz(-0.66682839) q[0];
sx q[0];
rz(2.4163321) q[0];
rz(-pi) q[1];
rz(1.4508171) q[2];
sx q[2];
rz(-1.2960805) q[2];
sx q[2];
rz(-3.0697696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7187319) q[1];
sx q[1];
rz(-1.2694468) q[1];
sx q[1];
rz(-1.919793) q[1];
rz(1.1317838) q[3];
sx q[3];
rz(-1.7903504) q[3];
sx q[3];
rz(-1.1846402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2008449) q[2];
sx q[2];
rz(-0.36483279) q[2];
sx q[2];
rz(2.9546837) q[2];
rz(2.9777891) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(-3.0189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(0.69766587) q[0];
rz(-0.54667073) q[1];
sx q[1];
rz(-0.29269871) q[1];
sx q[1];
rz(-1.1950511) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.228926) q[0];
sx q[0];
rz(-2.1536835) q[0];
sx q[0];
rz(2.7253662) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61364737) q[2];
sx q[2];
rz(-0.23229182) q[2];
sx q[2];
rz(-0.42343806) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4711485) q[1];
sx q[1];
rz(-0.19757825) q[1];
sx q[1];
rz(1.0001282) q[1];
rz(0.073506876) q[3];
sx q[3];
rz(-1.8544844) q[3];
sx q[3];
rz(-2.7105041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4635072) q[2];
sx q[2];
rz(-1.7184075) q[2];
sx q[2];
rz(-0.19742337) q[2];
rz(3.0366963) q[3];
sx q[3];
rz(-1.1095122) q[3];
sx q[3];
rz(-3.0535898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5474434) q[0];
sx q[0];
rz(-0.53590411) q[0];
sx q[0];
rz(1.0092258) q[0];
rz(0.028845305) q[1];
sx q[1];
rz(-1.4776769) q[1];
sx q[1];
rz(0.65863329) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.490059) q[0];
sx q[0];
rz(-2.5680827) q[0];
sx q[0];
rz(-2.308368) q[0];
x q[1];
rz(2.9372327) q[2];
sx q[2];
rz(-1.0931226) q[2];
sx q[2];
rz(-2.7088239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4411021) q[1];
sx q[1];
rz(-1.4785188) q[1];
sx q[1];
rz(0.16903318) q[1];
rz(-pi) q[2];
rz(1.8010718) q[3];
sx q[3];
rz(-0.88876969) q[3];
sx q[3];
rz(-1.3547699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.37904) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(-0.23400447) q[2];
rz(1.0903357) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(1.0696629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2284018) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(1.2983904) q[0];
rz(-2.3550854) q[1];
sx q[1];
rz(-0.85580099) q[1];
sx q[1];
rz(1.3589121) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8131065) q[0];
sx q[0];
rz(-2.4755602) q[0];
sx q[0];
rz(0.12909992) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4491389) q[2];
sx q[2];
rz(-0.34231774) q[2];
sx q[2];
rz(-0.76215012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.003215) q[1];
sx q[1];
rz(-1.5724465) q[1];
sx q[1];
rz(1.5720075) q[1];
x q[2];
rz(1.2511061) q[3];
sx q[3];
rz(-2.0123008) q[3];
sx q[3];
rz(0.46605817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7274373) q[2];
sx q[2];
rz(-1.6748019) q[2];
sx q[2];
rz(2.9690913) q[2];
rz(2.6532459) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(1.1138227) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1682424) q[0];
sx q[0];
rz(-0.872648) q[0];
sx q[0];
rz(-0.33997047) q[0];
rz(-0.32644692) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(-1.9065769) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1025196) q[0];
sx q[0];
rz(-0.93329859) q[0];
sx q[0];
rz(1.1028642) q[0];
x q[1];
rz(0.65151524) q[2];
sx q[2];
rz(-2.5861202) q[2];
sx q[2];
rz(-3.0998203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81081796) q[1];
sx q[1];
rz(-1.0913117) q[1];
sx q[1];
rz(-1.564784) q[1];
x q[2];
rz(-0.49576393) q[3];
sx q[3];
rz(-2.7003482) q[3];
sx q[3];
rz(-0.10291162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1116703) q[2];
sx q[2];
rz(-0.80283529) q[2];
sx q[2];
rz(-2.5862582) q[2];
rz(2.9739144) q[3];
sx q[3];
rz(-1.6043112) q[3];
sx q[3];
rz(-0.47784561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(2.4324206) q[0];
sx q[0];
rz(-0.92996159) q[0];
sx q[0];
rz(-2.0322556) q[0];
rz(-2.0093911) q[1];
sx q[1];
rz(-2.7448476) q[1];
sx q[1];
rz(-2.1999377) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8245715) q[0];
sx q[0];
rz(-2.5618658) q[0];
sx q[0];
rz(1.8648022) q[0];
rz(-1.2852816) q[2];
sx q[2];
rz(-1.4219577) q[2];
sx q[2];
rz(0.51644737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1600765) q[1];
sx q[1];
rz(-1.4443411) q[1];
sx q[1];
rz(-0.92595788) q[1];
x q[2];
rz(0.38576084) q[3];
sx q[3];
rz(-1.7683791) q[3];
sx q[3];
rz(2.945874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38810101) q[2];
sx q[2];
rz(-1.0270303) q[2];
sx q[2];
rz(-1.5388185) q[2];
rz(1.3711035) q[3];
sx q[3];
rz(-1.95581) q[3];
sx q[3];
rz(0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082315363) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(2.1737461) q[0];
rz(1.2713426) q[1];
sx q[1];
rz(-1.8495193) q[1];
sx q[1];
rz(0.06180067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7907352) q[0];
sx q[0];
rz(-0.684832) q[0];
sx q[0];
rz(0.63836309) q[0];
rz(2.4872753) q[2];
sx q[2];
rz(-2.8192602) q[2];
sx q[2];
rz(0.92609012) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95434953) q[1];
sx q[1];
rz(-2.4496485) q[1];
sx q[1];
rz(-2.991597) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8093407) q[3];
sx q[3];
rz(-1.8470754) q[3];
sx q[3];
rz(1.075923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1900078) q[2];
sx q[2];
rz(-0.20415674) q[2];
sx q[2];
rz(3.0657213) q[2];
rz(-1.9742981) q[3];
sx q[3];
rz(-2.0397525) q[3];
sx q[3];
rz(-0.30437881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17700125) q[0];
sx q[0];
rz(-0.27934203) q[0];
sx q[0];
rz(2.8569073) q[0];
rz(3.0668861) q[1];
sx q[1];
rz(-1.9120522) q[1];
sx q[1];
rz(-1.7237192) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79399949) q[0];
sx q[0];
rz(-0.82836223) q[0];
sx q[0];
rz(1.30778) q[0];
rz(-1.8080797) q[2];
sx q[2];
rz(-2.5717989) q[2];
sx q[2];
rz(-0.89536086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98445856) q[1];
sx q[1];
rz(-2.1878831) q[1];
sx q[1];
rz(2.0871215) q[1];
rz(-pi) q[2];
rz(3.10793) q[3];
sx q[3];
rz(-2.2164248) q[3];
sx q[3];
rz(-2.1098304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3384) q[2];
sx q[2];
rz(-1.1407547) q[2];
sx q[2];
rz(1.8624064) q[2];
rz(0.98215669) q[3];
sx q[3];
rz(-1.6833865) q[3];
sx q[3];
rz(0.88472432) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2411156) q[0];
sx q[0];
rz(-1.3237088) q[0];
sx q[0];
rz(1.5644912) q[0];
rz(1.0152394) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(-0.10812689) q[2];
sx q[2];
rz(-1.3095289) q[2];
sx q[2];
rz(1.1332569) q[2];
rz(-2.1463263) q[3];
sx q[3];
rz(-2.223816) q[3];
sx q[3];
rz(1.73903) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
