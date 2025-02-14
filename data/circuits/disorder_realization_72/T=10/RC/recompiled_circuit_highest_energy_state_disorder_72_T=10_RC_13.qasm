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
rz(1.3067955) q[0];
sx q[0];
rz(-0.84796325) q[0];
sx q[0];
rz(-0.037394878) q[0];
rz(-1.0132064) q[1];
sx q[1];
rz(-0.4816882) q[1];
sx q[1];
rz(1.4451292) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73611063) q[0];
sx q[0];
rz(-2.6713712) q[0];
sx q[0];
rz(-0.41023631) q[0];
rz(-2.1712473) q[2];
sx q[2];
rz(-3.0171347) q[2];
sx q[2];
rz(-2.4336634) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5132222) q[1];
sx q[1];
rz(-1.0890102) q[1];
sx q[1];
rz(-2.1867153) q[1];
rz(-pi) q[2];
rz(-1.0057428) q[3];
sx q[3];
rz(-1.4742719) q[3];
sx q[3];
rz(2.9526476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41363132) q[2];
sx q[2];
rz(-3.0484338) q[2];
sx q[2];
rz(1.6348582) q[2];
rz(-0.19351752) q[3];
sx q[3];
rz(-2.3573124) q[3];
sx q[3];
rz(2.1957652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042260878) q[0];
sx q[0];
rz(-1.7758545) q[0];
sx q[0];
rz(0.17901626) q[0];
rz(1.5366813) q[1];
sx q[1];
rz(-2.8004526) q[1];
sx q[1];
rz(1.5515597) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5478599) q[0];
sx q[0];
rz(-1.6681801) q[0];
sx q[0];
rz(-1.2800526) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70945246) q[2];
sx q[2];
rz(-1.3724899) q[2];
sx q[2];
rz(2.545216) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4578456) q[1];
sx q[1];
rz(-0.48626712) q[1];
sx q[1];
rz(2.1621428) q[1];
rz(-pi) q[2];
rz(0.87249248) q[3];
sx q[3];
rz(-0.66965646) q[3];
sx q[3];
rz(0.39530259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4701074) q[2];
sx q[2];
rz(-1.476373) q[2];
sx q[2];
rz(-0.022484953) q[2];
rz(2.2454028) q[3];
sx q[3];
rz(-2.360207) q[3];
sx q[3];
rz(-1.5145068) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3317868) q[0];
sx q[0];
rz(-0.2066732) q[0];
sx q[0];
rz(1.608954) q[0];
rz(2.0017852) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(2.2679813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5748567) q[0];
sx q[0];
rz(-2.2169161) q[0];
sx q[0];
rz(2.2227051) q[0];
rz(-pi) q[1];
rz(2.1223091) q[2];
sx q[2];
rz(-0.9992558) q[2];
sx q[2];
rz(-0.71046605) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1908299) q[1];
sx q[1];
rz(-1.0894945) q[1];
sx q[1];
rz(-1.7442622) q[1];
rz(1.4954044) q[3];
sx q[3];
rz(-1.2986548) q[3];
sx q[3];
rz(-2.4484602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5736299) q[2];
sx q[2];
rz(-0.45845389) q[2];
sx q[2];
rz(-2.6652794) q[2];
rz(-1.025398) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(-0.29388139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3697701) q[0];
sx q[0];
rz(-0.85404587) q[0];
sx q[0];
rz(-0.59637946) q[0];
rz(0.75309938) q[1];
sx q[1];
rz(-1.785708) q[1];
sx q[1];
rz(-2.7900043) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.387334) q[0];
sx q[0];
rz(-2.0947347) q[0];
sx q[0];
rz(1.41904) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6591107) q[2];
sx q[2];
rz(-2.1957779) q[2];
sx q[2];
rz(-0.04908726) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.233577) q[1];
sx q[1];
rz(-1.963841) q[1];
sx q[1];
rz(0.24891757) q[1];
rz(-pi) q[2];
rz(0.91686317) q[3];
sx q[3];
rz(-2.0987978) q[3];
sx q[3];
rz(-0.39378402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4190462) q[2];
sx q[2];
rz(-1.498035) q[2];
sx q[2];
rz(2.2124115) q[2];
rz(-1.9123214) q[3];
sx q[3];
rz(-0.036673948) q[3];
sx q[3];
rz(1.2693955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96255985) q[0];
sx q[0];
rz(-0.61229175) q[0];
sx q[0];
rz(-2.6137733) q[0];
rz(-2.5714696) q[1];
sx q[1];
rz(-2.5716883) q[1];
sx q[1];
rz(-0.39400563) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1264982) q[0];
sx q[0];
rz(-1.0929035) q[0];
sx q[0];
rz(-1.3518349) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2038241) q[2];
sx q[2];
rz(-1.0093371) q[2];
sx q[2];
rz(-2.9740564) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3574419) q[1];
sx q[1];
rz(-1.3489745) q[1];
sx q[1];
rz(0.12089575) q[1];
rz(-2.3287613) q[3];
sx q[3];
rz(-2.1413229) q[3];
sx q[3];
rz(-1.3798483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3923308) q[2];
sx q[2];
rz(-1.8529842) q[2];
sx q[2];
rz(-2.0070845) q[2];
rz(-0.025731651) q[3];
sx q[3];
rz(-2.0870356) q[3];
sx q[3];
rz(-0.64190763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5662017) q[0];
sx q[0];
rz(-1.3302777) q[0];
sx q[0];
rz(2.6824685) q[0];
rz(-0.8210012) q[1];
sx q[1];
rz(-2.2586925) q[1];
sx q[1];
rz(-2.6301036) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9022512) q[0];
sx q[0];
rz(-2.5731239) q[0];
sx q[0];
rz(1.0064378) q[0];
x q[1];
rz(1.3364001) q[2];
sx q[2];
rz(-2.2210741) q[2];
sx q[2];
rz(2.5087207) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3977511) q[1];
sx q[1];
rz(-1.4201179) q[1];
sx q[1];
rz(0.57445261) q[1];
rz(-pi) q[2];
rz(-1.4019764) q[3];
sx q[3];
rz(-0.69069117) q[3];
sx q[3];
rz(-2.9258779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56200999) q[2];
sx q[2];
rz(-2.0333813) q[2];
sx q[2];
rz(0.75135922) q[2];
rz(2.5747418) q[3];
sx q[3];
rz(-1.5375429) q[3];
sx q[3];
rz(-1.3273299) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.421627) q[0];
sx q[0];
rz(-2.9599074) q[0];
sx q[0];
rz(3.1340461) q[0];
rz(-0.94002062) q[1];
sx q[1];
rz(-0.66497856) q[1];
sx q[1];
rz(3.0090581) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36245334) q[0];
sx q[0];
rz(-2.4041688) q[0];
sx q[0];
rz(-2.2064184) q[0];
rz(-pi) q[1];
rz(-0.33905115) q[2];
sx q[2];
rz(-2.7586852) q[2];
sx q[2];
rz(0.47557451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.346437) q[1];
sx q[1];
rz(-2.7043685) q[1];
sx q[1];
rz(-0.45262419) q[1];
x q[2];
rz(2.0293983) q[3];
sx q[3];
rz(-2.762714) q[3];
sx q[3];
rz(1.0058114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14564766) q[2];
sx q[2];
rz(-1.6463248) q[2];
sx q[2];
rz(-0.10124595) q[2];
rz(0.64949399) q[3];
sx q[3];
rz(-2.5443304) q[3];
sx q[3];
rz(0.35972843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7262909) q[0];
sx q[0];
rz(-1.8446209) q[0];
sx q[0];
rz(1.2671965) q[0];
rz(1.8750809) q[1];
sx q[1];
rz(-2.9837065) q[1];
sx q[1];
rz(-0.74323851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948461) q[0];
sx q[0];
rz(-2.136793) q[0];
sx q[0];
rz(1.8106242) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2879804) q[2];
sx q[2];
rz(-0.79999812) q[2];
sx q[2];
rz(-1.204042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.083588138) q[1];
sx q[1];
rz(-0.28540719) q[1];
sx q[1];
rz(0.98085515) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2970639) q[3];
sx q[3];
rz(-2.2457079) q[3];
sx q[3];
rz(2.9425987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1242975) q[2];
sx q[2];
rz(-0.20169078) q[2];
sx q[2];
rz(-1.7151493) q[2];
rz(2.1258449) q[3];
sx q[3];
rz(-1.4371212) q[3];
sx q[3];
rz(0.82971853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5468686) q[0];
sx q[0];
rz(-2.4877553) q[0];
sx q[0];
rz(0.015901707) q[0];
rz(-1.6390027) q[1];
sx q[1];
rz(-0.67633164) q[1];
sx q[1];
rz(-0.99448386) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976006) q[0];
sx q[0];
rz(-0.74403896) q[0];
sx q[0];
rz(0.97980325) q[0];
rz(0.15764641) q[2];
sx q[2];
rz(-0.38262832) q[2];
sx q[2];
rz(-1.1747109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82518903) q[1];
sx q[1];
rz(-1.087552) q[1];
sx q[1];
rz(-2.2903328) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9904305) q[3];
sx q[3];
rz(-1.0619319) q[3];
sx q[3];
rz(2.9921032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9859621) q[2];
sx q[2];
rz(-1.7113643) q[2];
sx q[2];
rz(-2.963781) q[2];
rz(1.6832247) q[3];
sx q[3];
rz(-2.7198313) q[3];
sx q[3];
rz(-1.873707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.47852248) q[0];
sx q[0];
rz(-1.2489742) q[0];
sx q[0];
rz(1.2836237) q[0];
rz(0.78370699) q[1];
sx q[1];
rz(-0.77373928) q[1];
sx q[1];
rz(1.3073889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147025) q[0];
sx q[0];
rz(-0.44825867) q[0];
sx q[0];
rz(2.9284555) q[0];
rz(3.1112017) q[2];
sx q[2];
rz(-1.9585397) q[2];
sx q[2];
rz(3.0770242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.02919217) q[1];
sx q[1];
rz(-1.629273) q[1];
sx q[1];
rz(-1.7450208) q[1];
rz(-2.1940007) q[3];
sx q[3];
rz(-0.64631185) q[3];
sx q[3];
rz(-2.050296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54138294) q[2];
sx q[2];
rz(-1.4098488) q[2];
sx q[2];
rz(2.6516338) q[2];
rz(2.0269035) q[3];
sx q[3];
rz(-1.1763108) q[3];
sx q[3];
rz(-0.32976845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212696) q[0];
sx q[0];
rz(-1.2280432) q[0];
sx q[0];
rz(-1.8228774) q[0];
rz(-1.8439138) q[1];
sx q[1];
rz(-2.4083125) q[1];
sx q[1];
rz(-0.42793035) q[1];
rz(0.25925706) q[2];
sx q[2];
rz(-1.8198063) q[2];
sx q[2];
rz(1.168269) q[2];
rz(1.8036203) q[3];
sx q[3];
rz(-0.719984) q[3];
sx q[3];
rz(0.58926982) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
