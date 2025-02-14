OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.46848133) q[0];
sx q[0];
rz(-0.10804478) q[0];
sx q[0];
rz(-2.151902) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(1.4563814) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0385168) q[0];
sx q[0];
rz(-1.4355735) q[0];
sx q[0];
rz(-1.9413906) q[0];
rz(-pi) q[1];
rz(-0.26881071) q[2];
sx q[2];
rz(-3.0033026) q[2];
sx q[2];
rz(0.88218305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8238358) q[1];
sx q[1];
rz(-1.2878622) q[1];
sx q[1];
rz(-2.1499277) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0979324) q[3];
sx q[3];
rz(-2.147445) q[3];
sx q[3];
rz(2.1917564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1504537) q[2];
sx q[2];
rz(-3.1296215) q[2];
sx q[2];
rz(-2.0963304) q[2];
rz(2.1446877) q[3];
sx q[3];
rz(-3.1361339) q[3];
sx q[3];
rz(-1.3158984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5528706) q[0];
sx q[0];
rz(-1.2405688) q[0];
sx q[0];
rz(-1.3528104) q[0];
rz(3.1006587) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(1.5418242) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7517325) q[0];
sx q[0];
rz(-2.9417188) q[0];
sx q[0];
rz(-0.78011192) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5920609) q[2];
sx q[2];
rz(-1.5545115) q[2];
sx q[2];
rz(0.88693888) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92506986) q[1];
sx q[1];
rz(-2.1669186) q[1];
sx q[1];
rz(-0.03742569) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0066693) q[3];
sx q[3];
rz(-0.65186497) q[3];
sx q[3];
rz(-1.8396371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6359977) q[2];
sx q[2];
rz(-0.032568585) q[2];
sx q[2];
rz(-0.51844281) q[2];
rz(-2.017766) q[3];
sx q[3];
rz(-2.3321407) q[3];
sx q[3];
rz(2.7037485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.6464624) q[0];
sx q[0];
rz(-0.050124425) q[0];
sx q[0];
rz(1.6019524) q[0];
rz(0.72499544) q[1];
sx q[1];
rz(-0.031818964) q[1];
sx q[1];
rz(2.4620788) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6294265) q[0];
sx q[0];
rz(-1.1376732) q[0];
sx q[0];
rz(-2.4999451) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9271487) q[2];
sx q[2];
rz(-1.5661998) q[2];
sx q[2];
rz(-1.4397804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3287769) q[1];
sx q[1];
rz(-0.41197244) q[1];
sx q[1];
rz(-2.2669492) q[1];
x q[2];
rz(-2.1587426) q[3];
sx q[3];
rz(-1.3508537) q[3];
sx q[3];
rz(-1.4877121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9069549) q[2];
sx q[2];
rz(-0.24181557) q[2];
sx q[2];
rz(2.6406636) q[2];
rz(0.45589724) q[3];
sx q[3];
rz(-0.026345043) q[3];
sx q[3];
rz(-2.946089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86922115) q[0];
sx q[0];
rz(-3.0932194) q[0];
sx q[0];
rz(-0.79917556) q[0];
rz(-0.29798206) q[1];
sx q[1];
rz(-2.8831392) q[1];
sx q[1];
rz(-2.2395649) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.018343) q[0];
sx q[0];
rz(-2.3756174) q[0];
sx q[0];
rz(-0.19874707) q[0];
rz(-pi) q[1];
rz(2.8935585) q[2];
sx q[2];
rz(-0.55540652) q[2];
sx q[2];
rz(-0.75128864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2711401) q[1];
sx q[1];
rz(-1.5803403) q[1];
sx q[1];
rz(1.4231667) q[1];
x q[2];
rz(3.0864363) q[3];
sx q[3];
rz(-1.4934908) q[3];
sx q[3];
rz(2.5173924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67348376) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(1.6710949) q[2];
rz(-2.9521613) q[3];
sx q[3];
rz(-0.048308689) q[3];
sx q[3];
rz(0.76907492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9162132) q[0];
sx q[0];
rz(-0.12752859) q[0];
sx q[0];
rz(-0.39176971) q[0];
rz(1.0852934) q[1];
sx q[1];
rz(-0.009805209) q[1];
sx q[1];
rz(-2.772803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1802496) q[0];
sx q[0];
rz(-1.215544) q[0];
sx q[0];
rz(0.60549462) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9874175) q[2];
sx q[2];
rz(-1.3236127) q[2];
sx q[2];
rz(0.32933035) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6435191) q[1];
sx q[1];
rz(-0.0066702492) q[1];
sx q[1];
rz(1.4590864) q[1];
x q[2];
rz(1.8043133) q[3];
sx q[3];
rz(-0.87461014) q[3];
sx q[3];
rz(0.96674572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5607249) q[2];
sx q[2];
rz(-0.11595011) q[2];
sx q[2];
rz(1.7725393) q[2];
rz(3.0154058) q[3];
sx q[3];
rz(-2.7044665) q[3];
sx q[3];
rz(-1.060846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.5996025) q[0];
sx q[0];
rz(-0.76102155) q[0];
sx q[0];
rz(-1.5498932) q[0];
rz(-2.1688993) q[1];
sx q[1];
rz(-0.21316554) q[1];
sx q[1];
rz(-2.1125643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9950819) q[0];
sx q[0];
rz(-1.8712988) q[0];
sx q[0];
rz(-2.5142558) q[0];
rz(-pi) q[1];
rz(-2.6112399) q[2];
sx q[2];
rz(-0.27410045) q[2];
sx q[2];
rz(-2.0567187) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29063836) q[1];
sx q[1];
rz(-0.096157638) q[1];
sx q[1];
rz(1.5753217) q[1];
x q[2];
rz(3.1054822) q[3];
sx q[3];
rz(-1.4034158) q[3];
sx q[3];
rz(-0.13425628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35919967) q[2];
sx q[2];
rz(-1.426037) q[2];
sx q[2];
rz(-0.4314118) q[2];
rz(-2.1443478) q[3];
sx q[3];
rz(-0.037171818) q[3];
sx q[3];
rz(0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3417086) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(-2.0836015) q[0];
rz(-2.3175088) q[1];
sx q[1];
rz(-1.7015142e-05) q[1];
sx q[1];
rz(-2.3242059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0426328) q[0];
sx q[0];
rz(-1.3606679) q[0];
sx q[0];
rz(-1.1636583) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36211966) q[2];
sx q[2];
rz(-3.1220925) q[2];
sx q[2];
rz(1.8060613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.73002023) q[1];
sx q[1];
rz(-1.8121094) q[1];
sx q[1];
rz(3.1012721) q[1];
rz(-pi) q[2];
rz(-1.041374) q[3];
sx q[3];
rz(-2.2499871) q[3];
sx q[3];
rz(0.27959945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99007964) q[2];
sx q[2];
rz(-1.238287) q[2];
sx q[2];
rz(1.5129169) q[2];
rz(-0.40142909) q[3];
sx q[3];
rz(-0.028086834) q[3];
sx q[3];
rz(-0.95429558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.7734739) q[0];
sx q[0];
rz(-0.064726949) q[0];
sx q[0];
rz(-1.7777959) q[0];
rz(-0.041944567) q[1];
sx q[1];
rz(-0.13101235) q[1];
sx q[1];
rz(-1.0036453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59597385) q[0];
sx q[0];
rz(-2.501308) q[0];
sx q[0];
rz(0.15592928) q[0];
rz(3.1334468) q[2];
sx q[2];
rz(-1.2334494) q[2];
sx q[2];
rz(-2.103919) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3253206) q[1];
sx q[1];
rz(-0.31041103) q[1];
sx q[1];
rz(2.2903633) q[1];
rz(-pi) q[2];
rz(-1.2673301) q[3];
sx q[3];
rz(-0.93765536) q[3];
sx q[3];
rz(1.2337402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7046788) q[2];
sx q[2];
rz(-0.059160058) q[2];
sx q[2];
rz(1.7436279) q[2];
rz(0.49020234) q[3];
sx q[3];
rz(-3.0981045) q[3];
sx q[3];
rz(-1.9628261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040084664) q[0];
sx q[0];
rz(-3.0136216) q[0];
sx q[0];
rz(0.21139938) q[0];
rz(-2.0231817) q[1];
sx q[1];
rz(-0.011592955) q[1];
sx q[1];
rz(1.6308019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64563066) q[0];
sx q[0];
rz(-1.2617153) q[0];
sx q[0];
rz(1.2675955) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5671179) q[2];
sx q[2];
rz(-0.7099289) q[2];
sx q[2];
rz(1.0021051) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.37247) q[1];
sx q[1];
rz(-1.587953) q[1];
sx q[1];
rz(-0.15461289) q[1];
rz(-pi) q[2];
rz(0.49064891) q[3];
sx q[3];
rz(-2.9534441) q[3];
sx q[3];
rz(2.8110281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7819034) q[2];
sx q[2];
rz(-0.017904559) q[2];
sx q[2];
rz(-0.27931279) q[2];
rz(-2.277788) q[3];
sx q[3];
rz(-3.1371208) q[3];
sx q[3];
rz(0.68252501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.913468) q[0];
sx q[0];
rz(-1.1044015) q[0];
sx q[0];
rz(1.3155235) q[0];
rz(-0.77519351) q[1];
sx q[1];
rz(-2.9029791) q[1];
sx q[1];
rz(1.388789) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086292656) q[0];
sx q[0];
rz(-0.32150966) q[0];
sx q[0];
rz(-2.6525524) q[0];
x q[1];
rz(-2.6493373) q[2];
sx q[2];
rz(-2.6894719) q[2];
sx q[2];
rz(1.4398354) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9541746) q[1];
sx q[1];
rz(-3.1404964) q[1];
sx q[1];
rz(0.61598678) q[1];
rz(-2.7241012) q[3];
sx q[3];
rz(-1.1840828) q[3];
sx q[3];
rz(-1.8698805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76643252) q[2];
sx q[2];
rz(-3.1258686) q[2];
sx q[2];
rz(1.6895705) q[2];
rz(0.25894138) q[3];
sx q[3];
rz(-2.9539234) q[3];
sx q[3];
rz(-1.9848721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6257085) q[0];
sx q[0];
rz(-2.4242171) q[0];
sx q[0];
rz(-1.7051359) q[0];
rz(1.6547849) q[1];
sx q[1];
rz(-2.868352) q[1];
sx q[1];
rz(0.19788338) q[1];
rz(2.941359) q[2];
sx q[2];
rz(-1.5774792) q[2];
sx q[2];
rz(1.7984628) q[2];
rz(1.4987569) q[3];
sx q[3];
rz(-2.7938953) q[3];
sx q[3];
rz(0.10216879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
