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
rz(0.19652551) q[0];
sx q[0];
rz(2.3152469) q[0];
sx q[0];
rz(8.5066975) q[0];
rz(-0.11153587) q[1];
sx q[1];
rz(4.9354878) q[1];
sx q[1];
rz(7.8505904) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2409256) q[0];
sx q[0];
rz(-1.5673715) q[0];
sx q[0];
rz(-1.5802556) q[0];
rz(-pi) q[1];
rz(-2.5794633) q[2];
sx q[2];
rz(-1.8888753) q[2];
sx q[2];
rz(-1.095301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5985377) q[1];
sx q[1];
rz(-0.32078136) q[1];
sx q[1];
rz(-1.2307274) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49310103) q[3];
sx q[3];
rz(-1.4639018) q[3];
sx q[3];
rz(-1.3809539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.427318) q[2];
sx q[2];
rz(-1.6552507) q[2];
sx q[2];
rz(-1.2662668) q[2];
rz(-1.0424987) q[3];
sx q[3];
rz(-1.6008585) q[3];
sx q[3];
rz(-1.7460167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.105044) q[0];
sx q[0];
rz(-0.62750134) q[0];
sx q[0];
rz(0.096916048) q[0];
rz(-0.80822432) q[1];
sx q[1];
rz(-2.7327635) q[1];
sx q[1];
rz(-1.2867297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0109804) q[0];
sx q[0];
rz(-1.535278) q[0];
sx q[0];
rz(2.0207441) q[0];
rz(-2.1551844) q[2];
sx q[2];
rz(-2.2041568) q[2];
sx q[2];
rz(-1.6478761) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3073439) q[1];
sx q[1];
rz(-1.1888224) q[1];
sx q[1];
rz(2.3267507) q[1];
x q[2];
rz(1.2169514) q[3];
sx q[3];
rz(-1.8895738) q[3];
sx q[3];
rz(0.48976605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5010995) q[2];
sx q[2];
rz(-2.2576136) q[2];
sx q[2];
rz(-2.7596149) q[2];
rz(0.52465087) q[3];
sx q[3];
rz(-1.7347521) q[3];
sx q[3];
rz(0.14373556) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3749738) q[0];
sx q[0];
rz(-1.5559649) q[0];
sx q[0];
rz(-2.4554456) q[0];
rz(2.0740017) q[1];
sx q[1];
rz(-0.99718863) q[1];
sx q[1];
rz(0.080726191) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19182675) q[0];
sx q[0];
rz(-1.6097665) q[0];
sx q[0];
rz(-1.7944337) q[0];
rz(2.3600134) q[2];
sx q[2];
rz(-2.321876) q[2];
sx q[2];
rz(2.9224861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0559433) q[1];
sx q[1];
rz(-2.9102444) q[1];
sx q[1];
rz(-0.10142188) q[1];
rz(-pi) q[2];
rz(-0.46960652) q[3];
sx q[3];
rz(-2.2966343) q[3];
sx q[3];
rz(-0.71401789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8778235) q[2];
sx q[2];
rz(-0.70085415) q[2];
sx q[2];
rz(-2.7027255) q[2];
rz(-1.8435439) q[3];
sx q[3];
rz(-0.38709199) q[3];
sx q[3];
rz(0.41518655) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2271659) q[0];
sx q[0];
rz(-0.54410797) q[0];
sx q[0];
rz(-0.67657226) q[0];
rz(-3.0034972) q[1];
sx q[1];
rz(-1.0905677) q[1];
sx q[1];
rz(2.9409883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4355136) q[0];
sx q[0];
rz(-1.9301751) q[0];
sx q[0];
rz(0.91889221) q[0];
rz(2.5121763) q[2];
sx q[2];
rz(-1.417629) q[2];
sx q[2];
rz(-1.6352194) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.54040376) q[1];
sx q[1];
rz(-2.5174383) q[1];
sx q[1];
rz(-2.9140364) q[1];
rz(1.7142606) q[3];
sx q[3];
rz(-2.1273566) q[3];
sx q[3];
rz(-2.0549162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5033919) q[2];
sx q[2];
rz(-1.8663422) q[2];
sx q[2];
rz(1.4208043) q[2];
rz(3.0270789) q[3];
sx q[3];
rz(-1.6679461) q[3];
sx q[3];
rz(-0.99944559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81922174) q[0];
sx q[0];
rz(-2.350816) q[0];
sx q[0];
rz(-2.1694699) q[0];
rz(0.32360336) q[1];
sx q[1];
rz(-1.4515896) q[1];
sx q[1];
rz(-1.1824898) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2401449) q[0];
sx q[0];
rz(-1.6109544) q[0];
sx q[0];
rz(1.3754649) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57020541) q[2];
sx q[2];
rz(-1.5657164) q[2];
sx q[2];
rz(-0.069381086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4515498) q[1];
sx q[1];
rz(-1.4937972) q[1];
sx q[1];
rz(-1.4119215) q[1];
rz(-1.0504405) q[3];
sx q[3];
rz(-2.0791868) q[3];
sx q[3];
rz(0.45724487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16367308) q[2];
sx q[2];
rz(-1.7719496) q[2];
sx q[2];
rz(-0.53691205) q[2];
rz(2.4162857) q[3];
sx q[3];
rz(-0.0497497) q[3];
sx q[3];
rz(2.3427826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86730114) q[0];
sx q[0];
rz(-2.0766356) q[0];
sx q[0];
rz(1.6987479) q[0];
rz(0.2105712) q[1];
sx q[1];
rz(-1.6981533) q[1];
sx q[1];
rz(0.64705667) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66523401) q[0];
sx q[0];
rz(-1.826735) q[0];
sx q[0];
rz(-0.38689918) q[0];
x q[1];
rz(0.17433106) q[2];
sx q[2];
rz(-2.5781879) q[2];
sx q[2];
rz(-2.8772815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81457227) q[1];
sx q[1];
rz(-1.2271735) q[1];
sx q[1];
rz(-2.6966641) q[1];
rz(-pi) q[2];
x q[2];
rz(1.296455) q[3];
sx q[3];
rz(-1.8963433) q[3];
sx q[3];
rz(0.59301585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1232542) q[2];
sx q[2];
rz(-1.0945357) q[2];
sx q[2];
rz(1.1798165) q[2];
rz(-1.2849464) q[3];
sx q[3];
rz(-0.40950567) q[3];
sx q[3];
rz(3.0783317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0354075) q[0];
sx q[0];
rz(-1.4466865) q[0];
sx q[0];
rz(-2.2440198) q[0];
rz(-1.3342185) q[1];
sx q[1];
rz(-1.7709657) q[1];
sx q[1];
rz(2.1468377) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332116) q[0];
sx q[0];
rz(-1.2612088) q[0];
sx q[0];
rz(2.3214843) q[0];
rz(2.1734235) q[2];
sx q[2];
rz(-1.1873086) q[2];
sx q[2];
rz(2.9619975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2487476) q[1];
sx q[1];
rz(-1.1214646) q[1];
sx q[1];
rz(-0.32545089) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2901102) q[3];
sx q[3];
rz(-1.5940021) q[3];
sx q[3];
rz(0.9965903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0576386) q[2];
sx q[2];
rz(-1.3292162) q[2];
sx q[2];
rz(-0.29339054) q[2];
rz(3.0883664) q[3];
sx q[3];
rz(-2.2727727) q[3];
sx q[3];
rz(-1.7085541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5720125) q[0];
sx q[0];
rz(-1.7520289) q[0];
sx q[0];
rz(-2.8386175) q[0];
rz(1.6527893) q[1];
sx q[1];
rz(-1.3158512) q[1];
sx q[1];
rz(0.66910076) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2888063) q[0];
sx q[0];
rz(-1.1256721) q[0];
sx q[0];
rz(-1.9381802) q[0];
rz(-pi) q[1];
rz(1.2190518) q[2];
sx q[2];
rz(-1.3430077) q[2];
sx q[2];
rz(0.57727376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.41324612) q[1];
sx q[1];
rz(-1.9921682) q[1];
sx q[1];
rz(-1.0131809) q[1];
x q[2];
rz(1.8101997) q[3];
sx q[3];
rz(-1.1576011) q[3];
sx q[3];
rz(-1.5904782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9992493) q[2];
sx q[2];
rz(-1.8084904) q[2];
sx q[2];
rz(2.2744501) q[2];
rz(-2.0182746) q[3];
sx q[3];
rz(-2.2893548) q[3];
sx q[3];
rz(1.1423906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2808696) q[0];
sx q[0];
rz(-2.6779802) q[0];
sx q[0];
rz(-3.0238357) q[0];
rz(2.2640758) q[1];
sx q[1];
rz(-1.0209457) q[1];
sx q[1];
rz(2.1868475) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60459899) q[0];
sx q[0];
rz(-2.1736896) q[0];
sx q[0];
rz(-2.0379809) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47378503) q[2];
sx q[2];
rz(-1.1228645) q[2];
sx q[2];
rz(0.21096274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6462209) q[1];
sx q[1];
rz(-2.2471161) q[1];
sx q[1];
rz(0.27797525) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68902632) q[3];
sx q[3];
rz(-2.195916) q[3];
sx q[3];
rz(2.2415015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.8942326) q[2];
sx q[2];
rz(-0.33984137) q[2];
sx q[2];
rz(1.6575238) q[2];
rz(1.5225211) q[3];
sx q[3];
rz(-1.7584453) q[3];
sx q[3];
rz(1.1744261) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5840983) q[0];
sx q[0];
rz(-1.7099986) q[0];
sx q[0];
rz(-0.75310055) q[0];
rz(1.9901216) q[1];
sx q[1];
rz(-2.617372) q[1];
sx q[1];
rz(-1.6023191) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1881018) q[0];
sx q[0];
rz(-1.3878146) q[0];
sx q[0];
rz(-1.1142) q[0];
rz(-pi) q[1];
rz(2.2686917) q[2];
sx q[2];
rz(-1.1660327) q[2];
sx q[2];
rz(1.0412316) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11048296) q[1];
sx q[1];
rz(-1.3785161) q[1];
sx q[1];
rz(0.39327217) q[1];
x q[2];
rz(-2.4866546) q[3];
sx q[3];
rz(-1.7823151) q[3];
sx q[3];
rz(1.7251273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3573542) q[2];
sx q[2];
rz(-0.18111649) q[2];
sx q[2];
rz(-0.46585807) q[2];
rz(-1.0957796) q[3];
sx q[3];
rz(-0.992479) q[3];
sx q[3];
rz(0.54212681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0453542) q[0];
sx q[0];
rz(-1.5753373) q[0];
sx q[0];
rz(-1.5684431) q[0];
rz(1.0678328) q[1];
sx q[1];
rz(-0.44034958) q[1];
sx q[1];
rz(-0.82028295) q[1];
rz(-1.3374511) q[2];
sx q[2];
rz(-1.3139616) q[2];
sx q[2];
rz(-1.4001605) q[2];
rz(-2.0910083) q[3];
sx q[3];
rz(-2.0512085) q[3];
sx q[3];
rz(-2.7762085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
