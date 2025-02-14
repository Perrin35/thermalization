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
rz(-2.5833997) q[0];
sx q[0];
rz(-0.71891958) q[0];
sx q[0];
rz(-0.67607003) q[0];
rz(-2.8683635) q[1];
sx q[1];
rz(-0.9571119) q[1];
sx q[1];
rz(-0.91125429) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.229743) q[0];
sx q[0];
rz(-1.5485244) q[0];
sx q[0];
rz(-3.1215206) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62556556) q[2];
sx q[2];
rz(-1.564933) q[2];
sx q[2];
rz(-1.0926343) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82750597) q[1];
sx q[1];
rz(-0.99775278) q[1];
sx q[1];
rz(1.5324462) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4449213) q[3];
sx q[3];
rz(-0.81842234) q[3];
sx q[3];
rz(-2.6361806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91997826) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(-2.3294219) q[2];
rz(-0.075856097) q[3];
sx q[3];
rz(-2.4935738) q[3];
sx q[3];
rz(-0.59504741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49038637) q[0];
sx q[0];
rz(-2.6631329) q[0];
sx q[0];
rz(1.867021) q[0];
rz(1.6365341) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(1.9777745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464435) q[0];
sx q[0];
rz(-2.4730485) q[0];
sx q[0];
rz(-0.77495452) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50588079) q[2];
sx q[2];
rz(-1.9569279) q[2];
sx q[2];
rz(-0.94368151) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5608985) q[1];
sx q[1];
rz(-2.1836062) q[1];
sx q[1];
rz(-1.3152907) q[1];
x q[2];
rz(0.31812654) q[3];
sx q[3];
rz(-1.9788747) q[3];
sx q[3];
rz(2.7216743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85094467) q[2];
sx q[2];
rz(-0.023357563) q[2];
sx q[2];
rz(-1.5604431) q[2];
rz(-0.23049878) q[3];
sx q[3];
rz(-1.0483619) q[3];
sx q[3];
rz(0.44052625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480963) q[0];
sx q[0];
rz(-0.31262147) q[0];
sx q[0];
rz(3.0189959) q[0];
rz(0.7971881) q[1];
sx q[1];
rz(-1.0120579) q[1];
sx q[1];
rz(3.0534993) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3159972) q[0];
sx q[0];
rz(-1.3154217) q[0];
sx q[0];
rz(-2.3488464) q[0];
rz(-1.2253907) q[2];
sx q[2];
rz(-0.24125762) q[2];
sx q[2];
rz(3.0221107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3077724) q[1];
sx q[1];
rz(-1.4188179) q[1];
sx q[1];
rz(-3.003158) q[1];
x q[2];
rz(-0.54666211) q[3];
sx q[3];
rz(-2.529749) q[3];
sx q[3];
rz(-2.8091968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88283551) q[2];
sx q[2];
rz(-1.5587403) q[2];
sx q[2];
rz(-1.3239512) q[2];
rz(0.49805182) q[3];
sx q[3];
rz(-2.1784454) q[3];
sx q[3];
rz(0.74670416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28904706) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(-0.14719851) q[0];
rz(-2.4920801) q[1];
sx q[1];
rz(-1.6308866) q[1];
sx q[1];
rz(1.46896) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8513059) q[0];
sx q[0];
rz(-1.3477579) q[0];
sx q[0];
rz(-3.1342491) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9584453) q[2];
sx q[2];
rz(-1.1467271) q[2];
sx q[2];
rz(-1.1893502) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4893579) q[1];
sx q[1];
rz(-1.918004) q[1];
sx q[1];
rz(2.2330797) q[1];
rz(-1.1280414) q[3];
sx q[3];
rz(-2.1078133) q[3];
sx q[3];
rz(2.2396416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9327675) q[2];
sx q[2];
rz(-0.494445) q[2];
sx q[2];
rz(0.24173582) q[2];
rz(0.73511165) q[3];
sx q[3];
rz(-1.3999516) q[3];
sx q[3];
rz(-0.66489768) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46297896) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(-3.0188634) q[0];
rz(-2.2562476) q[1];
sx q[1];
rz(-1.0095936) q[1];
sx q[1];
rz(2.5811894) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3353676) q[0];
sx q[0];
rz(-2.0809552) q[0];
sx q[0];
rz(0.2038184) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3483062) q[2];
sx q[2];
rz(-1.8610536) q[2];
sx q[2];
rz(0.09387108) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6736629) q[1];
sx q[1];
rz(-1.1129675) q[1];
sx q[1];
rz(0.28854847) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0521096) q[3];
sx q[3];
rz(-1.9718687) q[3];
sx q[3];
rz(0.29386753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.6178599) q[2];
sx q[2];
rz(-2.3475519) q[2];
sx q[2];
rz(2.9368371) q[2];
rz(0.93920416) q[3];
sx q[3];
rz(-1.771628) q[3];
sx q[3];
rz(-3.052616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27696779) q[0];
sx q[0];
rz(-3.0245916) q[0];
sx q[0];
rz(-2.4846039) q[0];
rz(0.14343801) q[1];
sx q[1];
rz(-1.5564432) q[1];
sx q[1];
rz(2.4030446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.030169372) q[0];
sx q[0];
rz(-2.2403754) q[0];
sx q[0];
rz(1.4233146) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22566585) q[2];
sx q[2];
rz(-2.0516925) q[2];
sx q[2];
rz(-1.7250329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1164897) q[1];
sx q[1];
rz(-1.8855576) q[1];
sx q[1];
rz(2.2584951) q[1];
rz(0.62885999) q[3];
sx q[3];
rz(-0.81506461) q[3];
sx q[3];
rz(-0.12303837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8831545) q[2];
sx q[2];
rz(-0.65885764) q[2];
sx q[2];
rz(-3.1385885) q[2];
rz(-2.352412) q[3];
sx q[3];
rz(-2.7572258) q[3];
sx q[3];
rz(1.9743617) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0810735) q[0];
sx q[0];
rz(-0.34927148) q[0];
sx q[0];
rz(2.9935167) q[0];
rz(-2.608346) q[1];
sx q[1];
rz(-2.8956469) q[1];
sx q[1];
rz(1.5124793) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5083744) q[0];
sx q[0];
rz(-2.0323159) q[0];
sx q[0];
rz(1.870023) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.327269) q[2];
sx q[2];
rz(-1.0934798) q[2];
sx q[2];
rz(0.60233357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.032669205) q[1];
sx q[1];
rz(-1.9100238) q[1];
sx q[1];
rz(-0.12813385) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.809172) q[3];
sx q[3];
rz(-0.73638232) q[3];
sx q[3];
rz(-0.045698085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.44275722) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(-0.093078144) q[2];
rz(-3.0794411) q[3];
sx q[3];
rz(-0.39669865) q[3];
sx q[3];
rz(1.2232346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.020697866) q[0];
sx q[0];
rz(-0.59497213) q[0];
sx q[0];
rz(0.66463941) q[0];
rz(-1.7918034) q[1];
sx q[1];
rz(-0.87769687) q[1];
sx q[1];
rz(-2.452449) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6751373) q[0];
sx q[0];
rz(-0.27091089) q[0];
sx q[0];
rz(-2.6938426) q[0];
rz(0.87920636) q[2];
sx q[2];
rz(-1.1779281) q[2];
sx q[2];
rz(-2.526432) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0943567) q[1];
sx q[1];
rz(-2.7942889) q[1];
sx q[1];
rz(-2.7014665) q[1];
rz(-pi) q[2];
rz(2.9767959) q[3];
sx q[3];
rz(-1.5110803) q[3];
sx q[3];
rz(2.7415558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0906715) q[2];
sx q[2];
rz(-2.4854269) q[2];
sx q[2];
rz(1.7655168) q[2];
rz(-1.9164267) q[3];
sx q[3];
rz(-0.71198946) q[3];
sx q[3];
rz(-0.041697748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10251481) q[0];
sx q[0];
rz(-3.097105) q[0];
sx q[0];
rz(-2.5503889) q[0];
rz(2.8996331) q[1];
sx q[1];
rz(-1.0958902) q[1];
sx q[1];
rz(-2.718149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71899881) q[0];
sx q[0];
rz(-0.43433493) q[0];
sx q[0];
rz(2.3584764) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9482163) q[2];
sx q[2];
rz(-0.49623734) q[2];
sx q[2];
rz(-1.9691182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8264173) q[1];
sx q[1];
rz(-1.8674382) q[1];
sx q[1];
rz(1.0147995) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1657646) q[3];
sx q[3];
rz(-0.99294239) q[3];
sx q[3];
rz(-0.33194968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.55687904) q[2];
sx q[2];
rz(-2.5355279) q[2];
sx q[2];
rz(1.1663743) q[2];
rz(-1.1276468) q[3];
sx q[3];
rz(-0.96846628) q[3];
sx q[3];
rz(-2.6166272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.128085) q[0];
sx q[0];
rz(-0.43376827) q[0];
sx q[0];
rz(0.67310131) q[0];
rz(-2.290944) q[1];
sx q[1];
rz(-1.3530082) q[1];
sx q[1];
rz(0.32352111) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.464768) q[0];
sx q[0];
rz(-2.6861405) q[0];
sx q[0];
rz(0.97525979) q[0];
rz(0.27574972) q[2];
sx q[2];
rz(-1.1080896) q[2];
sx q[2];
rz(1.5327061) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0464152) q[1];
sx q[1];
rz(-2.2942465) q[1];
sx q[1];
rz(-2.0207538) q[1];
rz(-pi) q[2];
rz(0.58384044) q[3];
sx q[3];
rz(-1.5642318) q[3];
sx q[3];
rz(-3.1232426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2181776) q[2];
sx q[2];
rz(-0.18109334) q[2];
sx q[2];
rz(3.0695445) q[2];
rz(-2.1273023) q[3];
sx q[3];
rz(-2.225596) q[3];
sx q[3];
rz(-2.5924276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2608248) q[0];
sx q[0];
rz(-1.7171971) q[0];
sx q[0];
rz(1.1203753) q[0];
rz(-0.33521677) q[1];
sx q[1];
rz(-2.4991279) q[1];
sx q[1];
rz(2.0229708) q[1];
rz(-3.0006888) q[2];
sx q[2];
rz(-2.0456516) q[2];
sx q[2];
rz(2.4608154) q[2];
rz(0.14306457) q[3];
sx q[3];
rz(-1.5456556) q[3];
sx q[3];
rz(1.7329334) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
