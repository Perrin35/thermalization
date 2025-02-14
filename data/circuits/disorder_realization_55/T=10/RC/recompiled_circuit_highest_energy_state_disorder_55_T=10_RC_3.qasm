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
rz(0.48756227) q[0];
sx q[0];
rz(-0.3806448) q[0];
sx q[0];
rz(0.08925499) q[0];
rz(2.6788977) q[1];
sx q[1];
rz(-0.2806288) q[1];
sx q[1];
rz(-1.7791003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93264047) q[0];
sx q[0];
rz(-0.37368837) q[0];
sx q[0];
rz(-0.057500827) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5342185) q[2];
sx q[2];
rz(-1.7192516) q[2];
sx q[2];
rz(2.0853993) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.465724) q[1];
sx q[1];
rz(-1.7935408) q[1];
sx q[1];
rz(-0.59091076) q[1];
rz(-pi) q[2];
rz(1.4815349) q[3];
sx q[3];
rz(-2.0253455) q[3];
sx q[3];
rz(0.7346357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.51916277) q[2];
sx q[2];
rz(-0.76202718) q[2];
sx q[2];
rz(-2.2117174) q[2];
rz(-0.30396384) q[3];
sx q[3];
rz(-1.3320987) q[3];
sx q[3];
rz(-2.9728594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836483) q[0];
sx q[0];
rz(-1.5644263) q[0];
sx q[0];
rz(1.2516578) q[0];
rz(-2.9257863) q[1];
sx q[1];
rz(-1.0390751) q[1];
sx q[1];
rz(0.13914093) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5128327) q[0];
sx q[0];
rz(-2.2248575) q[0];
sx q[0];
rz(0.20166986) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49824841) q[2];
sx q[2];
rz(-1.8115269) q[2];
sx q[2];
rz(1.5125076) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5522668) q[1];
sx q[1];
rz(-1.3331474) q[1];
sx q[1];
rz(2.4219805) q[1];
rz(-0.71897935) q[3];
sx q[3];
rz(-2.1785979) q[3];
sx q[3];
rz(-1.833235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4591878) q[2];
sx q[2];
rz(-0.46899691) q[2];
sx q[2];
rz(-2.0501308) q[2];
rz(1.5952716) q[3];
sx q[3];
rz(-2.0681486) q[3];
sx q[3];
rz(-2.6923164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062155398) q[0];
sx q[0];
rz(-1.2261483) q[0];
sx q[0];
rz(-3.1297744) q[0];
rz(-2.2305409) q[1];
sx q[1];
rz(-1.3900737) q[1];
sx q[1];
rz(-1.8131088) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62579836) q[0];
sx q[0];
rz(-2.5949202) q[0];
sx q[0];
rz(0.31577073) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6027868) q[2];
sx q[2];
rz(-2.9201815) q[2];
sx q[2];
rz(0.087269727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40483002) q[1];
sx q[1];
rz(-1.1859845) q[1];
sx q[1];
rz(2.3865108) q[1];
rz(-1.6823009) q[3];
sx q[3];
rz(-1.498025) q[3];
sx q[3];
rz(-2.7804216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.38982424) q[2];
sx q[2];
rz(-1.6168892) q[2];
sx q[2];
rz(-0.41979182) q[2];
rz(1.3587562) q[3];
sx q[3];
rz(-0.32521453) q[3];
sx q[3];
rz(1.9512008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6841986) q[0];
sx q[0];
rz(-2.0036819) q[0];
sx q[0];
rz(-0.5140636) q[0];
rz(-0.38814107) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(1.911389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83004763) q[0];
sx q[0];
rz(-1.8732949) q[0];
sx q[0];
rz(1.8414458) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72429652) q[2];
sx q[2];
rz(-1.8593374) q[2];
sx q[2];
rz(0.1741102) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.14452) q[1];
sx q[1];
rz(-1.9341661) q[1];
sx q[1];
rz(2.4097869) q[1];
rz(1.3421273) q[3];
sx q[3];
rz(-2.2284751) q[3];
sx q[3];
rz(0.87008892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.5019044) q[2];
sx q[2];
rz(-2.1268667) q[2];
sx q[2];
rz(0.91628966) q[2];
rz(2.6472951) q[3];
sx q[3];
rz(-1.1980134) q[3];
sx q[3];
rz(0.77427197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-2.2275527) q[0];
sx q[0];
rz(-0.95714772) q[0];
sx q[0];
rz(0.48602948) q[0];
rz(-0.60274974) q[1];
sx q[1];
rz(-1.175368) q[1];
sx q[1];
rz(2.026162) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19768342) q[0];
sx q[0];
rz(-0.75431529) q[0];
sx q[0];
rz(-0.59342845) q[0];
rz(-pi) q[1];
rz(0.79430842) q[2];
sx q[2];
rz(-1.2758848) q[2];
sx q[2];
rz(-2.8821534) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72157114) q[1];
sx q[1];
rz(-0.90769115) q[1];
sx q[1];
rz(1.0350735) q[1];
rz(-2.0511469) q[3];
sx q[3];
rz(-1.3355888) q[3];
sx q[3];
rz(3.1344828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9976161) q[2];
sx q[2];
rz(-0.3207427) q[2];
sx q[2];
rz(-1.0616659) q[2];
rz(-1.6210506) q[3];
sx q[3];
rz(-1.942626) q[3];
sx q[3];
rz(-2.2475713) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068037085) q[0];
sx q[0];
rz(-0.03980045) q[0];
sx q[0];
rz(1.918248) q[0];
rz(-2.6262737) q[1];
sx q[1];
rz(-0.66910187) q[1];
sx q[1];
rz(-1.3656778) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66636234) q[0];
sx q[0];
rz(-1.7710553) q[0];
sx q[0];
rz(-0.11563511) q[0];
rz(1.7530982) q[2];
sx q[2];
rz(-1.6627835) q[2];
sx q[2];
rz(-2.0192041) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9881043) q[1];
sx q[1];
rz(-2.0447013) q[1];
sx q[1];
rz(-1.7405645) q[1];
x q[2];
rz(-2.8441098) q[3];
sx q[3];
rz(-1.6053891) q[3];
sx q[3];
rz(2.5643947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.078865163) q[2];
sx q[2];
rz(-1.8026423) q[2];
sx q[2];
rz(1.8143181) q[2];
rz(-1.4158538) q[3];
sx q[3];
rz(-3.0684107) q[3];
sx q[3];
rz(2.2590526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.13407229) q[0];
sx q[0];
rz(-2.6044758) q[0];
sx q[0];
rz(-0.19675955) q[0];
rz(2.9255731) q[1];
sx q[1];
rz(-2.1327503) q[1];
sx q[1];
rz(0.78741995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8926933) q[0];
sx q[0];
rz(-1.8149879) q[0];
sx q[0];
rz(-0.1949744) q[0];
rz(-pi) q[1];
rz(2.6419053) q[2];
sx q[2];
rz(-2.5737107) q[2];
sx q[2];
rz(-0.5918006) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0328788) q[1];
sx q[1];
rz(-1.5252264) q[1];
sx q[1];
rz(0.5918064) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18936746) q[3];
sx q[3];
rz(-1.8872617) q[3];
sx q[3];
rz(2.6031983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.54062033) q[2];
sx q[2];
rz(-1.5121907) q[2];
sx q[2];
rz(-2.7579894) q[2];
rz(2.2331734) q[3];
sx q[3];
rz(-2.3288265) q[3];
sx q[3];
rz(2.9878591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8613794) q[0];
sx q[0];
rz(-2.6160243) q[0];
sx q[0];
rz(-2.5049765) q[0];
rz(-2.5909297) q[1];
sx q[1];
rz(-1.5148342) q[1];
sx q[1];
rz(2.6689463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7028684) q[0];
sx q[0];
rz(-1.4819549) q[0];
sx q[0];
rz(0.6779365) q[0];
x q[1];
rz(-0.51690312) q[2];
sx q[2];
rz(-0.10714018) q[2];
sx q[2];
rz(0.60984367) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9437406) q[1];
sx q[1];
rz(-1.6216369) q[1];
sx q[1];
rz(-0.13266854) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2718464) q[3];
sx q[3];
rz(-2.1808694) q[3];
sx q[3];
rz(1.7536193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7134646) q[2];
sx q[2];
rz(-1.8992004) q[2];
sx q[2];
rz(0.098527519) q[2];
rz(-3.1108917) q[3];
sx q[3];
rz(-1.5714329) q[3];
sx q[3];
rz(-2.9876685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.6105662) q[0];
sx q[0];
rz(-0.95777804) q[0];
sx q[0];
rz(0.73262334) q[0];
rz(0.0064119617) q[1];
sx q[1];
rz(-1.0259722) q[1];
sx q[1];
rz(1.4220672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32860562) q[0];
sx q[0];
rz(-1.3805132) q[0];
sx q[0];
rz(1.8237795) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49894519) q[2];
sx q[2];
rz(-1.8932668) q[2];
sx q[2];
rz(-2.6889963) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.040557794) q[1];
sx q[1];
rz(-0.15753105) q[1];
sx q[1];
rz(0.21871241) q[1];
rz(1.645817) q[3];
sx q[3];
rz(-1.0524233) q[3];
sx q[3];
rz(3.1278267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0473359) q[2];
sx q[2];
rz(-1.7862659) q[2];
sx q[2];
rz(0.52200851) q[2];
rz(3.1212854) q[3];
sx q[3];
rz(-0.69965196) q[3];
sx q[3];
rz(-2.844753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8727528) q[0];
sx q[0];
rz(-1.6167384) q[0];
sx q[0];
rz(1.1312477) q[0];
rz(0.77992431) q[1];
sx q[1];
rz(-1.9386539) q[1];
sx q[1];
rz(0.47526971) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4161454) q[0];
sx q[0];
rz(-2.3088881) q[0];
sx q[0];
rz(-2.0717784) q[0];
x q[1];
rz(0.90311994) q[2];
sx q[2];
rz(-1.756487) q[2];
sx q[2];
rz(0.3214432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.77716) q[1];
sx q[1];
rz(-1.8231943) q[1];
sx q[1];
rz(-0.22478687) q[1];
x q[2];
rz(-1.4251208) q[3];
sx q[3];
rz(-2.2620107) q[3];
sx q[3];
rz(-0.6941877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7705226) q[2];
sx q[2];
rz(-1.5757685) q[2];
sx q[2];
rz(1.2062262) q[2];
rz(-1.3953588) q[3];
sx q[3];
rz(-2.5985056) q[3];
sx q[3];
rz(-0.079308184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4974925) q[0];
sx q[0];
rz(-2.946749) q[0];
sx q[0];
rz(-0.83644833) q[0];
rz(-2.1438228) q[1];
sx q[1];
rz(-1.0918959) q[1];
sx q[1];
rz(0.087654884) q[1];
rz(1.3523819) q[2];
sx q[2];
rz(-1.6383258) q[2];
sx q[2];
rz(0.93371423) q[2];
rz(-3.0042778) q[3];
sx q[3];
rz(-0.68344322) q[3];
sx q[3];
rz(-2.4333011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
