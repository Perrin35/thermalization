OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(-0.45698693) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2109283) q[0];
sx q[0];
rz(-1.352248) q[0];
sx q[0];
rz(1.4320089) q[0];
rz(-0.692042) q[2];
sx q[2];
rz(-3.0243052) q[2];
sx q[2];
rz(2.8534944) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(0.061518016) q[1];
rz(-pi) q[2];
rz(2.071322) q[3];
sx q[3];
rz(-2.4300017) q[3];
sx q[3];
rz(2.4430032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(-2.3993649) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052658) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(-0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(-0.80274686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.734056) q[0];
sx q[0];
rz(-1.3450087) q[0];
sx q[0];
rz(1.1733823) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0741834) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(0.78546055) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3992577) q[1];
sx q[1];
rz(-2.5207673) q[1];
sx q[1];
rz(1.4820815) q[1];
x q[2];
rz(-0.12985274) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(0.065404281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(-1.7155898) q[2];
rz(0.60570335) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-3.0595996) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95214343) q[0];
sx q[0];
rz(-1.9316297) q[0];
sx q[0];
rz(2.8564575) q[0];
rz(2.6248698) q[1];
sx q[1];
rz(-1.4915833) q[1];
sx q[1];
rz(-1.3396938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0403255) q[0];
sx q[0];
rz(-2.2650532) q[0];
sx q[0];
rz(-0.56157748) q[0];
rz(-0.3091829) q[2];
sx q[2];
rz(-0.90484607) q[2];
sx q[2];
rz(0.14112976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(-1.3819441) q[1];
x q[2];
rz(0.6891059) q[3];
sx q[3];
rz(-1.2561241) q[3];
sx q[3];
rz(0.93022197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.7335256) q[2];
rz(0.18313289) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(2.9273422) q[0];
rz(0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(2.4750211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19569451) q[0];
sx q[0];
rz(-1.7576522) q[0];
sx q[0];
rz(0.20901434) q[0];
x q[1];
rz(-0.33118403) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(-0.77510288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.22073332) q[1];
sx q[1];
rz(-2.1174866) q[1];
sx q[1];
rz(-2.2842555) q[1];
rz(0.14931071) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(-3.1130476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8551222) q[2];
sx q[2];
rz(-2.9292967) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(-1.1874229) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51263556) q[0];
sx q[0];
rz(-1.2197138) q[0];
sx q[0];
rz(0.83706013) q[0];
rz(-2.191026) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(-1.1594835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5975208) q[0];
sx q[0];
rz(-2.263875) q[0];
sx q[0];
rz(1.1580628) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97255623) q[2];
sx q[2];
rz(-1.4417357) q[2];
sx q[2];
rz(-1.8051403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0598462) q[1];
sx q[1];
rz(-1.0330327) q[1];
sx q[1];
rz(0.50077011) q[1];
rz(-pi) q[2];
rz(-2.5043143) q[3];
sx q[3];
rz(-2.0162458) q[3];
sx q[3];
rz(2.3973936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7065457) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(-2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(-0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(-2.0264453) q[0];
rz(0.043958157) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(-0.10087068) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8773008) q[0];
sx q[0];
rz(-1.722324) q[0];
sx q[0];
rz(-2.507188) q[0];
rz(0.48730536) q[2];
sx q[2];
rz(-1.3107745) q[2];
sx q[2];
rz(-0.11815182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3979891) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-0.15921758) q[1];
x q[2];
rz(0.40525873) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(2.1586771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-2.2163088) q[2];
rz(0.18151367) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(-1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7589384) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(0.96310258) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(-0.68774736) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6531658) q[0];
sx q[0];
rz(-1.7894206) q[0];
sx q[0];
rz(-1.2290918) q[0];
x q[1];
rz(-2.9155832) q[2];
sx q[2];
rz(-1.2603716) q[2];
sx q[2];
rz(1.7058536) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75377611) q[1];
sx q[1];
rz(-3.0101335) q[1];
sx q[1];
rz(-2.9629575) q[1];
x q[2];
rz(2.1391159) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(-0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-0.95454803) q[0];
sx q[0];
rz(3.0623073) q[0];
rz(1.1212564) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(-1.1987196) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90889701) q[0];
sx q[0];
rz(-0.43982738) q[0];
sx q[0];
rz(-0.16172414) q[0];
rz(-pi) q[1];
rz(-1.347581) q[2];
sx q[2];
rz(-0.37149059) q[2];
sx q[2];
rz(-1.5944634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8504387) q[1];
sx q[1];
rz(-1.6206695) q[1];
sx q[1];
rz(-0.50317851) q[1];
rz(-1.8317322) q[3];
sx q[3];
rz(-1.121215) q[3];
sx q[3];
rz(0.98228067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85598677) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(-2.3525227) q[2];
rz(-2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(2.5532706) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(2.2081597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3073472) q[0];
sx q[0];
rz(-1.630162) q[0];
sx q[0];
rz(-3.0349915) q[0];
rz(-pi) q[1];
rz(-2.80745) q[2];
sx q[2];
rz(-1.9311937) q[2];
sx q[2];
rz(0.44597086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0837299) q[1];
sx q[1];
rz(-1.607556) q[1];
sx q[1];
rz(-0.68395331) q[1];
rz(-1.9282354) q[3];
sx q[3];
rz(-0.94944984) q[3];
sx q[3];
rz(2.3693717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0788706) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(-2.771647) q[2];
rz(2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(-0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(2.9577589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7659347) q[0];
sx q[0];
rz(-1.6871534) q[0];
sx q[0];
rz(3.0202306) q[0];
rz(-0.76639558) q[2];
sx q[2];
rz(-2.1470214) q[2];
sx q[2];
rz(0.11608427) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2059584) q[1];
sx q[1];
rz(-1.7365672) q[1];
sx q[1];
rz(0.88840719) q[1];
rz(3.1179908) q[3];
sx q[3];
rz(-2.5581103) q[3];
sx q[3];
rz(-0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68207537) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(0.36322414) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(2.6828962) q[2];
sx q[2];
rz(-1.9445322) q[2];
sx q[2];
rz(-1.4945488) q[2];
rz(-2.1609859) q[3];
sx q[3];
rz(-1.7173613) q[3];
sx q[3];
rz(2.2955256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
