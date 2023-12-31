OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6025699) q[0];
sx q[0];
rz(-0.56350001) q[0];
sx q[0];
rz(-2.6846057) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(-2.4609202) q[1];
sx q[1];
rz(-1.2759804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5029992) q[0];
sx q[0];
rz(-0.25829691) q[0];
sx q[0];
rz(-0.55708416) q[0];
x q[1];
rz(-0.090473526) q[2];
sx q[2];
rz(-1.6455368) q[2];
sx q[2];
rz(-2.5475516) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(3.0800746) q[1];
rz(-0.92313487) q[3];
sx q[3];
rz(-1.8895518) q[3];
sx q[3];
rz(1.2649328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(2.9764552) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5052658) q[0];
sx q[0];
rz(-0.16997448) q[0];
sx q[0];
rz(0.045036137) q[0];
rz(2.8024407) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(-0.80274686) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40753663) q[0];
sx q[0];
rz(-1.3450087) q[0];
sx q[0];
rz(-1.9682103) q[0];
rz(-pi) q[1];
rz(-1.3834125) q[2];
sx q[2];
rz(-2.6307081) q[2];
sx q[2];
rz(0.62142205) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85125605) q[1];
sx q[1];
rz(-0.95278059) q[1];
sx q[1];
rz(3.0783155) q[1];
x q[2];
rz(2.2400168) q[3];
sx q[3];
rz(-1.672861) q[3];
sx q[3];
rz(-1.7166184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95214343) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(2.6248698) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.3396938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8782608) q[0];
sx q[0];
rz(-0.86255951) q[0];
sx q[0];
rz(-2.1398628) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2604733) q[2];
sx q[2];
rz(-1.8124049) q[2];
sx q[2];
rz(1.9067665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(1.7596485) q[1];
x q[2];
rz(2.6684746) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(-2.8603922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.7335256) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6614439) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(-2.9273422) q[0];
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
rz(-2.9458981) q[0];
sx q[0];
rz(-1.3839405) q[0];
sx q[0];
rz(-2.9325783) q[0];
x q[1];
rz(0.71422691) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(2.5941338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2506634) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(-2.3198747) q[1];
rz(-pi) q[2];
x q[2];
rz(0.781226) q[3];
sx q[3];
rz(-2.9328212) q[3];
sx q[3];
rz(-2.3695932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2864705) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(-0.66037035) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.9821092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75393049) q[0];
sx q[0];
rz(-1.2571063) q[0];
sx q[0];
rz(-2.4051106) q[0];
x q[1];
rz(1.7972838) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(-0.047705334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0817464) q[1];
sx q[1];
rz(-2.1085599) q[1];
sx q[1];
rz(-0.50077011) q[1];
x q[2];
rz(-0.63727832) q[3];
sx q[3];
rz(-2.0162458) q[3];
sx q[3];
rz(-2.3973936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.435047) q[2];
sx q[2];
rz(-2.0303625) q[2];
sx q[2];
rz(-1.2716028) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.2579066) q[3];
sx q[3];
rz(3.0767373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(0.043958157) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(-3.040722) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80412241) q[0];
sx q[0];
rz(-0.94479783) q[0];
sx q[0];
rz(1.3834329) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6542873) q[2];
sx q[2];
rz(-1.8308182) q[2];
sx q[2];
rz(-3.0234408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74360352) q[1];
sx q[1];
rz(-1.3685703) q[1];
sx q[1];
rz(-2.9823751) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3376065) q[3];
sx q[3];
rz(-1.1753852) q[3];
sx q[3];
rz(0.49664859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(0.18151367) q[3];
sx q[3];
rz(-0.74129024) q[3];
sx q[3];
rz(1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(-0.96310258) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(-2.4538453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46524099) q[0];
sx q[0];
rz(-0.40333336) q[0];
sx q[0];
rz(-2.1562735) q[0];
x q[1];
rz(-2.9155832) q[2];
sx q[2];
rz(-1.8812211) q[2];
sx q[2];
rz(-1.7058536) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75377611) q[1];
sx q[1];
rz(-3.0101335) q[1];
sx q[1];
rz(-2.9629575) q[1];
rz(-pi) q[2];
rz(1.0024768) q[3];
sx q[3];
rz(-1.5573313) q[3];
sx q[3];
rz(-0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8035651) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(-0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-0.079285346) q[0];
rz(2.0203363) q[1];
sx q[1];
rz(-2.2032578) q[1];
sx q[1];
rz(-1.1987196) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0872902) q[0];
sx q[0];
rz(-2.0044921) q[0];
sx q[0];
rz(-1.646423) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.933868) q[2];
sx q[2];
rz(-1.6512401) q[2];
sx q[2];
rz(2.9094839) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5115576) q[1];
sx q[1];
rz(-2.6361598) q[1];
sx q[1];
rz(0.10314421) q[1];
rz(-pi) q[2];
rz(2.6783887) q[3];
sx q[3];
rz(-1.8052881) q[3];
sx q[3];
rz(-0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(-2.7834535) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.566074) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(0.026780216) q[1];
sx q[1];
rz(-2.2356922) q[1];
sx q[1];
rz(2.2081597) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844922) q[0];
sx q[0];
rz(-1.677209) q[0];
sx q[0];
rz(-1.5110925) q[0];
x q[1];
rz(2.80745) q[2];
sx q[2];
rz(-1.2103989) q[2];
sx q[2];
rz(-2.6956218) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48298207) q[1];
sx q[1];
rz(-0.88739363) q[1];
sx q[1];
rz(1.618209) q[1];
rz(0.65255717) q[3];
sx q[3];
rz(-1.2823294) q[3];
sx q[3];
rz(1.0126589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0788706) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(-2.771647) q[2];
rz(-0.89921078) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(2.1883011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(-2.9577589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659347) q[0];
sx q[0];
rz(-1.4544393) q[0];
sx q[0];
rz(-0.12136202) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3046938) q[2];
sx q[2];
rz(-0.95015929) q[2];
sx q[2];
rz(-1.2037954) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7067691) q[1];
sx q[1];
rz(-2.4425047) q[1];
sx q[1];
rz(1.8301151) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1179908) q[3];
sx q[3];
rz(-0.58348237) q[3];
sx q[3];
rz(-0.49707801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(-1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68207537) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(0.45869641) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(1.6470439) q[2];
rz(-1.8300874) q[3];
sx q[3];
rz(-2.5355831) q[3];
sx q[3];
rz(-2.2021962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
