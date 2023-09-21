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
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(-1.8656123) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514483) q[0];
sx q[0];
rz(-1.4353308) q[0];
sx q[0];
rz(-0.22060237) q[0];
rz(-pi) q[1];
x q[1];
rz(0.090473526) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(0.59404101) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0253804) q[1];
sx q[1];
rz(-0.21157163) q[1];
sx q[1];
rz(1.2799353) q[1];
x q[2];
rz(2.2184578) q[3];
sx q[3];
rz(-1.2520408) q[3];
sx q[3];
rz(1.8766599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1093381) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(2.3993649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(2.8024407) q[1];
sx q[1];
rz(-1.76666) q[1];
sx q[1];
rz(0.80274686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.734056) q[0];
sx q[0];
rz(-1.3450087) q[0];
sx q[0];
rz(1.9682103) q[0];
x q[1];
rz(1.3834125) q[2];
sx q[2];
rz(-2.6307081) q[2];
sx q[2];
rz(2.5201706) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85125605) q[1];
sx q[1];
rz(-0.95278059) q[1];
sx q[1];
rz(0.063277146) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12985274) q[3];
sx q[3];
rz(-0.90568554) q[3];
sx q[3];
rz(-3.0761884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.4260028) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.4915833) q[1];
sx q[1];
rz(1.8018988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1012672) q[0];
sx q[0];
rz(-0.87653941) q[0];
sx q[0];
rz(2.5800152) q[0];
rz(-pi) q[1];
rz(2.8324098) q[2];
sx q[2];
rz(-2.2367466) q[2];
sx q[2];
rz(-0.14112976) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36741751) q[1];
sx q[1];
rz(-1.7304825) q[1];
sx q[1];
rz(1.3819441) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6684746) q[3];
sx q[3];
rz(-0.74672532) q[3];
sx q[3];
rz(0.28120041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67119917) q[2];
sx q[2];
rz(-0.931804) q[2];
sx q[2];
rz(-1.7335256) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-2.1268842) q[3];
sx q[3];
rz(-3.0696707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4801487) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(-0.42901531) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(-2.4750211) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4857805) q[0];
sx q[0];
rz(-0.27944767) q[0];
sx q[0];
rz(-0.73894545) q[0];
x q[1];
rz(0.33118403) q[2];
sx q[2];
rz(-0.74197717) q[2];
sx q[2];
rz(2.3664898) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9208593) q[1];
sx q[1];
rz(-2.1174866) q[1];
sx q[1];
rz(2.2842555) q[1];
rz(2.9922819) q[3];
sx q[3];
rz(-1.4243323) q[3];
sx q[3];
rz(3.1130476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8551222) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(2.4812223) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(-0.56046265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289571) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(2.3045325) q[0];
rz(-0.95056668) q[1];
sx q[1];
rz(-2.4730813) q[1];
sx q[1];
rz(1.1594835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449074) q[0];
sx q[0];
rz(-0.7888182) q[0];
sx q[0];
rz(2.6916654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3443089) q[2];
sx q[2];
rz(-0.61033568) q[2];
sx q[2];
rz(-0.047705334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0817464) q[1];
sx q[1];
rz(-1.0330327) q[1];
sx q[1];
rz(-0.50077011) q[1];
x q[2];
rz(-1.0347576) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-1.1151474) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(-0.10087068) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4911762) q[0];
sx q[0];
rz(-0.6498148) q[0];
sx q[0];
rz(-2.8894436) q[0];
x q[1];
rz(0.51668824) q[2];
sx q[2];
rz(-0.54737216) q[2];
sx q[2];
rz(-1.237243) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74360352) q[1];
sx q[1];
rz(-1.7730224) q[1];
sx q[1];
rz(2.9823751) q[1];
x q[2];
rz(1.8039861) q[3];
sx q[3];
rz(-1.9662074) q[3];
sx q[3];
rz(0.49664859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4042525) q[2];
sx q[2];
rz(-1.5161783) q[2];
sx q[2];
rz(2.2163088) q[2];
rz(0.18151367) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7589384) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(-1.5902279) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(2.4538453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46524099) q[0];
sx q[0];
rz(-2.7382593) q[0];
sx q[0];
rz(0.98531918) q[0];
rz(-2.9155832) q[2];
sx q[2];
rz(-1.8812211) q[2];
sx q[2];
rz(-1.7058536) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75377611) q[1];
sx q[1];
rz(-0.13145914) q[1];
sx q[1];
rz(-0.17863518) q[1];
rz(1.5958105) q[3];
sx q[3];
rz(-2.5731312) q[3];
sx q[3];
rz(0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(0.079285346) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.942873) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2326956) q[0];
sx q[0];
rz(-2.7017653) q[0];
sx q[0];
rz(-0.16172414) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2077246) q[2];
sx q[2];
rz(-1.4903526) q[2];
sx q[2];
rz(0.23210873) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.7477914) q[1];
sx q[1];
rz(-2.0732905) q[1];
sx q[1];
rz(-1.5138813) q[1];
rz(-pi) q[2];
rz(2.6783887) q[3];
sx q[3];
rz(-1.3363046) q[3];
sx q[3];
rz(0.47298688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2856059) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(-2.3525227) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5877682) q[3];
sx q[3];
rz(-1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.566074) q[0];
sx q[0];
rz(-1.7552567) q[0];
sx q[0];
rz(0.58832204) q[0];
rz(3.1148124) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(-0.9334329) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.371884) q[0];
sx q[0];
rz(-3.0196307) q[0];
sx q[0];
rz(2.6321649) q[0];
x q[1];
rz(0.8546631) q[2];
sx q[2];
rz(-2.6551506) q[2];
sx q[2];
rz(-1.9180627) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.583608) q[1];
sx q[1];
rz(-2.4568111) q[1];
sx q[1];
rz(0.058137356) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45456072) q[3];
sx q[3];
rz(-2.436736) q[3];
sx q[3];
rz(-2.9398033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.062722) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(-0.36994568) q[2];
rz(-2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(-2.1883011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(0.6533587) q[0];
rz(1.6409142) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18098022) q[0];
sx q[0];
rz(-1.4502589) q[0];
sx q[0];
rz(1.4535849) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3751971) q[2];
sx q[2];
rz(-2.1470214) q[2];
sx q[2];
rz(3.0255084) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.373133) q[1];
sx q[1];
rz(-0.8994973) q[1];
sx q[1];
rz(-2.9292604) q[1];
x q[2];
rz(1.5863745) q[3];
sx q[3];
rz(-0.9874978) q[3];
sx q[3];
rz(-2.616236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(2.9392021) q[2];
rz(1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595173) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(-2.6828962) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(1.6470439) q[2];
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