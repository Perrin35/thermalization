OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0527394) q[0];
sx q[0];
rz(-2.7272447) q[0];
sx q[0];
rz(0.9847087) q[0];
rz(-0.86496487) q[1];
sx q[1];
rz(-0.84671658) q[1];
sx q[1];
rz(0.72189483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013896) q[0];
sx q[0];
rz(-0.63562515) q[0];
sx q[0];
rz(0.74936015) q[0];
rz(-1.2217916) q[2];
sx q[2];
rz(-0.86599444) q[2];
sx q[2];
rz(-2.8794131) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1221296) q[1];
sx q[1];
rz(-1.5999428) q[1];
sx q[1];
rz(-1.4997838) q[1];
rz(1.2087565) q[3];
sx q[3];
rz(-2.7305121) q[3];
sx q[3];
rz(0.53888881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7559173) q[2];
sx q[2];
rz(-0.52272457) q[2];
sx q[2];
rz(0.44504607) q[2];
rz(-1.2612777) q[3];
sx q[3];
rz(-1.6877561) q[3];
sx q[3];
rz(1.6311579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234914) q[0];
sx q[0];
rz(-2.0913251) q[0];
sx q[0];
rz(2.4871248) q[0];
rz(-1.0819613) q[1];
sx q[1];
rz(-1.315821) q[1];
sx q[1];
rz(-1.6771603) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0102756) q[0];
sx q[0];
rz(-0.47476381) q[0];
sx q[0];
rz(1.6807589) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0686321) q[2];
sx q[2];
rz(-0.96809298) q[2];
sx q[2];
rz(2.631244) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1480746) q[1];
sx q[1];
rz(-1.6119527) q[1];
sx q[1];
rz(-0.74034526) q[1];
rz(1.9212747) q[3];
sx q[3];
rz(-2.0168024) q[3];
sx q[3];
rz(-1.1498017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4945041) q[2];
sx q[2];
rz(-0.096179811) q[2];
sx q[2];
rz(0.75867009) q[2];
rz(1.8132973) q[3];
sx q[3];
rz(-2.0601065) q[3];
sx q[3];
rz(1.7648511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7741622) q[0];
sx q[0];
rz(-1.4752911) q[0];
sx q[0];
rz(0.037516315) q[0];
rz(2.1027749) q[1];
sx q[1];
rz(-0.33375868) q[1];
sx q[1];
rz(2.2822101) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0276703) q[0];
sx q[0];
rz(-1.3827208) q[0];
sx q[0];
rz(-2.5187224) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6994292) q[2];
sx q[2];
rz(-2.3273902) q[2];
sx q[2];
rz(1.278879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1990314) q[1];
sx q[1];
rz(-1.1479401) q[1];
sx q[1];
rz(-2.7884952) q[1];
rz(-1.0802876) q[3];
sx q[3];
rz(-2.1119364) q[3];
sx q[3];
rz(-3.0681572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77217707) q[2];
sx q[2];
rz(-0.93537664) q[2];
sx q[2];
rz(0.76425648) q[2];
rz(-2.1956826) q[3];
sx q[3];
rz(-1.2063682) q[3];
sx q[3];
rz(-2.2435772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94364828) q[0];
sx q[0];
rz(-2.2479489) q[0];
sx q[0];
rz(1.8290895) q[0];
rz(-2.8370044) q[1];
sx q[1];
rz(-0.99921387) q[1];
sx q[1];
rz(-0.33590683) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2755796) q[0];
sx q[0];
rz(-0.073750138) q[0];
sx q[0];
rz(1.3763675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44563771) q[2];
sx q[2];
rz(-1.012371) q[2];
sx q[2];
rz(1.4134917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.057317352) q[1];
sx q[1];
rz(-1.0773939) q[1];
sx q[1];
rz(-1.1019023) q[1];
rz(-1.5139989) q[3];
sx q[3];
rz(-1.3496163) q[3];
sx q[3];
rz(-0.086139679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37561068) q[2];
sx q[2];
rz(-0.42372647) q[2];
sx q[2];
rz(-1.83164) q[2];
rz(-1.8464108) q[3];
sx q[3];
rz(-0.83596197) q[3];
sx q[3];
rz(2.0415993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64422166) q[0];
sx q[0];
rz(-0.27232429) q[0];
sx q[0];
rz(-2.3261133) q[0];
rz(2.4244335) q[1];
sx q[1];
rz(-1.6970789) q[1];
sx q[1];
rz(1.1517634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4916954) q[0];
sx q[0];
rz(-2.5743352) q[0];
sx q[0];
rz(-2.7277566) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4732501) q[2];
sx q[2];
rz(-2.7615504) q[2];
sx q[2];
rz(-0.10157346) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95204853) q[1];
sx q[1];
rz(-0.77922339) q[1];
sx q[1];
rz(-2.0019931) q[1];
rz(-0.11687704) q[3];
sx q[3];
rz(-2.8377551) q[3];
sx q[3];
rz(-0.60910329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0723476) q[2];
sx q[2];
rz(-0.77584156) q[2];
sx q[2];
rz(1.5577215) q[2];
rz(-0.97258687) q[3];
sx q[3];
rz(-1.2310622) q[3];
sx q[3];
rz(-1.4371654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69777456) q[0];
sx q[0];
rz(-2.8688718) q[0];
sx q[0];
rz(-0.66993237) q[0];
rz(-2.2386235) q[1];
sx q[1];
rz(-0.65330708) q[1];
sx q[1];
rz(3.0526551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.406946) q[0];
sx q[0];
rz(-1.7714689) q[0];
sx q[0];
rz(-2.1964551) q[0];
rz(-pi) q[1];
rz(-1.3207664) q[2];
sx q[2];
rz(-0.93968117) q[2];
sx q[2];
rz(-2.4970412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0873333) q[1];
sx q[1];
rz(-1.3218193) q[1];
sx q[1];
rz(-0.072467828) q[1];
rz(-1.813671) q[3];
sx q[3];
rz(-2.6481934) q[3];
sx q[3];
rz(2.3034434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6637806) q[2];
sx q[2];
rz(-2.0974443) q[2];
sx q[2];
rz(0.021154724) q[2];
rz(-2.2145005) q[3];
sx q[3];
rz(-1.6095716) q[3];
sx q[3];
rz(-2.482614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0878736) q[0];
sx q[0];
rz(-1.7807732) q[0];
sx q[0];
rz(-1.1370283) q[0];
rz(-0.46788767) q[1];
sx q[1];
rz(-1.7158022) q[1];
sx q[1];
rz(-2.5877924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62821913) q[0];
sx q[0];
rz(-2.1495753) q[0];
sx q[0];
rz(2.7952162) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8973068) q[2];
sx q[2];
rz(-2.6061432) q[2];
sx q[2];
rz(-2.698512) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1688331) q[1];
sx q[1];
rz(-2.085003) q[1];
sx q[1];
rz(2.7943816) q[1];
rz(-pi) q[2];
rz(-0.99618995) q[3];
sx q[3];
rz(-1.6227727) q[3];
sx q[3];
rz(1.0135723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0566473) q[2];
sx q[2];
rz(-2.2979996) q[2];
sx q[2];
rz(0.6960558) q[2];
rz(2.7737235) q[3];
sx q[3];
rz(-1.1080247) q[3];
sx q[3];
rz(0.28611103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2680161) q[0];
sx q[0];
rz(-1.3883256) q[0];
sx q[0];
rz(2.3263113) q[0];
rz(-2.6878327) q[1];
sx q[1];
rz(-0.90180698) q[1];
sx q[1];
rz(-2.544983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19102272) q[0];
sx q[0];
rz(-2.2248984) q[0];
sx q[0];
rz(-2.4838537) q[0];
rz(-pi) q[1];
rz(3.0456226) q[2];
sx q[2];
rz(-2.308613) q[2];
sx q[2];
rz(-2.677409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.836005) q[1];
sx q[1];
rz(-2.506425) q[1];
sx q[1];
rz(-2.0912384) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43858975) q[3];
sx q[3];
rz(-1.0209173) q[3];
sx q[3];
rz(-0.97766961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10670184) q[2];
sx q[2];
rz(-1.7009578) q[2];
sx q[2];
rz(1.2197257) q[2];
rz(1.8969511) q[3];
sx q[3];
rz(-1.605426) q[3];
sx q[3];
rz(0.19967782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91584665) q[0];
sx q[0];
rz(-0.93996489) q[0];
sx q[0];
rz(-1.4858656) q[0];
rz(2.0303717) q[1];
sx q[1];
rz(-1.8467555) q[1];
sx q[1];
rz(1.750754) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0593587) q[0];
sx q[0];
rz(-1.9080722) q[0];
sx q[0];
rz(-1.5683334) q[0];
rz(-1.6860854) q[2];
sx q[2];
rz(-1.2284245) q[2];
sx q[2];
rz(-1.1631249) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0842485) q[1];
sx q[1];
rz(-2.930713) q[1];
sx q[1];
rz(1.4046304) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3284055) q[3];
sx q[3];
rz(-2.6549746) q[3];
sx q[3];
rz(-0.77017654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2755462) q[2];
sx q[2];
rz(-1.486843) q[2];
sx q[2];
rz(-2.9545968) q[2];
rz(-2.9343119) q[3];
sx q[3];
rz(-2.5076301) q[3];
sx q[3];
rz(-2.0980339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.2279376) q[0];
sx q[0];
rz(-0.74988237) q[0];
sx q[0];
rz(0.0023181152) q[0];
rz(-1.9193316) q[1];
sx q[1];
rz(-1.5536676) q[1];
sx q[1];
rz(1.4171756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9681536) q[0];
sx q[0];
rz(-1.7002956) q[0];
sx q[0];
rz(2.9347675) q[0];
rz(-1.8882636) q[2];
sx q[2];
rz(-1.0523044) q[2];
sx q[2];
rz(-0.41050342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.05310381) q[1];
sx q[1];
rz(-1.3657059) q[1];
sx q[1];
rz(-3.0916832) q[1];
rz(1.8270666) q[3];
sx q[3];
rz(-0.60171222) q[3];
sx q[3];
rz(-0.0046241143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2554539) q[2];
sx q[2];
rz(-1.5172493) q[2];
sx q[2];
rz(-0.21190602) q[2];
rz(-0.40031561) q[3];
sx q[3];
rz(-2.2369657) q[3];
sx q[3];
rz(1.3306085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048007456) q[0];
sx q[0];
rz(-1.891991) q[0];
sx q[0];
rz(1.4461507) q[0];
rz(-0.49225898) q[1];
sx q[1];
rz(-1.4795563) q[1];
sx q[1];
rz(1.5835887) q[1];
rz(-1.7854431) q[2];
sx q[2];
rz(-1.3673906) q[2];
sx q[2];
rz(-2.2833952) q[2];
rz(1.0633123) q[3];
sx q[3];
rz(-2.0290658) q[3];
sx q[3];
rz(3.1334044) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
