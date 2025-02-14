OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(-0.25622955) q[0];
rz(3.5186634) q[1];
sx q[1];
rz(4.9405603) q[1];
sx q[1];
rz(8.3648051) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.357517) q[0];
sx q[0];
rz(-1.4526443) q[0];
sx q[0];
rz(1.4489001) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69251142) q[2];
sx q[2];
rz(-1.5438809) q[2];
sx q[2];
rz(1.7153502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6457451) q[1];
sx q[1];
rz(-1.331067) q[1];
sx q[1];
rz(1.354753) q[1];
rz(-pi) q[2];
rz(-2.5422402) q[3];
sx q[3];
rz(-0.17440344) q[3];
sx q[3];
rz(-2.6389183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9609191) q[2];
sx q[2];
rz(-0.91327614) q[2];
sx q[2];
rz(-1.712435) q[2];
rz(0.6116496) q[3];
sx q[3];
rz(-1.6983756) q[3];
sx q[3];
rz(3.070224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83577689) q[0];
sx q[0];
rz(-0.80128765) q[0];
sx q[0];
rz(2.3961156) q[0];
rz(-1.2019134) q[1];
sx q[1];
rz(-1.8169836) q[1];
sx q[1];
rz(-2.1048996) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7734402) q[0];
sx q[0];
rz(-0.37640171) q[0];
sx q[0];
rz(1.3163811) q[0];
x q[1];
rz(1.5444504) q[2];
sx q[2];
rz(-0.80412946) q[2];
sx q[2];
rz(0.94601099) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17192852) q[1];
sx q[1];
rz(-1.8542037) q[1];
sx q[1];
rz(0.71782459) q[1];
rz(-pi) q[2];
rz(-0.58425957) q[3];
sx q[3];
rz(-1.3887032) q[3];
sx q[3];
rz(2.3626773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1408954) q[2];
sx q[2];
rz(-1.2709728) q[2];
sx q[2];
rz(-0.97266436) q[2];
rz(-0.85727143) q[3];
sx q[3];
rz(-2.075115) q[3];
sx q[3];
rz(1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.5001517) q[0];
sx q[0];
rz(-2.2730136) q[0];
sx q[0];
rz(2.3386173) q[0];
rz(2.4909486) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(-0.21779901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686435) q[0];
sx q[0];
rz(-1.9470707) q[0];
sx q[0];
rz(2.2479288) q[0];
x q[1];
rz(1.3399383) q[2];
sx q[2];
rz(-2.1862324) q[2];
sx q[2];
rz(0.99986693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4068702) q[1];
sx q[1];
rz(-2.0592505) q[1];
sx q[1];
rz(-0.10481701) q[1];
rz(-pi) q[2];
rz(2.5434615) q[3];
sx q[3];
rz(-0.39970988) q[3];
sx q[3];
rz(0.78228116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.190072) q[2];
sx q[2];
rz(-2.0900487) q[2];
sx q[2];
rz(1.6311084) q[2];
rz(-1.914628) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25201061) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(2.4776283) q[0];
rz(0.12531677) q[1];
sx q[1];
rz(-1.434451) q[1];
sx q[1];
rz(1.0391611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0200203) q[0];
sx q[0];
rz(-1.5735605) q[0];
sx q[0];
rz(-3.1317668) q[0];
rz(-pi) q[1];
rz(-1.9661994) q[2];
sx q[2];
rz(-1.8186453) q[2];
sx q[2];
rz(-1.445766) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3287828) q[1];
sx q[1];
rz(-2.610958) q[1];
sx q[1];
rz(-1.1483436) q[1];
rz(-pi) q[2];
rz(0.19840513) q[3];
sx q[3];
rz(-0.45061776) q[3];
sx q[3];
rz(3.0991447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8404954) q[2];
sx q[2];
rz(-1.7258464) q[2];
sx q[2];
rz(0.077979716) q[2];
rz(2.0237427) q[3];
sx q[3];
rz(-2.100914) q[3];
sx q[3];
rz(-1.0361205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2074821) q[0];
sx q[0];
rz(-0.20620646) q[0];
sx q[0];
rz(-2.7284486) q[0];
rz(-1.235599) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(-1.8280169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.006047) q[0];
sx q[0];
rz(-2.1925547) q[0];
sx q[0];
rz(0.22613392) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35982015) q[2];
sx q[2];
rz(-2.0096547) q[2];
sx q[2];
rz(3.1391337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6550555) q[1];
sx q[1];
rz(-1.948272) q[1];
sx q[1];
rz(2.0546163) q[1];
x q[2];
rz(0.78727874) q[3];
sx q[3];
rz(-1.9036999) q[3];
sx q[3];
rz(1.6407788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95973394) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(0.44446298) q[2];
rz(-1.2005165) q[3];
sx q[3];
rz(-1.2512755) q[3];
sx q[3];
rz(3.0357231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7926517) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(0.077202395) q[0];
rz(2.1791747) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(-0.033800689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8706521) q[0];
sx q[0];
rz(-2.0217289) q[0];
sx q[0];
rz(-1.1876039) q[0];
x q[1];
rz(-0.044174657) q[2];
sx q[2];
rz(-2.5934873) q[2];
sx q[2];
rz(1.7063315) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9056184) q[1];
sx q[1];
rz(-1.8840569) q[1];
sx q[1];
rz(-2.226023) q[1];
rz(-pi) q[2];
rz(0.77671364) q[3];
sx q[3];
rz(-1.4956105) q[3];
sx q[3];
rz(-3.1326339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17322156) q[2];
sx q[2];
rz(-1.5366448) q[2];
sx q[2];
rz(0.25924337) q[2];
rz(2.1974468) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060870085) q[0];
sx q[0];
rz(-2.935077) q[0];
sx q[0];
rz(-2.5282705) q[0];
rz(-2.2382286) q[1];
sx q[1];
rz(-2.3009243) q[1];
sx q[1];
rz(2.9936252) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6981993) q[0];
sx q[0];
rz(-0.86981378) q[0];
sx q[0];
rz(0.92595358) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5388707) q[2];
sx q[2];
rz(-2.0096001) q[2];
sx q[2];
rz(-1.9500458) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.24616775) q[1];
sx q[1];
rz(-2.2895799) q[1];
sx q[1];
rz(0.5753998) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0056173) q[3];
sx q[3];
rz(-1.5508754) q[3];
sx q[3];
rz(1.3764149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.938574) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(-0.99986783) q[2];
rz(1.3153007) q[3];
sx q[3];
rz(-2.857693) q[3];
sx q[3];
rz(2.4711173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95241791) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(-2.1536105) q[0];
rz(1.8723764) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(0.16407897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7994499) q[0];
sx q[0];
rz(-0.14524325) q[0];
sx q[0];
rz(2.5119315) q[0];
x q[1];
rz(-1.0751749) q[2];
sx q[2];
rz(-1.7640503) q[2];
sx q[2];
rz(2.9444864) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8089701) q[1];
sx q[1];
rz(-1.5676985) q[1];
sx q[1];
rz(-2.0634406) q[1];
rz(2.9470575) q[3];
sx q[3];
rz(-2.2119868) q[3];
sx q[3];
rz(1.3536782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8998731) q[2];
sx q[2];
rz(-0.49979979) q[2];
sx q[2];
rz(2.0797753) q[2];
rz(-0.1344943) q[3];
sx q[3];
rz(-2.2420292) q[3];
sx q[3];
rz(0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6580842) q[0];
sx q[0];
rz(-2.4299419) q[0];
sx q[0];
rz(0.2051951) q[0];
rz(0.23458734) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(-2.9451784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.043606) q[0];
sx q[0];
rz(-0.65904407) q[0];
sx q[0];
rz(-1.0511257) q[0];
x q[1];
rz(-1.4142198) q[2];
sx q[2];
rz(-2.0897667) q[2];
sx q[2];
rz(1.9760946) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82746303) q[1];
sx q[1];
rz(-1.2478831) q[1];
sx q[1];
rz(2.9548378) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78077448) q[3];
sx q[3];
rz(-1.9650243) q[3];
sx q[3];
rz(-3.0633283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3966169) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(1.9632001) q[2];
rz(-2.7922503) q[3];
sx q[3];
rz(-1.1290461) q[3];
sx q[3];
rz(1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1402533) q[0];
sx q[0];
rz(-1.1118735) q[0];
sx q[0];
rz(0.70571357) q[0];
rz(0.69752518) q[1];
sx q[1];
rz(-1.768528) q[1];
sx q[1];
rz(-2.2360905) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5030497) q[0];
sx q[0];
rz(-2.7075504) q[0];
sx q[0];
rz(0.040157138) q[0];
rz(-2.7780813) q[2];
sx q[2];
rz(-1.264252) q[2];
sx q[2];
rz(-0.53378045) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21852979) q[1];
sx q[1];
rz(-2.8273812) q[1];
sx q[1];
rz(-2.988222) q[1];
x q[2];
rz(-1.2071183) q[3];
sx q[3];
rz(-1.8810026) q[3];
sx q[3];
rz(2.9077173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2424348) q[2];
sx q[2];
rz(-2.0826715) q[2];
sx q[2];
rz(2.4230797) q[2];
rz(0.68273035) q[3];
sx q[3];
rz(-1.9971137) q[3];
sx q[3];
rz(3.0375286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.8231507) q[0];
sx q[0];
rz(-2.5704076) q[0];
sx q[0];
rz(3.0388863) q[0];
rz(1.5009343) q[1];
sx q[1];
rz(-1.2702912) q[1];
sx q[1];
rz(2.695695) q[1];
rz(0.14328843) q[2];
sx q[2];
rz(-1.454151) q[2];
sx q[2];
rz(3.1232338) q[2];
rz(0.80397687) q[3];
sx q[3];
rz(-1.0729651) q[3];
sx q[3];
rz(3.1333522) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
