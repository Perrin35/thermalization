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
rz(2.8853631) q[0];
rz(0.37707075) q[1];
sx q[1];
rz(-1.7989676) q[1];
sx q[1];
rz(-2.0816198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9794344) q[0];
sx q[0];
rz(-0.16955626) q[0];
sx q[0];
rz(2.3441699) q[0];
x q[1];
rz(0.04214123) q[2];
sx q[2];
rz(-2.4486447) q[2];
sx q[2];
rz(2.9646089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6457451) q[1];
sx q[1];
rz(-1.331067) q[1];
sx q[1];
rz(1.354753) q[1];
rz(-pi) q[2];
rz(0.59935244) q[3];
sx q[3];
rz(-2.9671892) q[3];
sx q[3];
rz(-0.50267437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18067351) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(1.712435) q[2];
rz(-0.6116496) q[3];
sx q[3];
rz(-1.6983756) q[3];
sx q[3];
rz(-3.070224) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83577689) q[0];
sx q[0];
rz(-0.80128765) q[0];
sx q[0];
rz(-2.3961156) q[0];
rz(1.2019134) q[1];
sx q[1];
rz(-1.8169836) q[1];
sx q[1];
rz(2.1048996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03464493) q[0];
sx q[0];
rz(-1.6634403) q[0];
sx q[0];
rz(-1.205446) q[0];
rz(0.76684029) q[2];
sx q[2];
rz(-1.5518223) q[2];
sx q[2];
rz(2.5350867) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9696641) q[1];
sx q[1];
rz(-1.2873889) q[1];
sx q[1];
rz(-0.71782459) q[1];
rz(-pi) q[2];
rz(2.8193982) q[3];
sx q[3];
rz(-2.5327842) q[3];
sx q[3];
rz(2.6170066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1408954) q[2];
sx q[2];
rz(-1.8706198) q[2];
sx q[2];
rz(2.1689283) q[2];
rz(0.85727143) q[3];
sx q[3];
rz(-2.075115) q[3];
sx q[3];
rz(1.2681786) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5001517) q[0];
sx q[0];
rz(-0.86857906) q[0];
sx q[0];
rz(0.80297536) q[0];
rz(-2.4909486) q[1];
sx q[1];
rz(-2.3425808) q[1];
sx q[1];
rz(0.21779901) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686435) q[0];
sx q[0];
rz(-1.1945219) q[0];
sx q[0];
rz(-2.2479288) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3399383) q[2];
sx q[2];
rz(-2.1862324) q[2];
sx q[2];
rz(2.1417257) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9552635) q[1];
sx q[1];
rz(-0.4986838) q[1];
sx q[1];
rz(-1.7651943) q[1];
rz(1.337255) q[3];
sx q[3];
rz(-1.8982049) q[3];
sx q[3];
rz(-1.722432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.190072) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(1.6311084) q[2];
rz(1.2269646) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(-1.0043043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.889582) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(0.66396436) q[0];
rz(0.12531677) q[1];
sx q[1];
rz(-1.434451) q[1];
sx q[1];
rz(1.0391611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17500278) q[0];
sx q[0];
rz(-0.010207264) q[0];
sx q[0];
rz(0.27423476) q[0];
x q[1];
rz(1.9661994) q[2];
sx q[2];
rz(-1.8186453) q[2];
sx q[2];
rz(1.445766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3323556) q[1];
sx q[1];
rz(-1.0910104) q[1];
sx q[1];
rz(-2.9055041) q[1];
x q[2];
rz(-0.19840513) q[3];
sx q[3];
rz(-2.6909749) q[3];
sx q[3];
rz(-0.042447986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8404954) q[2];
sx q[2];
rz(-1.4157462) q[2];
sx q[2];
rz(3.0636129) q[2];
rz(-1.1178499) q[3];
sx q[3];
rz(-2.100914) q[3];
sx q[3];
rz(2.1054721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9341105) q[0];
sx q[0];
rz(-0.20620646) q[0];
sx q[0];
rz(0.41314405) q[0];
rz(1.235599) q[1];
sx q[1];
rz(-1.8159591) q[1];
sx q[1];
rz(1.8280169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13554561) q[0];
sx q[0];
rz(-0.94903799) q[0];
sx q[0];
rz(-2.9154587) q[0];
rz(-pi) q[1];
rz(1.1059472) q[2];
sx q[2];
rz(-1.8951956) q[2];
sx q[2];
rz(1.731763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8660248) q[1];
sx q[1];
rz(-1.1236262) q[1];
sx q[1];
rz(-0.42110301) q[1];
rz(0.4540654) q[3];
sx q[3];
rz(-0.84058529) q[3];
sx q[3];
rz(0.24475748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.95973394) q[2];
sx q[2];
rz(-1.9071969) q[2];
sx q[2];
rz(-2.6971297) q[2];
rz(1.2005165) q[3];
sx q[3];
rz(-1.8903172) q[3];
sx q[3];
rz(-0.1058696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34894094) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(0.077202395) q[0];
rz(-2.1791747) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(-3.107792) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8706521) q[0];
sx q[0];
rz(-1.1198638) q[0];
sx q[0];
rz(-1.1876039) q[0];
x q[1];
rz(0.54767139) q[2];
sx q[2];
rz(-1.5477836) q[2];
sx q[2];
rz(0.097824899) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2359743) q[1];
sx q[1];
rz(-1.8840569) q[1];
sx q[1];
rz(2.226023) q[1];
x q[2];
rz(0.77671364) q[3];
sx q[3];
rz(-1.4956105) q[3];
sx q[3];
rz(0.0089587072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17322156) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(0.25924337) q[2];
rz(-0.94414583) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(1.3313782) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060870085) q[0];
sx q[0];
rz(-2.935077) q[0];
sx q[0];
rz(-2.5282705) q[0];
rz(2.2382286) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(-0.14796743) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3053235) q[0];
sx q[0];
rz(-0.9137872) q[0];
sx q[0];
rz(0.61886529) q[0];
x q[1];
rz(1.5388707) q[2];
sx q[2];
rz(-1.1319926) q[2];
sx q[2];
rz(-1.9500458) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1175673) q[1];
sx q[1];
rz(-0.88741517) q[1];
sx q[1];
rz(2.1272117) q[1];
x q[2];
rz(-2.1359753) q[3];
sx q[3];
rz(-1.5907173) q[3];
sx q[3];
rz(-1.7651778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2030187) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(-0.99986783) q[2];
rz(-1.8262919) q[3];
sx q[3];
rz(-2.857693) q[3];
sx q[3];
rz(2.4711173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.1891747) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(-2.1536105) q[0];
rz(-1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(0.16407897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9768499) q[0];
sx q[0];
rz(-1.4535507) q[0];
sx q[0];
rz(-1.4848723) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21890004) q[2];
sx q[2];
rz(-1.0852119) q[2];
sx q[2];
rz(-1.6644434) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8976479) q[1];
sx q[1];
rz(-2.6489394) q[1];
sx q[1];
rz(1.5773462) q[1];
rz(-pi) q[2];
rz(-2.9470575) q[3];
sx q[3];
rz(-2.2119868) q[3];
sx q[3];
rz(-1.3536782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8998731) q[2];
sx q[2];
rz(-0.49979979) q[2];
sx q[2];
rz(-2.0797753) q[2];
rz(-3.0070983) q[3];
sx q[3];
rz(-2.2420292) q[3];
sx q[3];
rz(-0.053248052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6580842) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(-2.9363976) q[0];
rz(-2.9070053) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(0.1964143) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0935555) q[0];
sx q[0];
rz(-1.8797848) q[0];
sx q[0];
rz(0.97889203) q[0];
rz(0.26650776) q[2];
sx q[2];
rz(-0.54000137) q[2];
sx q[2];
rz(-0.85734474) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80322335) q[1];
sx q[1];
rz(-1.3937989) q[1];
sx q[1];
rz(-1.8990252) q[1];
rz(-pi) q[2];
rz(2.3608182) q[3];
sx q[3];
rz(-1.1765683) q[3];
sx q[3];
rz(-3.0633283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3966169) q[2];
sx q[2];
rz(-1.9605512) q[2];
sx q[2];
rz(1.9632001) q[2];
rz(-2.7922503) q[3];
sx q[3];
rz(-2.0125466) q[3];
sx q[3];
rz(-1.0497302) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1402533) q[0];
sx q[0];
rz(-1.1118735) q[0];
sx q[0];
rz(0.70571357) q[0];
rz(-0.69752518) q[1];
sx q[1];
rz(-1.768528) q[1];
sx q[1];
rz(2.2360905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6827992) q[0];
sx q[0];
rz(-1.1371277) q[0];
sx q[0];
rz(-1.5894029) q[0];
x q[1];
rz(0.72737965) q[2];
sx q[2];
rz(-0.47107163) q[2];
sx q[2];
rz(-1.4339652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6433555) q[1];
sx q[1];
rz(-1.61803) q[1];
sx q[1];
rz(0.31075732) q[1];
rz(-pi) q[2];
rz(1.9344744) q[3];
sx q[3];
rz(-1.26059) q[3];
sx q[3];
rz(0.23387533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8991578) q[2];
sx q[2];
rz(-2.0826715) q[2];
sx q[2];
rz(2.4230797) q[2];
rz(-0.68273035) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-0.1040641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.318442) q[0];
sx q[0];
rz(-0.57118509) q[0];
sx q[0];
rz(-0.1027064) q[0];
rz(-1.5009343) q[1];
sx q[1];
rz(-1.8713015) q[1];
sx q[1];
rz(-0.4458977) q[1];
rz(-0.14328843) q[2];
sx q[2];
rz(-1.6874416) q[2];
sx q[2];
rz(-0.018358827) q[2];
rz(0.90632306) q[3];
sx q[3];
rz(-0.88574468) q[3];
sx q[3];
rz(1.1024324) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
