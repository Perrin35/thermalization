OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2849046) q[0];
sx q[0];
rz(-1.8104799) q[0];
sx q[0];
rz(2.3321505) q[0];
rz(-2.5419905) q[1];
sx q[1];
rz(2.0166346) q[1];
sx q[1];
rz(12.044608) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3974575) q[0];
sx q[0];
rz(-1.43072) q[0];
sx q[0];
rz(3.1404353) q[0];
rz(-pi) q[1];
rz(2.6894301) q[2];
sx q[2];
rz(-2.5322862) q[2];
sx q[2];
rz(0.094560187) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.45568902) q[1];
sx q[1];
rz(-2.0200811) q[1];
sx q[1];
rz(-2.8390769) q[1];
rz(-pi) q[2];
rz(-1.6572324) q[3];
sx q[3];
rz(-1.0364136) q[3];
sx q[3];
rz(2.5680281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1521505) q[2];
sx q[2];
rz(-0.70789727) q[2];
sx q[2];
rz(-1.36261) q[2];
rz(1.4710434) q[3];
sx q[3];
rz(-1.5971284) q[3];
sx q[3];
rz(2.7267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3926369) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(2.0718527) q[0];
rz(-2.3906129) q[1];
sx q[1];
rz(-2.2396478) q[1];
sx q[1];
rz(0.78838563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716924) q[0];
sx q[0];
rz(-0.65061749) q[0];
sx q[0];
rz(0.31221892) q[0];
rz(-0.66547243) q[2];
sx q[2];
rz(-1.5147527) q[2];
sx q[2];
rz(1.9201345) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0505817) q[1];
sx q[1];
rz(-1.3772703) q[1];
sx q[1];
rz(2.0975608) q[1];
x q[2];
rz(-2.0072924) q[3];
sx q[3];
rz(-1.1009163) q[3];
sx q[3];
rz(1.3677471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.02701935) q[2];
sx q[2];
rz(-2.7832649) q[2];
sx q[2];
rz(-0.9066073) q[2];
rz(-1.00057) q[3];
sx q[3];
rz(-2.0228736) q[3];
sx q[3];
rz(2.6226543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24666102) q[0];
sx q[0];
rz(-0.38235679) q[0];
sx q[0];
rz(1.3038127) q[0];
rz(1.1600102) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(1.4132285) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055186836) q[0];
sx q[0];
rz(-2.9759183) q[0];
sx q[0];
rz(2.590986) q[0];
rz(2.0264852) q[2];
sx q[2];
rz(-2.1564061) q[2];
sx q[2];
rz(2.2024296) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0121514) q[1];
sx q[1];
rz(-2.0291162) q[1];
sx q[1];
rz(1.767551) q[1];
rz(-2.2858544) q[3];
sx q[3];
rz(-1.6659587) q[3];
sx q[3];
rz(2.9449938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7580737) q[2];
sx q[2];
rz(-1.5531837) q[2];
sx q[2];
rz(-0.12913945) q[2];
rz(0.50390759) q[3];
sx q[3];
rz(-1.7241071) q[3];
sx q[3];
rz(-2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1103766) q[0];
sx q[0];
rz(-1.9289368) q[0];
sx q[0];
rz(0.84621286) q[0];
rz(0.57202488) q[1];
sx q[1];
rz(-1.6366448) q[1];
sx q[1];
rz(1.9073073) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9548428) q[0];
sx q[0];
rz(-0.34663793) q[0];
sx q[0];
rz(0.28916547) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0364785) q[2];
sx q[2];
rz(-2.6200175) q[2];
sx q[2];
rz(-1.2273481) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14951359) q[1];
sx q[1];
rz(-1.7798335) q[1];
sx q[1];
rz(-1.4884218) q[1];
rz(-pi) q[2];
rz(-2.0999072) q[3];
sx q[3];
rz(-2.2941781) q[3];
sx q[3];
rz(0.46284562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7531551) q[2];
sx q[2];
rz(-0.93418241) q[2];
sx q[2];
rz(-1.6969121) q[2];
rz(-1.8709024) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-1.1893893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63950771) q[0];
sx q[0];
rz(-1.6725699) q[0];
sx q[0];
rz(1.0953267) q[0];
rz(0.0630088) q[1];
sx q[1];
rz(-1.4728225) q[1];
sx q[1];
rz(-0.94620401) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5216833) q[0];
sx q[0];
rz(-0.54542002) q[0];
sx q[0];
rz(-1.5167461) q[0];
rz(-pi) q[1];
rz(0.49503742) q[2];
sx q[2];
rz(-0.69827852) q[2];
sx q[2];
rz(-0.78294047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.19095382) q[1];
sx q[1];
rz(-2.6304448) q[1];
sx q[1];
rz(-0.010015476) q[1];
rz(-1.1796451) q[3];
sx q[3];
rz(-1.3538651) q[3];
sx q[3];
rz(2.2892078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63518628) q[2];
sx q[2];
rz(-1.0794285) q[2];
sx q[2];
rz(1.5080473) q[2];
rz(-0.15277282) q[3];
sx q[3];
rz(-1.2341876) q[3];
sx q[3];
rz(0.93301409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2372811) q[0];
sx q[0];
rz(-2.5560684) q[0];
sx q[0];
rz(-0.10865077) q[0];
rz(0.12734224) q[1];
sx q[1];
rz(-1.2510977) q[1];
sx q[1];
rz(-2.2551575) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3868635) q[0];
sx q[0];
rz(-2.6219256) q[0];
sx q[0];
rz(-3.1082252) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8511925) q[2];
sx q[2];
rz(-0.2125769) q[2];
sx q[2];
rz(-2.5857298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.34102124) q[1];
sx q[1];
rz(-1.6280975) q[1];
sx q[1];
rz(-1.351746) q[1];
x q[2];
rz(2.5648067) q[3];
sx q[3];
rz(-2.9401719) q[3];
sx q[3];
rz(-2.6164428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.844937) q[2];
sx q[2];
rz(-2.3462494) q[2];
sx q[2];
rz(-0.33357683) q[2];
rz(-2.2070456) q[3];
sx q[3];
rz(-2.3881674) q[3];
sx q[3];
rz(0.88753382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4535256) q[0];
sx q[0];
rz(-1.6418566) q[0];
sx q[0];
rz(0.30853477) q[0];
rz(0.60641369) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(-0.44953129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21643695) q[0];
sx q[0];
rz(-1.0569658) q[0];
sx q[0];
rz(1.4975182) q[0];
rz(0.78747998) q[2];
sx q[2];
rz(-1.228294) q[2];
sx q[2];
rz(2.1975348) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14913371) q[1];
sx q[1];
rz(-2.1861052) q[1];
sx q[1];
rz(-3.0165218) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1949439) q[3];
sx q[3];
rz(-0.80393857) q[3];
sx q[3];
rz(-0.80292668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0799847) q[2];
sx q[2];
rz(-1.4915165) q[2];
sx q[2];
rz(1.1770581) q[2];
rz(-2.4401149) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(-0.85566068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.6845282) q[0];
sx q[0];
rz(-0.51790154) q[0];
sx q[0];
rz(-1.49217) q[0];
rz(1.0555142) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(-2.7141056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9857835) q[0];
sx q[0];
rz(-1.0571348) q[0];
sx q[0];
rz(0.86648057) q[0];
x q[1];
rz(2.7107377) q[2];
sx q[2];
rz(-1.4240032) q[2];
sx q[2];
rz(1.0041179) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52054469) q[1];
sx q[1];
rz(-1.2579903) q[1];
sx q[1];
rz(-1.5473844) q[1];
x q[2];
rz(-0.31510809) q[3];
sx q[3];
rz(-1.8297394) q[3];
sx q[3];
rz(-1.068598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.18275729) q[2];
sx q[2];
rz(-1.8104825) q[2];
sx q[2];
rz(-1.9367564) q[2];
rz(-0.13293535) q[3];
sx q[3];
rz(-0.098630579) q[3];
sx q[3];
rz(-1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8590915) q[0];
sx q[0];
rz(-1.4097255) q[0];
sx q[0];
rz(2.6891563) q[0];
rz(0.85894194) q[1];
sx q[1];
rz(-2.5622538) q[1];
sx q[1];
rz(1.4525684) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4389324) q[0];
sx q[0];
rz(-2.0468674) q[0];
sx q[0];
rz(2.4804546) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8088687) q[2];
sx q[2];
rz(-1.0411658) q[2];
sx q[2];
rz(-0.67654787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7048129) q[1];
sx q[1];
rz(-1.5680934) q[1];
sx q[1];
rz(1.0648492) q[1];
rz(-2.7730016) q[3];
sx q[3];
rz(-1.2342216) q[3];
sx q[3];
rz(-0.14813885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9026044) q[2];
sx q[2];
rz(-1.8391823) q[2];
sx q[2];
rz(2.7461309) q[2];
rz(2.0579193) q[3];
sx q[3];
rz(-0.62267059) q[3];
sx q[3];
rz(-1.4130939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.85635066) q[0];
sx q[0];
rz(-0.1150035) q[0];
sx q[0];
rz(-1.850199) q[0];
rz(-2.3639288) q[1];
sx q[1];
rz(-1.5592557) q[1];
sx q[1];
rz(1.2819598) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7972231) q[0];
sx q[0];
rz(-0.31569052) q[0];
sx q[0];
rz(-0.80202054) q[0];
x q[1];
rz(0.28996946) q[2];
sx q[2];
rz(-2.4735138) q[2];
sx q[2];
rz(1.5529322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29867783) q[1];
sx q[1];
rz(-2.5576572) q[1];
sx q[1];
rz(0.42386492) q[1];
rz(-pi) q[2];
rz(3.0324803) q[3];
sx q[3];
rz(-2.0874573) q[3];
sx q[3];
rz(-1.9369879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.97734863) q[2];
sx q[2];
rz(-0.35173309) q[2];
sx q[2];
rz(-2.742761) q[2];
rz(-2.2073958) q[3];
sx q[3];
rz(-0.7312921) q[3];
sx q[3];
rz(3.0493693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4451404) q[0];
sx q[0];
rz(-2.2284989) q[0];
sx q[0];
rz(-3.1365119) q[0];
rz(2.483881) q[1];
sx q[1];
rz(-1.7291768) q[1];
sx q[1];
rz(-1.9531858) q[1];
rz(0.80398139) q[2];
sx q[2];
rz(-1.8009427) q[2];
sx q[2];
rz(1.367205) q[2];
rz(1.6647958) q[3];
sx q[3];
rz(-0.58782676) q[3];
sx q[3];
rz(1.6997433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
