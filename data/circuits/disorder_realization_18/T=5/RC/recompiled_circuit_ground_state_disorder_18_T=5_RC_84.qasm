OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9590149) q[0];
sx q[0];
rz(0.20354095) q[0];
sx q[0];
rz(8.1770906) q[0];
rz(-1.7893451) q[1];
sx q[1];
rz(-2.4440553) q[1];
sx q[1];
rz(0.69812671) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20839918) q[0];
sx q[0];
rz(-2.0740447) q[0];
sx q[0];
rz(2.0397951) q[0];
x q[1];
rz(-0.013555992) q[2];
sx q[2];
rz(-2.6883311) q[2];
sx q[2];
rz(-2.1238985) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.082939741) q[1];
sx q[1];
rz(-1.2695244) q[1];
sx q[1];
rz(0.79853398) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1269253) q[3];
sx q[3];
rz(-1.3157744) q[3];
sx q[3];
rz(1.6153112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84833604) q[2];
sx q[2];
rz(-1.2245327) q[2];
sx q[2];
rz(-0.87948925) q[2];
rz(0.73421156) q[3];
sx q[3];
rz(-1.1314393) q[3];
sx q[3];
rz(0.82936275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64049983) q[0];
sx q[0];
rz(-1.9556029) q[0];
sx q[0];
rz(-0.99047852) q[0];
rz(-2.1462006) q[1];
sx q[1];
rz(-1.6973015) q[1];
sx q[1];
rz(-2.676414) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.746691) q[0];
sx q[0];
rz(-1.1550265) q[0];
sx q[0];
rz(-1.7703992) q[0];
rz(-0.37613531) q[2];
sx q[2];
rz(-2.1636164) q[2];
sx q[2];
rz(1.3786292) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0439321) q[1];
sx q[1];
rz(-2.2273835) q[1];
sx q[1];
rz(-2.2022579) q[1];
rz(-pi) q[2];
rz(0.2618913) q[3];
sx q[3];
rz(-1.3757816) q[3];
sx q[3];
rz(3.0660925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5140932) q[2];
sx q[2];
rz(-1.2343312) q[2];
sx q[2];
rz(-2.9631183) q[2];
rz(-0.56525362) q[3];
sx q[3];
rz(-1.3924799) q[3];
sx q[3];
rz(0.55725151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8239215) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(0.70096651) q[0];
rz(0.40755454) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(2.0661381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385745) q[0];
sx q[0];
rz(-2.10963) q[0];
sx q[0];
rz(-1.1807022) q[0];
rz(-pi) q[1];
rz(1.2811866) q[2];
sx q[2];
rz(-2.4649005) q[2];
sx q[2];
rz(-1.8829568) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0044864) q[1];
sx q[1];
rz(-1.1710694) q[1];
sx q[1];
rz(-2.9997679) q[1];
x q[2];
rz(-1.7549137) q[3];
sx q[3];
rz(-0.26911165) q[3];
sx q[3];
rz(-1.5003913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0740697) q[2];
sx q[2];
rz(-1.3347551) q[2];
sx q[2];
rz(-1.4074116) q[2];
rz(-0.47809005) q[3];
sx q[3];
rz(-2.7985088) q[3];
sx q[3];
rz(2.7966255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5673229) q[0];
sx q[0];
rz(-1.747921) q[0];
sx q[0];
rz(2.0950914) q[0];
rz(2.8872755) q[1];
sx q[1];
rz(-0.81552243) q[1];
sx q[1];
rz(0.3124803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3182718) q[0];
sx q[0];
rz(-1.4411363) q[0];
sx q[0];
rz(1.983485) q[0];
rz(-pi) q[1];
rz(-2.0489659) q[2];
sx q[2];
rz(-0.96726438) q[2];
sx q[2];
rz(-0.1665394) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29280832) q[1];
sx q[1];
rz(-0.40221805) q[1];
sx q[1];
rz(-2.4821698) q[1];
rz(-pi) q[2];
rz(-2.813597) q[3];
sx q[3];
rz(-2.061279) q[3];
sx q[3];
rz(-2.1759263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8326524) q[2];
sx q[2];
rz(-2.4231484) q[2];
sx q[2];
rz(0.18860513) q[2];
rz(-0.23379937) q[3];
sx q[3];
rz(-0.89582396) q[3];
sx q[3];
rz(-0.58376592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4466062) q[0];
sx q[0];
rz(-1.3124895) q[0];
sx q[0];
rz(0.0083228668) q[0];
rz(-2.7024929) q[1];
sx q[1];
rz(-1.3689901) q[1];
sx q[1];
rz(1.8956553) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3942892) q[0];
sx q[0];
rz(-2.6458594) q[0];
sx q[0];
rz(0.41252211) q[0];
rz(-pi) q[1];
rz(-1.8077085) q[2];
sx q[2];
rz(-1.7157752) q[2];
sx q[2];
rz(-2.0857834) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8716395) q[1];
sx q[1];
rz(-1.9325745) q[1];
sx q[1];
rz(-0.33906083) q[1];
x q[2];
rz(0.28723081) q[3];
sx q[3];
rz(-0.93614791) q[3];
sx q[3];
rz(-1.3040445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40450725) q[2];
sx q[2];
rz(-0.046701996) q[2];
sx q[2];
rz(1.6339634) q[2];
rz(2.5493933) q[3];
sx q[3];
rz(-1.7630354) q[3];
sx q[3];
rz(-0.73970214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020788766) q[0];
sx q[0];
rz(-1.4494267) q[0];
sx q[0];
rz(2.8756496) q[0];
rz(1.9256915) q[1];
sx q[1];
rz(-2.3561056) q[1];
sx q[1];
rz(-1.9507834) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5191089) q[0];
sx q[0];
rz(-1.5734696) q[0];
sx q[0];
rz(2.3038505) q[0];
x q[1];
rz(-1.6295678) q[2];
sx q[2];
rz(-2.588039) q[2];
sx q[2];
rz(-2.2433081) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9126351) q[1];
sx q[1];
rz(-1.5630782) q[1];
sx q[1];
rz(-0.75386366) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5134638) q[3];
sx q[3];
rz(-1.2964037) q[3];
sx q[3];
rz(-0.17255863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8094981) q[2];
sx q[2];
rz(-1.3846493) q[2];
sx q[2];
rz(0.14796999) q[2];
rz(-0.45251265) q[3];
sx q[3];
rz(-1.0871525) q[3];
sx q[3];
rz(-0.40542671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0026523503) q[0];
sx q[0];
rz(-1.8754706) q[0];
sx q[0];
rz(1.7953605) q[0];
rz(3.1345308) q[1];
sx q[1];
rz(-0.66771737) q[1];
sx q[1];
rz(0.5074358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1088828) q[0];
sx q[0];
rz(-1.0531972) q[0];
sx q[0];
rz(0.93982307) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9865737) q[2];
sx q[2];
rz(-1.2298202) q[2];
sx q[2];
rz(0.91779681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.145157) q[1];
sx q[1];
rz(-1.5801589) q[1];
sx q[1];
rz(-2.7476176) q[1];
rz(-3.1321615) q[3];
sx q[3];
rz(-1.2985029) q[3];
sx q[3];
rz(1.997153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4355882) q[2];
sx q[2];
rz(-2.0861349) q[2];
sx q[2];
rz(-2.9662507) q[2];
rz(-0.38241479) q[3];
sx q[3];
rz(-1.1924815) q[3];
sx q[3];
rz(-0.46149883) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14313702) q[0];
sx q[0];
rz(-2.1303506) q[0];
sx q[0];
rz(1.3125516) q[0];
rz(-3.0328499) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(0.67799062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4512166) q[0];
sx q[0];
rz(-1.548598) q[0];
sx q[0];
rz(2.064384) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5505232) q[2];
sx q[2];
rz(-0.82189362) q[2];
sx q[2];
rz(-0.66416364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99876632) q[1];
sx q[1];
rz(-1.2825535) q[1];
sx q[1];
rz(-1.0267797) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1021712) q[3];
sx q[3];
rz(-2.2153004) q[3];
sx q[3];
rz(2.927305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8678191) q[2];
sx q[2];
rz(-1.5742233) q[2];
sx q[2];
rz(1.1567953) q[2];
rz(-0.50446883) q[3];
sx q[3];
rz(-1.0320832) q[3];
sx q[3];
rz(0.57197905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2061763) q[0];
sx q[0];
rz(-1.6125866) q[0];
sx q[0];
rz(2.5231498) q[0];
rz(0.84856021) q[1];
sx q[1];
rz(-2.6237374) q[1];
sx q[1];
rz(2.3458164) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2413809) q[0];
sx q[0];
rz(-2.6707895) q[0];
sx q[0];
rz(-1.6241252) q[0];
x q[1];
rz(1.6219989) q[2];
sx q[2];
rz(-2.3948811) q[2];
sx q[2];
rz(-1.8987617) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91589815) q[1];
sx q[1];
rz(-1.2978329) q[1];
sx q[1];
rz(1.2429487) q[1];
x q[2];
rz(1.2451781) q[3];
sx q[3];
rz(-0.44105704) q[3];
sx q[3];
rz(2.8056322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15647469) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(1.1953243) q[2];
rz(0.34504238) q[3];
sx q[3];
rz(-2.2796567) q[3];
sx q[3];
rz(-2.2752458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7095551) q[0];
sx q[0];
rz(-0.23831743) q[0];
sx q[0];
rz(-3.0639783) q[0];
rz(1.9084825) q[1];
sx q[1];
rz(-2.1014919) q[1];
sx q[1];
rz(-0.79201039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38538489) q[0];
sx q[0];
rz(-2.1127708) q[0];
sx q[0];
rz(-2.3535924) q[0];
rz(-pi) q[1];
rz(1.7260688) q[2];
sx q[2];
rz(-1.7285992) q[2];
sx q[2];
rz(-1.4580245) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.92160901) q[1];
sx q[1];
rz(-1.7243544) q[1];
sx q[1];
rz(1.1034225) q[1];
x q[2];
rz(-0.11802063) q[3];
sx q[3];
rz(-0.80720085) q[3];
sx q[3];
rz(-1.0420052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21976694) q[2];
sx q[2];
rz(-2.2302088) q[2];
sx q[2];
rz(-2.0683973) q[2];
rz(2.8940767) q[3];
sx q[3];
rz(-1.9297587) q[3];
sx q[3];
rz(3.1246576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76219227) q[0];
sx q[0];
rz(-1.8397377) q[0];
sx q[0];
rz(2.4158438) q[0];
rz(-0.87860592) q[1];
sx q[1];
rz(-0.85522063) q[1];
sx q[1];
rz(-1.9488889) q[1];
rz(1.0868308) q[2];
sx q[2];
rz(-1.7499583) q[2];
sx q[2];
rz(-1.3365895) q[2];
rz(3.0847753) q[3];
sx q[3];
rz(-0.84470261) q[3];
sx q[3];
rz(-0.62538996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
