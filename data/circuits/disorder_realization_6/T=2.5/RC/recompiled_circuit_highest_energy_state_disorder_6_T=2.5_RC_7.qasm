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
rz(1.2441128) q[0];
sx q[0];
rz(-2.8219858) q[0];
sx q[0];
rz(0.50080103) q[0];
rz(0.29419857) q[1];
sx q[1];
rz(3.8837641) q[1];
sx q[1];
rz(9.972218) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1337903) q[0];
sx q[0];
rz(-2.7977075) q[0];
sx q[0];
rz(-2.6965428) q[0];
rz(-pi) q[1];
rz(2.6260761) q[2];
sx q[2];
rz(-0.52243865) q[2];
sx q[2];
rz(1.9439657) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3708122) q[1];
sx q[1];
rz(-2.6264084) q[1];
sx q[1];
rz(0.18411247) q[1];
x q[2];
rz(-2.1087113) q[3];
sx q[3];
rz(-1.8807372) q[3];
sx q[3];
rz(-0.021837318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1529634) q[2];
sx q[2];
rz(-2.2965778) q[2];
sx q[2];
rz(2.5353234) q[2];
rz(1.9340949) q[3];
sx q[3];
rz(-1.2946125) q[3];
sx q[3];
rz(1.5906364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9015053) q[0];
sx q[0];
rz(-1.5272239) q[0];
sx q[0];
rz(0.69183451) q[0];
rz(1.8531307) q[1];
sx q[1];
rz(-1.9970147) q[1];
sx q[1];
rz(2.8537234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8765204) q[0];
sx q[0];
rz(-1.8567883) q[0];
sx q[0];
rz(-3.0080481) q[0];
rz(-pi) q[1];
rz(-0.98357304) q[2];
sx q[2];
rz(-2.1701239) q[2];
sx q[2];
rz(2.9119208) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1296245) q[1];
sx q[1];
rz(-1.5245364) q[1];
sx q[1];
rz(-2.6955626) q[1];
rz(-pi) q[2];
rz(-0.13707844) q[3];
sx q[3];
rz(-1.6774639) q[3];
sx q[3];
rz(2.9003473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79491478) q[2];
sx q[2];
rz(-1.4823464) q[2];
sx q[2];
rz(-2.1686926) q[2];
rz(1.8074624) q[3];
sx q[3];
rz(-2.3502246) q[3];
sx q[3];
rz(2.5202675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7767104) q[0];
sx q[0];
rz(-2.8948687) q[0];
sx q[0];
rz(0.2488939) q[0];
rz(-1.4502672) q[1];
sx q[1];
rz(-1.9906809) q[1];
sx q[1];
rz(-2.2379564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8758276) q[0];
sx q[0];
rz(-0.96185124) q[0];
sx q[0];
rz(-1.5399571) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9157639) q[2];
sx q[2];
rz(-1.6716701) q[2];
sx q[2];
rz(0.28988923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1978479) q[1];
sx q[1];
rz(-1.1863803) q[1];
sx q[1];
rz(0.77898394) q[1];
x q[2];
rz(2.5375705) q[3];
sx q[3];
rz(-1.9088749) q[3];
sx q[3];
rz(-0.018751831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9768208) q[2];
sx q[2];
rz(-0.73126078) q[2];
sx q[2];
rz(2.8110647) q[2];
rz(-2.9183689) q[3];
sx q[3];
rz(-2.0285172) q[3];
sx q[3];
rz(2.7631675) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6966003) q[0];
sx q[0];
rz(-1.7855423) q[0];
sx q[0];
rz(-0.14183216) q[0];
rz(-0.089135535) q[1];
sx q[1];
rz(-0.24239692) q[1];
sx q[1];
rz(-2.1884122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0713816) q[0];
sx q[0];
rz(-1.1888904) q[0];
sx q[0];
rz(2.3863263) q[0];
x q[1];
rz(-0.29229477) q[2];
sx q[2];
rz(-2.7461548) q[2];
sx q[2];
rz(-1.866445) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2948896) q[1];
sx q[1];
rz(-1.0422455) q[1];
sx q[1];
rz(-2.5336086) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7715376) q[3];
sx q[3];
rz(-1.0827409) q[3];
sx q[3];
rz(-2.1303915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3700833) q[2];
sx q[2];
rz(-0.13088317) q[2];
sx q[2];
rz(2.6123602) q[2];
rz(-1.0445163) q[3];
sx q[3];
rz(-1.7525201) q[3];
sx q[3];
rz(2.9862826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9411377) q[0];
sx q[0];
rz(-1.2417685) q[0];
sx q[0];
rz(0.45503765) q[0];
rz(3.1269238) q[1];
sx q[1];
rz(-2.585482) q[1];
sx q[1];
rz(-2.1591149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.608484) q[0];
sx q[0];
rz(-0.31705515) q[0];
sx q[0];
rz(-0.32755537) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.234794) q[2];
sx q[2];
rz(-2.325146) q[2];
sx q[2];
rz(1.0524257) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8807672) q[1];
sx q[1];
rz(-1.2952571) q[1];
sx q[1];
rz(1.6795621) q[1];
rz(1.807688) q[3];
sx q[3];
rz(-1.8474527) q[3];
sx q[3];
rz(-1.93678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5386397) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(0.61720103) q[2];
rz(1.7547102) q[3];
sx q[3];
rz(-0.91104561) q[3];
sx q[3];
rz(-0.174218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0844326) q[0];
sx q[0];
rz(-0.76157695) q[0];
sx q[0];
rz(-2.1630951) q[0];
rz(-1.017336) q[1];
sx q[1];
rz(-0.2778191) q[1];
sx q[1];
rz(0.4496347) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89116865) q[0];
sx q[0];
rz(-2.7319067) q[0];
sx q[0];
rz(-0.72575672) q[0];
x q[1];
rz(-2.041018) q[2];
sx q[2];
rz(-2.3210249) q[2];
sx q[2];
rz(2.5448397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75737917) q[1];
sx q[1];
rz(-1.5259201) q[1];
sx q[1];
rz(2.0424538) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68296229) q[3];
sx q[3];
rz(-1.0691081) q[3];
sx q[3];
rz(1.3987808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6189239) q[2];
sx q[2];
rz(-1.7974682) q[2];
sx q[2];
rz(-1.3783003) q[2];
rz(1.1528692) q[3];
sx q[3];
rz(-1.9882354) q[3];
sx q[3];
rz(2.8809179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.5744837) q[0];
sx q[0];
rz(-0.672988) q[0];
sx q[0];
rz(-2.8048977) q[0];
rz(2.4163272) q[1];
sx q[1];
rz(-1.5545574) q[1];
sx q[1];
rz(2.7560962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1565026) q[0];
sx q[0];
rz(-0.96961951) q[0];
sx q[0];
rz(1.1416092) q[0];
rz(2.9988937) q[2];
sx q[2];
rz(-1.8178504) q[2];
sx q[2];
rz(2.3076535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25076097) q[1];
sx q[1];
rz(-2.4422993) q[1];
sx q[1];
rz(-1.1108524) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9970421) q[3];
sx q[3];
rz(-0.9916456) q[3];
sx q[3];
rz(-0.338011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54726279) q[2];
sx q[2];
rz(-1.2196093) q[2];
sx q[2];
rz(0.29898137) q[2];
rz(0.79214054) q[3];
sx q[3];
rz(-2.7463089) q[3];
sx q[3];
rz(1.3778752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91656536) q[0];
sx q[0];
rz(-1.7074317) q[0];
sx q[0];
rz(2.3634971) q[0];
rz(-1.3968702) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(0.088833749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9561712) q[0];
sx q[0];
rz(-1.9156792) q[0];
sx q[0];
rz(-1.6817143) q[0];
x q[1];
rz(-2.2065032) q[2];
sx q[2];
rz(-0.82391058) q[2];
sx q[2];
rz(1.3694135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1792887) q[1];
sx q[1];
rz(-1.8306762) q[1];
sx q[1];
rz(1.13729) q[1];
x q[2];
rz(3.1109654) q[3];
sx q[3];
rz(-1.2032685) q[3];
sx q[3];
rz(-0.54545975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5953411) q[2];
sx q[2];
rz(-0.83117861) q[2];
sx q[2];
rz(2.0383932) q[2];
rz(2.5698419) q[3];
sx q[3];
rz(-0.56643707) q[3];
sx q[3];
rz(-1.6167238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53200805) q[0];
sx q[0];
rz(-3.0201865) q[0];
sx q[0];
rz(-0.65015656) q[0];
rz(-1.0981759) q[1];
sx q[1];
rz(-1.8868586) q[1];
sx q[1];
rz(0.065902725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9377334) q[0];
sx q[0];
rz(-0.90918505) q[0];
sx q[0];
rz(2.1801722) q[0];
rz(-pi) q[1];
rz(-2.8464855) q[2];
sx q[2];
rz(-2.9431651) q[2];
sx q[2];
rz(2.0917542) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3014631) q[1];
sx q[1];
rz(-2.3345229) q[1];
sx q[1];
rz(2.927344) q[1];
rz(0.66146694) q[3];
sx q[3];
rz(-1.3693878) q[3];
sx q[3];
rz(3.0440743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1419475) q[2];
sx q[2];
rz(-1.1885175) q[2];
sx q[2];
rz(0.56335062) q[2];
rz(0.79936409) q[3];
sx q[3];
rz(-1.0146217) q[3];
sx q[3];
rz(-2.8857901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946741) q[0];
sx q[0];
rz(-0.020337157) q[0];
sx q[0];
rz(0.49816966) q[0];
rz(2.8794471) q[1];
sx q[1];
rz(-1.405895) q[1];
sx q[1];
rz(-1.795305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.41723) q[0];
sx q[0];
rz(-2.2688393) q[0];
sx q[0];
rz(2.5265187) q[0];
rz(-pi) q[1];
rz(-1.5602697) q[2];
sx q[2];
rz(-1.5731647) q[2];
sx q[2];
rz(0.51711581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.87835364) q[1];
sx q[1];
rz(-1.2290956) q[1];
sx q[1];
rz(0.71818761) q[1];
rz(-pi) q[2];
rz(-0.22110053) q[3];
sx q[3];
rz(-2.3490864) q[3];
sx q[3];
rz(-1.4421876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6654309) q[2];
sx q[2];
rz(-0.68887812) q[2];
sx q[2];
rz(0.8821744) q[2];
rz(0.65794182) q[3];
sx q[3];
rz(-1.123817) q[3];
sx q[3];
rz(2.1420124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.5049725) q[0];
sx q[0];
rz(-1.5047147) q[0];
sx q[0];
rz(2.104105) q[0];
rz(0.50719117) q[1];
sx q[1];
rz(-0.7829983) q[1];
sx q[1];
rz(-1.0601039) q[1];
rz(0.29303356) q[2];
sx q[2];
rz(-1.3923981) q[2];
sx q[2];
rz(-0.57033718) q[2];
rz(-2.2118123) q[3];
sx q[3];
rz(-2.2140183) q[3];
sx q[3];
rz(2.1341013) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
