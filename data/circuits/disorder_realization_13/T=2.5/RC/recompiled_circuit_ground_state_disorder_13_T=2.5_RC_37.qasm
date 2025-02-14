OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(-0.14578851) q[0];
sx q[0];
rz(1.2989651) q[0];
rz(-0.19620148) q[1];
sx q[1];
rz(-1.3407522) q[1];
sx q[1];
rz(3.0412716) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.318905) q[0];
sx q[0];
rz(-1.6226557) q[0];
sx q[0];
rz(2.479631) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6703628) q[2];
sx q[2];
rz(-1.8719851) q[2];
sx q[2];
rz(0.1567086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4397706) q[1];
sx q[1];
rz(-1.3622829) q[1];
sx q[1];
rz(-2.8609402) q[1];
rz(-pi) q[2];
rz(0.59009976) q[3];
sx q[3];
rz(-1.3076772) q[3];
sx q[3];
rz(2.2424169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4002865) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(0.8068223) q[2];
rz(0.72871366) q[3];
sx q[3];
rz(-0.7586793) q[3];
sx q[3];
rz(-1.2046643) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5173986) q[0];
sx q[0];
rz(-0.26573467) q[0];
sx q[0];
rz(-0.30701315) q[0];
rz(-1.864805) q[1];
sx q[1];
rz(-2.0025608) q[1];
sx q[1];
rz(-1.101864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3392838) q[0];
sx q[0];
rz(-1.7442987) q[0];
sx q[0];
rz(1.0109148) q[0];
rz(-0.33874933) q[2];
sx q[2];
rz(-0.18998665) q[2];
sx q[2];
rz(2.7471586) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9872322) q[1];
sx q[1];
rz(-1.4353961) q[1];
sx q[1];
rz(-1.7495278) q[1];
x q[2];
rz(-2.4073007) q[3];
sx q[3];
rz(-0.9356539) q[3];
sx q[3];
rz(-1.4189625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.028194204) q[2];
sx q[2];
rz(-2.5598309) q[2];
sx q[2];
rz(-3.0451575) q[2];
rz(0.15245572) q[3];
sx q[3];
rz(-1.5072482) q[3];
sx q[3];
rz(0.69795394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089712791) q[0];
sx q[0];
rz(-1.1002325) q[0];
sx q[0];
rz(2.7413947) q[0];
rz(1.6290889) q[1];
sx q[1];
rz(-2.9632603) q[1];
sx q[1];
rz(2.8834744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54329693) q[0];
sx q[0];
rz(-0.27320293) q[0];
sx q[0];
rz(-1.0323204) q[0];
rz(-pi) q[1];
rz(-0.67419184) q[2];
sx q[2];
rz(-0.12756824) q[2];
sx q[2];
rz(-1.3889165) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8363179) q[1];
sx q[1];
rz(-1.5658251) q[1];
sx q[1];
rz(-0.012955091) q[1];
rz(2.9003224) q[3];
sx q[3];
rz(-1.7381958) q[3];
sx q[3];
rz(0.77360717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1199946) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(-0.032133948) q[2];
rz(0.26767996) q[3];
sx q[3];
rz(-1.5175502) q[3];
sx q[3];
rz(-1.5816429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5716008) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(-0.084240325) q[0];
rz(-0.035942297) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(0.34119225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223136) q[0];
sx q[0];
rz(-1.3094762) q[0];
sx q[0];
rz(0.65066353) q[0];
x q[1];
rz(2.555055) q[2];
sx q[2];
rz(-1.7866724) q[2];
sx q[2];
rz(-3.0160273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0964716) q[1];
sx q[1];
rz(-1.1323788) q[1];
sx q[1];
rz(2.3271266) q[1];
rz(-pi) q[2];
rz(-2.0218418) q[3];
sx q[3];
rz(-1.0195135) q[3];
sx q[3];
rz(-1.0663084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42698947) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(1.6061456) q[2];
rz(2.8793907) q[3];
sx q[3];
rz(-1.6180792) q[3];
sx q[3];
rz(2.1867627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3839805) q[0];
sx q[0];
rz(-2.7181427) q[0];
sx q[0];
rz(-1.0773995) q[0];
rz(-2.7491838) q[1];
sx q[1];
rz(-0.078374021) q[1];
sx q[1];
rz(-1.0307301) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2868256) q[0];
sx q[0];
rz(-2.1041901) q[0];
sx q[0];
rz(2.0598165) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0692798) q[2];
sx q[2];
rz(-1.7297812) q[2];
sx q[2];
rz(-2.0778401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2490152) q[1];
sx q[1];
rz(-1.4059781) q[1];
sx q[1];
rz(-1.7958162) q[1];
x q[2];
rz(1.0593653) q[3];
sx q[3];
rz(-1.8728731) q[3];
sx q[3];
rz(1.5010639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.82421676) q[2];
sx q[2];
rz(-2.4874096) q[2];
sx q[2];
rz(-2.3023494) q[2];
rz(0.9134891) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(2.998108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6827253) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(1.6656026) q[0];
rz(0.39235517) q[1];
sx q[1];
rz(-1.0959492) q[1];
sx q[1];
rz(-2.5501693) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05706035) q[0];
sx q[0];
rz(-2.1520237) q[0];
sx q[0];
rz(2.6956062) q[0];
x q[1];
rz(2.5673836) q[2];
sx q[2];
rz(-0.32445766) q[2];
sx q[2];
rz(-0.41947075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6063664) q[1];
sx q[1];
rz(-2.6312345) q[1];
sx q[1];
rz(-1.7061069) q[1];
x q[2];
rz(-0.37844946) q[3];
sx q[3];
rz(-2.2317118) q[3];
sx q[3];
rz(-1.4405516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8412987) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(2.5692614) q[2];
rz(0.22219292) q[3];
sx q[3];
rz(-0.43155813) q[3];
sx q[3];
rz(-0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(0.38729024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(2.733316) q[0];
rz(-2.4093742) q[1];
sx q[1];
rz(-3.0156942) q[1];
sx q[1];
rz(-2.8439723) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8430804) q[0];
sx q[0];
rz(-1.9924) q[0];
sx q[0];
rz(-1.9527736) q[0];
rz(2.0655259) q[2];
sx q[2];
rz(-1.552993) q[2];
sx q[2];
rz(-0.36843637) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.054033259) q[1];
sx q[1];
rz(-2.7440024) q[1];
sx q[1];
rz(1.041574) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4206395) q[3];
sx q[3];
rz(-0.87267002) q[3];
sx q[3];
rz(1.0696749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.97062651) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(0.69620281) q[2];
rz(2.2512186) q[3];
sx q[3];
rz(-1.9647157) q[3];
sx q[3];
rz(-1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17906976) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(-2.9288375) q[0];
rz(2.6720324) q[1];
sx q[1];
rz(-2.1766365) q[1];
sx q[1];
rz(-2.3874217) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.045005) q[0];
sx q[0];
rz(-1.2822598) q[0];
sx q[0];
rz(0.19591422) q[0];
rz(2.8163359) q[2];
sx q[2];
rz(-1.3768679) q[2];
sx q[2];
rz(-0.12530279) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.27867815) q[1];
sx q[1];
rz(-0.92611852) q[1];
sx q[1];
rz(-1.7134604) q[1];
rz(-pi) q[2];
rz(-0.92408224) q[3];
sx q[3];
rz(-1.0093186) q[3];
sx q[3];
rz(-1.1738861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0065877) q[2];
sx q[2];
rz(-0.90234119) q[2];
sx q[2];
rz(-2.3593486) q[2];
rz(1.4462645) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049147216) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(2.1771722) q[0];
rz(1.2760705) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(1.5020348) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14012155) q[0];
sx q[0];
rz(-1.3706511) q[0];
sx q[0];
rz(-1.8076623) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4110498) q[2];
sx q[2];
rz(-1.9078095) q[2];
sx q[2];
rz(-0.72357976) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6148541) q[1];
sx q[1];
rz(-0.51797082) q[1];
sx q[1];
rz(1.6171842) q[1];
rz(-pi) q[2];
rz(1.4235017) q[3];
sx q[3];
rz(-1.5039267) q[3];
sx q[3];
rz(1.221958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2532578) q[2];
sx q[2];
rz(-1.8586681) q[2];
sx q[2];
rz(2.1962732) q[2];
rz(-2.3156598) q[3];
sx q[3];
rz(-1.632894) q[3];
sx q[3];
rz(-0.53502423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0655521) q[0];
sx q[0];
rz(-1.9726418) q[0];
sx q[0];
rz(-2.3384576) q[0];
rz(-1.5777292) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(-2.8520083) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0286547) q[0];
sx q[0];
rz(-1.583033) q[0];
sx q[0];
rz(-1.4652246) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8990191) q[2];
sx q[2];
rz(-1.8716836) q[2];
sx q[2];
rz(-2.2102578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8378505) q[1];
sx q[1];
rz(-1.2134064) q[1];
sx q[1];
rz(0.2405432) q[1];
rz(-pi) q[2];
rz(2.2749994) q[3];
sx q[3];
rz(-2.7150053) q[3];
sx q[3];
rz(2.085919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.76558602) q[2];
sx q[2];
rz(-0.12066081) q[2];
sx q[2];
rz(0.96013367) q[2];
rz(-2.5586186) q[3];
sx q[3];
rz(-2.4834902) q[3];
sx q[3];
rz(-0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4509907) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(0.040738978) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(1.1004943) q[2];
sx q[2];
rz(-1.294877) q[2];
sx q[2];
rz(-0.69653947) q[2];
rz(0.95016877) q[3];
sx q[3];
rz(-3.0541854) q[3];
sx q[3];
rz(-1.0704609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
