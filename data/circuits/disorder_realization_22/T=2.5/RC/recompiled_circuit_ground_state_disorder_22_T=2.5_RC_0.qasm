OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.206267) q[0];
sx q[0];
rz(-1.6214108) q[0];
sx q[0];
rz(1.4186463) q[0];
rz(6.3267646) q[1];
sx q[1];
rz(5.0001231) q[1];
sx q[1];
rz(12.725886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.537764) q[0];
sx q[0];
rz(-1.7366752) q[0];
sx q[0];
rz(0.028847522) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8279499) q[2];
sx q[2];
rz(-0.97460496) q[2];
sx q[2];
rz(1.8108846) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7456513) q[1];
sx q[1];
rz(-2.6439754) q[1];
sx q[1];
rz(1.2785589) q[1];
rz(-pi) q[2];
rz(-0.84955022) q[3];
sx q[3];
rz(-1.8462984) q[3];
sx q[3];
rz(1.4412137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45304766) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(-3.0513406) q[2];
rz(3.0669751) q[3];
sx q[3];
rz(-1.4679694) q[3];
sx q[3];
rz(0.41372764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8892882) q[0];
sx q[0];
rz(-1.9828718) q[0];
sx q[0];
rz(1.4537551) q[0];
rz(0.74268913) q[1];
sx q[1];
rz(-1.366726) q[1];
sx q[1];
rz(2.3765366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9738282) q[0];
sx q[0];
rz(-2.544738) q[0];
sx q[0];
rz(0.11886998) q[0];
rz(-pi) q[1];
rz(0.80185955) q[2];
sx q[2];
rz(-0.86881402) q[2];
sx q[2];
rz(3.1131256) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6366406) q[1];
sx q[1];
rz(-0.8209629) q[1];
sx q[1];
rz(-2.0269462) q[1];
rz(-0.0043930769) q[3];
sx q[3];
rz(-0.54800225) q[3];
sx q[3];
rz(-2.9145384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.226977) q[2];
sx q[2];
rz(-1.0116297) q[2];
sx q[2];
rz(1.6870618) q[2];
rz(-2.4948273) q[3];
sx q[3];
rz(-0.55964595) q[3];
sx q[3];
rz(-1.2796992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890559) q[0];
sx q[0];
rz(-1.6707957) q[0];
sx q[0];
rz(-3.0964417) q[0];
rz(-1.2308944) q[1];
sx q[1];
rz(-1.8321593) q[1];
sx q[1];
rz(-2.0848354) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8674744) q[0];
sx q[0];
rz(-1.2454659) q[0];
sx q[0];
rz(2.1528457) q[0];
rz(-2.5185539) q[2];
sx q[2];
rz(-2.3941561) q[2];
sx q[2];
rz(1.7926271) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9581117) q[1];
sx q[1];
rz(-0.90881244) q[1];
sx q[1];
rz(0.32420154) q[1];
rz(0.23356861) q[3];
sx q[3];
rz(-2.2456777) q[3];
sx q[3];
rz(0.44665953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4061654) q[2];
sx q[2];
rz(-0.96144095) q[2];
sx q[2];
rz(2.4889634) q[2];
rz(1.2930219) q[3];
sx q[3];
rz(-2.3725464) q[3];
sx q[3];
rz(3.0434216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80766135) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(1.3941143) q[0];
rz(1.8380503) q[1];
sx q[1];
rz(-1.9775672) q[1];
sx q[1];
rz(-2.0533144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63193479) q[0];
sx q[0];
rz(-0.83197537) q[0];
sx q[0];
rz(-0.37507071) q[0];
x q[1];
rz(0.59594229) q[2];
sx q[2];
rz(-0.94178994) q[2];
sx q[2];
rz(-3.0160835) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7409121) q[1];
sx q[1];
rz(-2.0517382) q[1];
sx q[1];
rz(-2.4467642) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6234947) q[3];
sx q[3];
rz(-1.0720338) q[3];
sx q[3];
rz(0.76585211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7465705) q[2];
sx q[2];
rz(-1.3777379) q[2];
sx q[2];
rz(2.6611888) q[2];
rz(0.26614842) q[3];
sx q[3];
rz(-1.1689309) q[3];
sx q[3];
rz(-0.84421414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0323623) q[0];
sx q[0];
rz(-2.5205595) q[0];
sx q[0];
rz(1.2048298) q[0];
rz(-0.042639848) q[1];
sx q[1];
rz(-1.3748906) q[1];
sx q[1];
rz(-2.1991594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.401817) q[0];
sx q[0];
rz(-1.7063921) q[0];
sx q[0];
rz(-1.1419161) q[0];
rz(-pi) q[1];
rz(-2.3187014) q[2];
sx q[2];
rz(-2.3510755) q[2];
sx q[2];
rz(-2.2728777) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1548295) q[1];
sx q[1];
rz(-1.0752605) q[1];
sx q[1];
rz(0.43191989) q[1];
rz(-pi) q[2];
rz(-2.5025401) q[3];
sx q[3];
rz(-1.6585232) q[3];
sx q[3];
rz(-1.7114342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11613906) q[2];
sx q[2];
rz(-2.4567273) q[2];
sx q[2];
rz(1.9174513) q[2];
rz(-2.4141566) q[3];
sx q[3];
rz(-1.6614611) q[3];
sx q[3];
rz(0.050189655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8023476) q[0];
sx q[0];
rz(-0.50528637) q[0];
sx q[0];
rz(-2.8832866) q[0];
rz(2.2282265) q[1];
sx q[1];
rz(-2.2759627) q[1];
sx q[1];
rz(0.9679274) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0845522) q[0];
sx q[0];
rz(-1.0121317) q[0];
sx q[0];
rz(2.0303805) q[0];
rz(1.1656233) q[2];
sx q[2];
rz(-0.76442598) q[2];
sx q[2];
rz(2.7583099) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.591394) q[1];
sx q[1];
rz(-1.0934783) q[1];
sx q[1];
rz(-1.6106907) q[1];
rz(-pi) q[2];
rz(-2.4658396) q[3];
sx q[3];
rz(-2.3738513) q[3];
sx q[3];
rz(2.2217285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4919058) q[2];
sx q[2];
rz(-2.457149) q[2];
sx q[2];
rz(0.41927949) q[2];
rz(2.3573917) q[3];
sx q[3];
rz(-1.0201642) q[3];
sx q[3];
rz(1.3721344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.075031) q[0];
sx q[0];
rz(-3.0176268) q[0];
sx q[0];
rz(0.01874622) q[0];
rz(3.0491414) q[1];
sx q[1];
rz(-1.8094614) q[1];
sx q[1];
rz(3.0466373) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088647) q[0];
sx q[0];
rz(-1.0148541) q[0];
sx q[0];
rz(-1.9535319) q[0];
rz(1.9906391) q[2];
sx q[2];
rz(-1.0165086) q[2];
sx q[2];
rz(1.3406164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9419721) q[1];
sx q[1];
rz(-2.6411861) q[1];
sx q[1];
rz(0.54897658) q[1];
x q[2];
rz(-2.9800013) q[3];
sx q[3];
rz(-1.1757903) q[3];
sx q[3];
rz(2.145951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9065325) q[2];
sx q[2];
rz(-1.1842714) q[2];
sx q[2];
rz(2.2825784) q[2];
rz(0.61339316) q[3];
sx q[3];
rz(-2.9412013) q[3];
sx q[3];
rz(1.8092417) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61619806) q[0];
sx q[0];
rz(-1.8680251) q[0];
sx q[0];
rz(1.6599576) q[0];
rz(1.067465) q[1];
sx q[1];
rz(-2.4165202) q[1];
sx q[1];
rz(-1.8225089) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5123915) q[0];
sx q[0];
rz(-0.90921697) q[0];
sx q[0];
rz(3.0223836) q[0];
x q[1];
rz(-0.51609938) q[2];
sx q[2];
rz(-1.8987499) q[2];
sx q[2];
rz(-2.536377) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7072152) q[1];
sx q[1];
rz(-2.5188252) q[1];
sx q[1];
rz(-1.9799443) q[1];
rz(-0.33101636) q[3];
sx q[3];
rz(-2.7661471) q[3];
sx q[3];
rz(2.7956243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8760425) q[2];
sx q[2];
rz(-2.7713113) q[2];
sx q[2];
rz(-0.041157095) q[2];
rz(1.1228849) q[3];
sx q[3];
rz(-1.4833781) q[3];
sx q[3];
rz(-0.52367228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6599643) q[0];
sx q[0];
rz(-2.1631503) q[0];
sx q[0];
rz(-1.6690669) q[0];
rz(2.4166079) q[1];
sx q[1];
rz(-1.9703194) q[1];
sx q[1];
rz(-2.5670126) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32550463) q[0];
sx q[0];
rz(-0.34127745) q[0];
sx q[0];
rz(-2.9983872) q[0];
rz(-pi) q[1];
rz(-2.8594744) q[2];
sx q[2];
rz(-1.9777672) q[2];
sx q[2];
rz(2.5438683) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7704627) q[1];
sx q[1];
rz(-1.5411301) q[1];
sx q[1];
rz(2.4231623) q[1];
rz(1.7506041) q[3];
sx q[3];
rz(-0.69560183) q[3];
sx q[3];
rz(-0.66015388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4943115) q[2];
sx q[2];
rz(-2.8936671) q[2];
sx q[2];
rz(2.1809123) q[2];
rz(0.48323768) q[3];
sx q[3];
rz(-1.5944696) q[3];
sx q[3];
rz(-1.5620935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305337) q[0];
sx q[0];
rz(-2.647825) q[0];
sx q[0];
rz(0.44179398) q[0];
rz(0.68253016) q[1];
sx q[1];
rz(-1.4492757) q[1];
sx q[1];
rz(1.0535047) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5614612) q[0];
sx q[0];
rz(-0.86505689) q[0];
sx q[0];
rz(-2.8012026) q[0];
x q[1];
rz(-1.9437213) q[2];
sx q[2];
rz(-1.6536384) q[2];
sx q[2];
rz(1.5035455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39643416) q[1];
sx q[1];
rz(-2.4188383) q[1];
sx q[1];
rz(2.0682343) q[1];
rz(-pi) q[2];
rz(-1.6445475) q[3];
sx q[3];
rz(-1.9637655) q[3];
sx q[3];
rz(2.949209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0423476) q[2];
sx q[2];
rz(-1.8727563) q[2];
sx q[2];
rz(0.22259268) q[2];
rz(-1.7706722) q[3];
sx q[3];
rz(-1.418908) q[3];
sx q[3];
rz(-2.5194949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.901578) q[0];
sx q[0];
rz(-2.1299025) q[0];
sx q[0];
rz(1.0255751) q[0];
rz(1.9758132) q[1];
sx q[1];
rz(-2.5318601) q[1];
sx q[1];
rz(-0.21808521) q[1];
rz(-2.0596621) q[2];
sx q[2];
rz(-1.0837383) q[2];
sx q[2];
rz(2.4080129) q[2];
rz(1.8924186) q[3];
sx q[3];
rz(-1.5703441) q[3];
sx q[3];
rz(2.3342568) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
