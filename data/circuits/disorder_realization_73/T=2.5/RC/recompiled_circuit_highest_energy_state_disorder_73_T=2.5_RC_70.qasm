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
rz(-1.7492548) q[0];
sx q[0];
rz(-0.37819401) q[0];
sx q[0];
rz(0.35051546) q[0];
rz(0.89138436) q[1];
sx q[1];
rz(3.9337629) q[1];
sx q[1];
rz(13.73929) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.590345) q[0];
sx q[0];
rz(-2.8714193) q[0];
sx q[0];
rz(-0.7858289) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95845285) q[2];
sx q[2];
rz(-2.8757189) q[2];
sx q[2];
rz(2.079351) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8280331) q[1];
sx q[1];
rz(-2.782722) q[1];
sx q[1];
rz(0.48017217) q[1];
rz(-1.7426026) q[3];
sx q[3];
rz(-2.0366002) q[3];
sx q[3];
rz(-1.727263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9709836) q[2];
sx q[2];
rz(-1.8656518) q[2];
sx q[2];
rz(2.9288536) q[2];
rz(0.085266026) q[3];
sx q[3];
rz(-1.8906967) q[3];
sx q[3];
rz(-1.8703478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79412115) q[0];
sx q[0];
rz(-2.329282) q[0];
sx q[0];
rz(2.0982657) q[0];
rz(-3.0087545) q[1];
sx q[1];
rz(-0.21505198) q[1];
sx q[1];
rz(-1.1588233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2272233) q[0];
sx q[0];
rz(-0.4242737) q[0];
sx q[0];
rz(2.7844564) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01387502) q[2];
sx q[2];
rz(-1.4094947) q[2];
sx q[2];
rz(-2.814295) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.093250153) q[1];
sx q[1];
rz(-0.56386098) q[1];
sx q[1];
rz(-2.6287931) q[1];
rz(-0.93166931) q[3];
sx q[3];
rz(-2.6216345) q[3];
sx q[3];
rz(-0.66321841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1966689) q[2];
sx q[2];
rz(-2.479574) q[2];
sx q[2];
rz(2.4169253) q[2];
rz(-2.479018) q[3];
sx q[3];
rz(-0.96810883) q[3];
sx q[3];
rz(2.841943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1835943) q[0];
sx q[0];
rz(-1.2677001) q[0];
sx q[0];
rz(0.9147574) q[0];
rz(2.5547408) q[1];
sx q[1];
rz(-1.4205168) q[1];
sx q[1];
rz(-0.14952001) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6059593) q[0];
sx q[0];
rz(-2.0894755) q[0];
sx q[0];
rz(-2.9796322) q[0];
x q[1];
rz(-0.36386602) q[2];
sx q[2];
rz(-1.9155028) q[2];
sx q[2];
rz(0.41512903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65626493) q[1];
sx q[1];
rz(-1.0228923) q[1];
sx q[1];
rz(2.5119971) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2633576) q[3];
sx q[3];
rz(-0.98537579) q[3];
sx q[3];
rz(0.4957605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84825039) q[2];
sx q[2];
rz(-1.205227) q[2];
sx q[2];
rz(-2.3411574) q[2];
rz(2.6750001) q[3];
sx q[3];
rz(-1.1060017) q[3];
sx q[3];
rz(0.65565562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7993497) q[0];
sx q[0];
rz(-1.8514587) q[0];
sx q[0];
rz(-2.2829862) q[0];
rz(-0.41060064) q[1];
sx q[1];
rz(-2.938439) q[1];
sx q[1];
rz(2.1536749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77495518) q[0];
sx q[0];
rz(-1.8995842) q[0];
sx q[0];
rz(-2.9820697) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3484389) q[2];
sx q[2];
rz(-2.6359865) q[2];
sx q[2];
rz(2.8067775) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6776442) q[1];
sx q[1];
rz(-1.8419187) q[1];
sx q[1];
rz(-1.2605002) q[1];
x q[2];
rz(-2.6060054) q[3];
sx q[3];
rz(-0.37112826) q[3];
sx q[3];
rz(0.27419004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9504488) q[2];
sx q[2];
rz(-1.7473651) q[2];
sx q[2];
rz(-2.6960755) q[2];
rz(-2.4281003) q[3];
sx q[3];
rz(-2.4092509) q[3];
sx q[3];
rz(-0.92002404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5934061) q[0];
sx q[0];
rz(-2.0596518) q[0];
sx q[0];
rz(-1.5997546) q[0];
rz(-1.2708739) q[1];
sx q[1];
rz(-1.4022695) q[1];
sx q[1];
rz(-0.52938968) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75225509) q[0];
sx q[0];
rz(-1.9175954) q[0];
sx q[0];
rz(2.8777468) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5240493) q[2];
sx q[2];
rz(-1.7360592) q[2];
sx q[2];
rz(-2.3338855) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8258851) q[1];
sx q[1];
rz(-1.3199908) q[1];
sx q[1];
rz(2.6205089) q[1];
rz(-2.9934817) q[3];
sx q[3];
rz(-1.221773) q[3];
sx q[3];
rz(-2.7452041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9993837) q[2];
sx q[2];
rz(-1.9500407) q[2];
sx q[2];
rz(1.3231529) q[2];
rz(2.7600944) q[3];
sx q[3];
rz(-1.1161209) q[3];
sx q[3];
rz(2.748446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2926272) q[0];
sx q[0];
rz(-2.3830074) q[0];
sx q[0];
rz(-1.0211771) q[0];
rz(-1.5317597) q[1];
sx q[1];
rz(-2.1949218) q[1];
sx q[1];
rz(0.80345947) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55554237) q[0];
sx q[0];
rz(-1.4110421) q[0];
sx q[0];
rz(-2.8303888) q[0];
rz(2.7293936) q[2];
sx q[2];
rz(-1.4934818) q[2];
sx q[2];
rz(-2.9207723) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.66113055) q[1];
sx q[1];
rz(-1.6906889) q[1];
sx q[1];
rz(2.5698623) q[1];
rz(2.5147076) q[3];
sx q[3];
rz(-1.6520368) q[3];
sx q[3];
rz(2.8089942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44567406) q[2];
sx q[2];
rz(-0.58738223) q[2];
sx q[2];
rz(0.43061259) q[2];
rz(2.5066091) q[3];
sx q[3];
rz(-0.018714232) q[3];
sx q[3];
rz(-1.8733321) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1465313) q[0];
sx q[0];
rz(-0.54556161) q[0];
sx q[0];
rz(1.5437641) q[0];
rz(0.68430463) q[1];
sx q[1];
rz(-1.4477718) q[1];
sx q[1];
rz(2.3421471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9198624) q[0];
sx q[0];
rz(-1.5654025) q[0];
sx q[0];
rz(-3.0649867) q[0];
rz(-1.6580901) q[2];
sx q[2];
rz(-0.71719522) q[2];
sx q[2];
rz(-2.4503477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23161665) q[1];
sx q[1];
rz(-1.5607395) q[1];
sx q[1];
rz(2.6811428) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43947222) q[3];
sx q[3];
rz(-2.2547561) q[3];
sx q[3];
rz(0.33392936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53576175) q[2];
sx q[2];
rz(-2.3039218) q[2];
sx q[2];
rz(-2.3568995) q[2];
rz(3.0253518) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(-1.2522662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5410974) q[0];
sx q[0];
rz(-1.9116115) q[0];
sx q[0];
rz(-2.1121209) q[0];
rz(0.12044278) q[1];
sx q[1];
rz(-1.8414626) q[1];
sx q[1];
rz(-0.57074237) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1627038) q[0];
sx q[0];
rz(-1.9523639) q[0];
sx q[0];
rz(-2.0636153) q[0];
rz(-pi) q[1];
rz(-2.4230401) q[2];
sx q[2];
rz(-1.4863803) q[2];
sx q[2];
rz(-0.29044232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9214258) q[1];
sx q[1];
rz(-1.1296185) q[1];
sx q[1];
rz(0.36018546) q[1];
rz(1.086497) q[3];
sx q[3];
rz(-2.7653793) q[3];
sx q[3];
rz(2.0813326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.268078) q[2];
sx q[2];
rz(-0.87468481) q[2];
sx q[2];
rz(-2.6178005) q[2];
rz(-1.1526456) q[3];
sx q[3];
rz(-0.58061424) q[3];
sx q[3];
rz(-1.4279648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6607587) q[0];
sx q[0];
rz(-1.2013712) q[0];
sx q[0];
rz(0.42385605) q[0];
rz(-2.2747874) q[1];
sx q[1];
rz(-2.1173756) q[1];
sx q[1];
rz(-2.9383235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5728666) q[0];
sx q[0];
rz(-0.6558658) q[0];
sx q[0];
rz(-3.0650861) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31000455) q[2];
sx q[2];
rz(-1.6223202) q[2];
sx q[2];
rz(-1.0088986) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79938526) q[1];
sx q[1];
rz(-0.3179271) q[1];
sx q[1];
rz(-2.6683776) q[1];
rz(-pi) q[2];
rz(-2.8839797) q[3];
sx q[3];
rz(-1.2748162) q[3];
sx q[3];
rz(-1.6483893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.086143494) q[2];
sx q[2];
rz(-1.92675) q[2];
sx q[2];
rz(-2.8622368) q[2];
rz(1.8481567) q[3];
sx q[3];
rz(-1.8297628) q[3];
sx q[3];
rz(-1.2156585) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16348895) q[0];
sx q[0];
rz(-0.80253974) q[0];
sx q[0];
rz(1.4168903) q[0];
rz(-2.2143927) q[1];
sx q[1];
rz(-1.5798774) q[1];
sx q[1];
rz(-1.4097811) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15148189) q[0];
sx q[0];
rz(-1.6674433) q[0];
sx q[0];
rz(2.2352909) q[0];
rz(-pi) q[1];
rz(-0.3566202) q[2];
sx q[2];
rz(-2.0276514) q[2];
sx q[2];
rz(2.8518554) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43391398) q[1];
sx q[1];
rz(-2.0053889) q[1];
sx q[1];
rz(2.4151925) q[1];
rz(0.042155592) q[3];
sx q[3];
rz(-1.32844) q[3];
sx q[3];
rz(2.4740296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37567821) q[2];
sx q[2];
rz(-1.3216852) q[2];
sx q[2];
rz(-2.1709757) q[2];
rz(-2.4961903) q[3];
sx q[3];
rz(-2.161945) q[3];
sx q[3];
rz(1.0055044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476625) q[0];
sx q[0];
rz(-1.6456589) q[0];
sx q[0];
rz(1.8346067) q[0];
rz(2.7597799) q[1];
sx q[1];
rz(-1.3419071) q[1];
sx q[1];
rz(0.76795427) q[1];
rz(-2.1099595) q[2];
sx q[2];
rz(-0.8349541) q[2];
sx q[2];
rz(-0.63944774) q[2];
rz(0.80841165) q[3];
sx q[3];
rz(-2.5403317) q[3];
sx q[3];
rz(-2.9797277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
