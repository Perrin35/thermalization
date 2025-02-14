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
rz(1.3923378) q[0];
sx q[0];
rz(-2.7633986) q[0];
sx q[0];
rz(2.7910772) q[0];
rz(-2.2502083) q[1];
sx q[1];
rz(-0.79217029) q[1];
sx q[1];
rz(-1.1729191) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25305155) q[0];
sx q[0];
rz(-1.3810147) q[0];
sx q[0];
rz(1.3773328) q[0];
rz(-pi) q[1];
rz(1.3515377) q[2];
sx q[2];
rz(-1.4191896) q[2];
sx q[2];
rz(1.104179) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3135595) q[1];
sx q[1];
rz(-2.782722) q[1];
sx q[1];
rz(0.48017217) q[1];
rz(2.6698106) q[3];
sx q[3];
rz(-1.4174479) q[3];
sx q[3];
rz(2.9073496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9709836) q[2];
sx q[2];
rz(-1.2759408) q[2];
sx q[2];
rz(-2.9288536) q[2];
rz(0.085266026) q[3];
sx q[3];
rz(-1.8906967) q[3];
sx q[3];
rz(1.2712449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3474715) q[0];
sx q[0];
rz(-0.81231064) q[0];
sx q[0];
rz(2.0982657) q[0];
rz(-0.13283816) q[1];
sx q[1];
rz(-0.21505198) q[1];
sx q[1];
rz(1.1588233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3029859) q[0];
sx q[0];
rz(-1.1748472) q[0];
sx q[0];
rz(1.4141757) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7321132) q[2];
sx q[2];
rz(-1.5571014) q[2];
sx q[2];
rz(-1.2412702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1082127) q[1];
sx q[1];
rz(-1.8361109) q[1];
sx q[1];
rz(0.50361522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93166931) q[3];
sx q[3];
rz(-0.51995819) q[3];
sx q[3];
rz(2.4783742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1966689) q[2];
sx q[2];
rz(-2.479574) q[2];
sx q[2];
rz(0.72466737) q[2];
rz(-2.479018) q[3];
sx q[3];
rz(-2.1734838) q[3];
sx q[3];
rz(0.29964963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1835943) q[0];
sx q[0];
rz(-1.2677001) q[0];
sx q[0];
rz(-0.9147574) q[0];
rz(2.5547408) q[1];
sx q[1];
rz(-1.7210759) q[1];
sx q[1];
rz(0.14952001) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6059593) q[0];
sx q[0];
rz(-1.0521172) q[0];
sx q[0];
rz(2.9796322) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36386602) q[2];
sx q[2];
rz(-1.9155028) q[2];
sx q[2];
rz(0.41512903) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.55184522) q[1];
sx q[1];
rz(-1.044163) q[1];
sx q[1];
rz(0.92411516) q[1];
rz(-0.87823509) q[3];
sx q[3];
rz(-2.1562169) q[3];
sx q[3];
rz(-2.6458322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2933423) q[2];
sx q[2];
rz(-1.9363656) q[2];
sx q[2];
rz(-2.3411574) q[2];
rz(-0.46659255) q[3];
sx q[3];
rz(-1.1060017) q[3];
sx q[3];
rz(-2.485937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34224299) q[0];
sx q[0];
rz(-1.8514587) q[0];
sx q[0];
rz(0.85860646) q[0];
rz(0.41060064) q[1];
sx q[1];
rz(-0.20315367) q[1];
sx q[1];
rz(-0.98791775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3666375) q[0];
sx q[0];
rz(-1.8995842) q[0];
sx q[0];
rz(0.15952296) q[0];
x q[1];
rz(-1.0756827) q[2];
sx q[2];
rz(-1.4637814) q[2];
sx q[2];
rz(1.040695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.021182755) q[1];
sx q[1];
rz(-1.8694009) q[1];
sx q[1];
rz(2.8575767) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7668523) q[3];
sx q[3];
rz(-1.88797) q[3];
sx q[3];
rz(-2.848742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9504488) q[2];
sx q[2];
rz(-1.3942275) q[2];
sx q[2];
rz(0.44551715) q[2];
rz(2.4281003) q[3];
sx q[3];
rz(-0.73234171) q[3];
sx q[3];
rz(2.2215686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5934061) q[0];
sx q[0];
rz(-1.0819409) q[0];
sx q[0];
rz(-1.5997546) q[0];
rz(1.2708739) q[1];
sx q[1];
rz(-1.4022695) q[1];
sx q[1];
rz(0.52938968) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7177757) q[0];
sx q[0];
rz(-2.7090906) q[0];
sx q[0];
rz(-2.1958817) q[0];
x q[1];
rz(-2.8611456) q[2];
sx q[2];
rz(-0.63647645) q[2];
sx q[2];
rz(-2.150879) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84698757) q[1];
sx q[1];
rz(-2.5683476) q[1];
sx q[1];
rz(-2.6662988) q[1];
rz(-pi) q[2];
rz(-1.2182203) q[3];
sx q[3];
rz(-1.7099172) q[3];
sx q[3];
rz(-1.2253882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9993837) q[2];
sx q[2];
rz(-1.1915519) q[2];
sx q[2];
rz(1.3231529) q[2];
rz(-0.38149825) q[3];
sx q[3];
rz(-2.0254717) q[3];
sx q[3];
rz(-2.748446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8489654) q[0];
sx q[0];
rz(-0.75858527) q[0];
sx q[0];
rz(-2.1204156) q[0];
rz(-1.5317597) q[1];
sx q[1];
rz(-0.94667089) q[1];
sx q[1];
rz(-0.80345947) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96413104) q[0];
sx q[0];
rz(-1.2636856) q[0];
sx q[0];
rz(-1.7384647) q[0];
rz(2.7293936) q[2];
sx q[2];
rz(-1.4934818) q[2];
sx q[2];
rz(-2.9207723) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66113055) q[1];
sx q[1];
rz(-1.6906889) q[1];
sx q[1];
rz(0.57173034) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4705974) q[3];
sx q[3];
rz(-0.9462983) q[3];
sx q[3];
rz(-1.1794832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6959186) q[2];
sx q[2];
rz(-2.5542104) q[2];
sx q[2];
rz(-0.43061259) q[2];
rz(-2.5066091) q[3];
sx q[3];
rz(-3.1228784) q[3];
sx q[3];
rz(1.2682605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9950614) q[0];
sx q[0];
rz(-0.54556161) q[0];
sx q[0];
rz(1.5437641) q[0];
rz(-2.457288) q[1];
sx q[1];
rz(-1.6938208) q[1];
sx q[1];
rz(0.79944557) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.862683) q[0];
sx q[0];
rz(-0.076795243) q[0];
sx q[0];
rz(3.0712295) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85549037) q[2];
sx q[2];
rz(-1.6281307) q[2];
sx q[2];
rz(-2.3279026) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7974248) q[1];
sx q[1];
rz(-1.1103715) q[1];
sx q[1];
rz(-1.5820222) q[1];
rz(-pi) q[2];
rz(-2.0518028) q[3];
sx q[3];
rz(-0.79339992) q[3];
sx q[3];
rz(-0.97360669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53576175) q[2];
sx q[2];
rz(-2.3039218) q[2];
sx q[2];
rz(2.3568995) q[2];
rz(0.11624087) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(-1.8893265) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5410974) q[0];
sx q[0];
rz(-1.2299812) q[0];
sx q[0];
rz(-1.0294718) q[0];
rz(0.12044278) q[1];
sx q[1];
rz(-1.8414626) q[1];
sx q[1];
rz(2.5708503) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3523173) q[0];
sx q[0];
rz(-2.0253775) q[0];
sx q[0];
rz(-0.42735703) q[0];
rz(-pi) q[1];
rz(1.6827356) q[2];
sx q[2];
rz(-2.2862391) q[2];
sx q[2];
rz(1.2067522) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2201669) q[1];
sx q[1];
rz(-2.0119742) q[1];
sx q[1];
rz(2.7814072) q[1];
x q[2];
rz(1.9071155) q[3];
sx q[3];
rz(-1.7426995) q[3];
sx q[3];
rz(3.086103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.268078) q[2];
sx q[2];
rz(-0.87468481) q[2];
sx q[2];
rz(0.52379215) q[2];
rz(1.1526456) q[3];
sx q[3];
rz(-2.5609784) q[3];
sx q[3];
rz(1.7136278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
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
rz(-1.024217) q[1];
sx q[1];
rz(2.9383235) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.476442) q[0];
sx q[0];
rz(-0.91717952) q[0];
sx q[0];
rz(1.6295432) q[0];
rz(-pi) q[1];
rz(0.16746232) q[2];
sx q[2];
rz(-0.31412087) q[2];
sx q[2];
rz(0.7213074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3422074) q[1];
sx q[1];
rz(-2.8236656) q[1];
sx q[1];
rz(-0.47321508) q[1];
rz(-pi) q[2];
rz(1.2653233) q[3];
sx q[3];
rz(-1.3246228) q[3];
sx q[3];
rz(-0.15428972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.086143494) q[2];
sx q[2];
rz(-1.2148427) q[2];
sx q[2];
rz(2.8622368) q[2];
rz(1.8481567) q[3];
sx q[3];
rz(-1.8297628) q[3];
sx q[3];
rz(1.9259341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16348895) q[0];
sx q[0];
rz(-2.3390529) q[0];
sx q[0];
rz(-1.7247024) q[0];
rz(-0.92719999) q[1];
sx q[1];
rz(-1.5798774) q[1];
sx q[1];
rz(-1.7318116) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3438691) q[0];
sx q[0];
rz(-2.2316405) q[0];
sx q[0];
rz(-3.0190574) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0538834) q[2];
sx q[2];
rz(-1.8894686) q[2];
sx q[2];
rz(2.0234194) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7076787) q[1];
sx q[1];
rz(-2.0053889) q[1];
sx q[1];
rz(0.7264002) q[1];
x q[2];
rz(1.8133598) q[3];
sx q[3];
rz(-1.6117192) q[3];
sx q[3];
rz(0.89311069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37567821) q[2];
sx q[2];
rz(-1.3216852) q[2];
sx q[2];
rz(-2.1709757) q[2];
rz(-2.4961903) q[3];
sx q[3];
rz(-0.97964764) q[3];
sx q[3];
rz(-1.0055044) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29393016) q[0];
sx q[0];
rz(-1.6456589) q[0];
sx q[0];
rz(1.8346067) q[0];
rz(-2.7597799) q[1];
sx q[1];
rz(-1.7996856) q[1];
sx q[1];
rz(-2.3736384) q[1];
rz(-2.3293498) q[2];
sx q[2];
rz(-1.1803738) q[2];
sx q[2];
rz(1.3132172) q[2];
rz(-0.44245023) q[3];
sx q[3];
rz(-1.1493324) q[3];
sx q[3];
rz(2.4449287) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
