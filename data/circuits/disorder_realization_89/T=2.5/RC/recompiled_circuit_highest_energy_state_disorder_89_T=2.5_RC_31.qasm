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
rz(0.97857082) q[0];
sx q[0];
rz(3.0244654) q[0];
sx q[0];
rz(9.2000118) q[0];
rz(-0.4966785) q[1];
sx q[1];
rz(-1.574312) q[1];
sx q[1];
rz(1.1800676) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3008904) q[0];
sx q[0];
rz(-1.7127174) q[0];
sx q[0];
rz(3.0034756) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44344896) q[2];
sx q[2];
rz(-1.5649619) q[2];
sx q[2];
rz(-1.3590517) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2121316) q[1];
sx q[1];
rz(-2.5358905) q[1];
sx q[1];
rz(-0.49537658) q[1];
rz(-pi) q[2];
rz(1.7430844) q[3];
sx q[3];
rz(-2.7040561) q[3];
sx q[3];
rz(0.019339081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4687389) q[2];
sx q[2];
rz(-1.8127952) q[2];
sx q[2];
rz(1.9358181) q[2];
rz(-0.85637158) q[3];
sx q[3];
rz(-2.2090293) q[3];
sx q[3];
rz(-0.87730733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641651) q[0];
sx q[0];
rz(-2.7441315) q[0];
sx q[0];
rz(-2.933266) q[0];
rz(1.9147929) q[1];
sx q[1];
rz(-2.0715163) q[1];
sx q[1];
rz(0.56455451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78075942) q[0];
sx q[0];
rz(-1.5281023) q[0];
sx q[0];
rz(-1.7303995) q[0];
rz(2.4103569) q[2];
sx q[2];
rz(-1.952716) q[2];
sx q[2];
rz(2.3637091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7617088) q[1];
sx q[1];
rz(-1.392388) q[1];
sx q[1];
rz(-0.70990035) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0309593) q[3];
sx q[3];
rz(-1.587365) q[3];
sx q[3];
rz(1.4318717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5047001) q[2];
sx q[2];
rz(-1.8919287) q[2];
sx q[2];
rz(1.4862109) q[2];
rz(0.49556035) q[3];
sx q[3];
rz(-1.7335745) q[3];
sx q[3];
rz(1.0068309) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946852) q[0];
sx q[0];
rz(-2.4201396) q[0];
sx q[0];
rz(-1.7578693) q[0];
rz(1.9231298) q[1];
sx q[1];
rz(-1.0885025) q[1];
sx q[1];
rz(-2.7168435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.527194) q[0];
sx q[0];
rz(-1.5081811) q[0];
sx q[0];
rz(-1.539143) q[0];
rz(2.3282503) q[2];
sx q[2];
rz(-2.7995297) q[2];
sx q[2];
rz(-1.4881211) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37990824) q[1];
sx q[1];
rz(-1.1947922) q[1];
sx q[1];
rz(-1.6523182) q[1];
x q[2];
rz(-0.056428595) q[3];
sx q[3];
rz(-2.5241969) q[3];
sx q[3];
rz(-1.4790725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4212627) q[2];
sx q[2];
rz(-2.0739372) q[2];
sx q[2];
rz(-3.1171411) q[2];
rz(1.1484185) q[3];
sx q[3];
rz(-2.0423856) q[3];
sx q[3];
rz(-2.059977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.9352683) q[0];
sx q[0];
rz(-0.22980389) q[0];
sx q[0];
rz(-2.9277053) q[0];
rz(-2.3230486) q[1];
sx q[1];
rz(-1.7844424) q[1];
sx q[1];
rz(-2.4549386) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59348561) q[0];
sx q[0];
rz(-2.1186376) q[0];
sx q[0];
rz(2.0440897) q[0];
rz(0.29149194) q[2];
sx q[2];
rz(-2.1003869) q[2];
sx q[2];
rz(-1.518371) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80685341) q[1];
sx q[1];
rz(-0.76419965) q[1];
sx q[1];
rz(-3.0055226) q[1];
x q[2];
rz(-2.2591001) q[3];
sx q[3];
rz(-0.63418856) q[3];
sx q[3];
rz(-0.89804441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0246747) q[2];
sx q[2];
rz(-1.4093829) q[2];
sx q[2];
rz(-2.4521258) q[2];
rz(1.3245448) q[3];
sx q[3];
rz(-1.2169415) q[3];
sx q[3];
rz(2.569258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784742) q[0];
sx q[0];
rz(-0.093570396) q[0];
sx q[0];
rz(2.1043188) q[0];
rz(-2.4689238) q[1];
sx q[1];
rz(-1.172352) q[1];
sx q[1];
rz(2.345828) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0507495) q[0];
sx q[0];
rz(-1.2552869) q[0];
sx q[0];
rz(2.1054224) q[0];
rz(-1.1188385) q[2];
sx q[2];
rz(-2.5198433) q[2];
sx q[2];
rz(0.62134075) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2356253) q[1];
sx q[1];
rz(-0.86938953) q[1];
sx q[1];
rz(-0.4404621) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2656059) q[3];
sx q[3];
rz(-0.92545907) q[3];
sx q[3];
rz(-2.2923645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98046389) q[2];
sx q[2];
rz(-1.0995862) q[2];
sx q[2];
rz(-0.83578342) q[2];
rz(-0.25660822) q[3];
sx q[3];
rz(-2.0308688) q[3];
sx q[3];
rz(-0.49032828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68818727) q[0];
sx q[0];
rz(-1.9904933) q[0];
sx q[0];
rz(1.2123464) q[0];
rz(1.0351828) q[1];
sx q[1];
rz(-1.4915219) q[1];
sx q[1];
rz(-2.129668) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2741756) q[0];
sx q[0];
rz(-1.7417481) q[0];
sx q[0];
rz(-0.11049185) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8746568) q[2];
sx q[2];
rz(-1.5812318) q[2];
sx q[2];
rz(0.70994678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95532596) q[1];
sx q[1];
rz(-1.9791934) q[1];
sx q[1];
rz(0.46640654) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7583241) q[3];
sx q[3];
rz(-1.82194) q[3];
sx q[3];
rz(2.0519837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47784352) q[2];
sx q[2];
rz(-2.9426844) q[2];
sx q[2];
rz(-0.94856962) q[2];
rz(2.6712724) q[3];
sx q[3];
rz(-2.0430653) q[3];
sx q[3];
rz(1.2547913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54742852) q[0];
sx q[0];
rz(-1.0110039) q[0];
sx q[0];
rz(1.8773361) q[0];
rz(0.7803548) q[1];
sx q[1];
rz(-0.64096132) q[1];
sx q[1];
rz(2.9497214) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9148798) q[0];
sx q[0];
rz(-2.5396945) q[0];
sx q[0];
rz(0.65897091) q[0];
rz(0.81249313) q[2];
sx q[2];
rz(-0.50069649) q[2];
sx q[2];
rz(2.814584) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.807232) q[1];
sx q[1];
rz(-1.8340115) q[1];
sx q[1];
rz(-0.46137793) q[1];
rz(0.34911134) q[3];
sx q[3];
rz(-2.1421332) q[3];
sx q[3];
rz(3.0219363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4912305) q[2];
sx q[2];
rz(-1.0543084) q[2];
sx q[2];
rz(2.3336156) q[2];
rz(-0.84336495) q[3];
sx q[3];
rz(-1.9873514) q[3];
sx q[3];
rz(-1.5850867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845487) q[0];
sx q[0];
rz(-3.1175685) q[0];
sx q[0];
rz(2.3727544) q[0];
rz(-2.8688042) q[1];
sx q[1];
rz(-1.5130679) q[1];
sx q[1];
rz(1.3974894) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834806) q[0];
sx q[0];
rz(-1.8761823) q[0];
sx q[0];
rz(2.7598513) q[0];
x q[1];
rz(-1.0553841) q[2];
sx q[2];
rz(-1.2920818) q[2];
sx q[2];
rz(2.9936522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2915512) q[1];
sx q[1];
rz(-1.5102036) q[1];
sx q[1];
rz(-1.4539976) q[1];
rz(-pi) q[2];
rz(2.4726953) q[3];
sx q[3];
rz(-1.9188768) q[3];
sx q[3];
rz(2.477248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93691319) q[2];
sx q[2];
rz(-2.0759373) q[2];
sx q[2];
rz(-1.5020471) q[2];
rz(1.3192734) q[3];
sx q[3];
rz(-2.2917512) q[3];
sx q[3];
rz(-1.5892861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-2.4424292) q[0];
sx q[0];
rz(-0.84446877) q[0];
sx q[0];
rz(2.1635639) q[0];
rz(-2.6507822) q[1];
sx q[1];
rz(-1.699479) q[1];
sx q[1];
rz(1.4960272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.967088) q[0];
sx q[0];
rz(-1.8068815) q[0];
sx q[0];
rz(1.4437136) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7285887) q[2];
sx q[2];
rz(-1.4133412) q[2];
sx q[2];
rz(-2.8812489) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89801187) q[1];
sx q[1];
rz(-1.812938) q[1];
sx q[1];
rz(-2.8088074) q[1];
rz(-pi) q[2];
rz(-2.6096212) q[3];
sx q[3];
rz(-1.6364179) q[3];
sx q[3];
rz(2.402879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8533123) q[2];
sx q[2];
rz(-0.86709443) q[2];
sx q[2];
rz(1.8193998) q[2];
rz(-2.7754122) q[3];
sx q[3];
rz(-1.7145232) q[3];
sx q[3];
rz(-2.0100994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2243097) q[0];
sx q[0];
rz(-1.6095105) q[0];
sx q[0];
rz(-3.068058) q[0];
rz(-1.8247983) q[1];
sx q[1];
rz(-1.8837181) q[1];
sx q[1];
rz(-1.8427461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4816581) q[0];
sx q[0];
rz(-1.3941688) q[0];
sx q[0];
rz(0.28089471) q[0];
rz(-0.13330524) q[2];
sx q[2];
rz(-2.0553737) q[2];
sx q[2];
rz(-0.43371782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6886814) q[1];
sx q[1];
rz(-1.1961403) q[1];
sx q[1];
rz(-2.8883347) q[1];
rz(0.73963005) q[3];
sx q[3];
rz(-1.4907903) q[3];
sx q[3];
rz(2.7440939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1765882) q[2];
sx q[2];
rz(-2.1829288) q[2];
sx q[2];
rz(2.5625572) q[2];
rz(1.6252919) q[3];
sx q[3];
rz(-0.54280353) q[3];
sx q[3];
rz(-1.496284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8662921) q[0];
sx q[0];
rz(-1.0466812) q[0];
sx q[0];
rz(-0.064591797) q[0];
rz(-2.2304089) q[1];
sx q[1];
rz(-1.8235527) q[1];
sx q[1];
rz(1.8240737) q[1];
rz(0.42413128) q[2];
sx q[2];
rz(-1.5739706) q[2];
sx q[2];
rz(-1.7375873) q[2];
rz(0.95355036) q[3];
sx q[3];
rz(-1.336477) q[3];
sx q[3];
rz(-0.82930641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
