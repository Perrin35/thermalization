OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.81808972) q[0];
sx q[0];
rz(-0.40677318) q[0];
sx q[0];
rz(0.50954252) q[0];
rz(-0.31515631) q[1];
sx q[1];
rz(6.4767467) q[1];
sx q[1];
rz(11.560796) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5307546) q[0];
sx q[0];
rz(-1.6456398) q[0];
sx q[0];
rz(-0.2006528) q[0];
rz(-pi) q[1];
x q[1];
rz(1.124473) q[2];
sx q[2];
rz(-0.69689489) q[2];
sx q[2];
rz(2.6817037) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79709681) q[1];
sx q[1];
rz(-0.90289799) q[1];
sx q[1];
rz(-2.383166) q[1];
rz(-0.77862424) q[3];
sx q[3];
rz(-1.1168241) q[3];
sx q[3];
rz(0.046128143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2931508) q[2];
sx q[2];
rz(-1.2520496) q[2];
sx q[2];
rz(-1.1085917) q[2];
rz(-2.0075924) q[3];
sx q[3];
rz(-2.0847376) q[3];
sx q[3];
rz(0.081324287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.14652458) q[0];
sx q[0];
rz(-2.0159371) q[0];
sx q[0];
rz(-2.8296237) q[0];
rz(-2.9106855) q[1];
sx q[1];
rz(-2.0377908) q[1];
sx q[1];
rz(-2.013496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474079) q[0];
sx q[0];
rz(-2.0387406) q[0];
sx q[0];
rz(0.2121966) q[0];
x q[1];
rz(-2.0064193) q[2];
sx q[2];
rz(-0.95016236) q[2];
sx q[2];
rz(2.1114608) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.384239) q[1];
sx q[1];
rz(-1.6815255) q[1];
sx q[1];
rz(-2.6925681) q[1];
x q[2];
rz(0.45421447) q[3];
sx q[3];
rz(-0.81791211) q[3];
sx q[3];
rz(2.1061768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6850623) q[2];
sx q[2];
rz(-1.637746) q[2];
sx q[2];
rz(-3.077363) q[2];
rz(-1.8188933) q[3];
sx q[3];
rz(-2.5048246) q[3];
sx q[3];
rz(0.88671154) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54881683) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(1.2491666) q[0];
rz(-0.33117548) q[1];
sx q[1];
rz(-2.0024029) q[1];
sx q[1];
rz(-1.9452728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4713269) q[0];
sx q[0];
rz(-1.571079) q[0];
sx q[0];
rz(3.1413881) q[0];
rz(-3.1205503) q[2];
sx q[2];
rz(-2.8940563) q[2];
sx q[2];
rz(0.68138441) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1332449) q[1];
sx q[1];
rz(-2.3143594) q[1];
sx q[1];
rz(-0.72669795) q[1];
rz(-0.75342859) q[3];
sx q[3];
rz(-1.2810858) q[3];
sx q[3];
rz(-3.0353607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.62200481) q[2];
sx q[2];
rz(-0.76498166) q[2];
sx q[2];
rz(-2.2288442) q[2];
rz(-1.4384455) q[3];
sx q[3];
rz(-1.6310952) q[3];
sx q[3];
rz(-1.335817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66252935) q[0];
sx q[0];
rz(-2.0344069) q[0];
sx q[0];
rz(2.4712439) q[0];
rz(-2.4743075) q[1];
sx q[1];
rz(-2.1332462) q[1];
sx q[1];
rz(-0.60428062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8025148) q[0];
sx q[0];
rz(-2.3065563) q[0];
sx q[0];
rz(1.5123532) q[0];
rz(-pi) q[1];
rz(-1.7261581) q[2];
sx q[2];
rz(-1.948602) q[2];
sx q[2];
rz(-0.31766674) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62897077) q[1];
sx q[1];
rz(-1.781946) q[1];
sx q[1];
rz(-1.3107915) q[1];
x q[2];
rz(-2.9168455) q[3];
sx q[3];
rz(-0.852727) q[3];
sx q[3];
rz(-1.3026733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0303354) q[2];
sx q[2];
rz(-1.5680485) q[2];
sx q[2];
rz(-2.7643909) q[2];
rz(3.0180569) q[3];
sx q[3];
rz(-1.7081407) q[3];
sx q[3];
rz(1.4089353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971624) q[0];
sx q[0];
rz(-0.57752174) q[0];
sx q[0];
rz(0.29931983) q[0];
rz(-1.1209283) q[1];
sx q[1];
rz(-1.2052373) q[1];
sx q[1];
rz(2.1790806) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4967921) q[0];
sx q[0];
rz(-1.776665) q[0];
sx q[0];
rz(-1.0394761) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91470529) q[2];
sx q[2];
rz(-1.9491299) q[2];
sx q[2];
rz(-1.8112184) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1735517) q[1];
sx q[1];
rz(-2.3690201) q[1];
sx q[1];
rz(0.79141683) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5425517) q[3];
sx q[3];
rz(-1.1223464) q[3];
sx q[3];
rz(-0.34235172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0943429) q[2];
sx q[2];
rz(-0.80545682) q[2];
sx q[2];
rz(-2.1979525) q[2];
rz(-1.0654248) q[3];
sx q[3];
rz(-1.7908955) q[3];
sx q[3];
rz(-1.0984727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0672673) q[0];
sx q[0];
rz(-1.704498) q[0];
sx q[0];
rz(-3.1332916) q[0];
rz(-0.51757327) q[1];
sx q[1];
rz(-0.65474302) q[1];
sx q[1];
rz(1.5932721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8529274) q[0];
sx q[0];
rz(-1.2199243) q[0];
sx q[0];
rz(-0.025625833) q[0];
rz(1.8528884) q[2];
sx q[2];
rz(-1.0824426) q[2];
sx q[2];
rz(2.2887678) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0071809) q[1];
sx q[1];
rz(-2.6406277) q[1];
sx q[1];
rz(2.3208614) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7260439) q[3];
sx q[3];
rz(-0.99667785) q[3];
sx q[3];
rz(2.1972063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3090618) q[2];
sx q[2];
rz(-1.6522202) q[2];
sx q[2];
rz(1.8050516) q[2];
rz(-3.0858223) q[3];
sx q[3];
rz(-0.99168188) q[3];
sx q[3];
rz(-0.43558863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5585153) q[0];
sx q[0];
rz(-0.4774839) q[0];
sx q[0];
rz(-2.4537295) q[0];
rz(1.4631924) q[1];
sx q[1];
rz(-0.72931591) q[1];
sx q[1];
rz(-2.8942143) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5106544) q[0];
sx q[0];
rz(-1.3888089) q[0];
sx q[0];
rz(-2.8923558) q[0];
rz(-pi) q[1];
rz(-2.192657) q[2];
sx q[2];
rz(-2.1315711) q[2];
sx q[2];
rz(-1.6552629) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.069889594) q[1];
sx q[1];
rz(-2.1462198) q[1];
sx q[1];
rz(-2.6399122) q[1];
rz(-1.803627) q[3];
sx q[3];
rz(-1.4484753) q[3];
sx q[3];
rz(-1.2965152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4443724) q[2];
sx q[2];
rz(-1.8070544) q[2];
sx q[2];
rz(-0.42416254) q[2];
rz(1.0423202) q[3];
sx q[3];
rz(-0.76516953) q[3];
sx q[3];
rz(-0.55646363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384077) q[0];
sx q[0];
rz(-0.98588949) q[0];
sx q[0];
rz(2.5307122) q[0];
rz(-0.25018397) q[1];
sx q[1];
rz(-1.4004204) q[1];
sx q[1];
rz(0.47666034) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31075031) q[0];
sx q[0];
rz(-2.0526969) q[0];
sx q[0];
rz(-3.0974814) q[0];
x q[1];
rz(-1.3872434) q[2];
sx q[2];
rz(-1.6590365) q[2];
sx q[2];
rz(-2.1004408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6256801) q[1];
sx q[1];
rz(-1.971174) q[1];
sx q[1];
rz(-2.6039586) q[1];
rz(-pi) q[2];
rz(2.1564756) q[3];
sx q[3];
rz(-1.9352311) q[3];
sx q[3];
rz(1.6062615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1499947) q[2];
sx q[2];
rz(-1.0341045) q[2];
sx q[2];
rz(-2.4617713) q[2];
rz(0.46197915) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(-1.6405039) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3221472) q[0];
sx q[0];
rz(-1.3752022) q[0];
sx q[0];
rz(-0.15400259) q[0];
rz(0.18140659) q[1];
sx q[1];
rz(-2.0712974) q[1];
sx q[1];
rz(-3.0214686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4331872) q[0];
sx q[0];
rz(-2.7338203) q[0];
sx q[0];
rz(1.7208517) q[0];
rz(-1.3072877) q[2];
sx q[2];
rz(-1.7611836) q[2];
sx q[2];
rz(0.15979494) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5654345) q[1];
sx q[1];
rz(-1.6685969) q[1];
sx q[1];
rz(-1.3937772) q[1];
x q[2];
rz(-0.24415827) q[3];
sx q[3];
rz(-2.6006581) q[3];
sx q[3];
rz(1.8997826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4274365) q[2];
sx q[2];
rz(-1.3849881) q[2];
sx q[2];
rz(-0.45905054) q[2];
rz(0.51042026) q[3];
sx q[3];
rz(-2.3000058) q[3];
sx q[3];
rz(-0.69971219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089461483) q[0];
sx q[0];
rz(-0.94877807) q[0];
sx q[0];
rz(-0.38598886) q[0];
rz(1.5486859) q[1];
sx q[1];
rz(-0.49647757) q[1];
sx q[1];
rz(-1.5492424) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95741877) q[0];
sx q[0];
rz(-2.9702219) q[0];
sx q[0];
rz(-1.7383854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73076325) q[2];
sx q[2];
rz(-0.77789069) q[2];
sx q[2];
rz(0.90532263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4262271) q[1];
sx q[1];
rz(-1.5439529) q[1];
sx q[1];
rz(-1.6728841) q[1];
x q[2];
rz(-0.066927197) q[3];
sx q[3];
rz(-1.4176344) q[3];
sx q[3];
rz(-2.7471971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3647032) q[2];
sx q[2];
rz(-2.2025755) q[2];
sx q[2];
rz(-1.6949863) q[2];
rz(-2.1052836) q[3];
sx q[3];
rz(-2.2615137) q[3];
sx q[3];
rz(-0.026329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.885289) q[0];
sx q[0];
rz(-2.1623609) q[0];
sx q[0];
rz(1.4872861) q[0];
rz(-0.22008315) q[1];
sx q[1];
rz(-1.1047803) q[1];
sx q[1];
rz(1.6688375) q[1];
rz(0.34192495) q[2];
sx q[2];
rz(-2.3734063) q[2];
sx q[2];
rz(-2.0362324) q[2];
rz(2.6406399) q[3];
sx q[3];
rz(-1.8617478) q[3];
sx q[3];
rz(-0.72248722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
