OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1825778) q[0];
sx q[0];
rz(-0.20354095) q[0];
sx q[0];
rz(1.8939053) q[0];
rz(1.3522476) q[1];
sx q[1];
rz(-0.69753733) q[1];
sx q[1];
rz(-0.69812671) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20839918) q[0];
sx q[0];
rz(-2.0740447) q[0];
sx q[0];
rz(1.1017975) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1280367) q[2];
sx q[2];
rz(-0.45326158) q[2];
sx q[2];
rz(-1.0176941) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2064849) q[1];
sx q[1];
rz(-2.3000731) q[1];
sx q[1];
rz(-0.40928276) q[1];
rz(-pi) q[2];
rz(1.5145958) q[3];
sx q[3];
rz(-2.8861585) q[3];
sx q[3];
rz(-1.5572302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.84833604) q[2];
sx q[2];
rz(-1.2245327) q[2];
sx q[2];
rz(0.87948925) q[2];
rz(-0.73421156) q[3];
sx q[3];
rz(-2.0101533) q[3];
sx q[3];
rz(-2.3122299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.5010928) q[0];
sx q[0];
rz(-1.9556029) q[0];
sx q[0];
rz(-2.1511141) q[0];
rz(0.99539202) q[1];
sx q[1];
rz(-1.4442911) q[1];
sx q[1];
rz(2.676414) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93056193) q[0];
sx q[0];
rz(-0.4586691) q[0];
sx q[0];
rz(-2.7194897) q[0];
x q[1];
rz(1.0715818) q[2];
sx q[2];
rz(-0.68977654) q[2];
sx q[2];
rz(-2.3782955) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8930186) q[1];
sx q[1];
rz(-1.0842399) q[1];
sx q[1];
rz(2.3793329) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4898759) q[3];
sx q[3];
rz(-0.32518018) q[3];
sx q[3];
rz(-1.0203433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62749949) q[2];
sx q[2];
rz(-1.9072615) q[2];
sx q[2];
rz(0.17847432) q[2];
rz(-0.56525362) q[3];
sx q[3];
rz(-1.7491128) q[3];
sx q[3];
rz(2.5843411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8239215) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(-0.70096651) q[0];
rz(0.40755454) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(2.0661381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0658808) q[0];
sx q[0];
rz(-1.9032818) q[0];
sx q[0];
rz(-2.5677469) q[0];
rz(2.2268217) q[2];
sx q[2];
rz(-1.7505976) q[2];
sx q[2];
rz(0.54043661) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6523931) q[1];
sx q[1];
rz(-1.7013738) q[1];
sx q[1];
rz(-1.9741535) q[1];
rz(-pi) q[2];
rz(1.7549137) q[3];
sx q[3];
rz(-2.872481) q[3];
sx q[3];
rz(-1.5003913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0740697) q[2];
sx q[2];
rz(-1.8068376) q[2];
sx q[2];
rz(1.4074116) q[2];
rz(2.6635026) q[3];
sx q[3];
rz(-0.34308386) q[3];
sx q[3];
rz(-2.7966255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.5673229) q[0];
sx q[0];
rz(-1.747921) q[0];
sx q[0];
rz(-2.0950914) q[0];
rz(-2.8872755) q[1];
sx q[1];
rz(-2.3260702) q[1];
sx q[1];
rz(0.3124803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3182718) q[0];
sx q[0];
rz(-1.4411363) q[0];
sx q[0];
rz(1.1581076) q[0];
x q[1];
rz(0.66019571) q[2];
sx q[2];
rz(-1.9593628) q[2];
sx q[2];
rz(-1.6903433) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8487843) q[1];
sx q[1];
rz(-0.40221805) q[1];
sx q[1];
rz(-2.4821698) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.11361) q[3];
sx q[3];
rz(-0.58251721) q[3];
sx q[3];
rz(1.5503255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8326524) q[2];
sx q[2];
rz(-2.4231484) q[2];
sx q[2];
rz(2.9529875) q[2];
rz(0.23379937) q[3];
sx q[3];
rz(-2.2457687) q[3];
sx q[3];
rz(2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4466062) q[0];
sx q[0];
rz(-1.3124895) q[0];
sx q[0];
rz(-0.0083228668) q[0];
rz(-0.43909973) q[1];
sx q[1];
rz(-1.3689901) q[1];
sx q[1];
rz(1.2459374) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2089691) q[0];
sx q[0];
rz(-1.1198988) q[0];
sx q[0];
rz(-1.3572925) q[0];
rz(-pi) q[1];
rz(1.8077085) q[2];
sx q[2];
rz(-1.4258175) q[2];
sx q[2];
rz(-2.0857834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.627927) q[1];
sx q[1];
rz(-0.49066077) q[1];
sx q[1];
rz(2.2918244) q[1];
rz(1.9381136) q[3];
sx q[3];
rz(-2.4532336) q[3];
sx q[3];
rz(1.3753152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40450725) q[2];
sx q[2];
rz(-3.0948907) q[2];
sx q[2];
rz(1.6339634) q[2];
rz(-0.59219939) q[3];
sx q[3];
rz(-1.7630354) q[3];
sx q[3];
rz(2.4018905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020788766) q[0];
sx q[0];
rz(-1.692166) q[0];
sx q[0];
rz(2.8756496) q[0];
rz(1.2159011) q[1];
sx q[1];
rz(-2.3561056) q[1];
sx q[1];
rz(1.9507834) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.196249) q[0];
sx q[0];
rz(-0.73305819) q[0];
sx q[0];
rz(1.5747914) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036284281) q[2];
sx q[2];
rz(-2.1232833) q[2];
sx q[2];
rz(-0.82922574) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.792508) q[1];
sx q[1];
rz(-2.324632) q[1];
sx q[1];
rz(1.5602099) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9407752) q[3];
sx q[3];
rz(-2.8614223) q[3];
sx q[3];
rz(-0.38129378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8094981) q[2];
sx q[2];
rz(-1.7569434) q[2];
sx q[2];
rz(-2.9936227) q[2];
rz(-0.45251265) q[3];
sx q[3];
rz(-1.0871525) q[3];
sx q[3];
rz(2.7361659) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1389403) q[0];
sx q[0];
rz(-1.2661221) q[0];
sx q[0];
rz(1.7953605) q[0];
rz(3.1345308) q[1];
sx q[1];
rz(-0.66771737) q[1];
sx q[1];
rz(0.5074358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.00822) q[0];
sx q[0];
rz(-2.3485561) q[0];
sx q[0];
rz(-0.80312463) q[0];
rz(1.2260004) q[2];
sx q[2];
rz(-1.424768) q[2];
sx q[2];
rz(-0.70521077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7384711) q[1];
sx q[1];
rz(-0.39408052) q[1];
sx q[1];
rz(-3.1172063) q[1];
x q[2];
rz(-1.5370338) q[3];
sx q[3];
rz(-0.27245263) q[3];
sx q[3];
rz(1.1794943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4355882) q[2];
sx q[2];
rz(-2.0861349) q[2];
sx q[2];
rz(0.17534193) q[2];
rz(0.38241479) q[3];
sx q[3];
rz(-1.1924815) q[3];
sx q[3];
rz(-2.6800938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(3.0328499) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(-0.67799062) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16081339) q[0];
sx q[0];
rz(-2.6475472) q[0];
sx q[0];
rz(1.6176226) q[0];
rz(-pi) q[1];
rz(1.5505232) q[2];
sx q[2];
rz(-2.319699) q[2];
sx q[2];
rz(0.66416364) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1428263) q[1];
sx q[1];
rz(-1.8590392) q[1];
sx q[1];
rz(2.114813) q[1];
rz(-1.0394215) q[3];
sx q[3];
rz(-0.9262923) q[3];
sx q[3];
rz(-2.927305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2737736) q[2];
sx q[2];
rz(-1.5742233) q[2];
sx q[2];
rz(1.9847974) q[2];
rz(0.50446883) q[3];
sx q[3];
rz(-2.1095095) q[3];
sx q[3];
rz(0.57197905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.93541637) q[0];
sx q[0];
rz(-1.6125866) q[0];
sx q[0];
rz(0.61844283) q[0];
rz(-0.84856021) q[1];
sx q[1];
rz(-0.51785523) q[1];
sx q[1];
rz(2.3458164) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37695161) q[0];
sx q[0];
rz(-1.5466154) q[0];
sx q[0];
rz(-2.0410246) q[0];
rz(-pi) q[1];
rz(0.047330476) q[2];
sx q[2];
rz(-0.82529699) q[2];
sx q[2];
rz(-1.8290486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3251212) q[1];
sx q[1];
rz(-2.7181648) q[1];
sx q[1];
rz(-2.2864443) q[1];
x q[2];
rz(-1.9913748) q[3];
sx q[3];
rz(-1.707786) q[3];
sx q[3];
rz(1.610422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.985118) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(-1.1953243) q[2];
rz(-0.34504238) q[3];
sx q[3];
rz(-0.86193591) q[3];
sx q[3];
rz(-2.2752458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7095551) q[0];
sx q[0];
rz(-0.23831743) q[0];
sx q[0];
rz(0.077614345) q[0];
rz(-1.2331102) q[1];
sx q[1];
rz(-2.1014919) q[1];
sx q[1];
rz(-0.79201039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4778276) q[0];
sx q[0];
rz(-2.223513) q[0];
sx q[0];
rz(-2.2774612) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77113232) q[2];
sx q[2];
rz(-2.9206616) q[2];
sx q[2];
rz(0.67459092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4153667) q[1];
sx q[1];
rz(-2.0322406) q[1];
sx q[1];
rz(0.17166328) q[1];
rz(-0.80371334) q[3];
sx q[3];
rz(-1.6559542) q[3];
sx q[3];
rz(-0.61060315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9218257) q[2];
sx q[2];
rz(-0.91138387) q[2];
sx q[2];
rz(1.0731953) q[2];
rz(-2.8940767) q[3];
sx q[3];
rz(-1.9297587) q[3];
sx q[3];
rz(0.016935067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76219227) q[0];
sx q[0];
rz(-1.8397377) q[0];
sx q[0];
rz(2.4158438) q[0];
rz(2.2629867) q[1];
sx q[1];
rz(-0.85522063) q[1];
sx q[1];
rz(-1.9488889) q[1];
rz(-2.939777) q[2];
sx q[2];
rz(-2.0463662) q[2];
sx q[2];
rz(-3.0007888) q[2];
rz(-0.056817354) q[3];
sx q[3];
rz(-0.84470261) q[3];
sx q[3];
rz(-0.62538996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
