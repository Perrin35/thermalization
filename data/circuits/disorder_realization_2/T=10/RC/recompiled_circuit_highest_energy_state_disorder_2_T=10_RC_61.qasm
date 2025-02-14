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
rz(-0.18345565) q[0];
sx q[0];
rz(-0.7440716) q[0];
sx q[0];
rz(-1.0519354) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(2.5084578) q[1];
sx q[1];
rz(8.8517744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7318037) q[0];
sx q[0];
rz(-2.5136726) q[0];
sx q[0];
rz(1.464792) q[0];
rz(-pi) q[1];
rz(2.0016167) q[2];
sx q[2];
rz(-2.2099021) q[2];
sx q[2];
rz(0.72102816) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0172559) q[1];
sx q[1];
rz(-1.6410429) q[1];
sx q[1];
rz(0.4976561) q[1];
x q[2];
rz(-0.12270452) q[3];
sx q[3];
rz(-2.4413681) q[3];
sx q[3];
rz(2.6168857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28213349) q[2];
sx q[2];
rz(-2.8318475) q[2];
sx q[2];
rz(2.0852883) q[2];
rz(0.022196444) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(-1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2317155) q[0];
sx q[0];
rz(-0.48119369) q[0];
sx q[0];
rz(-0.43014446) q[0];
rz(-0.12769708) q[1];
sx q[1];
rz(-1.1558497) q[1];
sx q[1];
rz(1.7040303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50629726) q[0];
sx q[0];
rz(-2.7397635) q[0];
sx q[0];
rz(0.70962972) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84354894) q[2];
sx q[2];
rz(-3.0650716) q[2];
sx q[2];
rz(0.25329548) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75152698) q[1];
sx q[1];
rz(-1.2428778) q[1];
sx q[1];
rz(-0.71050127) q[1];
x q[2];
rz(1.7875769) q[3];
sx q[3];
rz(-1.0415029) q[3];
sx q[3];
rz(2.2640219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11640707) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(-2.6961668) q[2];
rz(2.7323501) q[3];
sx q[3];
rz(-2.0425115) q[3];
sx q[3];
rz(2.2621431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55260783) q[0];
sx q[0];
rz(-3.0486076) q[0];
sx q[0];
rz(-0.71480042) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-0.42066586) q[1];
sx q[1];
rz(-0.32726273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48178534) q[0];
sx q[0];
rz(-2.5707173) q[0];
sx q[0];
rz(-1.8828189) q[0];
rz(-2.1398323) q[2];
sx q[2];
rz(-1.1998917) q[2];
sx q[2];
rz(-0.47874622) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5048557) q[1];
sx q[1];
rz(-1.550192) q[1];
sx q[1];
rz(1.7223174) q[1];
rz(1.2233569) q[3];
sx q[3];
rz(-1.8104324) q[3];
sx q[3];
rz(-2.6289866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0032234) q[2];
sx q[2];
rz(-2.1681163) q[2];
sx q[2];
rz(1.6376015) q[2];
rz(-1.9231632) q[3];
sx q[3];
rz(-0.4156433) q[3];
sx q[3];
rz(2.159582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0729436) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(-2.9856227) q[0];
rz(0.34863696) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(-1.2329996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1240163) q[0];
sx q[0];
rz(-1.6917545) q[0];
sx q[0];
rz(1.8328387) q[0];
rz(-0.92278752) q[2];
sx q[2];
rz(-1.7402667) q[2];
sx q[2];
rz(-0.10709912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0291042) q[1];
sx q[1];
rz(-1.8524516) q[1];
sx q[1];
rz(-0.92093484) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8901222) q[3];
sx q[3];
rz(-1.5565539) q[3];
sx q[3];
rz(-2.3492658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4431346) q[2];
sx q[2];
rz(-3.105574) q[2];
sx q[2];
rz(1.6301463) q[2];
rz(2.002142) q[3];
sx q[3];
rz(-1.8099433) q[3];
sx q[3];
rz(0.99036923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.776942) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(1.1530217) q[0];
rz(-0.54368377) q[1];
sx q[1];
rz(-1.8749571) q[1];
sx q[1];
rz(-0.26184729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1516929) q[0];
sx q[0];
rz(-1.6542871) q[0];
sx q[0];
rz(-3.1079453) q[0];
rz(0.98723094) q[2];
sx q[2];
rz(-1.2886184) q[2];
sx q[2];
rz(-0.63167494) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8229554) q[1];
sx q[1];
rz(-1.3483817) q[1];
sx q[1];
rz(-2.5435124) q[1];
rz(-pi) q[2];
rz(2.35942) q[3];
sx q[3];
rz(-1.7165136) q[3];
sx q[3];
rz(2.4268933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5220945) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(-1.7363133) q[2];
rz(2.3960579) q[3];
sx q[3];
rz(-1.8757952) q[3];
sx q[3];
rz(-0.94329992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(0.0010506823) q[0];
sx q[0];
rz(-3.0241835) q[0];
sx q[0];
rz(0.35181272) q[0];
rz(-2.70128) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(2.9702759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85799828) q[0];
sx q[0];
rz(-2.4263315) q[0];
sx q[0];
rz(-2.7659228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75320737) q[2];
sx q[2];
rz(-1.79617) q[2];
sx q[2];
rz(-1.6373375) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.80444634) q[1];
sx q[1];
rz(-2.4427572) q[1];
sx q[1];
rz(-0.032445907) q[1];
rz(-pi) q[2];
rz(2.3790509) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(0.97555893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.064284023) q[2];
sx q[2];
rz(-1.7694387) q[2];
sx q[2];
rz(2.1997931) q[2];
rz(-1.1549548) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(1.340516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6963541) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(-3.136694) q[0];
rz(-2.694963) q[1];
sx q[1];
rz(-2.547867) q[1];
sx q[1];
rz(-0.68797025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5410091) q[0];
sx q[0];
rz(-2.2760411) q[0];
sx q[0];
rz(1.7974822) q[0];
x q[1];
rz(-1.2254459) q[2];
sx q[2];
rz(-0.98042578) q[2];
sx q[2];
rz(1.8777868) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0258642) q[1];
sx q[1];
rz(-2.0021571) q[1];
sx q[1];
rz(-0.23484767) q[1];
x q[2];
rz(-2.4554208) q[3];
sx q[3];
rz(-2.1451575) q[3];
sx q[3];
rz(-0.8818834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.530431) q[2];
sx q[2];
rz(-2.7335584) q[2];
sx q[2];
rz(2.2116057) q[2];
rz(1.9200578) q[3];
sx q[3];
rz(-1.113021) q[3];
sx q[3];
rz(-2.5482224) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68194836) q[0];
sx q[0];
rz(-1.3197897) q[0];
sx q[0];
rz(-0.090106877) q[0];
rz(1.4247165) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(-2.7065014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51564901) q[0];
sx q[0];
rz(-1.0727779) q[0];
sx q[0];
rz(1.7367212) q[0];
rz(1.3899562) q[2];
sx q[2];
rz(-0.9107843) q[2];
sx q[2];
rz(-1.7969799) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4415978) q[1];
sx q[1];
rz(-2.1837103) q[1];
sx q[1];
rz(-0.15721486) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1396542) q[3];
sx q[3];
rz(-2.7172305) q[3];
sx q[3];
rz(-0.75943179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4789751) q[2];
sx q[2];
rz(-1.0650029) q[2];
sx q[2];
rz(-0.62620658) q[2];
rz(0.19691697) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(1.7530493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597252) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(3.0873121) q[0];
rz(-0.42298969) q[1];
sx q[1];
rz(-1.7661679) q[1];
sx q[1];
rz(1.8064226) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8322382) q[0];
sx q[0];
rz(-2.1119499) q[0];
sx q[0];
rz(0.65017976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35690709) q[2];
sx q[2];
rz(-1.0400598) q[2];
sx q[2];
rz(2.4414947) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21653342) q[1];
sx q[1];
rz(-1.0581985) q[1];
sx q[1];
rz(-0.61116265) q[1];
x q[2];
rz(1.4328721) q[3];
sx q[3];
rz(-1.63878) q[3];
sx q[3];
rz(1.8658416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6173031) q[2];
sx q[2];
rz(-1.4681939) q[2];
sx q[2];
rz(-0.35150251) q[2];
rz(-1.7678123) q[3];
sx q[3];
rz(-0.51353729) q[3];
sx q[3];
rz(-2.8721749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7518625) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(-0.75916284) q[0];
rz(2.0579386) q[1];
sx q[1];
rz(-1.6108797) q[1];
sx q[1];
rz(-1.0677451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6408975) q[0];
sx q[0];
rz(-1.0984549) q[0];
sx q[0];
rz(2.2233836) q[0];
x q[1];
rz(-1.8482089) q[2];
sx q[2];
rz(-1.2870064) q[2];
sx q[2];
rz(2.0744155) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0800277) q[1];
sx q[1];
rz(-1.0288218) q[1];
sx q[1];
rz(0.47596495) q[1];
x q[2];
rz(2.8078733) q[3];
sx q[3];
rz(-0.90147831) q[3];
sx q[3];
rz(0.02726083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7119673) q[2];
sx q[2];
rz(-1.2099313) q[2];
sx q[2];
rz(2.9313226) q[2];
rz(0.78768864) q[3];
sx q[3];
rz(-1.382788) q[3];
sx q[3];
rz(-0.53019607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6982211) q[0];
sx q[0];
rz(-2.5826695) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(1.5086077) q[1];
sx q[1];
rz(-1.2651545) q[1];
sx q[1];
rz(-1.9427585) q[1];
rz(2.2589113) q[2];
sx q[2];
rz(-1.8919049) q[2];
sx q[2];
rz(1.6406825) q[2];
rz(0.15565025) q[3];
sx q[3];
rz(-1.2990824) q[3];
sx q[3];
rz(-1.6344447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
