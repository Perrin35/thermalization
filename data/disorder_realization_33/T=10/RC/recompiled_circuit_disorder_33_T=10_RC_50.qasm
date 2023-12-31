OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(1.8338058) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21068621) q[0];
sx q[0];
rz(-1.2331729) q[0];
sx q[0];
rz(0.36436413) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3762796) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(-0.050616654) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1085514) q[1];
sx q[1];
rz(-1.3290977) q[1];
sx q[1];
rz(-0.34717314) q[1];
rz(-pi) q[2];
rz(1.0899815) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(-2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(-2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(0.18584132) q[0];
rz(-2.5813685) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(-2.9247608) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502055) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(-2.015381) q[0];
x q[1];
rz(-0.47302834) q[2];
sx q[2];
rz(-2.0748667) q[2];
sx q[2];
rz(-0.31993983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.195897) q[1];
sx q[1];
rz(-0.90598124) q[1];
sx q[1];
rz(-2.871454) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5012245) q[3];
sx q[3];
rz(-0.98207563) q[3];
sx q[3];
rz(2.0224188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.310114) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.2878093) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(2.2089829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5843825) q[0];
sx q[0];
rz(-2.819448) q[0];
sx q[0];
rz(1.4003217) q[0];
rz(-pi) q[1];
rz(-1.7980174) q[2];
sx q[2];
rz(-1.7482687) q[2];
sx q[2];
rz(0.15572671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.052913594) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(-1.710379) q[1];
rz(-pi) q[2];
rz(-1.014939) q[3];
sx q[3];
rz(-1.9564637) q[3];
sx q[3];
rz(-2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.147826) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(1.0926584) q[2];
rz(0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-0.96737635) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(-2.6413667) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(-1.4979699) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.362975) q[0];
sx q[0];
rz(-0.59016363) q[0];
sx q[0];
rz(1.0233364) q[0];
rz(-1.7613212) q[2];
sx q[2];
rz(-1.5152144) q[2];
sx q[2];
rz(1.3560825) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8151617) q[1];
sx q[1];
rz(-1.8306499) q[1];
sx q[1];
rz(-1.3304779) q[1];
x q[2];
rz(-0.76969947) q[3];
sx q[3];
rz(-2.5315428) q[3];
sx q[3];
rz(1.1806012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.74636373) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(-0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(1.048208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883023) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(0.81272965) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1254966) q[2];
sx q[2];
rz(-0.65158366) q[2];
sx q[2];
rz(-0.96166699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2506927) q[1];
sx q[1];
rz(-2.7437468) q[1];
sx q[1];
rz(2.489151) q[1];
x q[2];
rz(0.95543315) q[3];
sx q[3];
rz(-1.2293929) q[3];
sx q[3];
rz(1.1036901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6158225) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(2.0416416) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(-0.88551372) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(-1.0166608) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(3.0117603) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31045612) q[0];
sx q[0];
rz(-1.9400915) q[0];
sx q[0];
rz(0.62311689) q[0];
x q[1];
rz(0.99545698) q[2];
sx q[2];
rz(-1.2577004) q[2];
sx q[2];
rz(-2.7196333) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0778724) q[1];
sx q[1];
rz(-1.2076326) q[1];
sx q[1];
rz(0.6086463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5927605) q[3];
sx q[3];
rz(-1.4026814) q[3];
sx q[3];
rz(-0.43743922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(2.9373346) q[2];
rz(1.2060818) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(-0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.7234574) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9335564) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(2.3963388) q[0];
rz(-pi) q[1];
rz(-2.5574066) q[2];
sx q[2];
rz(-2.2224732) q[2];
sx q[2];
rz(-0.78648957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64993091) q[1];
sx q[1];
rz(-1.4964536) q[1];
sx q[1];
rz(-1.7564303) q[1];
x q[2];
rz(2.7191914) q[3];
sx q[3];
rz(-1.8654612) q[3];
sx q[3];
rz(0.15448031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(0.0017722842) q[2];
rz(2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(1.6954533) q[0];
rz(2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(1.6400281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1743463) q[0];
sx q[0];
rz(-1.5902728) q[0];
sx q[0];
rz(-1.0780225) q[0];
x q[1];
rz(1.1649706) q[2];
sx q[2];
rz(-0.067194447) q[2];
sx q[2];
rz(-2.7213328) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48756632) q[1];
sx q[1];
rz(-1.8974202) q[1];
sx q[1];
rz(1.4766272) q[1];
rz(2.1498508) q[3];
sx q[3];
rz(-2.1175044) q[3];
sx q[3];
rz(-0.21608298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7897196) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(1.2119279) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(-0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(2.7899172) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6707014) q[0];
sx q[0];
rz(-1.1180709) q[0];
sx q[0];
rz(2.7265413) q[0];
rz(-0.36231626) q[2];
sx q[2];
rz(-10*pi/13) q[2];
sx q[2];
rz(-0.16513261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.069236) q[1];
sx q[1];
rz(-1.623739) q[1];
sx q[1];
rz(-2.3318021) q[1];
x q[2];
rz(3.0319801) q[3];
sx q[3];
rz(-1.622756) q[3];
sx q[3];
rz(-0.51652858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4984109) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(-0.19432755) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-2.5732645) q[1];
sx q[1];
rz(2.1077572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4170096) q[0];
sx q[0];
rz(-1.7010744) q[0];
sx q[0];
rz(2.2307322) q[0];
rz(-pi) q[1];
rz(0.98722234) q[2];
sx q[2];
rz(-0.7910896) q[2];
sx q[2];
rz(2.1949878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.918805) q[1];
sx q[1];
rz(-0.58090392) q[1];
sx q[1];
rz(1.3851628) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4284027) q[3];
sx q[3];
rz(-1.7367559) q[3];
sx q[3];
rz(0.99456577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0620492) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-1.6171261) q[2];
sx q[2];
rz(-0.6033069) q[2];
sx q[2];
rz(2.6848007) q[2];
rz(1.2407606) q[3];
sx q[3];
rz(-1.6374554) q[3];
sx q[3];
rz(-1.0367254) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
