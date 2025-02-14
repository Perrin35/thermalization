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
rz(4.4938402) q[1];
sx q[1];
rz(0.69753733) q[1];
sx q[1];
rz(11.868244) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1230948) q[0];
sx q[0];
rz(-2.4677489) q[0];
sx q[0];
rz(0.68742623) q[0];
x q[1];
rz(1.5641937) q[2];
sx q[2];
rz(-1.1175795) q[2];
sx q[2];
rz(-2.1389769) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.082939741) q[1];
sx q[1];
rz(-1.8720683) q[1];
sx q[1];
rz(0.79853398) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5145958) q[3];
sx q[3];
rz(-0.25543419) q[3];
sx q[3];
rz(-1.5572302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2932566) q[2];
sx q[2];
rz(-1.2245327) q[2];
sx q[2];
rz(-2.2621034) q[2];
rz(2.4073811) q[3];
sx q[3];
rz(-1.1314393) q[3];
sx q[3];
rz(2.3122299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5010928) q[0];
sx q[0];
rz(-1.9556029) q[0];
sx q[0];
rz(-2.1511141) q[0];
rz(-2.1462006) q[1];
sx q[1];
rz(-1.6973015) q[1];
sx q[1];
rz(-2.676414) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93056193) q[0];
sx q[0];
rz(-2.6829236) q[0];
sx q[0];
rz(0.42210292) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1976333) q[2];
sx q[2];
rz(-1.2612245) q[2];
sx q[2];
rz(-2.7322526) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.097660573) q[1];
sx q[1];
rz(-2.2273835) q[1];
sx q[1];
rz(0.93933479) q[1];
rz(-1.3690794) q[3];
sx q[3];
rz(-1.8276102) q[3];
sx q[3];
rz(-1.6981924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62749949) q[2];
sx q[2];
rz(-1.9072615) q[2];
sx q[2];
rz(0.17847432) q[2];
rz(2.576339) q[3];
sx q[3];
rz(-1.3924799) q[3];
sx q[3];
rz(0.55725151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31767118) q[0];
sx q[0];
rz(-1.3452106) q[0];
sx q[0];
rz(-0.70096651) q[0];
rz(2.7340381) q[1];
sx q[1];
rz(-1.2113672) q[1];
sx q[1];
rz(-2.0661381) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1140954) q[0];
sx q[0];
rz(-2.4878708) q[0];
sx q[0];
rz(2.5750676) q[0];
rz(-1.8604061) q[2];
sx q[2];
rz(-2.4649005) q[2];
sx q[2];
rz(-1.8829568) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6523931) q[1];
sx q[1];
rz(-1.4402188) q[1];
sx q[1];
rz(1.9741535) q[1];
x q[2];
rz(3.091142) q[3];
sx q[3];
rz(-1.3063432) q[3];
sx q[3];
rz(-1.8320302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.067523) q[2];
sx q[2];
rz(-1.8068376) q[2];
sx q[2];
rz(1.734181) q[2];
rz(0.47809005) q[3];
sx q[3];
rz(-0.34308386) q[3];
sx q[3];
rz(-0.34496719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57426977) q[0];
sx q[0];
rz(-1.747921) q[0];
sx q[0];
rz(-1.0465013) q[0];
rz(-0.25431713) q[1];
sx q[1];
rz(-2.3260702) q[1];
sx q[1];
rz(-0.3124803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69092551) q[0];
sx q[0];
rz(-1.1617799) q[0];
sx q[0];
rz(0.14139195) q[0];
rz(-pi) q[1];
rz(0.58860345) q[2];
sx q[2];
rz(-2.3905675) q[2];
sx q[2];
rz(2.5680755) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.65838012) q[1];
sx q[1];
rz(-1.3286031) q[1];
sx q[1];
rz(-2.8172452) q[1];
rz(-1.0572079) q[3];
sx q[3];
rz(-1.8589338) q[3];
sx q[3];
rz(0.44618928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3089402) q[2];
sx q[2];
rz(-0.71844429) q[2];
sx q[2];
rz(-0.18860513) q[2];
rz(-2.9077933) q[3];
sx q[3];
rz(-2.2457687) q[3];
sx q[3];
rz(2.5578267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69498649) q[0];
sx q[0];
rz(-1.3124895) q[0];
sx q[0];
rz(0.0083228668) q[0];
rz(-0.43909973) q[1];
sx q[1];
rz(-1.7726026) q[1];
sx q[1];
rz(-1.2459374) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9326235) q[0];
sx q[0];
rz(-1.1198988) q[0];
sx q[0];
rz(1.7843002) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1272922) q[2];
sx q[2];
rz(-2.864553) q[2];
sx q[2];
rz(0.02422457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51366561) q[1];
sx q[1];
rz(-0.49066077) q[1];
sx q[1];
rz(2.2918244) q[1];
x q[2];
rz(-2.8543618) q[3];
sx q[3];
rz(-2.2054447) q[3];
sx q[3];
rz(-1.8375481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7370854) q[2];
sx q[2];
rz(-0.046701996) q[2];
sx q[2];
rz(-1.6339634) q[2];
rz(-2.5493933) q[3];
sx q[3];
rz(-1.7630354) q[3];
sx q[3];
rz(0.73970214) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1208039) q[0];
sx q[0];
rz(-1.692166) q[0];
sx q[0];
rz(-0.26594308) q[0];
rz(-1.2159011) q[1];
sx q[1];
rz(-2.3561056) q[1];
sx q[1];
rz(1.1908092) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.196249) q[0];
sx q[0];
rz(-0.73305819) q[0];
sx q[0];
rz(-1.5747914) q[0];
rz(-1.6295678) q[2];
sx q[2];
rz(-2.588039) q[2];
sx q[2];
rz(-2.2433081) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3336181) q[1];
sx q[1];
rz(-0.75389538) q[1];
sx q[1];
rz(-3.1303166) q[1];
rz(2.9407752) q[3];
sx q[3];
rz(-2.8614223) q[3];
sx q[3];
rz(0.38129378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33209458) q[2];
sx q[2];
rz(-1.7569434) q[2];
sx q[2];
rz(-0.14796999) q[2];
rz(-2.68908) q[3];
sx q[3];
rz(-2.0544402) q[3];
sx q[3];
rz(-0.40542671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1389403) q[0];
sx q[0];
rz(-1.8754706) q[0];
sx q[0];
rz(-1.7953605) q[0];
rz(3.1345308) q[1];
sx q[1];
rz(-2.4738753) q[1];
sx q[1];
rz(-0.5074358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1333727) q[0];
sx q[0];
rz(-0.79303654) q[0];
sx q[0];
rz(-2.338468) q[0];
rz(-pi) q[1];
rz(1.981214) q[2];
sx q[2];
rz(-2.7682891) q[2];
sx q[2];
rz(0.48063916) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40312156) q[1];
sx q[1];
rz(-2.7475121) q[1];
sx q[1];
rz(0.02438633) q[1];
rz(0.0094311992) q[3];
sx q[3];
rz(-1.2985029) q[3];
sx q[3];
rz(-1.1444397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70600447) q[2];
sx q[2];
rz(-2.0861349) q[2];
sx q[2];
rz(-2.9662507) q[2];
rz(-0.38241479) q[3];
sx q[3];
rz(-1.1924815) q[3];
sx q[3];
rz(-0.46149883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14313702) q[0];
sx q[0];
rz(-2.1303506) q[0];
sx q[0];
rz(-1.829041) q[0];
rz(-0.10874272) q[1];
sx q[1];
rz(-2.5332632) q[1];
sx q[1];
rz(-0.67799062) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4512166) q[0];
sx q[0];
rz(-1.548598) q[0];
sx q[0];
rz(1.0772086) q[0];
rz(0.021804734) q[2];
sx q[2];
rz(-2.3924689) q[2];
sx q[2];
rz(-2.5072012) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7398518) q[1];
sx q[1];
rz(-1.051552) q[1];
sx q[1];
rz(-2.8080163) q[1];
rz(1.0394215) q[3];
sx q[3];
rz(-0.9262923) q[3];
sx q[3];
rz(-0.21428767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8678191) q[2];
sx q[2];
rz(-1.5742233) q[2];
sx q[2];
rz(-1.9847974) q[2];
rz(2.6371238) q[3];
sx q[3];
rz(-2.1095095) q[3];
sx q[3];
rz(2.5696136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93541637) q[0];
sx q[0];
rz(-1.529006) q[0];
sx q[0];
rz(-2.5231498) q[0];
rz(-0.84856021) q[1];
sx q[1];
rz(-2.6237374) q[1];
sx q[1];
rz(0.79577622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9600362) q[0];
sx q[0];
rz(-2.040876) q[0];
sx q[0];
rz(-3.114469) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82473849) q[2];
sx q[2];
rz(-1.5360263) q[2];
sx q[2];
rz(-0.29037133) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.578131) q[1];
sx q[1];
rz(-1.2555148) q[1];
sx q[1];
rz(2.8540846) q[1];
rz(1.9913748) q[3];
sx q[3];
rz(-1.4338067) q[3];
sx q[3];
rz(-1.5311706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15647469) q[2];
sx q[2];
rz(-1.5072344) q[2];
sx q[2];
rz(-1.9462684) q[2];
rz(-0.34504238) q[3];
sx q[3];
rz(-0.86193591) q[3];
sx q[3];
rz(-2.2752458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4320375) q[0];
sx q[0];
rz(-0.23831743) q[0];
sx q[0];
rz(-0.077614345) q[0];
rz(-1.2331102) q[1];
sx q[1];
rz(-1.0401007) q[1];
sx q[1];
rz(0.79201039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302942) q[0];
sx q[0];
rz(-2.2195507) q[0];
sx q[0];
rz(-2.4374967) q[0];
x q[1];
rz(0.77113232) q[2];
sx q[2];
rz(-2.9206616) q[2];
sx q[2];
rz(0.67459092) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35490552) q[1];
sx q[1];
rz(-0.49017683) q[1];
sx q[1];
rz(-1.9016983) q[1];
rz(-pi) q[2];
rz(2.3378793) q[3];
sx q[3];
rz(-1.4856385) q[3];
sx q[3];
rz(-2.5309895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9218257) q[2];
sx q[2];
rz(-0.91138387) q[2];
sx q[2];
rz(2.0683973) q[2];
rz(2.8940767) q[3];
sx q[3];
rz(-1.211834) q[3];
sx q[3];
rz(-3.1246576) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3794004) q[0];
sx q[0];
rz(-1.301855) q[0];
sx q[0];
rz(-0.72574885) q[0];
rz(0.87860592) q[1];
sx q[1];
rz(-2.286372) q[1];
sx q[1];
rz(1.1927037) q[1];
rz(-1.1996126) q[2];
sx q[2];
rz(-0.51358583) q[2];
sx q[2];
rz(-2.5805342) q[2];
rz(2.2976919) q[3];
sx q[3];
rz(-1.5283199) q[3];
sx q[3];
rz(0.98315317) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
